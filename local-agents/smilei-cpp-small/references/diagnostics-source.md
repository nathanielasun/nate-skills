# Diagnostics — source layer (small tier)

How diagnostic blocks are implemented in C++, HDF5 schema, happi integration, adding a new diagnostic.

## `Diagnostic` class hierarchy

`src/Diagnostic/`. Base:

```cpp
class Diagnostic {
    int every, flush_every;
    std::string filename;
    virtual void run(int itime, Patch* patch, ...) = 0;
    virtual void writeAll(int itime, ...) = 0;
    virtual void closeFile() = 0;
};
```

Concrete:

| Class | Namelist block |
|---|---|
| `DiagnosticScalar` | `DiagScalar` |
| `DiagnosticFields` | `DiagFields` |
| `DiagnosticProbe` | `DiagProbe` |
| `DiagnosticParticleBinning` | `DiagParticleBinning` |
| `DiagnosticTrackParticles` | `DiagTrackParticles` |
| `DiagnosticScreen` | `DiagScreen` |
| `DiagnosticRadiationSpectrum` | `DiagRadiationSpectrum` |
| `DiagnosticNewParticles` | `DiagNewParticles` |
| `DiagnosticPerformances` | `DiagPerformances` |

Each namelist block creates one instance, stored in `VectorPatch::globalDiags` (global-state diagnostics like Scalar) or `Patch::localDiags` (per-patch diagnostics like Fields, Probe, ParticleBinning).

## Per-step orchestration

`VectorPatch::runAllDiags(int itime, ...)`:

```cpp
void VectorPatch::runAllDiags(int itime, Params& params, SmileiMPI* smpi) {
    // Per-patch diagnostics
    #pragma omp for
    for (size_t ipatch = 0; ipatch < patches_.size(); ++ipatch) {
        for (auto* diag : patches_[ipatch]->localDiags) {
            if (itime % diag->every == 0) diag->run(itime, patches_[ipatch], ...);
        }
    }

    // Global diagnostics (MPI collectives)
    #pragma omp single
    for (auto* diag : globalDiags) {
        if (itime % diag->every == 0) diag->runAll(itime, smpi, ...);
    }

    // File flushing
    #pragma omp single
    for (auto* diag : globalDiags) if (itime % diag->flush_every == 0) diag->flush();
}
```

**Per-patch diagnostics** (Fields, Probe, ParticleBinning) run inside OMP-parallel loop. Each patch writes its slice independently via parallel HDF5.

**Global diagnostics** (Scalar) require cross-patch + cross-rank reduction → `omp single` + MPI collectives.

## `DiagnosticScalar`

```cpp
void DiagnosticScalar::runAll(int itime, SmileiMPI* smpi, ...) {
    double Ukin_rank = 0.0, Uelm_rank = 0.0;
    for (each patch) {
        for (each species) for (each particle i) Ukin_rank += w[i] * (gamma[i] - 1.0) * m_e;
        Uelm_rank += integrate_EM_energy(patch->EMfields);
    }
    double Ukin_total, Uelm_total;
    MPI_Reduce(&Ukin_rank, &Ukin_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Uelm_rank, &Uelm_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (smpi->isMaster()) scalars_file << itime << "  " << Ukin_total << "  " << Uelm_total << "\n";
}
```

Output: flat text `Scalars.txt` — row per output time, column per scalar. Cheap.

Adding a scalar: edit `runAll` to compute, add to column header, update happi's `Scalar` accessor.

## `DiagnosticFields`

```cpp
void DiagnosticFields::run(int itime, Patch* patch, ...) {
    for (auto& field_name : fields_to_output) {
        Field* field = patch->EMfields->getField(field_name);
        if (time_average > 1) {
            accumulate_into_buffer(field, this->buffer);
            if (itime + time_average/2 == output_itime) {
                write_buffer_to_hdf5(this->buffer);
                clear_buffer();
            }
        } else {
            write_field_to_hdf5(field, ...);
        }
    }
}
```

Parallel HDF5 — each patch writes its hyperslab. One file per `DiagFields` block (`Fields0.h5`, `Fields1.h5`, ...). One dataset per field per output timestep.

**Subgrid**: sub-slices before write. Reduces volume.
**time-averaging**: per-cell accumulator per patch, reset after each output.

Most expensive diagnostic. File size: 2D 1000×500 × 4 fields × every=200 × 1000 frames = 16 GB.

## `DiagnosticProbe`

Samples fields at user-defined geometric points. Cheap (output ~`n_points × n_fields × n_times`).

```cpp
void DiagnosticProbe::run(int itime, Patch* patch, ...) {
    for (size_t ipoint : probe_points_in_this_patch) {
        double x = probe_position[0][ipoint];
        double y = probe_position[1][ipoint];
        for (auto& field_name : probed_fields) {
            double value = interpolate_field(patch->EMfields, field_name, x, y);
            output_buffer[ipoint][field_name] = value;
        }
    }
}
```

Probe points statically distributed at startup — each point assigned to the patch containing it, never moves.

Geometry: 0/1/2/3 corners → 0D/1D/2D/3D probe (line, plane, box).

## `DiagnosticParticleBinning`

Flexible distribution projection.

```cpp
void DiagnosticParticleBinning::run(int itime, Patch* patch, ...) {
    Species* species = patch->vecSpecies[species_index];
    for (size_t ip = 0; ip < species->particles->size(); ++ip) {
        if (filter && !filter(species->particles, ip)) continue;

        int bin_indices[MAX_AXES];
        for (int a = 0; a < n_axes; ++a) {
            double value_a = compute_axis_value(axes[a], species->particles, ip);
            bin_indices[a] = (value_a - axes[a].min) / axes[a].width * axes[a].n_bins;
            if (bin_indices[a] < 0 || bin_indices[a] >= axes[a].n_bins) continue;
        }
        int flat_index = flatten(bin_indices, axes);
        double deposit = compute_deposit(deposited_quantity, species->particles, ip);

        #pragma omp atomic
        histogram[flat_index] += deposit;
    }
}
```

After all patches contribute, `MPI_Allreduce` aggregates across ranks, master writes HDF5.

**`deposited_quantity`** strings (`"weight"`, `"weight_ekin"`, `"weight_charge"`, ...) parsed at construction into enum; inner loop dispatches on enum.

**Axes**: `(name, min, max, n_bins, logscale)`. Name parsed to function pointer at construction.

**Filter callable**: Python lambda from namelist via C-Python bridge. Cheap if vectorized; expensive if per-particle.

## Other diagnostics

- **`DiagnosticTrackParticles`**: per-particle time series. Filter selects which. Memory and disk grow with tracked count × attributes × output times. Always narrow filter.
- **`DiagnosticScreen`**: time-integrated flux through virtual surface. Per-step accumulation of crossings.
- **`DiagnosticRadiationSpectrum`**: synchrotron-like emission per relativistic electron. Requires `reference_angular_frequency_SI`.
- **`DiagnosticNewParticles`**: every particle-creation event (ionization, pair generation, injection). High-volume.
- **`DiagnosticPerformances`**: per-patch timing/particle-count/memory. Essential for load-imbalance diagnosis.

## HDF5 layout

Each block produces files indexed by block order:

| Block | File |
|---|---|
| `DiagScalar` | `Scalars.txt` (text!) |
| `DiagFields[i]` | `Fields<i>.h5` |
| `DiagProbe[i]` | `Probes<i>.h5` |
| `DiagParticleBinning[i]` | `ParticleBinning<i>.h5` |
| `DiagTrackParticles[species]` | `TrackParticles_<species>.h5` |

Schema versioned via root-level attribute `SchemaVersion = "1.5"`. happi checks on open and refuses incompatible. Bump version when changing layout.

When adding a field to an existing HDF5 layout, prefer adding a new dataset rather than renaming/restructuring — tools (happi, external analyzers) depend on existing layout.

## happi integration

`happi/_Diagnostics/{Scalar.py, Fields.py, Probe.py, ...}` — one accessor class per diagnostic type. Each knows HDF5 layout, exposes uniform Python interface (`.getData()`, `.plot()`, `.slide()`).

Adding new diagnostic: create matching `happi/_Diagnostics/MyNew.py`, expose `S.MyNew(...)` in `happi/Happi.py`.

Pint integration: optional via `happi.Open(..., pint=True)`. Each diagnostic class knows its raw normalized units, converts to user-requested SI.

## Adding a new diagnostic

1. Create `src/Diagnostic/DiagnosticMyNew.{h,cpp}` deriving from `Diagnostic`.
2. Implement constructor (parses namelist args), `run`, `writeAll`, `createFile`.
3. Register in `DiagnosticFactory.h`.
4. Parse namelist block in `Params.cpp::parseDiagnosticBlocks(...)`.
5. Add `happi/_Diagnostics/MyNew.py` accessor.
6. Bump HDF5 schema version.
7. Document and test in `validation/`.

For niche analyses unlikely to be merged upstream, alternative: use existing diagnostics (`DiagParticleBinning` with custom axes, `DiagNewParticles` with custom filters) + post-process in Python.

## Performance considerations

Diagnostics typically ~5–10% of step time, but can blow up:
- `every = 1` on `DiagFields` — disk-filler.
- `time_average > 1` — per-cell buffers grow.
- `DiagTrackParticles` without filter — terabytes.
- `DiagNewParticles` in heavy ionization runs — millions of events per step.
- `DiagParticleBinning` with too many axes — `100³` = 1M cells (fine); 4-axis exceeds memory.

HDF5 I/O scales with output volume × write frequency. Parallel HDF5 collective writes more efficient than per-rank serial; Smilei uses collective by default.

`flush_every >> every` reduces I/O overhead via buffering.

Profile via `DiagPerformances`. Heavy diagnostic → reduce `every` or add `subgrid`.

## Common pitfalls

- **Schema drift**: new dataset on C++ side without happi update → user sees field in HDF5 but happi can't read it.
- **Missing OMP atomic on histogram**: silent off-by-some errors in `DiagnosticParticleBinning`.
- **Per-particle Python filter**: lambda `lambda p: p.x[ip] > 20` (non-vectorized) is millions of Python calls per step. Use vectorized `lambda p: p.x > 20`.
- **HDF5 non-collective I/O**: some builds default to independent — much slower than collective. Smilei sets collective explicitly; new I/O code should match.
- **`every = 0` silently disables**: helpful warning at parse would be improvement.
- **Forgetting `time_average` reset**: per-patch buffers must clear after output. Else accumulator runs forever.
- **Per-rank vs global state confusion**: `DiagnosticScalar` aggregates via `MPI_Reduce`. New scalar that doesn't reduce → wrong per-rank-master output.
- **Per-patch state lost on load balance**: when patches migrate, diagnostic state must migrate with them. Check `DiagnosticFields::serialize`/`deserialize` patterns.
- **Probe points outside box**: silently dropped at parse. Not runtime error; diagnostic has fewer points than requested.
- **HDF5 dtype mismatch**: write float, read double → garbage. Match types between writer and reader.
- **Memory growth in `DiagTrackParticles`**: each tracked particle per output buffered before flush. With `every=100`, `flush_every=100000`, 1000 frames buffered — substantial RAM.
