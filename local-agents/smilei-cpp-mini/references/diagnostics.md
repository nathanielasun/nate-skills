# Diagnostics — source layer (mini tier)

How diagnostic blocks work in C++, HDF5 schema, adding new diagnostics.

## `Diagnostic` class hierarchy

`src/Diagnostic/`. Base:

```cpp
class Diagnostic {
    int every, flush_every;
    std::string filename;
    virtual void run(int itime, Patch* patch, ...) = 0;
    virtual void writeAll(int itime, ...) = 0;
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

Each namelist block creates one instance, stored in `VectorPatch::globalDiags` (for `Scalar`) or `Patch::localDiags` (per-patch ones like `Fields`, `Probe`, `ParticleBinning`).

## Per-step orchestration

`VectorPatch::runAllDiags(int itime, ...)`:

```cpp
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
```

Per-patch: runs inside OMP-parallel loop, parallel HDF5 writes its slab. Global: requires cross-rank reduction (MPI), in `omp single`.

## `DiagnosticScalar`

Per-rank reduction → `MPI_Reduce` → master writes flat text `Scalars.txt`:

```cpp
double Ukin_rank = 0.0, Uelm_rank = 0.0;
for (each patch on rank) {
    for (each species) for (each particle) Ukin_rank += w * (gamma - 1) * m_e;
    Uelm_rank += integrate_EM_energy(patch->EMfields);
}
MPI_Reduce(&Ukin_rank, &Ukin_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
if (smpi->isMaster()) scalars_file << itime << "  " << Ukin_total << "\n";
```

Adding scalar: edit `runAll`, add column header, update happi's `Scalar` accessor.

## `DiagnosticFields`

Per-patch HDF5 hyperslab write per output time:

```cpp
for (auto& field_name : fields_to_output) {
    Field* field = patch->EMfields->getField(field_name);
    if (time_average > 1) {
        accumulate_into_buffer(field, this->buffer);
        if (itime + time_average/2 == output_itime) {
            write_buffer_to_hdf5(this->buffer);
            clear_buffer();
        }
    } else write_field_to_hdf5(field, ...);
}
```

One file per block (`Fields0.h5`, `Fields1.h5`, ...). Most expensive diagnostic — 2D 1000×500 × 4 fields × every=200 × 1000 frames = 16 GB.

## `DiagnosticProbe`

Samples fields at user-defined geometric points (0D/1D/2D/3D):

```cpp
for (size_t ipoint : probe_points_in_this_patch) {
    double x = probe_position[0][ipoint], y = probe_position[1][ipoint];
    for (auto& field_name : probed_fields) {
        double value = interpolate_field(patch->EMfields, field_name, x, y);
        output_buffer[ipoint][field_name] = value;
    }
}
```

Probe points statically distributed at startup (each assigned to containing patch, never moves). Cheap — output ~`n_points × n_fields × n_times`.

## `DiagnosticParticleBinning`

Histogram with `#pragma omp atomic`:

```cpp
for (size_t ip = 0; ip < species->particles->size(); ++ip) {
    if (filter && !filter(species->particles, ip)) continue;

    int bin_indices[MAX_AXES];
    for (int a = 0; a < n_axes; ++a) {
        double value = compute_axis_value(axes[a], species->particles, ip);
        bin_indices[a] = (value - axes[a].min) / axes[a].width * axes[a].n_bins;
        if (bin_indices[a] < 0 || bin_indices[a] >= axes[a].n_bins) continue;
    }
    int flat_index = flatten(bin_indices, axes);
    double deposit = compute_deposit(deposited_quantity, species->particles, ip);

    #pragma omp atomic
    histogram[flat_index] += deposit;
}
```

Then `MPI_Allreduce` across ranks. Master writes HDF5.

`deposited_quantity` strings parsed to enum at construction. Axes: `name → function pointer` at construction. Filter is Python lambda via C-Python bridge — cheap if vectorized, expensive if per-particle.

## Other diagnostics

- `DiagnosticTrackParticles` — per-particle time series, filter selects which. Memory/disk = tracked × attributes × times. **Always narrow filter.**
- `DiagnosticScreen` — time-integrated flux through virtual surface.
- `DiagnosticRadiationSpectrum` — synchrotron emission per relativistic electron. Requires `reference_angular_frequency_SI`.
- `DiagnosticNewParticles` — every particle-creation event (ionization, pair, injection). High-volume.
- `DiagnosticPerformances` — per-patch timing. Essential for load-imbalance diagnosis.

## HDF5 layout

| Block | File |
|---|---|
| `DiagScalar` | `Scalars.txt` (text!) |
| `DiagFields[i]` | `Fields<i>.h5` |
| `DiagProbe[i]` | `Probes<i>.h5` |
| `DiagParticleBinning[i]` | `ParticleBinning<i>.h5` |
| `DiagTrackParticles[species]` | `TrackParticles_<species>.h5` |

Schema versioned: root attribute `SchemaVersion = "1.5"`. happi checks on open. **Bump version when changing layout.**

Prefer adding new datasets over renaming/restructuring — tools depend on existing layout.

## happi integration

`happi/_Diagnostics/{Scalar,Fields,Probe,...}.py` — one accessor class per diagnostic type. Each knows HDF5 layout, exposes `.getData()`, `.plot()`, `.slide()`.

Adding new diagnostic: create matching `happi/_Diagnostics/MyNew.py`, expose `S.MyNew(...)` in `happi/Happi.py`.

Pint integration: optional via `happi.Open(..., pint=True)`. Each class knows normalized units, converts to SI.

## Adding a new diagnostic

1. Create `src/Diagnostic/DiagnosticMyNew.{h,cpp}` deriving from `Diagnostic`.
2. Implement constructor (parse namelist args), `run`, `writeAll`, `createFile`.
3. Register in `DiagnosticFactory.h`.
4. Parse namelist block in `Params.cpp::parseDiagnosticBlocks(...)`.
5. Add `happi/_Diagnostics/MyNew.py`.
6. Bump HDF5 schema version.
7. Document and test.

For niche analyses unlikely to be merged upstream: use existing diagnostics (`DiagParticleBinning` with custom axes, `DiagNewParticles` with filters) + post-process in Python.

## Performance considerations

Diagnostics typically 5–10% of step time; can blow up if:
- `every = 1` on `DiagFields` — disk-filler.
- `time_average > 1` — per-cell buffers grow.
- `DiagTrackParticles` without filter — terabytes.
- `DiagNewParticles` in heavy ionization — millions of events/step.
- `DiagParticleBinning` with 4-axis fine bins — exceeds memory.

`flush_every >> every` reduces I/O overhead via buffering. Parallel HDF5 collective writes more efficient than per-rank serial; Smilei uses collective by default.

Profile via `DiagPerformances`. Heavy diagnostic → reduce `every` or add `subgrid`.

## Common pitfalls

- **Schema drift** — new C++ dataset without happi update → user sees data, can't read.
- **Missing OMP atomic** on histogram in `DiagnosticParticleBinning` → silent off-by-some.
- **Per-particle Python filter** — `lambda p: p.x[ip] > 20` (non-vectorized) is millions of Python calls/step. Vectorize: `lambda p: p.x > 20`.
- **HDF5 non-collective I/O** — much slower than collective. Smilei defaults to collective; new I/O should match.
- **`every = 0`** silently disables — helpful warning at parse would be improvement.
- **Forgetting time_average reset** — accumulator runs forever, outputs nonsense.
- **Per-rank vs global state** — `Scalar` aggregates via `MPI_Reduce`. New scalar without reduction → wrong output.
- **Per-patch state lost on load balance** — diagnostic state must migrate with patch. Check `serialize`/`deserialize`.
- **Probe points outside box** — silently dropped at parse, fewer points than requested.
- **HDF5 dtype mismatch** — write float, read double → garbage. Match types.
- **`DiagTrackParticles` memory growth** — each tracked particle per output buffered. Match `flush_every` to RAM.
