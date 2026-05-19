# Diagnostics — source layer

How diagnostic blocks are implemented in the C++ source, the HDF5 output schema, integration with happi, and how to add a new diagnostic. The output-half of the codebase.

## Table of contents

1. `Diagnostic` class hierarchy
2. Per-step orchestration — `runAllDiags`
3. `DiagnosticScalar` — implementation
4. `DiagnosticFields` — implementation
5. `DiagnosticProbe` — implementation
6. `DiagnosticParticleBinning` — implementation
7. `DiagnosticTrackParticles`, `DiagnosticScreen`, others
8. HDF5 layout and schema versioning
9. happi integration
10. Adding a new diagnostic
11. Performance considerations
12. Common pitfalls

---

## 1. `Diagnostic` class hierarchy

`src/Diagnostic/` contains the diagnostic implementation. Base class:

```cpp
class Diagnostic {
public:
    int every;                                   // output cadence (0 = disabled)
    int flush_every;                             // I/O buffer flush cadence
    std::string filename;

    virtual void run(int itime, Patch* patch, ...) = 0;    // per-patch per-step work
    virtual void writeAll(int itime, ...) = 0;             // global write
    virtual void closeFile() = 0;
};
```

Concrete derivations:

| File | Class | Namelist block |
|---|---|---|
| `DiagnosticScalar.cpp` | `DiagnosticScalar` | `DiagScalar` |
| `DiagnosticFields.cpp` | `DiagnosticFields` | `DiagFields` |
| `DiagnosticProbe.cpp` | `DiagnosticProbe` | `DiagProbe` |
| `DiagnosticParticleBinning.cpp` | `DiagnosticParticleBinning` | `DiagParticleBinning` |
| `DiagnosticTrackParticles.cpp` | `DiagnosticTrackParticles` | `DiagTrackParticles` |
| `DiagnosticScreen.cpp` | `DiagnosticScreen` | `DiagScreen` |
| `DiagnosticRadiationSpectrum.cpp` | `DiagnosticRadiationSpectrum` | `DiagRadiationSpectrum` |
| `DiagnosticNewParticles.cpp` | `DiagnosticNewParticles` | `DiagNewParticles` |
| `DiagnosticPerformances.cpp` | `DiagnosticPerformances` | `DiagPerformances` |

Each namelist `Diag*` block creates one instance, stored in `VectorPatch::globalDiags` (for global-state diagnostics like Scalar) or `Patch::localDiags` (for per-patch diagnostics like Fields, Probe, ParticleBinning).

## 2. Per-step orchestration — `runAllDiags`

`VectorPatch::runAllDiags(int itime, ...)` is called once per timestep from the main loop:

```cpp
void VectorPatch::runAllDiags(int itime, Params& params, SmileiMPI* smpi) {
    // Local (per-patch) diagnostics
    #pragma omp for
    for (size_t ipatch = 0; ipatch < patches_.size(); ++ipatch) {
        for (auto* diag : patches_[ipatch]->localDiags) {
            if (itime % diag->every == 0) {
                diag->run(itime, patches_[ipatch], ...);
            }
        }
    }

    // Global diagnostics (aggregate across patches and ranks)
    #pragma omp single
    {
        for (auto* diag : globalDiags) {
            if (itime % diag->every == 0) {
                diag->runAll(itime, smpi, ...);     // collective MPI ops here
            }
        }
    }

    // File flushing
    #pragma omp single
    {
        for (auto* diag : globalDiags) {
            if (itime % diag->flush_every == 0) diag->flush();
        }
    }
}
```

**Key distinction:**

- **Per-patch diagnostics** (`DiagnosticFields`, `DiagnosticProbe`, `DiagnosticParticleBinning`) run inside the OMP-parallel-over-patches loop. Each patch writes its own slice independently. HDF5 parallel I/O assembles the global view.

- **Global diagnostics** (`DiagnosticScalar`) require cross-patch + cross-rank reduction. They run in `#pragma omp single` regions with MPI collectives.

## 3. `DiagnosticScalar` — implementation

`DiagnosticScalar.cpp` computes global integrated quantities at each output time:

```cpp
void DiagnosticScalar::runAll(int itime, SmileiMPI* smpi, ...) {
    // 1. Per-rank reduction
    double Ukin_rank = 0.0;
    double Uelm_rank = 0.0;
    for (each patch on this rank) {
        for (each species on this patch) {
            for (each particle in this species) {
                Ukin_rank += w[i] * (gamma[i] - 1.0) * m_e;
            }
        }
        Uelm_rank += integrate_EM_energy(patch->EMfields);
    }

    // 2. MPI reduction
    double Ukin_total, Uelm_total;
    MPI_Reduce(&Ukin_rank, &Ukin_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Uelm_rank, &Uelm_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // 3. Master rank writes to Scalars.txt
    if (smpi->isMaster()) {
        scalars_file << itime << "  " << Ukin_total << "  " << Uelm_total << "  ...\n";
    }
}
```

Output is a flat text file (`Scalars.txt`) — one row per output time, columns per scalar. Cheap to compute, cheap to write. **Always enabled** in practice.

**Adding a new scalar**: edit `DiagnosticScalar::runAll` to compute the new quantity, add it to the column header in the file, update happi's `Scalar` accessor to know the new name.

## 4. `DiagnosticFields` — implementation

`DiagnosticFields.cpp` writes full-grid field snapshots:

```cpp
void DiagnosticFields::run(int itime, Patch* patch, ...) {
    // Each patch writes its slab into a parallel HDF5 file
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

**Parallel HDF5**: each patch writes to its hyperslab of the global dataset. The HDF5 file is one per `DiagFields` block (e.g., `Fields0.h5`), with one dataset per field per output timestep.

**Subgrid**: if `subgrid = [slice(...), slice(...)]` is set, the patch's slab is sub-sliced before write. Reduces output volume.

**Time-averaging**: averages over `time_average` timesteps centered on the output time. Each patch holds a per-cell accumulator; reset after each output.

**File size grows quickly**: 2D 1000×500 with 4 fields, every=200, 1000 frames = 16 GB. Most expensive diagnostic.

## 5. `DiagnosticProbe` — implementation

`DiagnosticProbe.cpp` samples fields at user-defined geometric points (0D, 1D, 2D, 3D shapes):

```cpp
void DiagnosticProbe::run(int itime, Patch* patch, ...) {
    for (size_t ipoint : probe_points_in_this_patch) {
        double x = probe_position[0][ipoint];
        double y = probe_position[1][ipoint];

        // Interpolate each requested field at (x, y)
        for (auto& field_name : probed_fields) {
            double value = interpolate_field(patch->EMfields, field_name, x, y);
            output_buffer[ipoint][field_name] = value;
        }
    }
    // Write only the points owned by this patch
}
```

Probe points are statically distributed at simulation start: each point is assigned to the patch containing it, never moves. The HDF5 output is one record per probe point per output time.

**Cheap diagnostic**: output volume is `n_points × n_fields × n_times`, much smaller than full-grid fields. Aggressive cadences (`every = 10`) are practical.

**Geometry**:

- 0 corners → 1 point at `origin`
- 1 corner → 1D line of `number[0]` points from `origin` to `corners[0]`
- 2 corners → 2D plane (parallelogram with `origin` + the two corners as edges), `number[0] × number[1]` points
- 3 corners → 3D box

## 6. `DiagnosticParticleBinning` — implementation

The flexible distribution-projection diagnostic. Per-patch:

```cpp
void DiagnosticParticleBinning::run(int itime, Patch* patch, ...) {
    Species* species = patch->vecSpecies[species_index];

    for (size_t ip = 0; ip < species->particles->size(); ++ip) {
        if (filter && !filter(species->particles, ip)) continue;

        // Compute bin indices for this particle
        int bin_indices[MAX_AXES];
        for (int a = 0; a < n_axes; ++a) {
            double value_a = compute_axis_value(axes[a], species->particles, ip);
            bin_indices[a] = (value_a - axes[a].min) / axes[a].width * axes[a].n_bins;
            if (bin_indices[a] < 0 || bin_indices[a] >= axes[a].n_bins) continue;
        }
        int flat_index = flatten(bin_indices, axes);

        // Compute deposited quantity
        double deposit = compute_deposit(deposited_quantity, species->particles, ip);

        // Atomic add for thread safety
        #pragma omp atomic
        histogram[flat_index] += deposit;
    }
}
```

After all patches contribute, `MPI_Allreduce` aggregates the histogram across ranks, and the master writes to HDF5.

**Filter callable**: a Python lambda from the namelist. Called via the C-Python bridge — cheap if vectorized (operates on arrays at once), expensive if per-particle. Smilei tries to detect vectorized lambdas; failure mode is per-particle Python calls.

**`deposited_quantity`** strings like `"weight"`, `"weight_charge"`, `"weight_ekin"` are parsed at construction into a small enum. The inner loop dispatches on the enum.

**Axes**: each axis has `(name, min, max, n_bins, logscale)`. The `name` is parsed into a function pointer at construction (e.g., `"x"` → `[](particles, ip) { return position[0][ip]; }`). The `logscale` flag rewires the binning to `log(value)`.

## 7. `DiagnosticTrackParticles`, `DiagnosticScreen`, others

**`DiagnosticTrackParticles`**: per-particle time series. Filter selects which particles to track; for each tracked particle, the state (position, momentum, etc.) is written each `every` steps. Memory and disk grow with `tracked_count × n_attributes × n_output_times`. Always narrow the filter.

**`DiagnosticScreen`**: time-integrated flux through a virtual surface. Implemented as per-step accumulation of particle crossings, weighted by the deposited quantity. Output schema is the same as `DiagnosticParticleBinning`.

**`DiagnosticRadiationSpectrum`**: per-relativistic-electron synchrotron-like emission spectrum. Computes the radiated power per frequency bin per output time. Requires `reference_angular_frequency_SI`.

**`DiagnosticNewParticles`**: records every particle-creation event (ionization, pair generation, injection). High-volume diagnostic; use sparingly.

**`DiagnosticPerformances`**: per-patch timing. Records compute time, particle count, memory per patch. Essential for diagnosing load imbalance on large clusters.

## 8. HDF5 layout and schema versioning

All diagnostics write HDF5 files in the simulation output directory. Each `Diag*` block produces one or more files indexed by block order:

| Block | File(s) |
|---|---|
| `DiagScalar` | `Scalars.txt` (text!) |
| `DiagFields[0]` | `Fields0.h5` |
| `DiagFields[1]` | `Fields1.h5` |
| `DiagProbe[0]` | `Probes0.h5` |
| `DiagParticleBinning[0]` | `ParticleBinning0.h5` |
| `DiagTrackParticles[species]` | `TrackParticles_<species>.h5` |

**HDF5 schema** is versioned. The schema version is stored as a root-level attribute:

```
/SchemaVersion = "1.5"
```

When the schema changes (new datasets, new metadata, restructured layouts), the version is bumped. happi checks the schema version on open and refuses incompatible versions.

**Backward compatibility**: when adding a field to an existing HDF5 layout, prefer adding a new dataset rather than renaming or restructuring. Tools (happi, external analyzers) parse the existing layout; restructuring breaks them.

**Schema sources**:
- `src/Diagnostic/DiagnosticScalar.cpp::createFile` — Scalar header
- `src/Diagnostic/DiagnosticFields.cpp::createFile` — Fields HDF5 structure
- `happi/_Diagnostics/*.py` — Python-side schema knowledge

If you change the schema, update both sides and bump the version. happi tests should catch incompatibilities.

## 9. happi integration

`happi/` is a Python module that wraps HDF5 file access. It is built separately via `make happi` and installed to user-site Python.

**Per-diagnostic accessor classes**: `happi/_Diagnostics/Scalar.py`, `Fields.py`, `Probe.py`, etc. Each class knows the HDF5 layout for its diagnostic type and exposes a uniform Python interface (`.getData()`, `.plot()`, `.slide()`, etc.).

**Adding a new accessor**: when adding a new diagnostic class on the C++ side, create a corresponding `happi/_Diagnostics/MyNewDiag.py` that reads the new HDF5 structure. Then expose `S.MyNewDiag(...)` in `happi/Happi.py`.

**Pint integration**: optional Pint-based unit conversion. Activated by `happi.Open(..., pint=True)`. Each diagnostic class knows its raw normalized units and converts to the user-requested SI representation.

## 10. Adding a new diagnostic

To add a new diagnostic (e.g., a per-cell histogram of mean particle energy):

1. **Create** `src/Diagnostic/DiagnosticMyNew.{h,cpp}` deriving from `Diagnostic`.

2. **Implement** at minimum:
   - `DiagnosticMyNew(params, block_data)` — constructor parsing namelist arguments.
   - `run(int itime, Patch* patch, ...)` — per-patch per-step work.
   - `writeAll(int itime, ...)` — per-rank or global write.
   - `createFile()` — HDF5 file creation with schema.

3. **Register** in `src/Diagnostic/DiagnosticFactory.h`:
   ```cpp
   else if (block_name == "DiagMyNew") {
       return new DiagnosticMyNew(params, block_data);
   }
   ```

4. **Parse** the namelist block in `src/Params/Params.cpp::parseDiagnosticBlocks(...)`.

5. **Add** the matching happi accessor in `happi/_Diagnostics/MyNew.py`.

6. **Bump** the HDF5 schema version.

7. **Document** in `doc/Sphinx/Source/namelist.rst`.

8. **Test** — `validation/tst_diag_mynew.py` exercising the new diagnostic.

For a niche analysis that won't be merged upstream, an alternative is to use existing diagnostics (`DiagParticleBinning` with custom axes, `DiagNewParticles` with custom filters) and post-process in Python.

## 11. Performance considerations

Diagnostics are typically a small fraction (~5–10%) of simulation cost — but can blow up if:

- **`every` too small**: `every = 1` on `DiagFields` is the canonical disk-filler.
- **`time_average` too large**: per-cell buffers grow proportionally.
- **`DiagTrackParticles` without filter**: tracks every particle, terabytes of output.
- **`DiagNewParticles` in heavy ionization runs**: every event recorded — millions per timestep possible.
- **`DiagParticleBinning` with too many axes**: 3-axis binning at `100 × 100 × 100` bins = 1M cells, fine; 4-axis or finer bins exceed memory.

**HDF5 I/O cost** scales with output volume and write frequency. Parallel HDF5 with collective writes is more efficient than per-rank serial writes; Smilei uses collective by default for `DiagFields`.

**Flush behavior**: `flush_every >> every` reduces I/O overhead by buffering multiple output steps before flushing.

**Profiling diagnostic cost**: enable `DiagPerformances` and check the per-patch breakdown. If a diagnostic dominates, reduce `every` or add `subgrid`.

## 12. Common pitfalls

**Schema drift**: adding a new dataset to `DiagFields` without updating happi's `Fields.py` → user sees the new field in the HDF5 file but happi can't read it.

**Race condition on histogram accumulation**: `histogram[flat_index] += deposit` from multiple threads requires `#pragma omp atomic` or per-thread histograms reduced after the parallel section. Omitting the atomic produces silent off-by-some errors.

**Filter callable with per-particle Python**: a non-vectorized lambda (`lambda p: p.x[ip] > 20`) is millions of Python calls per step. Always use vectorized form (`lambda p: p.x > 20`).

**HDF5 parallel I/O without collective mode**: some HDF5 builds default to independent I/O, which is much slower than collective. Smilei sets collective explicitly; if you add new I/O, follow the pattern.

**`every` set to 0 silently disables**: a user-facing edge case — if a parsed integer evaluates to 0, the diagnostic is silently disabled. Helpful warning would be to emit at namelist-parse time.

**Forgetting `time_average` reset**: per-patch buffers must be cleared after each output. Failure to clear → accumulator runs forever and outputs nonsense.

**Per-rank vs global state confusion**: `DiagnosticScalar` aggregates across ranks via MPI_Reduce. Adding a new scalar that doesn't reduce produces per-rank-master-output that's wrong.

**Per-patch diagnostic accumulators surviving load balance**: when patches migrate ranks, their diagnostic state must migrate with them. Check `DiagnosticFields::serialize`/`deserialize` patterns when adding new patch-local state.

**Probe points outside the box**: silently dropped at parse time. Not a runtime error; the diagnostic has fewer points than requested. happi displays the actual count.

**`getDataset` from HDF5 with wrong dtype**: HDF5 won't auto-cast. If you write `float` but read as `double`, you get garbage. Match types between writer and reader.

**Memory growth in `DiagTrackParticles`**: each tracked particle's state per output time accumulates in memory before flush. With `every = 100` and `flush_every = 100000`, that's 1000 frames buffered — substantial RAM. Match `flush_every` to available memory.
