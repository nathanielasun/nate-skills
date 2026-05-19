# Smilei architecture (small tier)

Patches, Hilbert ordering, MPI/OpenMP/GPU layers, build system, debugging. Compressed from frontier.

## Patches — the parallelism unit

A **patch** is the unit of parallel decomposition. Each patch holds:
- A contiguous block of grid cells (the interior)
- A halo of ghost cells for stencil operations
- All particles whose positions are inside the interior
- All field arrays for that subdomain

**Patches are not subdomains of MPI ranks.** Many patches per rank; the assignment is dynamic (load balancing). Per-rank patch list = OpenMP parallelism unit.

Patches uniform in size, set by `Main.number_of_patches` and `Main.grid_length`. Total patches: power of 2 per direction, ≥ 4 × MPI ranks. Typical size: 16–64 cells per direction.

## `Patch` and `VectorPatch`

`src/Patch/`:

- **`Patch`** (`Patch.h`): per-patch state — `ElectroMagn* EMfields`, `std::vector<Species*> vecSpecies`, particle BCs, local diagnostics, per-patch RNG. Pure data + per-patch methods.
- **`VectorPatch`** (`VectorPatch.h`): owns the patches of the current rank (`std::vector<Patch*>`). Main-loop methods (`dynamics`, `solveMaxwell`, `runAllDiags`) live here, dispatching OMP-parallel over `patch->...`.

```cpp
void VectorPatch::dynamics(Params& params, SmileiMPI* smpi, ...) {
    #pragma omp for schedule(runtime)
    for (size_t ipatch = 0; ipatch < patches_.size(); ipatch++) {
        patches_[ipatch]->dynamics(params, smpi, ...);
    }
    sumDensitiesAndCurrents(...);              // inter-patch communication
}
```

`schedule(runtime)` lets users tune via `OMP_SCHEDULE` — important for load imbalance.

## Hilbert curve and load balancing

Patches ordered along a Hilbert curve (`SmileiMPI::computeHilbertCurve()`), not row-major. Keeps spatially nearby patches close in the linear ordering — preserves locality during rank assignment.

`LoadBalancing` namelist block triggers redistribution:

1. Each patch reports load (`cell_load × n_cells + frequency_load × sum_particle_counts`).
2. Hilbert curve traversed to assign patches such that per-rank loads balance.
3. Migrating patches move with their particles, fields, diagnostic state via MPI point-to-point.

Migration is per-patch (no per-cell). Cost grows with particle count per patch. Tune `LoadBalancing(every=N)` so migration cost amortizes.

## MPI layer

`src/SmileiMPI/SmileiMPI.h`. One object per rank:
- The `MPI_Comm`
- Rank index, total ranks
- Hilbert table
- Send/recv buffers

Key collectives:
- `exchangeFields(...)` — ghost-cell exchange across patches
- `exchangeParticles(...)` — particles that crossed patch boundaries
- `sumDensities(...)` — current densities summed across patch boundaries

Non-blocking + waitall pattern. `MPI_THREAD_FUNNELED` — only master OMP thread calls MPI. **Don't add MPI inside `#pragma omp parallel`** without explicit `single`/`master` guards.

## OpenMP layer

Outer `#pragma omp parallel` opens once per simulation. `#pragma omp for` inside `VectorPatch::*` methods reuses it. MPI calls in `#pragma omp single` regions.

**Don't `omp_set_num_threads()` inside physics modules** — conflicts with outer region.

Patch state is mostly local. A few shared structures (patch-to-rank mapping, cumulative diagnostic buffers) need `#pragma omp critical` or `single`. When adding shared state, check for races.

## GPU layer

Smilei v5.1: OpenMP target (`#pragma omp target`) or OpenACC (`#pragma acc kernels`), chosen at compile time.

One MPI rank per GPU. Particle data and fields pinned on GPU. Patch concept maps cleanly — each GPU kernel processes one patch or a block.

`_GPU_KERNEL_` macros (`src/Tools/SmileiPP.h`) bracket GPU-offloaded code:

```cpp
void PusherBoris::operator()(Particles& particles, ...) {
    #ifdef _GPU
    #pragma omp target teams distribute parallel for ...
    #else
    #pragma omp simd
    #endif
    for (int ip = ipart_min; ip < ipart_max; ip++) { /* Boris push */ }
}
```

**Data movement**: particles/fields stay on GPU between kernels. Diagnostic dumps need host copy — frequent diagnostics (`every=1`) destroy GPU performance.

**Feature coverage** (v5.1): Cartesian fully supported; AM cylindrical partial; standard pushers yes; collisions partial (some CPU-only); `from_rate` ionization can't run on GPU (Python callable); envelope partial.

When adding a physics module, CPU first. GPU later with TODO acceptable.

## The main loop

`src/Smilei/Smilei.cpp::main()`. Simplified:

```cpp
int main(int argc, char* argv[]) {
    SmileiMPI smpi(&argc, &argv);
    Params params(&smpi, argv[1]);
    VectorPatch vecPatches(params, &smpi);
    Checkpoints checkpoints(params, &smpi);
    if (params.restart) checkpoints.restartAll(vecPatches, ...);

    #pragma omp parallel
    for (int itime = params.itime_start; itime < params.n_time; ++itime) {
        vecPatches.dynamics(params, &smpi, itime, ...);     // push + deposit
        vecPatches.sumDensities(params, &smpi, ...);
        vecPatches.solveMaxwell(params, &smpi, ...);
        vecPatches.applyBoundaryConditions(...);
        vecPatches.runAllDiags(params, &smpi, itime, ...);
        if (loadBalance.should(itime)) vecPatches.loadBalance(...);
        if (checkpoints.dump(itime)) checkpoints.dumpAll(...);
    }
    return 0;
}
```

## Build system

GNU Make. `Makefile` in repo root.

```bash
make                       # default build
make -j 16                 # parallel
make config=debug          # symbols + assertions
make config=debug_openmp   # + OMPT support
make machine=intel-skylake # apply machine file
make happi                 # Python module
```

Machine files: `scripts/compile_tools/machine/{intel-skylake, amd-rome, ibm-power9, nvidia-a100, ...}`. Each sets CXX, CXXFLAGS (with -march/SIMD), LDFLAGS, GPU flags.

Adding a machine: copy an existing file, edit, `make machine=mycluster`. Right machine file = 2–4× faster than generic.

## Compile-time options

`src/Tools/SmileiPP.h`:
- `_GPU_KERNEL_`, `_GPU`
- `SMILEI_USE_NUMPY`
- `SMILEI_ACCELERATOR_GPU_OMP`, `SMILEI_ACCELERATOR_GPU_OACC`

Active set in `defines.h` (build-generated). Use these macros to gate paths — don't use `__CUDA_ARCH__` or vendor-specific macros.

## Debugging

```bash
make config=debug                        # -O0 -g -DSMILEI_DEBUG
mpirun -np 1 ./smilei namelist.py        # single-rank first
./smilei --verbose namelist.py           # LOG() trace
```

In code:
```cpp
SMILEI_ASSERT(condition, "message");     // debug-only
ERROR("...");                             // hard error, MPI_Abort
WARNING("...");                           // soft warning
```

Sanitizers work on CPU builds (`-fsanitize=address/undefined/thread`). TSan + OpenMP has known false positives — use Intel Inspector for serious thread debugging.

GPU debugging: NVIDIA `compute-sanitizer`, AMD `rocgdb`.

## Performance

`src/Tools/Timers.h` — per-region timing. Active by default; results in `performance.txt`.

```cpp
Timer my_timer("my_physics_module");
my_timer.restart();
// ... expensive code ...
my_timer.update();
```

Profiling: NVIDIA Nsight Systems, Intel VTune, HPCToolkit. For load imbalance, enable `DiagPerformances`.

## Hard truths

- Changes to `Patch` or `VectorPatch` are high-stakes. Run the full validation suite.
- Hilbert ordering changes with patch count. Diagnostics referring to "patch 0" not portable across configurations.
- MPI/OpenMP boundaries are explicit. Don't add MPI inside `omp parallel` regions without `single`/`master`.
- Always rebuild with right machine file when benchmarking.
