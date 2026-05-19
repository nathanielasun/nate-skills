# Smilei architecture (mini tier)

Patches, MPI/OpenMP/GPU layers, build system. Compressed.

## Patch — the parallelism unit

A patch holds: contiguous block of cells (interior), halo of ghost cells, all particles in interior, all field arrays. Patches are NOT per-rank subdomains — many patches per rank, assignment is dynamic via load balancing.

`Main.number_of_patches` must be power-of-2 per direction. Total ≥ 4 × MPI ranks. Typical: 16–64 cells per direction.

## `Patch` and `VectorPatch`

`src/Patch/`:
- **`Patch`** (`Patch.h`): per-patch state. `ElectroMagn* EMfields`, `std::vector<Species*> vecSpecies`, particle BCs, local diagnostics, per-patch RNG.
- **`VectorPatch`** (`VectorPatch.h`): owns patches of the current MPI rank. Main-loop methods (`dynamics`, `solveMaxwell`, `runAllDiags`) live here, dispatching OMP-parallel over patches.

```cpp
void VectorPatch::dynamics(...) {
    #pragma omp for schedule(runtime)
    for (size_t ipatch = 0; ipatch < patches_.size(); ++ipatch) {
        patches_[ipatch]->dynamics(...);
    }
    sumDensitiesAndCurrents(...);              // inter-patch communication
}
```

## Hilbert curve + load balancing

Patches Hilbert-ordered (`SmileiMPI::computeHilbertCurve()`), not row-major. Keeps spatially nearby patches close linearly — preserves locality during rank assignment.

`LoadBalancing` redistributes periodically: each patch reports load, Hilbert curve traversed to balance, migrating patches move with particles/fields/diag state via MPI point-to-point. Per-patch migration (no per-cell).

Tune `LoadBalancing(every=N)` so migration cost amortizes.

## MPI layer

`src/SmileiMPI/SmileiMPI.h`. One object per rank. Key collectives:
- `exchangeFields(...)` — ghost cells across patches
- `exchangeParticles(...)` — particles crossing patch boundaries
- `sumDensities(...)` — current densities summed across patch boundaries

Non-blocking + waitall pattern. `MPI_THREAD_FUNNELED` — only master OMP thread calls MPI. **Don't add MPI inside `#pragma omp parallel`** without `single`/`master` guards.

## OpenMP layer

Outer `#pragma omp parallel` opens once per simulation. `#pragma omp for` inside `VectorPatch::*` methods reuses it. MPI calls in `#pragma omp single` regions.

Don't `omp_set_num_threads()` inside physics modules.

Patch state mostly local. Shared structures (rank mapping, cumulative diag buffers) need `critical` or `single`. Check races when adding shared state.

## GPU layer

v5.1: OpenMP target or OpenACC (compile-time). One MPI rank per GPU. Particles/fields stay on GPU between kernels. `_GPU_KERNEL_` macros bracket GPU code:

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

**Frequent diagnostics destroy GPU performance** (host copy overhead).

Coverage (v5.1): Cartesian full; AM cylindrical partial; standard pushers yes; collisions partial; `from_rate` ionization CPU-only (Python on GPU = no); envelope partial.

CPU first when adding modules. GPU later with TODO acceptable.

## Build

```bash
make                       # default
make -j 16                 # parallel
make config=debug          # symbols + assertions
make machine=intel-skylake # machine file
make happi                 # Python module
```

Machine files in `scripts/compile_tools/machine/`. Set CXX, CXXFLAGS (SIMD/-march), LDFLAGS, GPU flags. Right machine file: 2–4× faster than generic.

## Compile-time macros (`src/Tools/SmileiPP.h`)

- `_GPU_KERNEL_`, `_GPU`
- `SMILEI_USE_NUMPY`
- `SMILEI_ACCELERATOR_GPU_OMP` / `SMILEI_ACCELERATOR_GPU_OACC`

Active set in `defines.h` (build-generated). Use these to gate paths — don't use `__CUDA_ARCH__` or vendor-specific.

## Debugging

```bash
make config=debug                    # -O0 -g, assertions
mpirun -np 1 ./smilei namelist.py    # single rank first
./smilei --verbose namelist.py       # LOG() trace
```

In code:
```cpp
SMILEI_ASSERT(condition, "message"); // debug-only
ERROR("...");                         // hard error, MPI_Abort
WARNING("...");                       // soft
```

Sanitizers (AddressSan, UBSan, ThreadSan) work on CPU builds. TSan + OpenMP has false positives — use Intel Inspector for serious thread debugging.

GPU: NVIDIA `compute-sanitizer`, AMD `rocgdb`.

## Performance instrumentation

`src/Tools/Timers.h` — per-region timing. Active by default; results in `performance.txt`.

```cpp
Timer my_timer("my_module");
my_timer.restart();
// expensive code
my_timer.update();
```

Profiling: NVIDIA Nsight Systems, Intel VTune, HPCToolkit. For load imbalance, enable `DiagPerformances` and inspect per-patch breakdown.

## Hard truths

- Changes to `Patch`/`VectorPatch` are high-stakes — full validation suite.
- Hilbert ordering changes with patch count → diagnostics referencing "patch 0" not portable.
- MPI/OpenMP boundaries are explicit and intentional.
- Always rebuild with right machine file when benchmarking.
