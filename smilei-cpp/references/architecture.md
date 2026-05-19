# Smilei architecture

Patches, Hilbert ordering, MPI/OpenMP/GPU layers, the main loop in detail, build system, machine files, debugging machinery. The foundational reference for any source modification.

## Table of contents

1. The patch — Smilei's parallelism unit
2. `Patch` and `VectorPatch`
3. Hilbert space-filling curve and load balancing
4. MPI layer — `SmileiMPI`
5. OpenMP layer — patch-level parallelism
6. GPU layer — offload model
7. The main loop in detail
8. Build system and machine files
9. Compile-time options and macros
10. Directory layout (deeper than SKILL.md)
11. Debugging machinery
12. Performance instrumentation

---

## 1. The patch — Smilei's parallelism unit

A **patch** is the unit of parallel decomposition in Smilei. The simulation domain is divided into patches; each patch holds:

- A contiguous block of grid cells (the **interior**)
- A halo of **ghost cells** for stencil operations
- All particles whose positions are inside the interior
- All field arrays for that subdomain

**Patches are not subdomains of MPI ranks.** Many patches live within an MPI rank; the rank's set of patches is determined dynamically by load balancing. Each rank's patch list is the unit of OpenMP parallelism — `#pragma omp parallel for` iterates over the patches owned by the rank.

The patch concept is similar to but lighter than AMReX's "box": patches are not refined hierarchically (Smilei is a uniform-grid PIC code, no AMR). Patch sizes are uniform — set by `Main.number_of_patches` and `Main.grid_length` in the namelist.

Typical patch size: 16–64 cells per direction. Total patches: a power of 2 per direction, at least 4× the MPI rank count.

## 2. `Patch` and `VectorPatch`

Defined in `src/Patch/`:

- **`Patch`** (`src/Patch/Patch.h`) — owns the state for one patch: `ElectroMagn* EMfields`, `std::vector<Species*> vecSpecies`, particle BCs, diagnostics local state, the per-patch RNG. Pure data + per-patch methods.

- **`VectorPatch`** (`src/Patch/VectorPatch.h`) — owns the patches of the current MPI rank (`std::vector<Patch*>`) and orchestrates per-rank operations. The main-loop methods (`dynamics`, `solveMaxwell`, `runAllDiags`) live here.

The pattern: `VectorPatch` is the collective interface for one MPI rank's patches; `Patch` is the per-patch state. Most modifications to physics modules touch both layers — `VectorPatch::dynamics()` dispatches OpenMP-parallel over `patch->dynamics(...)`, which then calls into pushers, projectors, etc.

```cpp
// VectorPatch.cpp (sketch)
void VectorPatch::dynamics(Params& params, SmileiMPI* smpi, ...) {
    #pragma omp for schedule(runtime)
    for (size_t ipatch = 0; ipatch < patches_.size(); ipatch++) {
        patches_[ipatch]->dynamics(params, smpi, ...);
    }
    // Then inter-patch communication for currents
    sumDensitiesAndCurrents(...);
}
```

The `schedule(runtime)` lets users tune via `OMP_SCHEDULE=dynamic`/`guided` — important for load-imbalanced runs.

## 3. Hilbert space-filling curve and load balancing

Patches are ordered along a Hilbert curve, not a row-major raster. The Hilbert ordering keeps spatially nearby patches close in the linear ordering, which preserves locality when patches are assigned to MPI ranks.

The ordering lives in `src/SmileiMPI/SmileiMPI.cpp::computeHilbertCurve()`. Given a 2D or 3D patch grid, the Hilbert index is computed bit-interleaving the patch coordinates.

**Load balancing:**

`LoadBalancing` namelist block triggers periodic redistribution. Each redistribution:

1. Each patch reports its current "load" (estimated compute cost = `cell_load × n_cells + frequency_load × sum_of_particle_counts`).
2. The Hilbert curve is traversed to assign patches to ranks such that per-rank loads are balanced.
3. Patches that change ranks are migrated (particle data + field data move with the patch).

Implementation lives in `src/SmileiMPI/SmileiMPI.cpp::balanceLoad()`. Migration uses MPI point-to-point send/recv. Migration cost grows with particle count per patch — for dense plasmas with high PPC, load balancing can take seconds per call. Tune `LoadBalancing(every=N)` so migration cost amortizes over `N` steps.

**Important invariant**: patches are migrated as a unit. A patch's particles, fields, and diagnostic state all move together. There is no per-cell migration.

## 4. MPI layer — `SmileiMPI`

`src/SmileiMPI/SmileiMPI.h` defines `SmileiMPI`, the single object per rank that owns:

- The `MPI_Comm` (typically `MPI_COMM_WORLD`)
- The rank's index, total rank count
- The Hilbert curve table
- Per-rank send/recv buffers for inter-patch communication

**Key collective operations:**

- `exchangeFields(...)` — ghost-cell exchange for `ElectroMagn` fields between neighboring patches (across rank boundaries).
- `exchangeParticles(...)` — particles that crossed patch boundaries are sent to the new owning patch.
- `sumDensities(...)` — current densities `Jx/Jy/Jz/Rho` are summed across patch boundaries (current deposition is per-patch; sums fix the boundaries).

Communications are typically **non-blocking + waitall** for overlap with compute. Look in `src/Patch/VectorPatchAM.cpp` (or geometry-specific variants) for the canonical pattern.

**MPI thread safety**: Smilei uses `MPI_THREAD_FUNNELED` — only the master OpenMP thread calls MPI. Don't add MPI calls inside `#pragma omp parallel` regions without explicit master/single guards.

## 5. OpenMP layer — patch-level parallelism

Each MPI rank parallelizes over patches via OpenMP. The pattern is `#pragma omp for` inside a `#pragma omp parallel` region that wraps the entire timestep.

```cpp
// Smilei.cpp (sketch)
#pragma omp parallel
{
    for (timestep = 0; ...) {
        vecPatches.dynamics(...);       // omp for over patches inside
        #pragma omp single                  // MPI calls only from master
        {
            smpi->exchangeParticles(...);
        }
        vecPatches.solveMaxwell(...);   // omp for over patches inside
        // ...
    }
}
```

The outer `#pragma omp parallel` opens once per simulation, not once per timestep — this avoids fork/join overhead. The dispatch inside each `VectorPatch::*` method uses `#pragma omp for`, sharing the existing parallel region.

**Critical invariant**: shared state across patches must be guarded. Most patch state is local, but a few global structures (the patch-to-rank mapping, the cumulative diagnostic buffers) require explicit `#pragma omp critical` or `#pragma omp single` regions. When adding new shared state, double-check race conditions.

**Thread count**: set via `OMP_NUM_THREADS`. Smilei queries this at startup and uses it for the parallel region. Inside physics modules, **do not call `omp_set_num_threads()`** — it conflicts with the outer parallel region.

## 6. GPU layer — offload model

Smilei v5.1 supports GPU offload via OpenMP target (`#pragma omp target`) or OpenACC (`#pragma acc kernels`), chosen at compile time via the machine file.

**Offload model**: each MPI rank owns one GPU (the typical configuration). Particle data and fields are mirrored on the GPU; the patch-level loop dispatches kernels rather than CPU code. The patch concept maps cleanly: each GPU kernel processes one patch (or a block of patches) at a time.

**Code paths**: in `.cpp` files, GPU-offloaded code is bracketed by `_GPU_KERNEL_` macros (defined in `src/Tools/SmileiPP.h`). For example, in `src/Pusher/PusherBoris.cpp`:

```cpp
void PusherBoris::operator()(Particles& particles, ...) {
    #ifdef _GPU
    #pragma omp target teams distribute parallel for ...
    #else
    #pragma omp simd
    #endif
    for (int ip = ipart_min; ip < ipart_max; ip++) {
        // Boris push for particle ip
    }
}
```

**Data movement**: particle and field data are pinned on the GPU between kernel invocations. Diagnostic dumps require a host-side copy — these are issued at diagnostic cadence. Frequent diagnostics (`every=1`) destroy GPU performance because of constant data movement.

**GPU debugging**: NVIDIA `compute-sanitizer`, AMD `rocgdb`. Standard tools work but are slow.

**Feature coverage** (as of v5.1):
- Cartesian (1D/2D/3D): supported
- AM cylindrical: partial (check release notes)
- Standard pushers: yes
- Collisions: partial — some intra-species paths are CPU-only
- Field ionization: partial — `tunnel` works, `from_rate` callable doesn't (Python callable can't run on GPU)
- Envelope: partial

When adding a new physics module, write a CPU version first. Add GPU later when the CPU path is validated; document that GPU support is pending with a TODO.

## 7. The main loop in detail

`src/Smilei/Smilei.cpp::main()` is the entry point. The simplified loop:

```cpp
int main(int argc, char* argv[]) {
    // Initialization
    SmileiMPI smpi(&argc, &argv);
    Params params(&smpi, argv[1]);                 // parses namelist
    VectorPatch vecPatches(params, &smpi);

    Checkpoints checkpoints(params, &smpi);
    if (params.restart) checkpoints.restartAll(vecPatches, ...);

    // Main loop
    #pragma omp parallel
    for (int itime = params.itime_start; itime < params.n_time; ++itime) {
        // 1. Move particles, deposit currents
        vecPatches.dynamics(params, &smpi, itime, ...);

        // 2. Communicate currents across patch boundaries
        vecPatches.sumDensities(params, &smpi, ...);

        // 3. Solve Maxwell
        vecPatches.solveMaxwell(params, &smpi, ...);

        // 4. Particle and field boundary conditions
        vecPatches.applyBoundaryConditions(...);

        // 5. Diagnostics
        vecPatches.runAllDiags(params, &smpi, itime, ...);

        // 6. Periodic load balancing
        if (loadBalance.should(itime)) {
            vecPatches.loadBalance(params, &smpi, ...);
        }

        // 7. Checkpoint
        if (checkpoints.dump(itime)) {
            checkpoints.dumpAll(vecPatches, ...);
        }
    }

    return 0;
}
```

Each step is documented in subsequent references:
- `dynamics` → `references/particles-and-species.md` (push) + `solvers-and-fields.md` (deposit)
- `solveMaxwell` → `solvers-and-fields.md`
- `runAllDiags` → `diagnostics-source.md`
- collision step (within `dynamics`) → `collisions-source.md`
- ionization step (within `dynamics`) → `ionization-source.md`

## 8. Build system and machine files

Smilei uses **GNU Make**, not CMake. The `Makefile` is in the repo root.

**Common targets:**

```bash
make                       # default build
make -j 16                 # parallel
make config=debug          # debug symbols + assertions
make config=debug_openmp   # + OMPT support
make machine=intel-skylake # use a machine file
make happi                 # build the Python happi module
make doc                   # build Sphinx docs
make clean                 # remove object files
```

**Machine files** in `scripts/compile_tools/machine/` set compiler flags per architecture. Examples shipped: `intel-skylake`, `amd-rome`, `ibm-power9`, `nvidia-a100`. Each file sets:

- `CXX` (compiler)
- `CXXFLAGS` (optimization, SIMD flags)
- `LDFLAGS`
- GPU-specific flags if applicable

**Adding a new machine**:

```bash
cp scripts/compile_tools/machine/intel-skylake scripts/compile_tools/machine/my-cluster
# Edit my-cluster: set CXX, SIMD flags (-march=native or specific), MPI wrapper
make machine=my-cluster
```

Performance with the right machine file can be 2–4× faster than the generic build.

## 9. Compile-time options and macros

In `src/Tools/SmileiPP.h`:

- `_GPU_KERNEL_` — bracket GPU-offload code
- `_GPU` — set at compile time when GPU support is on
- `SMILEI_USE_NUMPY` — build with NumPy support in namelist callables
- `SMILEI_AMR` — adaptive mesh refinement (experimental)
- `SMILEI_ACCELERATOR_GPU_OMP` — GPU via OpenMP target
- `SMILEI_ACCELERATOR_GPU_OACC` — GPU via OpenACC

These are set by the Makefile based on machine file and `config=` arguments. Check `defines.h` (generated at build time) for the active set.

When writing portable code, use these macros to gate GPU/CPU paths. Don't use `__CUDA_ARCH__` or vendor-specific macros directly — they break across GPU vendors.

## 10. Directory layout (deeper)

Beyond the top-level layout in SKILL.md:

```
src/
├── Smilei/Smilei.cpp                       # main entry point
├── Patch/
│   ├── Patch.{h,cpp}                       # per-patch state
│   ├── VectorPatch.{h,cpp}                 # collective operations
│   └── VectorPatchAM.cpp                   # AM cylindrical specializations
├── SmileiMPI/
│   ├── SmileiMPI.{h,cpp}                   # MPI orchestration, Hilbert curve
│   └── LoadBalancing.{h,cpp}               # load-balancing algorithm
├── Params/
│   ├── Params.{h,cpp}                      # parsed namelist parameters
│   └── PicParams.{h,cpp}                   # core PIC parameters
├── Particles/
│   ├── Particles.{h,cpp}                   # SoA particle storage
│   ├── ParticleData.h                      # struct-of-vectors layout
│   └── Species.{h,cpp}                     # per-species container
├── ElectroMagn/
│   ├── ElectroMagn.{h,cpp}                 # base class
│   ├── ElectroMagn{1,2,3}D.{h,cpp}         # geometry specializations
│   ├── ElectroMagnAM.{h,cpp}               # cylindrical
│   └── MaxwellSolver/                      # Yee, Lehe, Cowan, Bouchard, ...
├── Field/
│   ├── Field.{h,cpp}                       # base field class
│   └── cField.{h,cpp}                      # complex field (for AM, envelope)
├── Pusher/
│   ├── Pusher.h                            # base class
│   ├── PusherBoris.{h,cpp}
│   ├── PusherVay.{h,cpp}
│   ├── PusherHigueraCary.{h,cpp}
│   ├── PusherBorisNR.{h,cpp}               # non-relativistic
│   └── PusherPhoton.{h,cpp}                # photons (mass=0)
├── Interpolator/                            # E and B → particle position
├── Projector/                               # particle → grid current deposition
├── Collisions/
│   ├── Collisions.{h,cpp}                  # binary collision driver
│   ├── CollisionalIonization.{h,cpp}       # impact ionization
│   ├── CollisionalNuclearReaction.{h,cpp}  # fusion reactions
│   └── BinaryProcesses.{h,cpp}             # wraps Coll + CollIon + NucReac
├── Ionization/
│   ├── Ionization.{h,cpp}                  # base
│   ├── IonizationTunnel.{h,cpp}            # ADK
│   ├── IonizationTunnelFullPPT.{h,cpp}     # PPT
│   ├── IonizationTunnelBSI.{h,cpp}         # ADK + BSI
│   ├── IonizationTunnelEnvelopeAveraged.{h,cpp}
│   └── IonizationFromRate.{h,cpp}          # user-defined
├── Diagnostic/                              # each diag type as a class
├── Laser/                                   # laser injection
├── ParticleBC/                              # particle boundary conditions
├── ElectroMagnBC/                           # field boundary conditions (incl. PML)
└── Tools/                                   # SmileiPP.h, Timers, HDF5 wrappers, etc.

happi/                                       # Python post-processing (built separately)
benchmarks/                                  # validation reference simulations
validation/                                  # test suite
scripts/compile_tools/                       # machine files
doc/                                         # Sphinx documentation source
```

## 11. Debugging machinery

**Debug build:**

```bash
make config=debug
```

Adds `-O0 -g -DSMILEI_DEBUG`. Many `SMILEI_ASSERT(...)` calls fire in debug mode.

**Asserting on a condition during development:**

```cpp
#include "Tools/Asserts.h"

SMILEI_ASSERT(particle_count_valid, "Particle count must be non-negative");
```

These are compiled out in release builds.

**ERROR macros:**

```cpp
ERROR("Detailed message including " << variable);    // hard error, MPI_Abort
WARNING("Detailed message");                          // soft warning
```

Errors trigger `MPI_Abort` — propagating cleanly across ranks. Use them for unrecoverable conditions.

**Verbose mode**:

```bash
./smilei --verbose namelist.py
```

Triggers `LOG()` calls — useful for tracing the main loop and initialization order.

**Sanitizer builds**:

```bash
make config=debug CXXFLAGS+="-fsanitize=address"
```

UBSan, ASan, TSan all work on CPU builds. Be aware: TSan and OpenMP have known false positives; use Intel Inspector for serious thread debugging.

## 12. Performance instrumentation

`src/Tools/Timers.h` provides per-region timing. Active by default; results dumped to `performance.txt` at simulation end.

Adding a custom timer:

```cpp
#include "Tools/Timers.h"

Timer my_timer("my_physics_module");
my_timer.restart();
// ... expensive code ...
my_timer.update();
```

The timer accumulates across the simulation and reports total time + count.

For profiler integration:
- **NVIDIA Nsight Systems**: `nsys profile ./smilei namelist.py`
- **Intel VTune**: `vtune -collect hotspots ./smilei namelist.py`
- **HPCToolkit**: works on most platforms

For per-rank load imbalance, enable `DiagPerformances` in the namelist (see `diagnostics-source.md` §10).

## Closing notes for architecture work

- Changes to `Patch` or `VectorPatch` are high-stakes — exercise extreme caution and run the full validation suite.
- The Hilbert ordering is stable across runs but **changes with patch count**. Diagnostics that look at "patch 0" are not portable across configurations.
- The MPI/OpenMP boundaries are explicit and intentional. Adding MPI calls inside parallel regions, or OpenMP parallelism outside the existing `#pragma omp parallel`, will break in subtle ways.
- Always rebuild with the right machine file when benchmarking. A `make -j 16` without `machine=...` can be 3× slower than the same code with proper flags.
