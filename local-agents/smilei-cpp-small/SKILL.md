---
name: smilei-cpp-small
description: Use when working in the Smilei C++ source code — modifying internals, adding physics modules (collision types, ionization models, solvers, pushers), debugging the BLAST-Smilei codebase, navigating its architecture (patches, Hilbert ordering, MPI/OpenMP/GPU layers), writing tests, or preparing PRs. Trigger on "Smilei C++", "BLAST-Smilei", `Patch`, `VectorPatch`, `ParticleContainer`, `ElectroMagn`, `Pusher`, `Interpolator`, `Projector`, `BinaryProcesses`, source paths under `src/`, or any non-trivial code modification. Optimized for low-temperature-plasma / excimer work but covers the full C++ codebase. For Python namelist work, use smilei-namelist-small instead.
---

# Smilei C++ Codebase Skill (small tier)

This skill covers work *inside* the Smilei C++ source — modifying the codebase, adding capability, debugging, contributing back. For namelist authoring see `smilei-namelist-small`.

Smilei is ~200k LOC of C++17. Unlike WarpX (built on AMReX), Smilei's data structures, parallelism, and I/O are all in-tree. Tested against v5.1.

## Hard rules for source modification

1. **Read `architecture.md` before touching the source.** Smilei's patch/Hilbert structure is the single most important invariant. Misunderstand it and changes deadlock or corrupt particle data.

2. **Match existing C++17 style.** No std::format, concepts, or ranges. Use `#pragma omp` for parallelism, not C++ atomics or `std::execution`. Read 3 existing files before writing.

3. **Particle data is structure-of-arrays.** `Particles` (`src/Particles/Particles.h`) holds `position[d]`, `momentum[d]`, `weight`, `charge`, `chi` as separate `std::vector<double>`. Never assume AoS.

4. **Patches are the parallelism unit.** MPI by patch (Hilbert-ordered); OpenMP over patches within a rank; cells are interior to patches.

5. **Ghost cells follow solver order.** Yee: 2 per face. Lehe: 3. Cowan: 3. Bouchard: 4–5. Solver code must respect the count.

6. **GPU paths via OpenMP target or OpenACC (compile-time).** `_GPU_KERNEL_` macros bracket GPU code. CPU and GPU paths exist in parallel — modifying one requires keeping the other in sync.

7. **HDF5 schema is versioned.** Adding/modifying diagnostics requires schema bump and happi update.

8. **Random streams are per-patch.** Use `patch->rand`, not `std::rand` or `<random>`.

9. **`reference_angular_frequency_SI` gates dimensional modules.** Collisions, ionization, radiation, MBW query `params.reference_angular_frequency_SI`. New SI-needing modules follow the same pattern.

10. **Tests are mandatory for new physics.** `validation/` has reference outputs. New module needs (a) a test namelist and (b) a reference output, within 1–2% of analytical solution where available.

## Voice for C++ work

- Prefer pointing to existing patterns: "Look at `src/Collisions/CollisionalIonization.cpp` — cleanest example of a per-species physics module hooked into the per-patch loop."
- State CPU/GPU side explicitly.
- Cite line ranges or function names, not abstractions.
- Don't invent class names. If unsure, confirm against current source.

## Routing

| Topic | Reference |
|---|---|
| Architecture, build system, MPI/OpenMP/GPU layers, patches, Hilbert curve, load balancing | `references/architecture.md` |
| `Particles` SoA, `Species`, pushers, particle BCs, injection, current deposition, Maxwell solvers, `ElectroMagn`, PML | `references/core-data-structures.md` |
| Pérez 2012 collision implementation, collisional/field ionization, ADK/PPT/BSI/envelope-averaged paths, adding a new model | `references/physics-modules.md` |
| Diagnostic implementation, HDF5 schema, adding new diagnostics, happi integration | `references/diagnostics-source.md` |

## Directory layout

```
src/
├── Smilei/Smilei.cpp              main entry, PIC loop
├── Patch/                          Patch, VectorPatch
├── SmileiMPI/                      MPI orchestration, Hilbert curve, LoadBalancing
├── Params/                         namelist parsing
├── Particles/                      Particles class (SoA), Species
├── ElectroMagn/                    EMFields, Maxwell solvers, geometry variants
├── ElectroMagnBC/                  EM boundary conditions (incl. PML)
├── Field/                          Field, cField (complex for AM/envelope)
├── Pusher/                         PusherBoris, Vay, HigueraCary, Photon, BorisNR
├── Interpolator/                   field → particle (Esirkepov)
├── Projector/                      particle → current grid (Esirkepov)
├── Collisions/                     Pérez 2012, CollisionalIonization, NuclearReaction
├── Ionization/                     field ionization (Tunnel, PPT, BSI, ...)
├── Laser/                          laser injection
├── Diagnostic/                     each Diag type as a class
├── Random/                         per-patch RNG
└── Tools/                          SmileiPP.h, Timers, HDF5 wrappers
happi/                              Python post-processing
benchmarks/                         validation reference simulations
validation/                         test suite
scripts/compile_tools/              machine files
```

## Main loop

`src/Smilei/Smilei.cpp` orchestrates:

```cpp
#pragma omp parallel
for (timestep = 0; timestep < n_timesteps; ++timestep) {
    vecPatches.dynamics(...);          // push + current deposition (per-patch, OMP)
    vecPatches.sumDensities(...);       // inter-patch current sums (MPI)
    vecPatches.solveMaxwell(...);       // Maxwell update (per-patch, OMP)
    vecPatches.applyBoundaryConditions(...);
    vecPatches.runAllDiags(...);
    if (loadBalance.should()) vecPatches.loadBalance(...);
    if (checkpoints.dump()) checkpoints.dumpAll(...);
}
```

OpenMP `parallel` opens once per simulation; per-method `#pragma omp for` reuses it. MPI calls happen in `#pragma omp single` regions outside the parallel-for loops.

## Workflow patterns

**Adding a new physics module:**
1. Read `architecture.md`.
2. Find the closest existing analog (e.g., `Collisions/CollisionalIonization.cpp` for a new collision type, `Pusher/PusherVay.cpp` for a new pusher).
3. Identify the namelist hook — find the string parse in `src/Params/PicParams.cpp` and the factory dispatch.
4. Implement with CPU and GPU paths (or TODO).
5. Add a `validation/tst_*.py` test.
6. Run full validation suite before submitting.

**Debugging:**
1. `make config=debug` for symbols + assertions.
2. `mpirun -np 1` first; many MPI hangs mask on single-rank.
3. For races, `make config=debug_openmp`.
4. GPU: `compute-sanitizer` (NVIDIA), `rocgdb` (AMD).

**PR:**
1. Branch from `develop`, not `master`.
2. Pass local validation suite.
3. Update `CHANGELOG.md`.
4. Match clang-format config.
5. For physics, include benchmark vs analytical or published result.

## Build

```bash
make                       # default
make -j 16                 # parallel
make config=debug          # symbols + assertions
make machine=intel-skylake # use machine file (in scripts/compile_tools/machine/)
make happi                 # build Python happi module
```

Machine files set CXX, CXXFLAGS (SIMD/-march), LDFLAGS, GPU flags. Adapting for your cluster: copy an existing file, edit, `make machine=mycluster`. Performance with the right machine file is 2–4× the generic build.

## Compile-time options

In `src/Tools/SmileiPP.h`:
- `_GPU_KERNEL_` — bracket GPU code
- `_GPU` — compile-time GPU enabled
- `SMILEI_USE_NUMPY` — NumPy in callables
- `SMILEI_ACCELERATOR_GPU_OMP` / `SMILEI_ACCELERATOR_GPU_OACC` — GPU backend

Set by Makefile based on machine file and `config=` flags. Check `defines.h` (build-generated) for the active set.

## Hard truths

- **Macros and `#ifdef`s shape the build.** Many code paths are compile-time-selected.
- **Headers do work.** Many definitions are inlined in `.h` for performance.
- **Geometry is explicit, not templated.** 1D/2D/3D/AM have separate `.cpp` files for many classes. Modifying one class usually means modifying all variants.
- **OpenMP and MPI are interleaved deliberately.** `#pragma omp single` guards MPI inside the parallel region.

## Closing notes

- For unfamiliar code, `grep -rn "Pattern" src/` before guessing.
- Build with `make -j 8 config=debug` during development.
- The validation suite is the ground truth.
- Maintainers (Grech, Pérez, Beck, Vinci) are responsive on GitHub Discussions.
