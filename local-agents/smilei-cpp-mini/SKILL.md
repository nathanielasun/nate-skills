---
name: smilei-cpp-mini
description: Use for working in the Smilei C++ source code — modifying internals, adding physics modules, debugging the BLAST-Smilei codebase, navigating its architecture (patches, MPI/OpenMP/GPU), writing tests, preparing PRs. Trigger on "Smilei C++", "BLAST-Smilei", `Patch`, `VectorPatch`, `Particles`, `Species`, `ElectroMagn`, `Pusher`, `Interpolator`, `Projector`, `BinaryProcesses`, source paths under `src/`, or any non-trivial source modification. Optimized for collision/ionization paths (excimer/LTP-relevant) but covers full C++ codebase. For Python namelist work, use smilei-namelist-mini instead.
---

# Smilei C++ Codebase Skill (mini tier)

Work *inside* the Smilei C++ source — modifying, adding capability, debugging. For namelist authoring see `smilei-namelist-mini`.

Smilei is ~200k LOC C++17. Unlike WarpX (built on AMReX), Smilei's data structures, parallelism, and I/O are all in-tree. Tested against v5.1.

## Hard rules — never violate

1. **Read `architecture.md` before touching the source.** Patch/Hilbert structure is the load-bearing invariant.

2. **Match existing C++17 style.** No std::format, concepts, ranges. `#pragma omp` for parallelism, not `std::execution`. Read 3 existing files before writing.

3. **Particle data is SoA.** `Particles` holds `position[d]`, `momentum[d]`, `weight`, `charge`, `chi` as separate `std::vector<double>`. Never AoS.

4. **Patches are the parallelism unit.** MPI by patch (Hilbert-ordered); OpenMP over patches within a rank. Cells are interior; no per-cell MPI.

5. **Ghost cells follow solver order.** Yee: 2. Lehe: 3. Cowan: 3. Bouchard: 4–5.

6. **GPU paths use `_GPU_KERNEL_` macros.** CPU and GPU exist in parallel — keep both in sync.

7. **HDF5 schema versioned.** Change diagnostics → bump version, update happi.

8. **Use `patch->rand`** for randomness, not `std::rand`.

9. **`reference_angular_frequency_SI` gates dimensional code paths** (collisions, ionization, radiation, MBW).

10. **Tests are mandatory for new physics.** `validation/` has reference outputs. New module needs test namelist + reference output, ideally within 1–2% of analytical solution.

## Directory layout

```
src/
├── Smilei/Smilei.cpp        main entry, PIC loop
├── Patch/                    Patch, VectorPatch
├── SmileiMPI/                MPI, Hilbert, LoadBalancing
├── Params/                   namelist parsing
├── Particles/                Particles SoA, Species
├── ElectroMagn/              EMFields, Maxwell solvers (per geometry)
├── ElectroMagnBC/            EM BC (Silver-Müller, PML, periodic, reflective)
├── Field/                    Field, cField (complex for AM/envelope)
├── Pusher/                   Boris, Vay, HigueraCary, Photon, BorisNR
├── Interpolator/             field → particle (Esirkepov)
├── Projector/                particle → current (Esirkepov, charge-conserving)
├── Collisions/               Pérez 2012 + CollisionalIonization + NuclearReaction
├── Ionization/               field ionization (Tunnel/PPT/BSI/envelope/from_rate)
├── Laser/                    laser injection
├── Diagnostic/               each diag type as class
├── Random/                   per-patch RNG
└── Tools/                    SmileiPP.h, Timers, HDF5 wrappers
happi/                        Python post-processing (built via `make happi`)
benchmarks/, validation/      reference simulations, test suite
scripts/compile_tools/        machine files
```

## Routing — read references when needed

| For details on… | Read |
|---|---|
| Architecture, patches, MPI/OMP/GPU layers, Hilbert curve, build system, debugging | `references/architecture.md` |
| `Particles` SoA, `Species`, pushers, projectors, interpolators, Maxwell solvers, fields, PML, collisions, ionization | `references/data-and-physics.md` |
| Diagnostic implementation, HDF5 schema, adding diagnostics, happi | `references/diagnostics.md` |

## Main loop

`src/Smilei/Smilei.cpp::main()`:

```cpp
SmileiMPI smpi(&argc, &argv);
Params params(&smpi, argv[1]);             // parses namelist
VectorPatch vecPatches(params, &smpi);

#pragma omp parallel
for (int itime = ...; itime < n_time; ++itime) {
    vecPatches.dynamics(...);              // push + deposit (per-patch, OMP)
    vecPatches.sumDensities(...);           // inter-patch current sums (MPI)
    vecPatches.solveMaxwell(...);           // Maxwell (per-patch, OMP)
    vecPatches.applyBoundaryConditions(...);
    vecPatches.runAllDiags(...);
    if (loadBalance.should()) vecPatches.loadBalance(...);
    if (checkpoints.dump()) checkpoints.dumpAll(...);
}
```

Outer `omp parallel` opens once per simulation. Per-method `omp for` reuses it. MPI calls in `omp single` regions.

## Per-patch dynamics

`Species::dynamics()`:

```cpp
if (Ionize) Ionize->apply(particles, EMfields);    // field ionization (before pusher)
Interp->fieldsAtParticle(EMfields, particles);     // E, B at particle positions
(*Push)(*particles, ...);                          // momentum + position
partBoundCond.apply(particles, ...);               // particle BCs (before projection)
Proj->depositCurrents(JEM, *particles);            // current deposition
```

Order matters: ionization → interpolation → push → BCs → projection.

## `Particles` SoA structure

```cpp
class Particles {
    std::vector<double> position[3];     // x,y,z arrays
    std::vector<double> momentum[3];     // px,py,pz arrays
    std::vector<double> weight;
    std::vector<short>  charge;          // per-particle (changes with ionization)
    std::vector<double> chi;             // radiation reaction quantum parameter
    std::vector<uint64_t> id;            // tracking

    size_t size() const { return weight.size(); }
    void resize(size_t n);
    void eraseParticle(size_t i);        // O(1) swap-with-last
};
```

All vectors same size. Order NOT stable across calls (collisions/sort/migration reorder).

## Build

```bash
make                       # default
make -j 16                 # parallel
make config=debug          # symbols + assertions
make machine=intel-skylake # machine file (in scripts/compile_tools/machine/)
make happi                 # Python module
```

Right machine file = 2–4× faster than generic. Copy existing file, edit CXX/CXXFLAGS, `make machine=mycluster`.

## Compile-time options (`src/Tools/SmileiPP.h`)

- `_GPU_KERNEL_` — bracket GPU code
- `_GPU` — compile-time GPU enabled
- `SMILEI_USE_NUMPY` — NumPy in callables
- `SMILEI_ACCELERATOR_GPU_OMP` / `SMILEI_ACCELERATOR_GPU_OACC` — GPU backend

Active set in `defines.h` (build-generated).

## Debugging

```bash
make config=debug                    # -O0 -g -DSMILEI_DEBUG
mpirun -np 1 ./smilei namelist.py    # single rank first
./smilei --verbose namelist.py
```

```cpp
SMILEI_ASSERT(condition, "message"); // debug-only
ERROR("...");                         // hard error, MPI_Abort
WARNING("...");                       // soft
```

Sanitizers work on CPU. TSan + OpenMP has false positives — use Intel Inspector. GPU: NVIDIA `compute-sanitizer`, AMD `rocgdb`.

Profiling: NVIDIA Nsight Systems, Intel VTune, HPCToolkit. For load imbalance, enable `DiagPerformances`.

## Workflow patterns

**Adding a new physics module:**
1. Read `architecture.md` for patch loop.
2. Find closest analog (`Collisions/CollisionalIonization.cpp` for collision-type, `Pusher/PusherVay.cpp` for pusher, `Ionization/IonizationTunnel.cpp` for ionization model).
3. Find the namelist string parse in `Params/PicParams.cpp` and the factory dispatch.
4. Implement with CPU + GPU paths (or TODO).
5. Add `validation/tst_*.py` with reference output.
6. Run full validation suite before submitting.

**PR:**
1. Branch from `develop`.
2. Pass validation suite.
3. Update `CHANGELOG.md`.
4. Match clang-format.
5. For physics, include benchmark vs analytical or published result.

## Common pitfalls

- **Order of particles not stable** across calls — don't cache indices.
- **CPU/GPU paths diverge** — when modifying one, update the other.
- **Per-particle vs species charge**: post-ionization use `particles.charge[ip]`, not `Species::charge`.
- **MPI inside `omp parallel`** is forbidden without `single`/`master` guards.
- **Allocations in hot loops**: pre-reserve, don't `new` in `dynamics`.
- **Ghost-cell stale reads**: refresh via `exchangeFields` first.
- **Wrong field location**: `Ex(i,j)` at `(i*dx + dx/2, j*dy)`, not `(i*dx, j*dy)`.
- **Patch smaller than ghost extent** → refuses startup.
- **GPU atomic missing in current deposition** → silent race.
- **`Field*` resize at runtime** is forbidden.
- **Missing `reference_angular_frequency_SI` with collisions/ionization**: silent NaN/order-of-magnitude error.
- **Test particles (`is_test`)** skip current deposition AND collisions — new modules must respect.
- **Schema drift**: change HDF5 layout on C++ side, forget happi → user sees data, can't read.

## Voice for C++ work

- Point to existing patterns. "Look at `src/Collisions/CollisionalIonization.cpp`, cleanest per-species hook into the patch loop."
- State CPU/GPU side explicitly.
- Cite line ranges or function names, not abstractions.
- Don't invent class names. If unsure, confirm against current source — API drifts across releases.

## Hard truths

- **Macros and `#ifdef`s** shape the build heavily.
- **Headers do work** — many definitions inlined in `.h` for performance.
- **Geometry is explicit**, not templated. 1D/2D/3D/AM have separate `.cpp` files for many classes.
- **OpenMP/MPI boundaries are explicit and intentional.**

## Closing

- `grep -rn "Pattern" src/` before guessing.
- Build `make -j 8 config=debug` during development.
- Validation suite is ground truth.
- Maintainers (Grech, Pérez, Beck, Vinci) responsive on GitHub Discussions.
- For low-temperature-plasma / excimer work (`Collisions`, `Ionization`, MCC paths), most relevant code lives in `src/Collisions/` and `src/Ionization/`.
