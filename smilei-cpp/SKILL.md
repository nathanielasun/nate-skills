---
name: smilei-cpp
description: Use when working in the Smilei C++ source code — modifying internals, adding physics modules (new collision types, new ionization models, new solvers), debugging the BLAST-Smilei codebase, navigating its architecture (patches, Hilbert ordering, MPI/OpenMP/GPU layers), writing tests, or preparing pull requests. Trigger on "Smilei C++", "BLAST-Smilei", `Patch`, `VectorPatch`, `Species` class (not the Python namelist block), `ParticleContainer`, `ElectroMagn`, `Pusher`, `Interpolator`, `Projector`, `BinaryProcesses`, `Ionization` source class, `Maxwell` solvers in code, source paths under `src/`, references to `picparams.h`, `simWindow`, or any non-trivial code modification. Optimized for low-temperature-plasma / excimer work (Collisions, MCC, Ionization paths) but covers the full C++ codebase. For Python namelist authoring or running Smilei from a user perspective, use smilei-namelist instead.
---

# Smilei C++ Codebase Skill (frontier tier)

This skill covers work *inside* the Smilei C++ source — modifying the codebase, adding capability, debugging, and contributing back. The audience is researchers extending Smilei, not end-users writing namelists. For namelist authoring see `smilei-namelist`.

Smilei is a maintained, well-engineered C++ codebase: about 200k LOC of header-rich, template-light C++17. It uses no third-party PIC framework (unlike WarpX which is built on AMReX) — Smilei's data structures, parallelism, and I/O are all in-tree. This means:

- **Smaller code surface to learn** than WarpX
- **Faster to navigate** — everything is in `src/`
- **More invariants are local** — fewer base-class hierarchies to chase
- **Less prior infrastructure** — if you want a new feature, you usually write it from scratch, not from a framework primitive

Tested against Smilei v5.1.

## Hard rules for source modification

1. **Read the architecture reference before touching the source.** Smilei's patch/Hilbert structure is the single most important thing to internalize; misunderstand it and your changes deadlock or corrupt particle data. See `references/architecture.md`.

2. **Match existing C++ style.** No std::format, no concepts, no ranges. The code targets C++17 across diverse compilers (GCC 7+, Intel, Clang, IBM). Use `#pragma omp` directives for parallelism (not C++ atomics, not std::execution). Read 3 existing files in the area you're modifying before writing new code.

3. **Particle data is structure-of-arrays.** `Particles` (`src/Particles/Particles.h`) holds `position[d]`, `momentum[d]`, `weight`, `charge`, `chi`, etc. as separate `std::vector<double>`. Never assume an AoS layout. See `particles-and-species.md`.

4. **Patches are the parallelism unit, not cells.** MPI domain decomposition is by patch (Hilbert-ordered). OpenMP parallelism is over patches within a rank. Cells are interior to a patch — there is no per-cell MPI structure. Modifying patch boundaries is a high-stakes operation.

5. **Ghost cells follow Maxwell-solver order.** Standard Yee uses 2 ghost cells per face; high-order solvers (Lehe, Cowan, Bouchard) use more. Solver code must respect the per-solver ghost-cell count. See `solvers-and-fields.md`.

6. **GPU paths use OpenMP target offload or OpenACC (compile-time choice).** GPU-aware code paths exist in parallel to CPU paths — modifying one requires keeping the other in sync. The `_GPU_KERNEL_` macros denote GPU-side code. See `architecture.md` §11.

7. **Diagnostics are written to HDF5 in a versioned schema.** Adding a new diagnostic or modifying an existing one requires bumping the schema and ensuring `happi` knows the new layout. See `diagnostics-source.md`.

8. **Random-number streams are per-patch.** Each patch carries its own `Random` instance (`src/Random/`). Don't call `std::rand()` or `<random>` generators from physics modules — use the patch's RNG to ensure reproducibility with `Main.random_seed`.

9. **`reference_angular_frequency_SI` controls dimensional code paths.** Collisions, ionization, radiation reaction, multiphoton Breit-Wheeler — these query `params.reference_angular_frequency_SI` (`src/Params/Params.h`). If you add a physics module needing SI conversion, follow the same pattern.

10. **Tests are mandatory for new physics modules.** Smilei has a `validation/` directory with reference outputs. A new physics module needs (a) a test namelist exercising it and (b) a reference output to compare against (often within 1–2% of an analytical solution where available).

## Voice and style for C++ work

When discussing or generating Smilei C++ code:

- **Prefer pointing the user to existing patterns.** "Look at `src/Collisions/CollisionalIonization.cpp` — it's the cleanest example of a per-species physics module hooked into the per-patch loop."
- **Be explicit about which side of the CPU/GPU split you're on.** Always state whether code lives in `_GPU_KERNEL_` paths or CPU-only paths.
- **Cite line ranges or function names, not abstractions.** "The Boris pusher lives in `src/Pusher/PusherBoris.cpp::operator()`, called from `VectorPatch::dynamics`" beats "the pusher is in the pusher class."
- **Never invent class names or method signatures.** If unsure, ask the user to confirm against the current source. Smilei's API has drifted across releases.

## When to use each reference

| If the question is about... | Read |
|---|---|
| Code architecture, build system (`Makefile`, machine files), directory layout, MPI/OpenMP/GPU layers, patch decomposition, Hilbert ordering, load balancing internals | `references/architecture.md` |
| `Particles` SoA layout, `Species` class, `ParticleContainer` semantics, particle pushers (Boris/Vay/HC), particle BC implementation, injection | `references/particles-and-species.md` |
| Maxwell solvers (Yee/Lehe/Cowan/Bouchard) internals, `ElectroMagn` / `ElectroMagnAM`, current deposition, field interpolation, PML implementation | `references/solvers-and-fields.md` |
| `Collisions` / `CollisionalIonization` source — the Pérez 2012 implementation, pairing, intra-species detection, the `BinaryProcesses` wrapper, screening, nuclear reactions | `references/collisions-source.md` |
| Field ionization source — ADK/PPT/BSI/envelope-averaged code paths, Monte-Carlo machinery, adding a new ionization model, the `from_rate` hook | `references/ionization-source.md` |
| Diagnostic implementation, HDF5 schema, adding a new diagnostic, integration with `happi`, `SmileiMPI::dump` patterns | `references/diagnostics-source.md` |

## Workflow patterns

**Adding a new physics module** (e.g., a new collision type, a new ionization model, a new pusher):

1. Read `references/architecture.md` to understand the per-patch loop.
2. Read the closest existing analog. For a new collision type, that's `src/Collisions/CollisionalIonization.cpp`. For a new pusher, `src/Pusher/PusherVay.cpp`.
3. Identify the namelist hook — most modules are selected by a string in the namelist. Find where the string is parsed in `src/Params/PicParams.cpp` and where the factory dispatches it.
4. Implement the module with both CPU and GPU code paths (or note that GPU is omitted with a TODO).
5. Add a test in `validation/` with a reference output.
6. Run the validation suite locally before submitting.

**Debugging the codebase:**

1. Compile with `make config=debug` for symbols and assertions.
2. Use `mpirun -np 1 ./smilei` first; many MPI hangs are masked by `-np 1`.
3. For race conditions, `make config=debug_openmp` adds OpenMP tools support.
4. For GPU debugging, `compute-sanitizer` (NVIDIA) or `roc-gdb` (AMD) are essential.

**Preparing a pull request:**

1. Fork and branch from `develop`, not `master`.
2. Run the local validation suite — passing the full test set is the bar.
3. Update `CHANGELOG.md`.
4. Match existing style (clang-format config in repo root).
5. For physics changes, include a benchmark against an analytical or published result.

## Directory layout at a glance

```
Smilei/
├── src/
│   ├── Diagnostic/      DiagScalar, DiagFields, DiagProbe, ..., HDF5 wrappers
│   ├── ElectroMagn/     EMFields, Maxwell solvers, PML
│   ├── ElectroMagnBC/   EM boundary conditions per geometry
│   ├── Field/           Field, cField, complex field classes
│   ├── Interpolator/    Field interpolation (Esirkepov, etc.)
│   ├── Ionization/      Field ionization models (Tunnel, PPT, BSI, ...)
│   ├── Collisions/      Binary collisions, CollisionalIonization, nuclear reactions
│   ├── Laser/           Laser injection (LaserGaussian2D, AM, envelope, ...)
│   ├── Merging/         Particle merging (for high-density regimes)
│   ├── MovWindow/       Moving simulation window
│   ├── Params/          Namelist parameter parsing
│   ├── Particles/       Particles class (SoA), Species, ParticleData
│   ├── Patch/           Patch, VectorPatch — the parallel decomposition
│   ├── Profiles/        Density / temperature profile evaluation
│   ├── ParticleBC/      Particle boundary conditions
│   ├── Projector/       Current deposition (Esirkepov)
│   ├── Pusher/          PusherBoris, PusherVay, PusherHigueraCary, PusherPhoton, ...
│   ├── Radiation/       Synchrotron radiation, photon emission
│   ├── Random/          Per-patch random number generators
│   ├── Smilei/          Top-level orchestration (main loop)
│   ├── SmileiMPI/       MPI domain decomposition, Hilbert ordering
│   └── Tools/           Utilities, timers, HDF5 wrappers
├── happi/               Python post-processing module
├── benchmarks/          Reference simulations for validation
├── validation/          Test suite (run via `python -m validation.run`)
├── scripts/
│   ├── compile_tools/   Machine files, compiler-specific flags
│   └── ...
├── doc/                 Sphinx documentation source
├── CHANGELOG.md
└── Makefile
```

## The main loop — orientation

`src/Smilei/Smilei.cpp` is the entry point. The PIC loop in pseudo-code:

```cpp
for (timestep = 0; timestep < n_timesteps; ++timestep) {
    vectorPatches.dynamics(...);       // particle push + current deposition (per-patch, OMP)
    vectorPatches.sumDensities(...);    // current densities communicated across patches
    vectorPatches.solveMaxwell(...);    // Maxwell solver (per-patch, OMP)
    vectorPatches.runAllDiags(...);    // diagnostics
    if (loadBalance.should()) {
        vectorPatches.loadBalance(...); // Hilbert-ordered redistribution
    }
}
```

Each `VectorPatch::*` method dispatches OpenMP-parallel over the patches owned by the current MPI rank. The patch-level methods then dispatch further (some are GPU kernels, others CPU-only).

When in doubt about where a piece of physics lives, look in `VectorPatch::dynamics` and follow the call chain.

## Hard truths about the codebase

- **Macros and `#ifdef`s shape the build.** GPU vs CPU, ZeroMQ vs MPI vs serial, 2D vs 3D vs AM — many code paths are compile-time-selected. Read `src/Tools/SmileiPP.h` and look for `_GPU_KERNEL_`, `SMILEI_USE_NUMPY`, `SMILEI_AMR`, etc.
- **Headers do a lot of work.** Many functions are inlined in `.h` files for performance. Don't assume `.cpp` is the only source of definitions.
- **Geometry is templated implicitly.** 1D/2D/3D and AM cylindrical have separate source files for many classes (e.g., `ElectroMagn1D.cpp`, `ElectroMagn2D.cpp`, etc.). When modifying a class, you usually need to modify all geometry variants.
- **OpenMP and MPI are interleaved.** OpenMP parallelizes over patches within a rank; MPI handles inter-rank communication. The interleaving is explicit (collective MPI calls outside OMP regions; OMP for-loops inside MPI sections).

## Closing notes

- For unfamiliar source files, run `grep -rn "PatternOfInterest" src/` before guessing.
- Build with `make -j 8 config=debug` for assertions during development.
- The validation suite is the ground truth — passing it is the contract.
- Smilei maintainers (Mickaël Grech, Frédéric Pérez, Arnaud Beck, Tommaso Vinci) are responsive on the GitHub Discussions page.
