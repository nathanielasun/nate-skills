---
name: warpx-cpp-mini
description: Use for any WarpX C++ codebase work — modifying internals, adding collision/scattering processes, particle pushers, field solvers, diagnostics; debugging C++ source under `Source/`; navigating the AMReX-based architecture; writing tests; preparing PRs. Trigger on "WarpX", "BLAST-WarpX", `WarpXParticleContainer`, `MultiParticleContainer`, `MultiFab` in WarpX context, MCC/DSMC, Boris/Vay/Higuera-Cary pushers, embedded boundaries, `Source/Particles/Collision`, electrostatic/electromagnetic/hybrid solvers, or paths under `Source/`. Optimized for low-temperature plasma. This is the ultra-compact variant tuned for small models; for PICMI Python scripts use warpx-python-mini.
---

# WarpX C++ (Mini)

For modifying WarpX C++ code. Not for PICMI input scripts (use `warpx-python-mini`).

WarpX is built on **AMReX**. The rules below are not stylistic — code that breaks them won't build, won't run on GPU, or will be rejected.

## ALWAYS verify APIs before using them

WarpX changes monthly. Class names, method signatures, and enum values shift. **Never invent an API.** Check:

1. **`AGENTS.md`** / **`CLAUDE.md`** at repo root
2. **Source code** in the user's checkout
3. **Doxygen**: `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/`

If you're not sure something exists, search the source or Doxygen first.

## The 8 hard rules

| # | Rule |
|---|---|
| 1 | **Field data**: `amrex::MultiFab` only. Access via `WarpX::m_fields` (the `MultiFabRegister`). Never add a raw `std::unique_ptr<amrex::MultiFab>` member to `WarpX`. |
| 2 | **Particle data**: `WarpXParticleContainer` (or subclass) only. Never write a custom container. |
| 3 | **Loops**: `MFIter` + `amrex::ParallelFor` for fields. `ParIter` + `ParallelFor` for particles. Never hand-roll `for (i, j, k)` in hot paths. |
| 4 | **Lambda captures**: copy member variables to local variables BEFORE `ParallelFor`. Never capture `this` into a device lambda. |
| 5 | **FP literals**: `0.0_rt` for `amrex::Real`; `0.0_prt` for `amrex::ParticleReal`. Never bare `0.0` or `0.0f`. |
| 6 | **Namespace**: full `amrex::` qualification. Never `using namespace amrex;` at file scope. |
| 7 | **New source files**: update BOTH `CMakeLists.txt` AND `Make.package` in the directory. Headers go in neither. |
| 8 | **Singletons**: avoid `WarpX::GetInstance()` in new code. Pass `WarpX&` or the specific container by reference. |

## What to read

| Task | Read |
|---|---|
| Writing C++ code — loops, captures, field/particle access | `references/patterns.md` |
| **MCC, DSMC, Coulomb, ionization (low-temp plasma)** | `references/collisions.md` |
| Building, running, testing | `references/build-and-run.md` |

## Code map (where things live)

| Subsystem | Location |
|---|---|
| Main class | `Source/WarpX.{H,cpp}` |
| PIC loop (Evolve, OneStep_nosub) | `Source/WarpXEvolve.cpp` |
| Field solvers | `Source/FieldSolver/{FiniteDifferenceSolver,SpectralSolver,ElectrostaticSolvers,MagnetostaticSolvers,HybridPICModel,ImplicitSolvers}/` |
| Particle containers | `Source/Particles/{WarpX,Multi,Physical,Laser,Photon}ParticleContainer.{H,cpp}` |
| Particle pushers | `Source/Particles/Pusher/` |
| Current/charge deposition | `Source/Particles/Deposition/` |
| Field gather | `Source/Particles/Gather/` |
| Collisions (MCC, DSMC, Coulomb, fusion) | `Source/Particles/Collision/` |
| Field/particle BCs | `Source/BoundaryConditions/`, `Source/Particles/ParticleBoundaries/` |
| Embedded boundaries | `Source/EmbeddedBoundary/` |
| Diagnostics | `Source/Diagnostics/` (reduced under `ReducedDiags/`) |
| Field registry | `Source/ablastr/fields/` |
| Initialization | `Source/Initialization/` |

## The PIC loop (per step)

`WarpX::OneStep_nosub` in `Source/WarpXEvolve.cpp`:

1. **Field gather** — from `Efield_aux`/`Bfield_aux` to particles
2. **Particle push** — Boris/Vay/Higuera-Cary
3. **Current deposition** — particles → J on grid
4. **Field solve** — `EvolveE`/`EvolveB` (FDTD/PSATD) or elliptic (ES/MS/hybrid)
5. **Collisions / reactions** — MCC, DSMC, Coulomb, ionization, fusion
6. **Boundaries / diagnostics / filtering / load-balancing**

## Critical AMReX facts

- **`amrex::MultiFab`**: holds all field data. Accessed via `WarpX::m_fields`.
- **`amrex::ParticleContainer`**: holds particles, organized per-Box per-tile.
- **`amrex::MFIter` and `amrex::ParIter`**: iterate over local Boxes/tiles. They walk the SAME decomposition — you can access fields and particles for the same Box together.
- **Patches**: `_fp` (fine, what solver updates), `_aux` (filtered, what particles gather from), `_cp` (coarse, MR), `_cax` (MR boundary). **Particles gather from `_aux`, not `_fp`.**

## Critical pitfalls

- **Capturing `this` in `ParallelFor`** copies the entire `WarpX` object to GPU. Copy members to locals first. The `m_` prefix on members exists to make this mistake visible.
- **`mfi.tilebox()` vs `mfi.fabbox()`**: use `tilebox()` for the valid region. `fabbox()` includes ghost cells.
- **After modifying particle positions**: call `pc.Redistribute()` to fix box ownership.
- **After updating a field**: call `FillBoundary()` to refresh ghosts.
- **Real vs ParticleReal**: field precision and particle precision are independent. Mixing them silently loses precision. Use `_rt` / `_prt` literals.

## Working pattern

1. **Find existing precedent.** Almost every new physics resembles existing physics. Copy and adapt.
2. **Read the relevant reference file.**
3. **Compile for one dimension first** (e.g., `WarpX_DIMS="3"`). Generalize later.
4. **Add a test.** WarpX requires tests for new features.
5. **Run `clang-tidy`** before pushing.
6. **One feature per PR.**

## Links

- Source: `https://github.com/BLAST-WarpX/warpx`
- Docs: `https://warpx.readthedocs.io/`
- Doxygen: `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/`
- AMReX: `https://amrex-codes.github.io/amrex/doxygen/`
- Discussions: `https://github.com/BLAST-WarpX/warpx/discussions`
- Cross sections: `https://github.com/BLAST-WarpX/warpx-data`
