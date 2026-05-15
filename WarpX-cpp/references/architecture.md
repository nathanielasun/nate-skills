# Architecture: AMReX foundation and WarpX class layout

This reference covers the code organization, the AMReX data structures WarpX is built on, and the principal WarpX classes. Read this when navigating the source tree, locating where a feature lives, or designing a new subsystem.

## Repository layout

All WarpX source code lives under `Source/` in the BLAST-WarpX repository. The subdirectory names are intentionally descriptive:

```
Source/
├── WarpX.{H,cpp}                      # The main class, owns everything
├── WarpXEvolve.cpp                    # PIC loop (Evolve, OneStep_nosub, OneStep_sub1)
├── Initialization/                    # Domain setup, MR levels, field/particle init
├── Parallelization/                   # MPI exchanges, guard cells, FillBoundary, SumBoundary
├── BoundaryConditions/                # Field BCs (PML, PEC, periodic, Silver-Mueller, ...)
├── EmbeddedBoundary/                  # EB geometry, ECT solver support
├── FieldSolver/
│   ├── FiniteDifferenceSolver/        # FDTD (Yee, CKC), templated by stencil algorithm
│   ├── SpectralSolver/                # PSATD (pseudo-spectral analytical time domain)
│   ├── ElectrostaticSolvers/          # Poisson-based ES (LabFrame, RelativisticExplicit, ...)
│   ├── MagnetostaticSolvers/          # Vector-potential magnetostatic
│   ├── HybridPICModel/                # Ohm's-law hybrid (kinetic ions, fluid electrons)
│   └── ImplicitSolvers/               # Theta-implicit, semi-implicit EM
├── Particles/
│   ├── WarpXParticleContainer.{H,cpp} # Base particle container
│   ├── MultiParticleContainer.{H,cpp} # Owner of all species containers (WarpX::mypc)
│   ├── PhysicalParticleContainer.{H,cpp}  # Plasma species (the main case)
│   ├── LaserParticleContainer.{H,cpp}     # Laser antenna particles
│   ├── PhotonParticleContainer.{H,cpp}    # Photons (QED)
│   ├── Pusher/                        # Boris, Vay, Higuera-Cary, copy-particle-to-buffer
│   ├── Deposition/                    # Current/charge deposition kernels
│   ├── Gather/                        # Field-gather from grid to particles
│   ├── Collision/                     # MCC, DSMC, Coulomb, fusion, nuclear (see collisions-mcc-dsmc.md)
│   ├── ParticleBoundaries/            # Reflection, absorption, periodic
│   ├── Sorting/                       # Per-cell sorting for collision modules
│   └── ParticleCreation/              # SmartCopy, SmartCreate, FilterCopyTransform
├── Diagnostics/                       # Full diagnostics, BTD, reduced diags
├── Utils/                             # Parser, IO, constants, profiler wrappers
├── Filter/                            # Bilinear, NCI-Godfrey filters
├── ablastr/                           # Shared library extracted from WarpX
│   └── fields/                        # MultiFabRegister — the field registry
└── Python/                            # pybind11 bindings (covered by warpx-python skill)
```

This is a guide, not exhaustive — subdirectories get added and renamed. Run `find Source -type d` on a current checkout for the truth.

**The PIC loop** is in `WarpX::Evolve` (file `Source/WarpXEvolve.cpp`). Its core iteration is `WarpX::OneStep_nosub` (no subcycling) or `WarpX::OneStep_sub1` (subcycled MR, method 1). When you want to understand "what happens at each time step," start in `OneStep_nosub` and follow the call graph.

## AMReX in 10 minutes

WarpX is built on AMReX. You cannot avoid AMReX abstractions, so internalize these:

### `amrex::Box`
A lightweight metadata class describing a rectangular region in *index space* (not physical space). It carries lower/upper integer indices in each dimension and a centering type (cell-centered, nodal, edge, face). Cheap to copy. Use `Box::contains()`, `Box::intersect()`, `Box::grow()`, `mfi.tilebox()` to obtain Boxes during loops.

### `amrex::BoxArray`
A collection of `Box`es covering a single AMR level. Combined with a `DistributionMapping` (which assigns each Box to an MPI rank), it defines the domain decomposition. You rarely touch `BoxArray` directly — you usually work with `MultiFab`s already built on top.

### `amrex::FArrayBox` (FAB)
A Fortran-ordered array of `amrex::Real` defined on a single `Box`, possibly with multiple components (`ncomp`). Represents the field data on one rank-owned subdomain. Don't loop over it manually — use the `Array4<Real>` view obtained via `fab.array()`.

### `amrex::MultiFab`
The workhorse field container. A `MultiFab` is a collection of FABs distributed across MPI ranks according to the `DistributionMapping`. It carries metadata about ghost cells (`nGrow`). All field data in WarpX (E, B, J, ρ, ϕ, ...) is stored as `MultiFab`s, typically owned by smart pointers in the `MultiFabRegister`.

### `amrex::ParticleContainer`
A collection of particles organized per-Box (and per-tile on CPU for cache efficiency, off on GPU). Particles within a tile are stored as:
- Struct-of-Arrays (SoA) of `amrex::ParticleReal` for floating-point per-particle data (positions, weight, momentum components, runtime real attributes).
- Struct-of-Arrays (SoA) of `int` for integer per-particle data (e.g., ionization level).
- A 64-bit unique ID per particle (40-bit id + 24-bit cpu sub-ID); used to track active/removed status.

The base AMReX template is `amrex::ParticleContainer_impl<SoAParticle<NReal, NInt>, NReal, NInt, Allocator, CellAssignor>`. WarpX wraps this in `WarpXParticleContainer`.

### Iterators
Two iterators are used pervasively:
- `amrex::MFIter` — iterates over local `FArrayBox`es in a `MultiFab` (with optional tiling).
- `amrex::ParIter` and its WarpX-derived `WarpXParIter` — iterates over local particle tiles in a `ParticleContainer`.

The crucial property: **MFIter and ParIter walk the same per-Box decomposition.** When you iterate over particles in a box, you can directly obtain the FABs for the field MultiFabs on that same box. This is how field-gather and current-deposition routines combine particle and grid data without communication.

### Mesh refinement structure
Every AMR level has a separate `BoxArray` and `DistributionMapping`. WarpX field MultiFabs are stored per-level in `amrex::Vector<...>` over `lev`. The patches:
- **Fine patch** (`_fp` suffix): the field at the actual resolution of level `lev`. Updated by the field solver.
- **Coarse patch** (`_cp`): same physical domain at the *coarser* resolution of `lev-1`. Bookkeeping for MR interpolation.
- **Auxiliary** (`_aux`): same resolution as `_fp`, but a copy that has been filtered/interpolated for field-gather. Particles gather from `_aux`, not `_fp`.
- **Coarse-aux** (`_cax`): only used at MR boundaries with buffers. Copy of the level-below `_aux` on the current level's grid.

For non-MR simulations (the common case in low-temperature-plasma work), only level 0 exists and only `_fp` and `_aux` are populated. Most code still iterates over `0 <= lev <= finest_level`, with `finest_level = 0` collapsing it to a single pass.

## WarpX class architecture

### `class WarpX`
Defined in `Source/WarpX.H` and `Source/WarpX.cpp`. The single owner of essentially all simulation state. Historically a Meyers singleton with `WarpX::GetInstance()`; new code is being weaned off this pattern by passing `WarpX&` or specific containers by reference.

Owns:
- All field MultiFabs (via the `MultiFabRegister`, accessible as `m_fields`).
- The `MultiParticleContainer` instance (`m_mypc` / historically `mypc`) — all particle species.
- The mesh geometry per level, `BoxArray`s, `DistributionMapping`s.
- All solver state (FDTD coefficients, PSATD spectral solver, Poisson solvers, ...).
- All collision handlers (via `MultiParticleContainer::m_collisionhandler`).
- Diagnostic objects, boundary condition state, moving window state, embedded boundary state.

Key entry points to remember:
- `WarpX::MakeWarpX()` — top-level setup (called once after construction).
- `WarpX::InitData()` — initialize fields and particles, MR levels, BCs.
- `WarpX::Evolve()` — main time-stepping loop.
- `WarpX::OneStep_nosub()` — one PIC step (no subcycling).
- `WarpX::AllocLevelMFs()` — allocate MultiFabs at a given level (called once per level).
- `WarpX::InitLevelData()` — initial values for fields at a level.
- `WarpX::EvolveE()` / `WarpX::EvolveB()` — field-solve advances (FDTD).

### `class MultiParticleContainer`
File: `Source/Particles/MultiParticleContainer.{H,cpp}`. Owned by `WarpX` as `m_mypc` (the legacy name `mypc` is still common). Holds a `Vector<std::unique_ptr<WarpXParticleContainer>>` called `allcontainers` — one per species.

Two kinds of methods:
- **Dispatchers**: iterate `allcontainers` and call the per-container function. Example: `MultiParticleContainer::Evolve` calls `WarpXParticleContainer::Evolve` for every species.
- **Multi-species coordinators**: handle physics that crosses species, like ionization product mapping (`mapSpeciesProduct`), QED Schwinger pair creation (`doQEDSchwinger`), reading input parameters that reference multiple species, and reaction product creation.

When you add a new physics module that creates particles of a different species than the parent particle (ionization, charge exchange, fusion), the cross-species bookkeeping lives in `MultiParticleContainer`, not in the per-species container.

### `class WarpXParticleContainer`
File: `Source/Particles/WarpXParticleContainer.{H,cpp}`. Base polymorphic class. Derives from `amrex::ParticleContainer_impl<SoAParticle<...>, ...>` (the exact template arguments have shifted across versions — check the current signature in `WarpXParticleContainer.H`).

Defines:
- Pure virtual functions (must be overridden): notably `Evolve()`.
- Concrete shared functionality: `DepositCurrent()`, `DepositCharge()`, particle initialization helpers, particle attribute management.

Subclassed by:
- `PhysicalParticleContainer` — the common case. Real plasma/beam/photon-target species. Has its own `Evolve` (field gather + push + deposit).
- `LaserParticleContainer` — non-physical antenna particles that emit a laser.
- `PhotonParticleContainer` — derives from `PhysicalParticleContainer`, specialized for photon push (no charge gather, different push).

### `class PhysicalParticleContainer`
The class you'll most often touch. Its `Evolve(...)` orchestrates the gather/push/deposit for a single species. Internally calls `PushPX(...)` to do the momentum update. The collision handlers live one level up (in `MultiParticleContainer`) and run *between* particle pushes.

The function signature for `Evolve` varies across versions and is long. Don't memorize it — look it up. The important thing is the *structure*: this is where field gather, current deposition, and particle push happen for one species per `MFIter`/`ParIter` step.

### `ablastr` library and `MultiFabRegister`
`ablastr` (Accelerator BLAST Library for Astro/accelerator Routines) is a code library that was extracted out of WarpX so it can be reused by other AMReX codes. Most importantly for field code: **field MultiFabs are now accessed through `ablastr::fields::MultiFabRegister`**, not as direct members of `WarpX`. The register is owned by `WarpX` and named `m_fields`.

To get a field MultiFab:

```cpp
// Get Ex on fine patch at level lev
amrex::MultiFab* Ex = m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{0}, lev);

// Get J as a 3-vector at level lev (returns ablastr::fields::VectorField)
auto J = m_fields.get_alldirs(warpx::fields::FieldType::current_fp, lev);
```

The `FieldType` enum is defined in the fields developer documentation and in the source under `ablastr::fields`. Common members: `Efield_fp`, `Bfield_fp`, `Efield_aux`, `Bfield_aux`, `Efield_cp`, `Bfield_cp`, `current_fp`, `rho_fp`, `phi_fp`, `F_fp` (divE cleaning), `G_fp` (divB cleaning), hybrid-Ohm fields (`hybrid_*`), embedded boundary fields (`distance_to_eb`, `edge_lengths`, `face_areas`).

When adding a new field type, add an enum value and register it through the standard registration flow — don't add a raw MultiFab member to `WarpX`.

## Dimensionality

WarpX supports five geometries, selected at compile time via `WarpX_DIMS`:
- `3` — 3D Cartesian.
- `2` — 2D Cartesian (typically `(x, z)` with `y` translationally invariant).
- `1` — 1D Cartesian (typically `z`).
- `RZ` — 2D cylindrical with azimuthal mode expansion (multi-mode).
- `RCYLINDER` — 1D cylindrical (radial only).
- `RSPHERE` — 1D spherical (radial only).

The selected dimension is exposed at compile time via the macro `AMREX_SPACEDIM` (and `WARPX_DIM_*` defines like `WARPX_DIM_RZ`, `WARPX_DIM_3D`). Code that needs to differ between dimensions uses `#if defined(WARPX_DIM_RZ)` style guards.

**Important:** WarpX produces *separate binaries per dimension*. The build system can compile multiple dimensions in one invocation (`WarpX_DIMS="1;2;3;RZ"`), but each binary is monolithic to one geometry. When you write dimension-dependent code, this is why the macros are compile-time rather than runtime.

The developer documentation page `Dimensionality` (currently sparse) lists which features support which geometries. Notably: embedded boundaries are *not* supported in RZ/RCYLINDER/RSPHERE; PSATD has restrictions in RZ; filters have restrictions in cylindrical geometries.

## How a new physics module slots in

If you're adding a new piece of physics, the question is: where does it sit relative to the canonical PIC loop?

| Physics | Typical location |
|---|---|
| New particle pusher | `Source/Particles/Pusher/`; selected by enum in `WarpXAlgorithmSelection.H` |
| New field solver | `Source/FieldSolver/<solver-family>/`; selected by `WarpX::electromagnetic_solver_id` |
| New collision/scattering process | `Source/Particles/Collision/...` (see `collisions-mcc-dsmc.md`) |
| New particle initialization profile | `Source/Initialization/` and `PlasmaInjector` machinery |
| New boundary condition (field) | `Source/BoundaryConditions/` |
| New boundary condition (particle) | `Source/Particles/ParticleBoundaries/` |
| New diagnostic | `Source/Diagnostics/`; reduced diags under `ReducedDiags/` |
| New filter for fields | `Source/Filter/` |
| New runtime particle attribute | `Source/Particles/` plus `MultiParticleContainer::ReadParameters` |

For each of these, find the closest existing example, copy and adapt. Trying to fit a new module without precedent is much more error-prone than mimicking an existing one. The team has spent years tuning the patterns; respect them.

## Some pitfalls to be aware of

- **Tile boundaries.** When iterating with `MFIter(*mf, TilingIfNotGPU())`, you get sub-Box tiles on CPU and full Boxes on GPU. Some assumptions break across tiles. Always use `mfi.tilebox()` (or its variant for the field's nodal type) inside the loop body, not `mfi.fabbox()` (which includes ghost cells).
- **Ghost cells.** Most field operations require valid data in some number of ghost cells. After updating a field, you typically need `FillBoundary` to refresh ghosts. Skipping this is a common bug.
- **Subcycling.** When subcycling MR is on, particle pushes and field solves on fine levels run multiple times per coarse step. New code that assumes "one step is one step" may be wrong on fine levels.
- **The `_aux` distinction.** Particles gather from `Efield_aux`/`Bfield_aux`, not directly from `Efield_fp`/`Bfield_fp`. The `aux` MultiFabs are kept synchronized in `WarpX::UpdateAuxilaryData`. Don't gather from `_fp` even if it "should be the same" — it may not be after filtering or interpolation.
- **Singletons under the hood.** Even though `WarpX::GetInstance()` is being deprecated, certain global state (random number generators, profiling, MPI) is still global. Resetting these in tests is non-trivial.
- **Real vs ParticleReal.** Field precision (`amrex::Real`) and particle precision (`amrex::ParticleReal`) are independently configurable. Mixing them in arithmetic can produce silent precision loss. The `_rt` and `_prt` literal suffixes exist to keep these straight.

When in doubt about how a class or function works, the Doxygen reference is authoritative. Run a search there before guessing, and check the actual source file before relying on remembered behavior.
