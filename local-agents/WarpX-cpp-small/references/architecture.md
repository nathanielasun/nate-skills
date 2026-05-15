# Architecture: AMReX, WarpX classes, code map

This is for navigating the source and understanding the type hierarchy. Compact reference — for deeper background, the canonical Doxygen and AMReX docs are linked at the bottom.

## Repository layout

All C++ source under `Source/`:

```
Source/
├── WarpX.{H,cpp}                  # Main class, owns everything
├── WarpXEvolve.cpp                # PIC loop (Evolve, OneStep_nosub, OneStep_sub1)
├── Initialization/                # Domain setup, MR levels, field/particle init
├── Parallelization/               # MPI exchanges, FillBoundary, SumBoundary
├── BoundaryConditions/            # PML, PEC, periodic, Silver-Mueller
├── EmbeddedBoundary/              # EB geometry, ECT solver support
├── FieldSolver/
│   ├── FiniteDifferenceSolver/    # FDTD (Yee, CKC)
│   ├── SpectralSolver/            # PSATD
│   ├── ElectrostaticSolvers/      # Poisson-based ES (LabFrame, RelativisticExplicit)
│   ├── MagnetostaticSolvers/      # Vector-potential MS
│   ├── HybridPICModel/            # Ohm's-law hybrid (kinetic ions, fluid electrons)
│   └── ImplicitSolvers/           # Theta-implicit, semi-implicit EM
├── Particles/
│   ├── WarpXParticleContainer.{H,cpp}      # Base particle container
│   ├── MultiParticleContainer.{H,cpp}      # Owner of all species (WarpX::mypc)
│   ├── PhysicalParticleContainer.{H,cpp}   # Plasma species
│   ├── LaserParticleContainer.{H,cpp}      # Laser antenna
│   ├── PhotonParticleContainer.{H,cpp}     # Photons (QED)
│   ├── Pusher/                    # Boris, Vay, Higuera-Cary, copy-to-buffer
│   ├── Deposition/                # Current/charge deposition kernels
│   ├── Gather/                    # Field-gather from grid to particles
│   ├── Collision/                 # MCC, DSMC, Coulomb, fusion (see collisions.md)
│   ├── ParticleBoundaries/        # Reflection, absorption, periodic
│   ├── Sorting/                   # Per-cell sorting for collision modules
│   └── ParticleCreation/          # SmartCopy, SmartCreate, FilterCopyTransform
├── Diagnostics/                   # Full diags, BTD, reduced diags
├── Utils/                         # Parser, IO, constants, profiler wrappers
├── Filter/                        # Bilinear, NCI-Godfrey filters
├── ablastr/fields/                # MultiFabRegister (field registry)
└── Python/                        # pybind11 bindings (covered by warpx-python)
```

Run `find Source -type d` on a current checkout for exact subdirectories.

## AMReX core types you must know

### `amrex::Box`
Lightweight rectangular index-space region. Has lower/upper integer indices and centering type (cell-centered, nodal, edge, face). Cheap to copy. Use `Box::contains()`, `Box::intersect()`, `Box::grow()`, `mfi.tilebox()`.

### `amrex::BoxArray` + `amrex::DistributionMapping`
The decomposition. `BoxArray` is the collection of Boxes; `DistributionMapping` assigns each Box to an MPI rank. You usually work with `MultiFab`s built on top, not these directly.

### `amrex::FArrayBox` (FAB)
Fortran-ordered array of `amrex::Real` on one `Box`, possibly multi-component (`ncomp`). The field data on one rank's subdomain. DO NOT loop manually — get `Array4<Real>` view via `fab.array()`.

### `amrex::MultiFab`
The workhorse field container. A `MultiFab` is a collection of FABs distributed across MPI ranks (per `DistributionMapping`). Carries ghost-cell metadata (`nGrow`). ALL field data in WarpX (E, B, J, ρ, ϕ, …) is `MultiFab`, owned by `MultiFabRegister`.

### `amrex::ParticleContainer`
Particles organized per-Box (and per-tile on CPU). Particles in a tile stored as:
- SoA of `amrex::ParticleReal` for floats (positions, weight, momentum, runtime real attributes)
- SoA of `int` for integer attributes (e.g., ionization level)
- 64-bit unique ID per particle (40-bit id + 24-bit cpu sub-ID); used for active/removed tracking

WarpX wraps the AMReX template in `WarpXParticleContainer`.

### Iterators
- `amrex::MFIter` — iterates local `FArrayBox`es in a `MultiFab` (with optional tiling)
- `amrex::ParIter` / `WarpXParIter` — iterates local particle tiles

**Critical property: MFIter and ParIter walk the same per-Box decomposition.** When iterating particles in a box, you can directly access FABs for fields on that same box. This is how gather/deposit avoid communication.

## Mesh refinement patches

Every AMR level has its own `BoxArray` + `DistributionMapping`. Field MultiFabs stored per-level in `amrex::Vector<...>` over `lev`. The suffixes:

| Suffix | Meaning |
|---|---|
| `_fp` (fine patch) | Field at actual resolution of level `lev`. Updated by solver. |
| `_cp` (coarse patch) | Same physical domain at `lev-1` resolution. MR bookkeeping. |
| `_aux` (auxiliary) | Same as `_fp` but filtered/interpolated for field-gather. **Particles gather from `_aux`, not `_fp`.** |
| `_cax` (coarse-aux) | Only at MR boundaries with buffers. |

Non-MR simulations (typical for low-temp plasma): only level 0, only `_fp` and `_aux` populated. Code still iterates `0 <= lev <= finest_level` — collapses to one pass.

## WarpX class hierarchy

### `class WarpX` (Source/WarpX.H)
The owner of everything. Historically a Meyers singleton (`WarpX::GetInstance()`); new code is moving away from this pattern.

Owns:
- All field MultiFabs (via `MultiFabRegister`, accessible as `m_fields`)
- `MultiParticleContainer` instance (`m_mypc` / legacy `mypc`)
- Mesh geometry per level, BoxArrays, DistributionMappings
- All solver state

Key entry points:
- `WarpX::MakeWarpX()` — top-level setup (once after construction)
- `WarpX::InitData()` — initialize fields, particles, MR levels, BCs
- `WarpX::Evolve()` — main time-stepping loop
- `WarpX::OneStep_nosub()` — one PIC step
- `WarpX::AllocLevelMFs()` — allocate MultiFabs at a level
- `WarpX::EvolveE()` / `WarpX::EvolveB()` — FDTD advances

### `class MultiParticleContainer` (Source/Particles/MultiParticleContainer.H)
Owned by `WarpX` as `m_mypc`. Holds `Vector<std::unique_ptr<WarpXParticleContainer>> allcontainers` — one per species.

Two roles:
- **Dispatcher**: iterates `allcontainers`, calls per-container function (`Evolve`)
- **Cross-species coordinator**: ionization product mapping, QED Schwinger, multi-species reactions

Cross-species bookkeeping (new ions from ionization, new neutrals from charge exchange) lives HERE, not in per-species containers.

### `class WarpXParticleContainer`
Base polymorphic class. Derives from `amrex::ParticleContainer_impl<SoAParticle<...>, ...>`.

Pure virtual: `Evolve()`. Concrete shared: `DepositCurrent()`, `DepositCharge()`, init helpers.

Subclassed by:
- `PhysicalParticleContainer` — real plasma species (the common case). Owns its own `Evolve` (gather + push + deposit).
- `LaserParticleContainer` — non-physical antenna particles emitting a laser.
- `PhotonParticleContainer` — derives from PhysicalParticleContainer; specialized for photons.

### `ablastr::fields::MultiFabRegister`
Field MultiFabs are accessed through `WarpX::m_fields`, not as direct members of `WarpX`.

```cpp
// Get Ex on fine patch at level lev
amrex::MultiFab* Ex = m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{0}, lev);

// Get J as a 3-vector at level lev
auto J = m_fields.get_alldirs(warpx::fields::FieldType::current_fp, lev);
```

Common `FieldType` values: `Efield_fp`, `Bfield_fp`, `Efield_aux`, `Bfield_aux`, `Efield_cp`, `Bfield_cp`, `current_fp`, `rho_fp`, `phi_fp`, `F_fp` (divE cleaning), `G_fp` (divB), hybrid-Ohm fields, EB fields (`distance_to_eb`, `edge_lengths`, `face_areas`).

**Adding a new field type**: add a `FieldType` enum value and register through standard flow. Don't add raw `MultiFab` members to `WarpX`.

## Dimensionality

Set at compile time via `WarpX_DIMS`:
- `3` — 3D Cartesian
- `2` — 2D Cartesian (typically (x, z) with y translationally invariant)
- `1` — 1D Cartesian
- `RZ` — 2D cylindrical with multi-mode azimuthal expansion
- `RCYLINDER` — 1D cylindrical
- `RSPHERE` — 1D spherical

Exposed at compile time as `AMREX_SPACEDIM` plus `WARPX_DIM_*` defines. Dimension-dependent code uses `#if defined(WARPX_DIM_RZ)` guards.

**Separate binaries per dimension.** Can build multiple in one invocation (`WarpX_DIMS="1;2;3;RZ"`), but each binary is monolithic. EB not supported in RZ/RCYLINDER/RSPHERE; PSATD has restrictions in RZ.

## Where new physics goes

| Adding... | Location |
|---|---|
| New particle pusher | `Source/Particles/Pusher/`; enum in `WarpXAlgorithmSelection.H` |
| New field solver | `Source/FieldSolver/<family>/`; `WarpX::electromagnetic_solver_id` |
| New collision/scattering | `Source/Particles/Collision/...` (see collisions.md) |
| New particle init profile | `Source/Initialization/` + `PlasmaInjector` |
| New field BC | `Source/BoundaryConditions/` |
| New particle BC | `Source/Particles/ParticleBoundaries/` |
| New diagnostic | `Source/Diagnostics/` (reduced under `ReducedDiags/`) |
| New filter | `Source/Filter/` |
| New runtime particle attribute | `Source/Particles/` + `MultiParticleContainer::ReadParameters` |

For each: find closest existing example, copy and adapt. Don't design new modules without precedent.

## Common pitfalls

- **Tile boundaries**: `MFIter(*mf, TilingIfNotGPU())` gives sub-Box tiles on CPU, full Boxes on GPU. Use `mfi.tilebox()` (not `mfi.fabbox()` which has ghost cells) inside loops.
- **Ghost cells**: most field ops need valid ghost data. After updating a field, call `FillBoundary`. Skipping this is a frequent bug.
- **Subcycling**: with subcycled MR, particle pushes and field solves on fine levels run multiple times per coarse step. Don't assume "one step is one step" on fine levels.
- **`_aux` distinction**: particles gather from `Efield_aux`/`Bfield_aux`. After filtering or interpolation `_fp` and `_aux` may differ. Don't gather from `_fp`.
- **Real vs ParticleReal**: `amrex::Real` (field) and `amrex::ParticleReal` (particle) are independently configurable. Mixing them in arithmetic risks silent precision loss. Use `_rt` / `_prt` literal suffixes.
- **Singletons**: `WarpX::GetInstance()` is deprecated, but RNG/profiling/MPI state is still global; tests can't fully reset.

## Canonical references

- WarpX source: `https://github.com/BLAST-WarpX/warpx`
- WarpX Doxygen: `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/`
- AMReX docs: `https://amrex-codes.github.io/amrex/docs_html/`
- AMReX Doxygen: `https://amrex-codes.github.io/amrex/doxygen/`

Verify any API against Doxygen before relying on it.
