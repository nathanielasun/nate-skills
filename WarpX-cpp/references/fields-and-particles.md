# Fields and particles: data structures, loops, and GPU-portable kernels

This reference covers the field and particle data structures in detail, plus the canonical iteration patterns that you must use for any hot-path code in WarpX. Read this when modifying field solvers, current/charge deposition, field-gather, particle pushers, or anything that loops over grid cells or particles.

The single most important lesson: **never hand-roll nested `for` loops over cells or particles**. Use the iterator + `amrex::ParallelFor` pattern documented here. Code that doesn't follow this pattern won't be accepted in review, won't vectorize on CPU, and won't run on GPU.

## Table of contents

1. [Fields: storage and access](#fields-storage-and-access)
2. [Field MultiFab iteration pattern](#field-multifab-iteration-pattern)
3. [Field solver implementation pattern](#field-solver-implementation-pattern)
4. [Particles: container hierarchy and storage layout](#particles-container-hierarchy-and-storage-layout)
5. [Particle iteration pattern](#particle-iteration-pattern)
6. [Field-gather and current-deposition: linking particles and fields](#field-gather-and-current-deposition-linking-particles-and-fields)
7. [Particle attributes: built-in and runtime](#particle-attributes-built-in-and-runtime)
8. [GPU portability: rules that aren't optional](#gpu-portability-rules-that-arent-optional)

## Fields: storage and access

Fields in WarpX are owned by the `MultiFabRegister` accessible from `WarpX` as `m_fields`. Each field is identified by a `warpx::fields::FieldType` enumerator and (for vector fields) a `Direction` index. The register replaces older direct-member access patterns; new code should use it exclusively.

### Common field types

The `warpx::fields::FieldType` enum is large but a handful of values cover most cases:

- **`Efield_fp`, `Bfield_fp`** — electromagnetic fields at the actual resolution of each level. Updated by the field solver each step.
- **`Efield_aux`, `Bfield_aux`** — fields the particles gather from. Synchronized from `_fp` (and `_cp` if MR is active) by `WarpX::UpdateAuxilaryData()`.
- **`current_fp`** — current density `J`, the source for the field solver. Deposited by particles each step.
- **`rho_fp`** — charge density `ρ`, used by electrostatic solvers and for current correction (Marder/divE cleaning).
- **`phi_fp`** — scalar potential, output of the Poisson solver in electrostatic mode.
- **`F_fp`, `G_fp`** — divergence-cleaning auxiliary fields for Maxwell's equations.
- **`Efield_fp_external`, `Bfield_fp_external`** — externally applied fields loaded from an openPMD file.
- **`distance_to_eb`, `edge_lengths`, `face_areas`** — embedded boundary geometry data.
- **`hybrid_*`** — fields used by the hybrid kinetic-fluid Ohm's law solver.
- **`pml_*_fp/cp`** — perfectly matched layer field components.

With mesh refinement active, each field type generally has both `_fp` and `_cp` variants (and `_aux`/`_cax` for `E` and `B`). Without MR, only `_fp` and `_aux` exist effectively.

### Access pattern

To get a `MultiFab*` for a scalar field or one component of a vector field:

```cpp
using namespace warpx::fields;
using ablastr::fields::Direction;

amrex::MultiFab* Ex = m_fields.get(FieldType::Efield_fp, Direction{0}, lev);
amrex::MultiFab* Ey = m_fields.get(FieldType::Efield_fp, Direction{1}, lev);
amrex::MultiFab* Ez = m_fields.get(FieldType::Efield_fp, Direction{2}, lev);
amrex::MultiFab* rho = m_fields.get(FieldType::rho_fp, lev);  // scalar
```

To get a vector field as a whole (returns `ablastr::fields::VectorField`, an alias for `std::array<MultiFab*, 3>`):

```cpp
auto J = m_fields.get_alldirs(FieldType::current_fp, lev);
amrex::MultiFab* Jx = J[0];
```

To get all levels:

```cpp
auto E_multilevel = m_fields.get_mr_levels_alldirs(FieldType::Efield_fp, finest_level);
```

### Adding a new field type

If your work introduces a new field (say, a new auxiliary multifab for a novel solver):

1. Add a new enumerator to `warpx::fields::FieldType` in `ablastr/fields/MultiFabRegister.H` (or the WarpX-side enum file, depending on the version's exact organization).
2. Document it in the enum's Doxygen comment block.
3. Register it during initialization with `m_fields.alloc_init(FieldType::YourNewField, ...)`. Look at how existing fields are allocated in `WarpX::AllocLevelMFs` (in `Source/Initialization/`) for the exact pattern in the current version.
4. Initialize it via `WarpX::InitLevelData` or a dedicated initialization function if it has nontrivial initial values.
5. Update guard-cell exchange routines if the field needs ghost cells synchronized across MPI ranks.

Don't add raw `std::unique_ptr<amrex::MultiFab>` members directly to `WarpX` for new fields. The register exists to keep field management uniform.

## Field MultiFab iteration pattern

The canonical pattern for any per-cell field operation:

```cpp
for (int lev = 0; lev <= finest_level; ++lev) {
    amrex::MultiFab* Ex = m_fields.get(FieldType::Efield_fp, Direction{0}, lev);

    // MFIter loops over local FABs (or tiles on CPU, full FABs on GPU)
    for (amrex::MFIter mfi(*Ex, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get the tile box. Use tilebox(nodal_flag) so the tile is correctly
        // sized for the field's centering (cell-centered, nodal, edge, face).
        const amrex::Box& tex = mfi.tilebox(Ex->ixType().toIntVect());

        // Extract a portable view to the FAB data
        amrex::Array4<amrex::Real> const& Ex_arr = Ex->array(mfi);
        amrex::Array4<amrex::Real const> const& Bz_arr = /* ... */;

        // Launch the kernel
        amrex::ParallelFor(tex,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Indexing: Ex_arr(i, j, k) for scalar, Ex_arr(i, j, k, comp) for multi-component
                Ex_arr(i, j, k) += /* ... */;
            });
    }
}
```

Key elements:

- **`amrex::TilingIfNotGPU()`** — enables tiling for cache efficiency on CPU, disables it on GPU (where each FAB is one launch). Use this for hot-path loops; for one-off operations, plain `MFIter(*mf)` is also acceptable.
- **`mfi.tilebox(nodal_flag)`** — the index range to iterate. Pass the nodal type of the MultiFab so the box is sized correctly. For cell-centered MultiFabs you can use `mfi.tilebox()` without arguments.
- **`Array4<T>`** — a thin view giving GPU-friendly access. `Array4<Real const>` for read-only fields. The view is captured by value in the lambda (cheap).
- **`amrex::ParallelFor`** — translates to a CUDA/HIP/SYCL kernel launch on GPU, OpenMP loop or `#pragma omp simd` on CPU. The `[=] AMREX_GPU_DEVICE` capture-by-value lambda is required.
- **No member access in the lambda.** Local variables only. See [GPU portability](#gpu-portability-rules-that-arent-optional) below.

For 2D, the lambda has signature `(int i, int j)`. For 1D, `(int i)`. There's also a 1D-flat form `ParallelFor(N, [=] AMREX_GPU_DEVICE (int n) { ... })` for index-based parallelism.

### Reductions

For per-cell reductions (sum, max, min) across a MultiFab:

```cpp
using ReduceOps = amrex::ReduceOps<amrex::ReduceOpSum>;
using ReduceData = amrex::ReduceData<amrex::Real>;
ReduceOps reduce_op;
ReduceData reduce_data(reduce_op);
using ReduceTuple = typename ReduceData::Type;

for (amrex::MFIter mfi(*Ex, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& tex = mfi.tilebox();
    amrex::Array4<amrex::Real const> const& Ex_arr = Ex->const_array(mfi);

    reduce_op.eval(tex, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return { Ex_arr(i, j, k) * Ex_arr(i, j, k) };
        });
}

ReduceTuple host_tuple = reduce_data.value(reduce_op);
amrex::Real result = amrex::get<0>(host_tuple);
amrex::ParallelDescriptor::ReduceRealSum(result);  // MPI reduce
```

Don't accumulate into a host-side variable inside a `ParallelFor` lambda — that's a data race on GPU and slow on CPU.

## Field solver implementation pattern

Field solvers in WarpX follow a templated pattern that lets the same outer loop code dispatch to multiple stencil algorithms. Take an FDTD electric-field update:

```cpp
template <typename T_Algo>
void EvolveECartesian (
    amrex::MultiFab* Ex, amrex::MultiFab* Ey, amrex::MultiFab* Ez,
    amrex::MultiFab const* Bx, amrex::MultiFab const* By, amrex::MultiFab const* Bz,
    amrex::MultiFab const* Jx, amrex::MultiFab const* Jy, amrex::MultiFab const* Jz,
    amrex::Real const dt)
{
    using namespace amrex::literals;

    for (amrex::MFIter mfi(*Ex, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& tex = mfi.tilebox(Ex->ixType().toIntVect());

        amrex::Array4<amrex::Real> const& Ex_arr = Ex->array(mfi);
        amrex::Array4<amrex::Real const> const& By_arr = By->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Bz_arr = Bz->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jx_arr = Jx->const_array(mfi);

        // Stencil coefficients (precomputed, captured by value)
        amrex::Real const c2dt = PhysConst::c2 * dt;
        amrex::Real const mu0c2dt = PhysConst::mu0 * c2dt;

        amrex::ParallelFor(tex,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Ex_arr(i, j, k) += c2dt * (
                    - T_Algo::DownwardDz(By_arr, /* coefs */, i, j, k)
                    + T_Algo::DownwardDy(Bz_arr, /* coefs */, i, j, k)
                    - mu0c2dt * Jx_arr(i, j, k)
                );
            });
        // Similar ParallelFor blocks for Ey, Ez
    }
}
```

The `T_Algo` template parameter is the stencil algorithm class (e.g., `CartesianYeeAlgorithm`, `CartesianCKCAlgorithm`). It provides static functions like `DownwardDz` that implement the discretized derivative. Adding a new field-solver algorithm typically means adding a new algorithm class with these static functions; the outer driver code stays the same.

See `Source/FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/` for existing algorithm classes.

### Communication after a field update

After modifying a field, you typically need to refresh its ghost cells:

```cpp
// Refresh ghosts for one MultiFab
Ex->FillBoundary(Geom(lev).periodicity());

// Or for a vector field, helper from WarpX:
WarpX::FillBoundaryE(lev, /* nghost */);
```

For current density, neighbors *add* (rather than overwrite) into ghost cells, since particles near a Box boundary deposit into the Box's ghost zones that overlap neighbors' valid zones:

```cpp
WarpX::SumBoundaryJ(/* ... */);
```

Communication functions for E/B/J/ρ are collected in `Source/Parallelization/`.

## Particles: container hierarchy and storage layout

### Hierarchy
```
amrex::ParticleContainer_impl<SoAParticle<NReal, NInt>, NReal, NInt, ...>
    └── WarpXParticleContainer (base, polymorphic)
        ├── PhysicalParticleContainer (plasma, beam, neutral species)
        │   └── PhotonParticleContainer (photons, QED)
        └── LaserParticleContainer (antenna particles)
```

`MultiParticleContainer` owns `Vector<unique_ptr<WarpXParticleContainer>> allcontainers` — one per species. Accessed from `WarpX` as `m_mypc` (some older code paths still call it `mypc`).

### Per-particle storage

Each particle is stored in SoA layout:
- **Real attributes** (compile-time + runtime): position (`x`, `y`, `z` in 3D), weight (`w`), momentum (`ux`, `uy`, `uz` — the *normalized* momentum γ·v in WarpX), plus runtime-added components like `prev_x`, `prev_y`, `prev_z`, `orig_x`, etc.
- **Integer attributes**: ionization level (when ionization is enabled), and a 64-bit unique ID (40-bit local index + 24-bit cpu sub-ID).

The compile-time real component indices are in `PIdx` (in `Source/Particles/Pusher/UpdateMomentumBoris.H` or related — exact location varies by version). Typical members: `PIdx::w`, `PIdx::ux`, `PIdx::uy`, `PIdx::uz`. Positions are accessed separately as `pti.GetPosition()` or via SoA real arrays.

A particle's status (active/marked-for-removal) is encoded in the high bit of the 64-bit ID. Use AMReX's helper macros (`AMREX_GPU_DEVICE bool is_valid(...)`) rather than checking the bit directly.

### Adding a new species

The conventional way: a new species is *declared* in the input file with `particles.species_names = ... new_species_name`, and at parse time `MultiParticleContainer::ReadParameters` instantiates an appropriate container. For most new physics, you don't add a new C++ container class — you configure the existing `PhysicalParticleContainer` with new parameters (mass, charge, initialization profile, collision processes).

Only add a new container class if the species has fundamentally different dynamics (no current deposition, special push, etc.). For example, `PhotonParticleContainer` overrides `PushPX` to use a photon-specific push (no Lorentz force).

## Particle iteration pattern

The canonical pattern:

```cpp
// 'pc' is a WarpXParticleContainer* or *this inside a member function
for (int lev = 0; lev <= finest_level; ++lev) {
    for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti) {
        // np = number of particles in this tile
        const long np = pti.numParticles();

        // Get SoA real attributes
        auto& attribs = pti.GetAttribs();
        amrex::ParticleReal* AMREX_RESTRICT w  = attribs[PIdx::w].dataPtr();
        amrex::ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
        amrex::ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
        amrex::ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

        // Get positions (helper, handles dimensionality)
        const auto GetPosition = GetParticlePosition(pti);
        const auto SetPosition = SetParticlePosition(pti);

        // Kernel
        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (long ip)
            {
                amrex::ParticleReal xp, yp, zp;
                GetPosition(ip, xp, yp, zp);

                // ... compute new ux, uy, uz, xp, yp, zp ...

                SetPosition(ip, xp, yp, zp);
            });
    }
}
```

Key elements:

- **`WarpXParIter`** — derives from `amrex::ParIter` and adds WarpX-specific conveniences. Each iteration corresponds to one tile (one Box on GPU).
- **`pti.GetAttribs()`** — returns a structure indexed by `PIdx`. Each entry has a `.dataPtr()` giving a raw `ParticleReal*` for GPU access.
- **`AMREX_RESTRICT`** — equivalent to `__restrict__`; promises non-aliasing to the compiler. Use on all pointers in hot loops.
- **`GetParticlePosition(pti)` / `SetParticlePosition(pti)`** — helper functors that abstract position storage across dimensions. Use these instead of accessing position arrays directly.
- **`ParallelFor(np, ...)`** — 1D index-based parallel loop over particles in the tile.
- **One particle per thread.** Don't try to share work between threads within a particle. Each particle is independent.

For tile-level information (the Box, the geometry):

```cpp
const amrex::Box& box = pti.validbox();
const auto plo = Geom(lev).ProbLoArray();    // physical lower corner of domain
const auto dxi = Geom(lev).InvCellSizeArray();  // 1/dx, 1/dy, 1/dz
```

### Adding and removing particles

To create particles inside a kernel, use the AMReX `filterCreateTransform` pattern:

```cpp
// This filters which cells produce particles, allocates the new particles
// in destination tiles, and applies a transform to set their attributes
amrex::Long count = filterCreateTransformFromFAB<...>(
    dst_tile, src_box, src_FABs, dst_index,
    /* filter functor */, /* create functor */, /* transform functor */
);
```

See `Source/Particles/ParticleCreation/` for the helpers and examples in QED Schwinger (`MultiParticleContainer::doQEDSchwinger`) and ionization. This pattern preserves GPU portability and handles the variable-output-per-input case (you don't know how many particles will be created until you check the filter).

To mark a particle for removal in a kernel, set the high bit of its ID. AMReX provides `amrex::ParticleIDWrapper` and the helper `id = amrex::SetParticleIDInvalid(id)` (exact name may vary). After the kernel, call `Redistribute()` on the container to actually purge invalid particles and rebalance.

Never directly remove particles from a vector during a parallel iteration — mark them and clean up afterward.

## Field-gather and current-deposition: linking particles and fields

Because particles and fields share the same per-Box decomposition, you can use the iteration to obtain matching FABs in the same loop:

```cpp
for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti) {
    // Get FABs for the fields on this Box
    const amrex::FArrayBox& Exfab = (*m_fields.get(FieldType::Efield_aux, Direction{0}, lev))[pti];
    const amrex::FArrayBox& Eyfab = (*m_fields.get(FieldType::Efield_aux, Direction{1}, lev))[pti];
    // ... and Ez, Bx, By, Bz

    amrex::FArrayBox& Jxfab = (*m_fields.get(FieldType::current_fp, Direction{0}, lev))[pti];
    // ... and Jy, Jz

    // Get particle data
    const long np = pti.numParticles();
    auto& attribs = pti.GetAttribs();
    // ...

    // The actual gather/push/deposit kernels go here. WarpX has these
    // factored into namespaces under Source/Particles/Gather/, /Pusher/, /Deposition/.
}
```

In production code, you don't write the gather and deposit kernels inline — you call into the established helpers in those subdirectories. For field-gather, look at `doGatherShapeN` (templated on shape order and gather algorithm). For current deposition, look at `doDepositionShapeN` (Esirkepov, Villasenor-Buneman, direct).

For mesh refinement, you need to handle buffer zones (particles near MR boundaries should gather/deposit on the coarser level for numerical reasons). The mechanism is in `PhysicalParticleContainer::PartitionParticlesInBuffers`, which uses bit masks (`gather_buffer_masks`, `current_buffer_masks`) stored in `WarpX`.

## Particle attributes: built-in and runtime

Compile-time particle attributes (always present):
- AoS-ish (now SoA in current versions): `position_x`, `position_y`, `position_z`, `cpu`, `id`.
- SoA: `PIdx::w` (weight), `PIdx::ux`, `PIdx::uy`, `PIdx::uz` (momentum).

Runtime attributes (added at runtime, conditionally on features):
- `ionizationLevel` (int) — added when impact or field ionization is enabled.
- `opticalDepthQSR`, `opticalDepthBW` (real) — added with PICSAR QED.
- `prev_x`/`prev_y`/`prev_z`, `orig_x`/`orig_y`/`orig_z` — user-requested in input file.
- User-defined real or int attributes — added via input file `addRealAttributes` / `addIntegerAttributes` with parser functions for initial values.

To add a runtime attribute programmatically:

```cpp
pc->AddRealComp("attr_name");  // or AddIntComp for integer
// ... and update PIdx-like enums or use particle_comps map to access it
auto attr_ptr = pti.GetAttribs(particle_comps["attr_name"]).dataPtr();
```

To use a user-defined attribute in input-file space, prefix any user-customizable component names with care — the convention is that vector components use underscores like `prev_x`, `prev_y`, `prev_z` (not `prev.x`).

## GPU portability: rules that aren't optional

The patterns above are not stylistic preferences — they are the GPU portability layer. Violating them yields code that compiles on CPU and silently fails or runs at 1% speed on GPU. The rules:

1. **All hot loops must use `amrex::ParallelFor`.** No manual `for (i = 0; i < N; ++i)` over particles or cells in performance-critical code.

2. **Lambdas must be `[=] AMREX_GPU_DEVICE`.** Capture-by-value, GPU-device-callable. Capture-by-reference (`[&]`) is silently wrong because the captured references point to host memory.

3. **Capture only local variables.** Member variables (`m_foo`) captured via implicit `this` will copy the entire enclosing object to GPU. The `m_` prefix exists specifically so that this mistake is visible at the call site:
   ```cpp
   // WRONG: m_dt captured via this
   amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
       ux[ip] += m_dt * Ex[ip];  // BAD: implicitly captures `this`
   });

   // RIGHT: capture local copy
   const amrex::Real dt_local = m_dt;
   amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
       ux[ip] += dt_local * Ex[ip];
   });
   ```

4. **No virtual calls on device.** Polymorphism does not work inside a `ParallelFor` lambda. If you need different behavior for different cases, dispatch *outside* the lambda (template, `if constexpr`, function pointer set up on host).

5. **No `std::vector`, no `std::map`, no STL containers on device.** Use AMReX's `amrex::Gpu::DeviceVector<T>` for device-resident dynamic arrays, and pass them in as raw pointers via `.dataPtr()`.

6. **Memory allocation patterns.** `amrex::PinnedArenaAllocator` for host memory that needs to be GPU-accessible. The `Arena` mechanism handles GPU device memory. Direct `new`/`malloc` on host memory is OK but won't be visible to GPU kernels.

7. **No `printf` debugging on device.** `amrex::Print()` on host only. For debug output from inside a GPU kernel, `AMREX_DEVICE_PRINTF` exists but it's slow and noisy.

8. **Random numbers on device.** Use `amrex::ParallelForRNG` or the per-thread `amrex::RandomEngine` to get a thread-local RNG state. Don't call host-side `rand()` from device code.

9. **MultiFab and ParticleContainer GPU placement.** These manage their own device storage; the data you access via `Array4<T>` or pointer is already on the device when running with GPU build. You don't need to do explicit transfers.

10. **Test on actual GPU.** Build with `WarpX_COMPUTE=CUDA` (or HIP/SYCL) and run a small simulation that exercises your new code path. CI runs GPU builds; uncommitted local code can break in non-obvious ways.

The single most common GPU mistake in new WarpX code is implicit `this` capture. If a CPU build works but the GPU build produces wrong numbers or crashes, suspect a member variable captured into a `ParallelFor` lambda first.

## A small concrete example: a Coulomb-friction term on a species

For illustration, here's what a small new particle-side operator might look like — a friction force `F = -ν * v` applied to a species, with `ν` a per-species coefficient:

```cpp
void PhysicalParticleContainer::ApplyFriction (amrex::Real dt)
{
    using namespace amrex::literals;

    const amrex::Real nu = m_friction_coefficient;  // local copy of member

    for (int lev = 0; lev <= finestLevel(); ++lev) {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            const long np = pti.numParticles();
            auto& attribs = pti.GetAttribs();
            amrex::ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            amrex::ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            amrex::ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

            const amrex::Real damping = std::exp(-nu * dt);

            amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE (long ip)
                {
                    ux[ip] *= damping;
                    uy[ip] *= damping;
                    uz[ip] *= damping;
                });
        }
    }
}
```

Note: local copy of `m_friction_coefficient`, restrict-qualified pointers, capture-by-value lambda, `_rt` not needed here because we're using `Real` directly and `damping` is computed on host, no member access in the kernel. This is the shape every new operator should have.
