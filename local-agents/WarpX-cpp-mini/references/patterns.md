# Code patterns

Copy-paste recipes for the most common operations. Replace placeholders, don't redesign the structure.

## 1. FP literals

```cpp
// BAD
amrex::Real x = 0.0;
amrex::Real y = 0.5f;
amrex::ParticleReal w = 1.0;

// GOOD
using namespace amrex::literals;     // inside function scope only
amrex::Real x = 0.0_rt;
amrex::Real y = 0.5_rt;
amrex::ParticleReal w = 1.0_prt;
```

## 2. Includes ordering

```cpp
// 1. Self header (.cpp only)
#include "MyClass.H"

// 2. WarpX internal headers (alphabetical)
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"

// 3. ablastr headers
#include <ablastr/fields/MultiFabRegister.H>

// 4. AMReX headers
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

// 5. Other libraries
#include <openPMD/openPMD.hpp>

// 6. C++ standard library
#include <memory>
#include <string>
```

Header files use `#ifndef WARPX_<PATH>_<NAME>_H_` include guards (not `#pragma once`).

## 3. Loop over field MultiFab + ParallelFor

```cpp
#include <AMReX_MultiFab.H>

void ApplyToField(amrex::MultiFab& mf, amrex::Real factor) {
    using namespace amrex::literals;
    const amrex::Real f = factor;                              // local copy

    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();                  // valid region (no ghosts)
        amrex::Array4<amrex::Real> const& arr = mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i, j, k) *= f;
        });
    }
}
```

For multi-component MultiFab:
```cpp
amrex::ParallelFor(bx, mf.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
    arr(i, j, k, n) *= f;
});
```

## 4. Loop over particles + ParallelFor

```cpp
#include <AMReX_ParticleContainer.H>
#include "Particles/WarpXParticleContainer.H"

void DampParticles(WarpXParticleContainer& pc, amrex::Real damp_factor) {
    using namespace amrex::literals;
    const amrex::ParticleReal damp = damp_factor;

    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
            auto& ptile = pti.GetParticleTile();
            auto& soa = ptile.GetStructOfArrays();
            const auto np = pti.numParticles();

            amrex::ParticleReal* AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            amrex::ParticleReal* AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            amrex::ParticleReal* AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) noexcept {
                ux[ip] *= damp;
                uy[ip] *= damp;
                uz[ip] *= damp;
            });
        }
    }
}
```

`PIdx` enum: standard particle indices (`w`, `ux`, `uy`, `uz`, etc.). Check `Source/Particles/WarpXParticleContainer.H` for the current list — it shifts.

## 5. Lambda captures — never `this` in device code

```cpp
// BAD — copies the entire WarpX object to GPU
class MyClass {
    amrex::Real m_dt;
    void Update(amrex::MultiFab& mf) {
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto arr = mf.array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=, this] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                arr(i, j, k) += m_dt;     // captures via `this`
            });
        }
    }
};

// GOOD — local copy before the lambda
class MyClass {
    amrex::Real m_dt;
    void Update(amrex::MultiFab& mf) {
        const amrex::Real dt = m_dt;      // copy to local
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto arr = mf.array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                arr(i, j, k) += dt;
            });
        }
    }
};
```

The `m_` prefix on members exists to make capture-by-`this` visually obvious. If you see `m_foo` inside a `ParallelFor` lambda, fix it.

## 6. Get a field from MultiFabRegister

```cpp
#include <ablastr/fields/MultiFabRegister.H>

void WarpX::DoSomething(int lev) {
    using ablastr::fields::Direction;
    namespace ft = warpx::fields;

    // Scalar fields (rho, phi, F, G):
    amrex::MultiFab* rho = m_fields.get(ft::FieldType::rho_fp, lev);

    // One component of a vector field:
    amrex::MultiFab* Ex = m_fields.get(ft::FieldType::Efield_fp, Direction{0}, lev);

    // All three components:
    auto E = m_fields.get_alldirs(ft::FieldType::Efield_fp, lev);
    // E[0], E[1], E[2] are MultiFab* for x, y, z

    // For field-gather, use _aux NOT _fp:
    amrex::MultiFab* Ex_aux = m_fields.get(ft::FieldType::Efield_aux, Direction{0}, lev);
}
```

Common `FieldType` values: `Efield_fp`, `Bfield_fp`, `Efield_aux`, `Bfield_aux`, `Efield_cp`, `Bfield_cp`, `current_fp`, `rho_fp`, `phi_fp`, `F_fp`, `G_fp`, hybrid-Ohm fields, EB fields.

**Never** add `std::unique_ptr<amrex::MultiFab>` as a new member of `WarpX`. Use the register.

## 7. Add a runtime particle attribute

```cpp
// In setup:
pc.AddRealComp("optical_depth_QSR");
pc.AddIntComp("ionization_level");

// Access later:
const int comp = pc.GetRealCompIndex("optical_depth_QSR");
amrex::ParticleReal* data = soa.GetRealData(comp).data();
```

## 8. Add a new source file

When adding `Source/Foo/Bar.cpp`:

**`Source/Foo/CMakeLists.txt`** — add `Bar.cpp`:
```cmake
target_sources(WarpX
    PRIVATE
        Bar.cpp
)
```

**`Source/Foo/Make.package`** — add `Bar.cpp`:
```
CEXE_sources += Bar.cpp
```

Headers (`.H`) are NOT listed in either file.

## 9. Mark a function device-callable

```cpp
AMREX_GPU_HOST_DEVICE
amrex::Real my_kernel_helper(amrex::Real x) {
    return x * x;
}
```

Without this annotation, calling from inside a `ParallelFor` lambda silently fails or won't compile on GPU.

## 10. GPU rules quick reference

| Don't | Do |
|---|---|
| `[&]` captures in device lambdas | `[=]` captures, copy values |
| `new` / `malloc` / STL containers in device code | Pre-allocate, use `amrex::Vector` on host |
| Recursion in device code | Iterative algorithms |
| `std::cout`, exceptions in device | `amrex::Print`, `AMREX_ASSERT_WITH_MESSAGE` |
| `std::rand`, `<random>` in device | `amrex::Random()`, `amrex::RandomNormal()` |
| Raw `+=` for shared accumulators | `amrex::Gpu::Atomic::Add()` |
| Calling unmarked functions from device | Mark with `AMREX_GPU_HOST_DEVICE` |

## 11. Use `WarpX&` instead of `GetInstance()`

```cpp
// BAD — legacy singleton pattern
void MyFunction() {
    auto& warpx = WarpX::GetInstance();
    auto* Ex = warpx.m_fields.get(...);
}

// GOOD — pass references
void MyFunction(WarpX& warpx, int lev) {
    auto* Ex = warpx.m_fields.get(...);
}

// EVEN BETTER — pass the specific container
void MyFunction(ablastr::fields::MultiFabRegister& fields, int lev) {
    auto* Ex = fields.get(...);
}
```

## 12. Refresh ghost cells after writes

```cpp
// After modifying a MultiFab, refresh ghost cells:
mf.FillBoundary(geom.periodicity());
```

Skipping this is a common bug — subsequent reads from ghost cells return stale data.

## 13. Redistribute particles after position changes

```cpp
// After modifying particle positions in any callback or routine:
pc.Redistribute();
```

Otherwise particles end up in the wrong box owner, and the next step's gather/deposit are wrong.

## 14. Namespace discipline

```cpp
// BAD
using namespace amrex;                  // never at file scope
namespace ax = amrex;                   // also avoid

// GOOD
amrex::Real x;
amrex::MultiFab* mf;
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { /* ... */ });

// Narrow function scope OK for literals:
void f() {
    using namespace amrex::literals;
    amrex::Real x = 0.5_rt;
}
```

## 15. Dimensionality guards

```cpp
#if defined(WARPX_DIM_3D)
    // 3D-specific code
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_2D)
    // 2D-specific
#elif defined(WARPX_DIM_RZ)
    // RZ-specific
#elif defined(WARPX_DIM_1D_Z) || defined(WARPX_DIM_1D)
    // 1D-specific
#endif
```

WarpX builds one binary per dimension. The macros are compile-time. `AMREX_SPACEDIM` also gives the dimension number.

## 16. Tile box vs FAB box

```cpp
// Use tilebox() — the valid region of this tile (no ghost cells)
const amrex::Box& bx = mfi.tilebox();

// Use fabbox() ONLY if you need ghost cells too — usually not
const amrex::Box& bx_with_ghosts = mfi.fabbox();

// For staggered fields, use the nodal variant:
const amrex::Box& bx_node = mfi.tilebox(amrex::IntVect::TheNodeVector());
const amrex::Box& bx_face = mfi.nodaltilebox(direction);
```

## Where new physics goes

| Adding | Location |
|---|---|
| Particle pusher | `Source/Particles/Pusher/` + enum in `WarpXAlgorithmSelection.H` |
| Field solver | `Source/FieldSolver/<family>/` + `WarpX::electromagnetic_solver_id` |
| Collision/scattering | `Source/Particles/Collision/...` (see `collisions.md`) |
| Particle init profile | `Source/Initialization/` + `PlasmaInjector` |
| Field BC | `Source/BoundaryConditions/` |
| Particle BC | `Source/Particles/ParticleBoundaries/` |
| Diagnostic | `Source/Diagnostics/` (reduced: `ReducedDiags/`) |
| Filter | `Source/Filter/` |
| Runtime particle attribute | `Source/Particles/` + `MultiParticleContainer::ReadParameters` |

For each: find the closest existing example, copy and adapt. Don't design new modules from scratch.

## Verify everything against Doxygen

`https://warpx.readthedocs.io/en/latest/_static/doxyhtml/` is the ground truth for class names, method signatures, and enum values. Always check there before relying on remembered APIs.
