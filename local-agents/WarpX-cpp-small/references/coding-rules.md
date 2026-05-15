# Coding rules with examples

These rules are mandatory. Code that breaks them won't pass review, may not build, or will fail on GPU. Examples are the fastest way to internalize the patterns — copy them.

## Table of contents

1. [FP literals: `_rt` and `_prt`](#fp-literals-_rt-and-_prt)
2. [Namespace discipline](#namespace-discipline)
3. [Includes ordering](#includes-ordering)
4. [Member variables: the `m_` prefix](#member-variables-the-m_-prefix)
5. [Lambda captures for ParallelFor](#lambda-captures-for-parallelfor)
6. [The MFIter + ParallelFor pattern (fields)](#the-mfiter--parallelfor-pattern-fields)
7. [The ParIter + ParallelFor pattern (particles)](#the-pariter--parallelfor-pattern-particles)
8. [Accessing fields through MultiFabRegister](#accessing-fields-through-multifabregister)
9. [Particle attributes: real and int components](#particle-attributes-real-and-int-components)
10. [Adding a new source file](#adding-a-new-source-file)
11. [GPU portability constraints](#gpu-portability-constraints)
12. [Avoid `WarpX::GetInstance()`](#avoid-warpxgetinstance)

## FP literals: `_rt` and `_prt`

Bare floating-point literals defeat WarpX's build-time precision configuration. Always use suffixes.

```cpp
// BAD
amrex::Real x = 0.0;           // wrong type if build is single-precision
amrex::Real y = 0.5f;          // explicit single — same problem
amrex::Real pi = 3.14159265;

// GOOD
using namespace amrex::literals;  // inside a function or narrow scope only
amrex::Real x = 0.0_rt;
amrex::Real y = 0.5_rt;
amrex::Real pi = 3.14159265_rt;

// Particle-precision values use _prt:
amrex::ParticleReal weight = 1.0_prt;
amrex::ParticleReal vx = 1.0e6_prt;
```

`_rt` = `amrex::Real` (field precision). `_prt` = `amrex::ParticleReal` (particle precision). These are independently configurable; mixing them invites silent precision loss.

## Namespace discipline

```cpp
// BAD
using namespace amrex;           // NEVER at file scope
namespace ax = amrex;            // aliases at file scope are also discouraged

// GOOD
amrex::Real x;
amrex::MultiFab* mf;
amrex::ParallelFor(box, [=] (int i, int j, int k) noexcept { /* ... */ });

// Narrow scope OK for literals only:
void my_function() {
    using namespace amrex::literals;
    amrex::Real x = 0.5_rt;
}
```

Always fully qualify `amrex::`. The only exception is `using namespace amrex::literals;` inside a function scope to enable `_rt`/`_prt` suffixes.

## Includes ordering

Includes in `.cpp` and `.H` files follow strict order:

```cpp
// 1. Self header (in .cpp only)
#include "MyClass.H"

// 2. WarpX internal headers (alphabetical within blocks)
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"

// 3. ablastr headers
#include <ablastr/fields/MultiFabRegister.H>

// 4. AMReX headers
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

// 5. Other library headers (e.g., openPMD, hdf5)
#include <openPMD/openPMD.hpp>

// 6. C++ standard library
#include <memory>
#include <string>
#include <vector>
```

Within each group: alphabetical. Forward-declaration headers (`_fwd.H`) are used widely to reduce compile time — when writing a new class that other headers reference by pointer, also write a `MyClass_fwd.H`.

Headers (`.H` files): use `#ifndef WARPX_<PATH>_<NAME>_H_` / `#define ... / #endif` include guards (not `#pragma once`).

## Member variables: the `m_` prefix

Member variables of WarpX classes are prefixed with `m_`. This is a deliberate visual marker — when you see `m_` you know capture-by-`this` would drag the whole object.

```cpp
class WarpX {
private:
    int m_n_step;
    amrex::Real m_dt;
    ablastr::fields::MultiFabRegister m_fields;
    std::unique_ptr<MultiParticleContainer> m_mypc;
};
```

You'll see legacy code without `m_` (e.g., `mypc` is still common). New code should follow the convention.

## Lambda captures for ParallelFor

This is the single most common source of GPU bugs. When passing a lambda to `amrex::ParallelFor`, the lambda runs on the device (GPU). Capturing `this` would copy the entire `WarpX` object to GPU. **Always copy needed members to local variables first.**

```cpp
// BAD — captures `this`; on GPU this copies the WarpX object
class MyClass {
    amrex::Real m_dt;
    void Update(amrex::MultiFab& mf) {
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto arr = mf.array(mfi);
            amrex::ParallelFor(bx, [=, this] (int i, int j, int k) noexcept {
                arr(i, j, k) += m_dt;     // m_dt accessed via captured `this`
            });
        }
    }
};

// GOOD — copy to local; lambda captures by value
class MyClass {
    amrex::Real m_dt;
    void Update(amrex::MultiFab& mf) {
        const amrex::Real dt = m_dt;       // local copy
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto arr = mf.array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) noexcept {
                arr(i, j, k) += dt;        // captures `dt` by value
            });
        }
    }
};
```

The `m_` prefix exists precisely to make this mistake visually obvious. If you see `m_foo` inside a ParallelFor lambda, you have a bug.

## The MFIter + ParallelFor pattern (fields)

The canonical loop over field data:

```cpp
#include <AMReX_MultiFab.H>

void ApplyToField(amrex::MultiFab& mf, amrex::Real factor) {
    using namespace amrex::literals;
    const amrex::Real f = factor;

    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // tilebox() returns the valid region of this tile (no ghost cells)
        const amrex::Box& bx = mfi.tilebox();
        // array() gives an Array4 view: arr(i, j, k, comp)
        amrex::Array4<amrex::Real> const& arr = mf.array(mfi);

        amrex::ParallelFor(bx, [=] (int i, int j, int k) noexcept {
            arr(i, j, k) *= f;
        });
    }
}
```

Key points:
- `amrex::TilingIfNotGPU()` enables CPU tiling (cache efficiency) but turns off on GPU (already parallel).
- `mfi.tilebox()` for the valid region. Use `mfi.fabbox()` only if you need ghost cells too — and you usually don't.
- `mf.array(mfi)` returns an `Array4<Real>` — index as `arr(i, j, k, comp)`.
- `ParallelFor` is non-blocking on GPU; AMReX handles synchronization at end-of-scope.

For multi-component MultiFabs:
```cpp
amrex::ParallelFor(bx, mf.nComp(), [=] (int i, int j, int k, int n) noexcept {
    arr(i, j, k, n) *= f;
});
```

For face-/edge-/node-centered fields: use the corresponding nodal-aware tile box (`mfi.nodaltilebox(dir)`, `mfi.tilebox(IntVect::TheNodeVector())`).

## The ParIter + ParallelFor pattern (particles)

Looping over particles:

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

            // Access SoA components by index (RealAoS::w, ux, uy, uz, etc.)
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

Key points:
- `WarpXParIter` walks particle tiles within a level.
- `pti.GetParticleTile()` gives the tile's data containers.
- `soa.GetRealData(PIdx::xxx).data()` returns a `ParticleReal*` to that attribute array.
- `PIdx` is an enum of standard particle indices: `w` (weight), `ux`, `uy`, `uz`, etc. Check `Source/Particles/WarpXParticleContainer.H` for the current list.
- `AMREX_RESTRICT` enables vectorization (equivalent to `__restrict__`).
- `AMREX_GPU_DEVICE` annotation is required for device-callable lambdas.

For runtime attributes (added after compile time):
```cpp
const int comp = pc.GetRealCompIndex("my_attribute");
amrex::ParticleReal* data = soa.GetRealData(comp).data();
```

## Accessing fields through MultiFabRegister

```cpp
#include <ablastr/fields/MultiFabRegister.H>

// In a WarpX member function:
void WarpX::DoSomething(int lev) {
    using ablastr::fields::Direction;
    // Get scalar field (rho, phi, F, G):
    amrex::MultiFab* rho = m_fields.get(warpx::fields::FieldType::rho_fp, lev);

    // Get one component of a vector field:
    amrex::MultiFab* Ex = m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{0}, lev);
    amrex::MultiFab* Ey = m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{1}, lev);
    amrex::MultiFab* Ez = m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{2}, lev);

    // Get all three components as a VectorField:
    auto E_vec = m_fields.get_alldirs(warpx::fields::FieldType::Efield_fp, lev);
    // E_vec[0], E_vec[1], E_vec[2] are MultiFab* for x, y, z

    // Particles gather from _aux, not _fp:
    amrex::MultiFab* Ex_aux = m_fields.get(warpx::fields::FieldType::Efield_aux, Direction{0}, lev);
}
```

NEVER add a `std::unique_ptr<amrex::MultiFab>` to `WarpX` as a new member. Use the register.

To register a new field type:
1. Add a value to `enum class FieldType` (in `Source/ablastr/fields/MultiFabRegister.H` or the WarpX-specific extension)
2. Call `m_fields.alloc_init(...)` during `WarpX::AllocLevelMFs` for that level
3. Access through `m_fields.get(...)`

## Particle attributes: real and int components

Compile-time attributes are listed in the `PIdx` enum (Source/Particles/WarpXParticleContainer.H). The standard ones for 3D include position internals, weight, ux/uy/uz, plus dimension-dependent extras. Look up the current enum — it shifts.

Runtime attributes (added per-species at runtime, e.g., ionization level, optical depth):

```cpp
// In a species setup function:
pc.AddRealComp("optical_depth_QSR");      // real attribute
pc.AddIntComp("ionization_level");         // int attribute

// Later, access by name:
const int comp = pc.GetRealCompIndex("optical_depth_QSR");
amrex::ParticleReal* od = soa.GetRealData(comp).data();
```

Runtime attribute setup happens in `MultiParticleContainer::ReadParameters` based on what physics is active for each species.

## Adding a new source file

When adding `Source/Foo/Bar.cpp`:

1. Add to **`Source/Foo/CMakeLists.txt`**:
   ```cmake
   target_sources(WarpX
       PRIVATE
           Bar.cpp
   )
   ```

2. Add to **`Source/Foo/Make.package`**:
   ```
   CEXE_sources += Bar.cpp
   ```

Headers (`.H` files) are NOT listed in either file. They're picked up automatically.

If you forget to update one, the file will be missing from one of the two build systems (CMake or GNU Make) and CI will fail.

## GPU portability constraints

Code that compiles fine on CPU can silently fail or produce wrong results on GPU:

- **Captures by reference `[&]`** inside `ParallelFor` are illegal on GPU. Always `[=]` and copy needed values.
- **Allocation** (`new`, `malloc`, STL containers) inside device kernels is forbidden. Use pre-allocated arrays.
- **Recursion** inside device code is forbidden (or constrained).
- **`std::cout`, `std::string`, exceptions** are not device-callable. Use `amrex::Print`, `amrex::PrintToFile`, `AMREX_ASSERT_WITH_MESSAGE`.
- **Function calls inside device kernels** must be `AMREX_GPU_HOST_DEVICE` (or `_DEVICE`). Calling a normal CPU function from a device lambda silently fails or refuses to compile.
- **Random numbers**: use `amrex::Random()` or `amrex::RandomNormal()` inside kernels. Don't use STL `<random>`.
- **Atomic operations**: use `amrex::Gpu::Atomic::Add()` etc. for atomics.

To mark a function device-callable:
```cpp
AMREX_GPU_HOST_DEVICE
amrex::Real my_kernel_helper(amrex::Real x) {
    return x * x;
}
```

Always test new kernels with `WarpX_COMPUTE=CUDA` (or HIP/SYCL) before assuming portability. Use `make CUDA=TRUE` for fast compile/test cycles on a GPU-equipped machine.

## Avoid `WarpX::GetInstance()`

Legacy code uses the singleton pattern:
```cpp
auto& warpx = WarpX::GetInstance();
auto* Ex = warpx.m_fields.get(...);
```

New code should accept `WarpX&` or the specific container by reference:
```cpp
// Better:
void MyNewFunction(WarpX& warpx, int lev) {
    auto* Ex = warpx.m_fields.get(warpx::fields::FieldType::Efield_fp, Direction{0}, lev);
}

// Or even better — pass the register directly:
void MyNewFunction(ablastr::fields::MultiFabRegister& fields, int lev) {
    auto* Ex = fields.get(warpx::fields::FieldType::Efield_fp, Direction{0}, lev);
}
```

This is being incrementally rolled out. New PRs that introduce more `GetInstance()` calls are likely to be asked to refactor.
