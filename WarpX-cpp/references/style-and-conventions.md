# Style and conventions

This reference covers WarpX's coding style — what the codebase looks like, what gets accepted in review, and what the conventions encode about the underlying architecture. Read this before writing new code or modifying existing code substantially.

WarpX style is mostly inherited from AMReX, with WarpX-specific tightening. Some conventions look pedantic but encode real constraints (GPU portability, build-time precision configuration, compile time). Don't relax them.

## Table of contents

1. [Formatting basics](#formatting-basics)
2. [Naming conventions](#naming-conventions)
3. [Member variables and the `m_` prefix](#member-variables-and-the-m_-prefix)
4. [Floating-point literals: `_rt` and `_prt`](#floating-point-literals-_rt-and-_prt)
5. [The `amrex::` namespace](#the-amrex-namespace)
6. [Include directives: order matters](#include-directives-order-matters)
7. [Forward-declaration headers (`*_fwd.H`)](#forward-declaration-headers-_fwdh)
8. [Function definitions vs calls](#function-definitions-vs-calls)
9. [Curly-brace conventions](#curly-brace-conventions)
10. [Templates and constexpr](#templates-and-constexpr)
11. [Error handling and assertions](#error-handling-and-assertions)
12. [GPU-portability rules in code](#gpu-portability-rules-in-code)
13. [Common review feedback](#common-review-feedback)

## Formatting basics

- **Indentation: four spaces, no tabs.** Configure your editor.
- **Line length: under 100 characters** for code. Use line breaks for long argument lists.
- **Documentation files** (`.rst`, `.md`): one sentence per line, ignore line length. This makes documentation diffs readable.
- **Trailing whitespace: remove.** Enable trailing-whitespace removal on save.
- **Final newline at end of file.** Standard.
- **Space around `=`** in assignments and initializations: `int x = 5;`, not `int x=5;`.

## Naming conventions

WarpX is *mostly* consistent on these, but you'll find old code that violates them. New code should follow:

| Entity | Convention | Example |
|---|---|---|
| File names | `CamelCase.{H,cpp}` | `MultiParticleContainer.cpp` |
| Class names | `CamelCase` | `class PhysicalParticleContainer` |
| Member variables | `m_snake_case` | `m_finest_level` |
| Static member variables | `m_snake_case` (some places use prefix `s_`) | `m_friction_coefficient` |
| Local variables | `snake_case` | `int particle_count;` |
| Function/method names | `MixedCase` (acceptable) or `snake_case` (newer code prefers) | `WarpX::EvolveE()`, `compute_collision_rate()` |
| Constants | `UPPER_SNAKE_CASE` | `MAX_PARTICLES_PER_TILE` |
| Macros | `UPPER_SNAKE_CASE`, prefixed `WARPX_` or `AMREX_` | `WARPX_DIM_3D` |
| Template parameters | `T_CamelCase` or `T_snake_case` | `template <typename T_Algo>` |
| Enums | `enum class CamelCase`; values `CamelCase` | `enum class PushType { Explicit, Implicit }` |

WarpX is migrating from older mixed-case function naming (`EvolveE`, `PushPX`) toward modern lowercase-snake names. Both styles exist in the current codebase; new functions can use either if it matches the surrounding context, but for newly created subsystems prefer `lower_snake_case` for free functions and `MixedCase` for class methods.

## Member variables and the `m_` prefix

**All new class member variables must be prefixed with `m_`.** This is documented in the contributing guide and is the single most-enforced naming rule.

The reason is not style — it's GPU safety. Lambda captures in `ParallelFor` look like:

```cpp
amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
    // ... what gets captured here? ...
});
```

Inside a class method, an unqualified member name implicitly resolves to `this->member`, which forces the lambda to capture `this`. On GPU, this either copies the entire class to device memory (slow) or produces a dangling host pointer (wrong). The `m_` prefix makes member access visible at the call site, so reviewers catch this immediately:

```cpp
// BAD — m_dt is captured via this
amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
    ux[ip] += m_dt * Ex[ip];  // implicit this-> here
});

// GOOD — local variable captured by value
const amrex::Real dt_local = m_dt;
amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
    ux[ip] += dt_local * Ex[ip];
});
```

Old WarpX code without the `m_` prefix exists. When you modify such code, applying the prefix is a welcome cleanup, but do it in a separate PR from any functional changes (per the contributing guide: "style changes are not included in the PR where new code is added").

## Floating-point literals: `_rt` and `_prt`

WarpX is configured at build time for single or double precision *independently* for fields (`amrex::Real`) and particles (`amrex::ParticleReal`). Bare floating-point literals like `0.0` or `3.14f` defeat this — they force a specific type.

Use the user-defined literals:

```cpp
amrex::Real x = 0.0_rt;         // amrex::Real, whatever precision Real is
amrex::Real y = 3.14_rt;
amrex::ParticleReal w = 1.0_prt; // amrex::ParticleReal

// In a function template parameterized on Real type:
T value = T(0.0);                // OK, explicit cast — but _rt is preferred when T == Real
```

To use these suffixes, you need `using namespace amrex::literals;` in the local scope. This is the *only* `using namespace` directive permitted in WarpX code:

```cpp
void EvolveECartesian (/* ... */) {
    using namespace amrex::literals;  // OK, narrow scope
    amrex::Real c2dt = PhysConst::c2 * dt;
    amrex::Real e_val = 0.0_rt;
    // ...
}
```

Integer literals (`0`, `1`) are fine — integers are exact and don't have precision configuration.

### Why this matters

WarpX's single-precision builds are used for fast prototyping and on GPUs where single precision is faster and uses less memory. If a single bare `0.5` literal sneaks in, it's a `double` that forces conversion at every use, defeating the precision configuration and slowing down the build. Worse, in single-precision builds, `3.14` (an `int`-then-implicit-conversion to `double`-then-narrowing to `float`) can introduce subtle precision-loss bugs.

Older code is *not* fully converted — you'll find raw `0.0`s in the codebase. Don't propagate them; use `_rt`/`_prt` in new code.

## The `amrex::` namespace

**Never `using namespace amrex;` at file scope.** Always write `amrex::MultiFab`, `amrex::Real`, `amrex::ParallelFor`. This is enforced.

The reason is that `amrex::` types have well-known names (`Real`, `Vector`, `Box`) that conflict with names in the standard library, third-party libraries, and WarpX itself. Implicit `amrex::` would create ambiguity that's painful to debug.

The *one* exception is `using namespace amrex::literals;`, which is needed for `_rt`/`_prt` suffixes and is unambiguous. Use it inside narrow scopes (functions, short blocks), never at file or namespace scope.

Inside a member function of an AMReX-derived class (e.g., `PhysicalParticleContainer` derives from `amrex::ParticleContainer_impl`), unqualified names of the base-class type members can be implicitly looked up via Argument-Dependent Lookup. Be explicit anyway:

```cpp
// In a method of PhysicalParticleContainer:
amrex::Long np = this->numParticles();   // Explicit, clear
auto& mfi = MFIter(...);                  // Bad — needs amrex::
```

## Include directives: order matters

The order of `#include` directives in WarpX is strict, both because it affects compile time (via forward declarations) and because it documents dependency direction. **The order, top to bottom, in a `.cpp` file:**

1. The corresponding header for this `.cpp` file: `#include "MyName.H"`
2. Other WarpX headers: `#include "Path/To/OtherHeader.H"`
3. WarpX forward-declaration headers: `#include "Path/To/OtherType_fwd.H"`
4. AMReX headers: `#include <AMReX_MultiFab.H>`
5. AMReX forward-declaration headers: `#include <AMReX_BaseFwd.H>`
6. PICSAR headers: `#include <picsar_qed/...>`
7. Other third-party: `#include <openPMD/...>`
8. Standard library: `#include <vector>`

In a `.H` header file, the order is the same except step 1 becomes "the corresponding `_fwd.H`":

1. `#include "MyName_fwd.H"` (if it exists)
2. WarpX headers
3. WarpX `_fwd.H` headers
4. AMReX headers
5. AMReX `Fwd.H` headers
6. PICSAR
7. Third-party
8. Standard library

**Within each group, sort alphabetically.** Insert blank lines between groups. This produces diffs that are easy to review and reduces merge conflicts.

### Quotation marks

- WarpX-local headers: `#include "Path/To/Header.H"` (with double quotes, path relative to `Source/`).
- External headers (AMReX, third-party, stdlib): `#include <Path/To/Header.H>` (angle brackets).

### Header inclusion philosophy

- **Forward-declare when possible.** If a header only needs to know that a type exists (e.g., it holds a pointer or reference to it), use the forward-declaration header. This dramatically reduces compile time across the project.
- **Include only what you use.** Don't include `WarpX.H` just because everyone else does. Include the minimum necessary.
- **In a `.cpp`, don't re-include what its `.H` already includes.** The compiler doesn't care, but it bloats the file and obscures dependencies.

If you don't know which header provides a given symbol, search the Doxygen reference, or `grep -r "class FooBar" Source/`.

## Forward-declaration headers (`*_fwd.H`)

For any non-trivial class declared in `MyClass.H`, there should be a matching `MyClass_fwd.H` containing just the forward declaration:

```cpp
// MyClass_fwd.H
#ifndef MY_CLASS_FWD_H
#define MY_CLASS_FWD_H

class MyClass;

// If MyClass is in a namespace, declare it there
namespace warpx::particles { class CollisionHandler; }

#endif // MY_CLASS_FWD_H
```

The full header `MyClass.H` then includes its own forward declaration:

```cpp
// MyClass.H
#ifndef MY_CLASS_H
#define MY_CLASS_H

#include "MyClass_fwd.H"  // own forward decl first
#include "OtherFullHeader.H"

class MyClass {
    void some_method();
};

#endif // MY_CLASS_H
```

This pattern lets headers that only need to know "`MyClass` is a type that exists" include `MyClass_fwd.H` instead of `MyClass.H` (which transitively pulls in everything `MyClass.H` includes). Compile time drops significantly when forward declarations are used consistently.

**Include guards are required** on `_fwd.H` files, just like full headers. Use `#define MY_CLASS_FWD_H` (the file name uppercased with underscores).

AMReX uses the convention `*Fwd.H` (no underscore, capital F) for its forward declarations: `AMReX_BaseFwd.H`, etc.

## Function definitions vs calls

WarpX uses a `git grep`-friendly convention to distinguish function definitions from calls:

- **Definitions**: space between function name and parenthesis:
  ```cpp
  void MyClass::do_thing ()
  {
      // ...
  }
  ```
- **Calls**: no space:
  ```cpp
  my_object.do_thing();
  ```

This lets `git grep "do_thing ()"` find only definitions, which is enormously useful in a large codebase. Don't put a space at the call site; don't omit the space at the definition site.

Same logic for classes:

- **Definitions**: `class MyClass` on the same line as the name. `git grep "class MyClass"` finds the definition.
- **Forward declarations**: also on one line: `class MyClass;`.

## Curly-brace conventions

Function and class definitions have the opening `{` on a *new* line:

```cpp
void some_function ()
{
    // body
}

class MyClass
{
public:
    void method ();
};
```

For control-flow blocks (`if`, `for`, `while`), the opening `{` is on the *same* line:

```cpp
if (condition) {
    do_thing();
} else {
    do_other_thing();
}

for (int n = 0; n < N; ++n) {
    process(n);
}
```

**Always use curly braces, even for single-statement blocks.** This is enforced:

```cpp
// GOOD
for (int n = 0; n < 10; ++n) {
    amrex::Print() << "Like this!";
}

// Also acceptable for very short bodies
for (int n = 0; n < 10; ++n) { amrex::Print() << "Or like this!"; }

// NOT acceptable
for (int n = 0; n < 10; ++n)
    amrex::Print() << "Not like this.";   // missing braces

// NOT acceptable
for (int n = 0; n < 10; ++n) amrex::Print() << "Nor like this.";  // missing braces
```

The braceless form is forbidden because it leads to bugs when adding statements ("add one more line to the loop body" → silently outside the loop because no braces).

## Templates and constexpr

WarpX uses templates extensively for compile-time algorithm dispatch. Field solvers are templated on stencil algorithm, particle pushers on push type, deposition kernels on shape order. The pattern:

```cpp
template <typename T_Algo>
void EvolveECartesian (/* ... */)
{
    // Use T_Algo::DownwardDz, T_Algo::UpwardDx, etc.
}

// Then dispatch at the call site:
if (algo == "yee") {
    EvolveECartesian<CartesianYeeAlgorithm>(/* ... */);
} else if (algo == "ckc") {
    EvolveECartesian<CartesianCKCAlgorithm>(/* ... */);
}
```

This pattern avoids per-cell virtual calls (which would prevent vectorization and GPU portability) at the cost of more compile work and binary size. Use it when:
- The choice between variants is fixed for the duration of a simulation (selected from input file once).
- The variants differ in a hot kernel.
- The variants share most of the surrounding driver code.

For one-off cases or rarely-called code, a regular `if`/`switch` is fine.

`constexpr` is used where possible for constants:

```cpp
constexpr amrex::Real PI = 3.14159265358979323846_rt;
constexpr int SHAPE_ORDER = 3;
```

For non-constant build-time values (single vs double precision), use the AMReX type machinery rather than `constexpr` literals.

## Error handling and assertions

WarpX uses a mix of AMReX-provided macros:

```cpp
// Always-checked assertion (release and debug)
WARPX_ALWAYS_ASSERT_WITH_MESSAGE(condition, "Helpful error message");

// Abort with a message
amrex::Abort("Something is wrong: " + descriptive_string);

// Print a message to all ranks (use sparingly)
amrex::Print() << "Diagnostic info\n";

// Print a message only from rank 0
amrex::Print() << "Always rank 0\n";
amrex::Print(amrex::Print::AllProcs) << "All ranks: " << my_value << "\n";

// Print to stderr
amrex::Warning("Non-fatal issue");
```

Use assertions liberally in setup code where the cost is one-time; use them sparingly inside per-particle or per-cell kernels. Inside GPU kernels, `AMREX_ASSERT` is the GPU-safe variant; standard `assert` doesn't work on device.

Exceptions are *not* idiomatic in WarpX kernels (they have overhead and don't work on GPUs). Use return codes or `amrex::Abort` for unrecoverable errors. For genuinely exceptional conditions in setup code (file not found, malformed input), throwing is acceptable.

### Input parsing

Parameters are parsed at startup. WarpX provides parser helpers:

```cpp
amrex::ParmParse pp("my_section");
int my_value = 0;
pp.query("my_value", my_value);   // optional, uses default if not present
pp.get("required_value", my_value); // required, errors if not present
```

The WarpX-extended versions (`utils::parser::queryWithParser`) accept arithmetic expressions in input values and should be preferred for parameters that could be expressions:

```cpp
utils::parser::queryWithParser(pp_warpx, "self_fields_required_precision", my_value);
```

## GPU-portability rules in code

Repeating from `fields-and-particles.md` because they belong here too:

### Hot loops use `amrex::ParallelFor`
Not hand-rolled. Not `std::for_each`. Not OpenMP pragmas directly (AMReX handles OpenMP under the hood).

### Lambdas: `[=] AMREX_GPU_DEVICE`
Capture by value, device-callable. Never `[&]` (captures by reference, which doesn't transfer to GPU).

### Local copies of member variables
The `m_` prefix makes the bug visible; the fix is to assign locally before the lambda:

```cpp
const amrex::Real dt = m_dt;
const int n_modes = m_num_modes;
amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip) {
    // use dt, n_modes
});
```

### No `this->` in kernels
Direct member access (which implicitly captures `this`) is the main GPU portability hazard. Always capture local copies.

### No STL containers on device
`std::vector`, `std::map`, `std::string` are host-only. Use `amrex::Gpu::DeviceVector<T>`, raw pointers, or `amrex::GpuArray<T, N>` (fixed-size, device-accessible).

### `AMREX_RESTRICT` on pointer parameters in hot kernels
Equivalent to `__restrict__`. Promises non-aliasing, allows the compiler to vectorize.

```cpp
amrex::ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
```

### `AMREX_GPU_DEVICE` annotation on helper functions
Functions called from inside a `ParallelFor` kernel need this annotation:

```cpp
[[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real my_helper (amrex::Real x)
{
    return x * x;
}
```

For functions used on both host and device:

```cpp
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real my_helper (amrex::Real x);
```

### Random numbers on device
Use `amrex::ParallelForRNG`:

```cpp
amrex::ParallelForRNG(np,
    [=] AMREX_GPU_DEVICE (long ip, amrex::RandomEngine const& engine) {
        amrex::Real rnd = amrex::Random(engine);
        // use rnd
    });
```

Never call host `rand()` or `std::rand()` from a kernel.

## Common review feedback

Things reviewers will flag, based on the contributing guide and patterns visible in recent PRs:

- Missing `m_` prefix on a new member variable.
- A bare `0.0` literal (should be `0.0_rt` or `0.0_prt`).
- `using namespace amrex;` at any scope outside literals.
- Hand-rolled `for` loop over cells or particles where `ParallelFor` should be used.
- Implicit `this` capture in a `ParallelFor` lambda (any unqualified member name).
- Include order out of spec or missing blank lines between groups.
- New `.cpp` file added to one of `CMakeLists.txt` or `Make.package` but not the other.
- New runtime parameter added without documentation in `parameters.rst`.
- New feature without an automated test in `Tests/`.
- Unused `using namespace` directives.
- Trailing whitespace, missing final newline, tab characters.
- Function defined without space before `(`; or called with space before `(`.

Most of these are caught by `clang-tidy`, which is run in CI. Running `clang-tidy` locally before opening a PR (see `references/build-and-test.md`) catches them in advance.
