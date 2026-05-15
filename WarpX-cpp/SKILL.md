---
name: warpx-cpp
description: Use whenever the user works on WarpX's C++ codebase (BLAST-WarpX particle-in-cell code) — modifying internals, adding physics modules, extending solvers, debugging C++ source, navigating the AMReX-based architecture, writing tests, or preparing PRs. Trigger on "WarpX", "BLAST-WarpX", `WarpXParticleContainer`, `MultiParticleContainer`, `MultiFab` in a WarpX context, MCC/DSMC modules, particle pushers (Boris/Vay/Higuera-Cary), embedded boundaries, `Source/Particles/Collision`, electrostatic/electromagnetic/hybrid solvers, or any path under `Source/` in the WarpX tree. Also trigger when asked to add a new collision type, scattering process, field solver, or diagnostic to WarpX, or to port C++ code between CPU and GPU. Optimized for low-temperature-plasma work (MCC/DSMC/ionization) but covers the full C++ codebase. For C++ work; for PICMI Python input scripts, defer to warpx-python skill if available.
---

# WarpX C++ Development Skill

This skill is for **modifying or extending WarpX's C++ codebase** — not for writing PICMI input decks (that's a separate concern). It is biased toward low-temperature plasma work (MCC, DSMC, ionization, atmospheric discharges) but covers the full C++ stack.

WarpX (BLAST-WarpX) is an exascale-class particle-in-cell code built on the AMReX adaptive mesh refinement framework. The codebase is large, performance-critical, and follows specific conventions that diverge from typical modern C++. Working with it productively requires understanding (a) the AMReX data-structure foundation, (b) the WarpX class hierarchy that sits on top, (c) the GPU-portable kernel pattern, and (d) the project's style and PR conventions.

## Before doing anything else: check for in-repo context

WarpX ships with its own LLM-facing context. If a local checkout of WarpX is available, read these *first*, before consulting this skill's reference files — they are authoritative and version-matched to the user's tree:

- **`AGENTS.md`** at the repo root (with `CLAUDE.md` as a symlink to it). Kept under 300 lines, contains build/test/style instructions specific to the current state of the codebase. Anything in there overrides the general guidance in this skill.
- **`.claude/skills/`** directory. WarpX defines its own per-task skills (currently `/warpx-answer-user-question` and `/warpx-new-paper-highlight`). If the user's task matches one of these workflows, prefer the in-repo skill.
- **Doxygen** at `https://warpx.readthedocs.io/en/latest/_static/doxyhtml` (or the in-tree generated copy). This is the ground-truth for C++ APIs — class hierarchies, member functions, signatures.

When in doubt about a specific API or symbol, the authoritative chain is: source code → Doxygen → readthedocs → this skill. Never invent an API; verify against one of the first three.

If MCP servers for Context7 are available, the WarpX, AMReX, pyAMReX, openPMD, PICSAR, and PICMI documentation are queryable through that channel. Prefer Context7 queries to web searches for these specific topics.

## What to read in this skill, when

This skill is intentionally guide-style. The SKILL.md body teaches structure and decision-making; depth lives in the `references/` directory. Read references as needed for the task at hand — don't pull all of them up front.

| If the task involves... | Read |
|---|---|
| Understanding code organization, navigating the source tree, finding where a feature lives | `references/architecture.md` |
| Touching field data (E/B/J/ρ), field solvers, MultiFabs, field-gather, current-deposition | `references/fields-and-particles.md` (fields section) |
| Touching particle data, particle pushers, species containers, runtime attributes, loops over particles | `references/fields-and-particles.md` (particles section) |
| **Adding/modifying MCC, DSMC, ionization, scattering, charge exchange, Coulomb collisions, fusion reactions** | `references/collisions-mcc-dsmc.md` |
| Writing new code: naming, includes, FP literals, GPU portability, lambda captures, member-variable conventions | `references/style-and-conventions.md` |
| Building locally, dependencies, CMake options, debugging compile failures, writing/running tests | `references/build-and-test.md` |
| Preparing a PR, git workflow, documentation requirements | `references/contributing-workflow.md` |

The references are designed so that **for any concrete task you can identify exactly one or two files to read**. If you find yourself needing all six, the task is probably actually multiple tasks — break it down.

## The PIC loop, briefly, as orientation

Every time something happens in WarpX, it happens in or around the PIC loop. Knowing where you are in the loop tells you which data structures are live and which assumptions hold.

The main loop is `WarpX::Evolve` in `Source/WarpXEvolve.cpp`. Its non-subcycled core is `WarpX::OneStep_nosub` (use `OneStep_sub1` for subcycled MR). At each time step the canonical PIC sequence is:

1. **Field gather** — interpolate `Efield_aux`/`Bfield_aux` from grid to particle positions. Per species in `PhysicalParticleContainer::Evolve` → `PushPX`.
2. **Particle push** — Boris/Vay/Higuera-Cary updates particle momentum and position. Selected by `algo.particle_pusher` (input file) and dispatched via `PushType`/`PushAlgo` enums.
3. **Current deposition** — particles deposit J onto the grid via `WarpXParticleContainer::DepositCurrent`. Charge-conserving variants exist (Esirkepov, Villasenor-Buneman). Charge deposition `DepositCharge` is separate.
4. **Field solve** — `WarpX::EvolveE` / `WarpX::EvolveB` advance Maxwell's equations (FDTD/PSATD/etc.) using J as source. Electrostatic/magnetostatic/hybrid-Ohm solvers replace this with elliptic solves.
5. **Collisions / reactions / other operators** — MCC, DSMC, Coulomb, ionization, fusion run at user-configurable intervals. See `Source/Particles/Collision/`.
6. **Boundary / diagnostics / filtering / load-balancing** — additional operators that run on a cadence.

Where in this sequence collisions run matters: as of late 2025 / early 2026, WarpX added a `collisions_split_position_push` option that places binary collisions in the middle of the position push for better energy conservation in electrostatic PIC. If you're touching the collision/push interface, check the current ordering in `OneStep_nosub` before making assumptions.

## Decision rules when writing new C++ code

A few high-value rules that come up over and over. The references expand each one.

**On data structures.** Never invent your own grid container or particle container. Use `amrex::MultiFab` for grid data and a `WarpXParticleContainer`-derived class for particles. Use `ablastr::fields::MultiFabRegister` (accessed via `WarpX::m_fields`) for new field types — don't add a raw `std::unique_ptr<amrex::MultiFab>` member to `WarpX` unless you have a strong reason. New field types go through the `warpx::fields::FieldType` enum.

**On loops.** Never write hand-rolled nested `for` loops over grid cells or particle indices in hot paths. Use the iterator + `amrex::ParallelFor` pattern documented in `references/fields-and-particles.md`. This is non-negotiable: hand-rolled loops won't vectorize, won't run on GPUs, and won't be accepted in review.

**On lambda captures.** Capture members by value into local variables *before* the `ParallelFor`. Member variables (prefixed `m_`) captured by `this` will copy the whole `WarpX` object to GPU. The `m_` prefix on members exists precisely to make this mistake visible.

**On floating-point literals.** Always `0.0_rt` (for `amrex::Real`) or `0.0_prt` (for `amrex::ParticleReal`), never `0.0` or `0.0f`. WarpX is configured at build time for single or double precision — bare literals defeat that.

**On includes.** Strict order, see `references/style-and-conventions.md`. Forward-declaration headers (`*_fwd.H`) are used aggressively to reduce compile time. When you write a new class, write a corresponding `_fwd.H` if other headers will reference it by pointer.

**On the `amrex::` namespace.** Never `using namespace amrex;` at file scope. Use full `amrex::` qualification. `using namespace amrex::literals;` is acceptable inside a narrow scope to enable the `_rt`/`_prt` suffixes.

**On adding new source files.** Both `CMakeLists.txt` *and* `Make.package` in the relevant subdirectory must be updated. Headers are not listed in either; only `.cpp` files.

**On member-variable access in functions exported to other compilation units.** Avoid `WarpX::GetInstance()` in new code where possible. Pass `WarpX&` or the specific containers you need by reference instead. The codebase is being gradually weaned off the singleton pattern.

## Low-temperature-plasma quick reference

If the user is doing low-temperature plasma work — neutral gas, MCC ionization, DSMC, atmospheric discharges, hybrid-PIC for ion kinetics — this is where the relevant code lives:

- **Background MCC** (collisions against a thermal neutral background): `Source/Particles/Collision/BackgroundMCC/`. Key classes: `BackgroundMCCCollision`, `MCCProcess`. Processes supported: elastic, back-scattering, charge exchange, excitation, ionization.
- **DSMC** (binary collisions between simulated particles): `Source/Particles/Collision/BinaryCollision/DSMC/`. Same scattering primitives as MCC, applied to pairs of simulated particles within a cell.
- **Coulomb collisions** (between simulated charged species): `Source/Particles/Collision/BinaryCollision/Coulomb/`. Perez et al. 2012 algorithm.
- **Nuclear fusion**: `Source/Particles/Collision/BinaryCollision/NuclearFusion/`. Several fusion reactions supported.
- **Particle utilities for scattering**: `ParticleUtils::getCollisionEnergy()` (relativistic), `ParticleUtils::doLorentzTransform()` (to/from COM frame), `ParticleUtils::RandomizeVelocity()` (isotropic scattering). Use these primitives for any new scattering process — don't reinvent them.
- **Cross-section data**: hosted in the separate `BLAST-WarpX/warpx-data` repository under `MCC_cross_sections/<species>/`. Tab-separated energy, cross-section format. Loaded at runtime via paths in the input file.
- **Hybrid kinetic-fluid (Ohm's law)**: `Source/FieldSolver/HybridPICModel/`. For ion-kinetic problems where electrons are a massless fluid.

For deeper guidance on extending or modifying any of these, see `references/collisions-mcc-dsmc.md`.

## Working pattern for a typical task

For most non-trivial tasks, this sequence works well:

1. **Restate the task** in terms of which C++ subsystems it touches (fields? particles? collisions? a new solver? a diagnostic?). If it spans more than two, propose a decomposition before writing code.
2. **Find existing precedent.** WarpX has been developed for years and most new physics resembles existing physics. For a new collision process, look at how existing MCC processes are structured. For a new pusher, look at Boris/Vay/Higuera-Cary. For a new field solver, look at FDTD/PSATD. Copying and adapting is overwhelmingly the right pattern; designing from scratch is overwhelmingly wrong.
3. **Read the relevant reference file(s).** Don't try to reason about AMReX patterns from C++ first principles — the conventions matter and they're documented.
4. **Write the code with a small first iteration** — get something that compiles for a single dimension and a single compute backend before generalizing. WarpX builds separate binaries per dimension, so a 2D-only first cut is a legitimate workflow.
5. **Add a test.** New features without tests will not be accepted in review. See `references/build-and-test.md` for how to write one.
6. **Run the linter and at least one local build.** `clang-tidy` is run in CI; running it locally catches embarrassments early.
7. **Add documentation.** New runtime parameters in `Docs/source/usage/parameters.rst`. New physics in the theory section. New algorithms in the implementation-details section if appropriate.
8. **Open a focused PR.** One feature per PR. Style changes separate from substance changes. Reference the issue or discussion if one exists.

## Honest caveats

- **The codebase changes weekly.** Cadenced monthly releases (e.g., 26.05, 26.06) and active development on `development` branch mean APIs evolve. Verify any API mentioned in this skill against the current source before relying on it.
- **The official "developer" docs are incomplete.** Several pages on readthedocs (`Initialization`, `Diagnostics`, `Portability`) are stubs as of early 2026. The references in this skill fill in some of those gaps but you should read the actual source for ground truth.
- **GPU portability has constraints not visible from C++ syntax.** Some patterns that compile fine on CPU silently produce wrong results or worse performance on GPU. Always check that new kernels run correctly with `WarpX_COMPUTE=CUDA` (or HIP/SYCL) before assuming portability.
- **This skill is for C++ work.** If the user is writing a PICMI input script, asks about Python bindings via pyAMReX, or wants to know how to *use* WarpX rather than modify it, point them to the (anticipated) `warpx-python` skill or to the user-facing documentation. Don't try to answer PICMI questions from this skill — the conventions are different.
- **LLM-generated code in WarpX has explicit community expectations.** The official `how_to_develop_with_llms.html` page is unambiguous: a careful, manual review by the human contributor is required before requesting maintainer review. Don't push code through this skill that the user has not read and understood.

## Where to find canonical references

- **WarpX GitHub**: `https://github.com/BLAST-WarpX/warpx`
- **Documentation**: `https://warpx.readthedocs.io/`
- **Doxygen**: `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/`
- **AMReX docs**: `https://amrex-codes.github.io/amrex/docs_html/`
- **AMReX Doxygen**: `https://amrex-codes.github.io/amrex/doxygen/`
- **Discussions** (for questions): `https://github.com/BLAST-WarpX/warpx/discussions`
- **warpx-data** (cross sections, large auxiliary data): `https://github.com/BLAST-WarpX/warpx-data`
- **AMReX repository**: `https://github.com/AMReX-Codes/amrex`
- **PICSAR (QED)**: `https://github.com/ECP-WarpX/picsar`

Start here when verifying APIs, looking for precedent, or answering questions you're not certain about. WarpX is sufficiently large and specialized that confident-sounding answers from training data are frequently wrong; verification against the actual source is the only reliable approach.
