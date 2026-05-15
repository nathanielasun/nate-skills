---
name: warpx-cpp-small
description: Use for any WarpX C++ codebase work — modifying internals, adding collision/scattering processes, particle pushers, field solvers, diagnostics; debugging C++ source under `Source/`; navigating the AMReX-based architecture; writing tests; preparing PRs. Trigger on "WarpX", "BLAST-WarpX", `WarpXParticleContainer`, `MultiParticleContainer`, `MultiFab` in WarpX context, MCC/DSMC, Boris/Vay/Higuera-Cary pushers, embedded boundaries, `Source/Particles/Collision`, electrostatic/electromagnetic/hybrid solvers, or any path under `Source/`. Optimized for low-temperature plasma (MCC, DSMC, ionization, atmospheric discharges). This is the compact variant tuned for mid-size models; for PICMI Python input scripts use the warpx-python skill.
---

# WarpX C++ Development (Compact)

This skill is for **modifying WarpX C++ code**. Not for PICMI input scripts — use `warpx-python` for those.

WarpX is built on **AMReX**. You cannot work productively without understanding AMReX patterns. The rules in this skill are not stylistic — code that breaks them won't build, won't run on GPUs, or will be rejected in review.

## Step 0: Check for in-repo context FIRST

Before using this skill, look for these in the user's checkout:

1. **`AGENTS.md`** (or `CLAUDE.md` symlink) at repo root — authoritative, version-matched.
2. **`.claude/skills/`** — in-repo workflow skills override general guidance.
3. **Doxygen** at `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/` — ground truth for C++ APIs.

**For any API: source → Doxygen → readthedocs → this skill.** Never invent an API. Always verify.

## What to read in this skill

| Your task involves... | Read |
|---|---|
| Navigating source, finding where things live, AMReX/WarpX class layout | `references/architecture.md` |
| Writing new C++: naming, includes, GPU portability, lambda captures | `references/coding-rules.md` |
| **MCC, DSMC, Coulomb, ionization, scattering processes (low-temp plasma)** | `references/collisions.md` |
| Building, testing, opening a PR | `references/build-test-pr.md` |

Read at most one or two per task.

## Hard rules (non-negotiable)

These appear repeatedly across all references. Internalize them now:

1. **Field data**: use `amrex::MultiFab`. Access via `WarpX::m_fields` (the `MultiFabRegister`). NEVER add a raw `std::unique_ptr<amrex::MultiFab>` member to `WarpX`.
2. **Particle data**: use `WarpXParticleContainer` (or a subclass). NEVER write your own particle container.
3. **Loops over cells or particles**: use `MFIter`/`ParIter` + `amrex::ParallelFor`. NEVER hand-roll nested `for` loops in hot paths — they won't vectorize and won't run on GPUs.
4. **Lambda captures**: copy members to local variables BEFORE the `ParallelFor`. The `m_` prefix on members exists to make capture-by-`this` bugs visible.
5. **FP literals**: always `0.0_rt` (for `amrex::Real`) or `0.0_prt` (for `amrex::ParticleReal`). NEVER bare `0.0` or `0.0f` — WarpX is precision-configurable at build.
6. **Namespace**: NEVER `using namespace amrex;` at file scope. Full `amrex::` qualification only. `using namespace amrex::literals;` inside a narrow scope is OK.
7. **New source files**: update BOTH `CMakeLists.txt` AND `Make.package` in the directory. Headers go in neither.
8. **Singletons**: avoid `WarpX::GetInstance()` in new code. Pass `WarpX&` or the specific container by reference.

## The PIC loop (orientation)

Every change happens in or around the PIC loop. The main loop is `WarpX::Evolve` in `Source/WarpXEvolve.cpp`. Per step (`WarpX::OneStep_nosub`):

1. **Field gather** — interpolate `Efield_aux`/`Bfield_aux` to particle positions (`PhysicalParticleContainer::Evolve` → `PushPX`)
2. **Particle push** — Boris/Vay/Higuera-Cary update of momentum and position
3. **Current deposition** — particles deposit J onto grid (`WarpXParticleContainer::DepositCurrent`)
4. **Field solve** — `WarpX::EvolveE` / `WarpX::EvolveB` advance Maxwell (FDTD/PSATD), or elliptic ES/MS/hybrid solve
5. **Collisions / reactions** — MCC, DSMC, Coulomb, ionization, fusion (`Source/Particles/Collision/`)
6. **Boundaries / diagnostics / filtering / load-balancing**

Note: as of late 2025+, the option `collisions_split_position_push` places binary collisions in the middle of the position push for energy conservation in ES PIC. If touching the collision/push interface, check current ordering in `OneStep_nosub`.

## Low-temp plasma quick reference

Code locations for MCC/DSMC/ionization work:

| Feature | Location |
|---|---|
| Background MCC (vs thermal neutral background) | `Source/Particles/Collision/BackgroundMCC/` |
| DSMC (binary collisions between simulated species) | `Source/Particles/Collision/BinaryCollision/DSMC/` |
| Coulomb collisions | `Source/Particles/Collision/BinaryCollision/Coulomb/` |
| Nuclear fusion | `Source/Particles/Collision/BinaryCollision/NuclearFusion/` |
| Hybrid kinetic-fluid (Ohm's law) | `Source/FieldSolver/HybridPICModel/` |
| Cross-section data | `BLAST-WarpX/warpx-data` repo, `MCC_cross_sections/<species>/` |
| Scattering primitives | `ParticleUtils::getCollisionEnergy`, `doLorentzTransform`, `RandomizeVelocity` |

For deeper coverage see `references/collisions.md`.

## Working pattern for a typical task

1. **Restate the task** in terms of subsystems touched (fields? particles? collisions?). If >2, propose a decomposition before writing code.
2. **Find existing precedent**. New physics almost always resembles existing physics. For a new collision process: look at existing MCC processes. For a new pusher: look at Boris/Vay. Copy and adapt — designing from scratch is almost always wrong.
3. **Read the relevant reference file** in this skill.
4. **Small first iteration**: get something that compiles for one dimension and one compute backend. Then generalize.
5. **Add a test** — WarpX won't accept new features without tests (see `references/build-test-pr.md`).
6. **Run `clang-tidy` locally** before pushing (it runs in CI).
7. **Add documentation**: runtime parameters → `Docs/source/usage/parameters.rst`; new physics → theory section.
8. **One feature per PR.** Style separate from substance.

## Critical caveats

- **The codebase changes monthly.** Verify any specific API against current source.
- **GPU portability has invisible constraints.** Code that compiles fine on CPU can silently fail on GPU. Test with `WarpX_COMPUTE=CUDA` (or HIP/SYCL) before assuming portability.
- **Developer docs are incomplete.** Some readthedocs pages are stubs; read actual source for ground truth.
- **LLM-generated code requires manual review** before requesting maintainer review. Don't push code you haven't read.

## Canonical links

- WarpX: `https://github.com/BLAST-WarpX/warpx`
- Docs: `https://warpx.readthedocs.io/`
- Doxygen: `https://warpx.readthedocs.io/en/latest/_static/doxyhtml/`
- AMReX: `https://amrex-codes.github.io/amrex/docs_html/` + `https://amrex-codes.github.io/amrex/doxygen/`
- warpx-data: `https://github.com/BLAST-WarpX/warpx-data`
- Discussions: `https://github.com/BLAST-WarpX/warpx/discussions`
