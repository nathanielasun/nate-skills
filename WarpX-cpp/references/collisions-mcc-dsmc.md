# Collisions, MCC, DSMC, and reactions

This reference covers the collision subsystem — the part of WarpX you'll spend the most time in for low-temperature-plasma work. It covers Background Monte Carlo Collisions (MCC), Direct Simulation Monte Carlo (DSMC), Coulomb collisions, nuclear fusion reactions, and the shared scattering primitives. Read this when adding or modifying any process that scatters, ionizes, recombines, or reacts particles.

## Table of contents

1. [Where collision code lives](#where-collision-code-lives)
2. [The collision-handler pattern](#the-collision-handler-pattern)
3. [Background MCC: against a thermal neutral background](#background-mcc-against-a-thermal-neutral-background)
4. [DSMC: binary collisions between simulated particles](#dsmc-binary-collisions-between-simulated-particles)
5. [Coulomb collisions](#coulomb-collisions)
6. [Nuclear fusion](#nuclear-fusion)
7. [Scattering primitives (shared)](#scattering-primitives-shared)
8. [Cross-section data: warpx-data and the LXCat format](#cross-section-data-warpx-data-and-the-lxcat-format)
9. [Particle creation during collisions](#particle-creation-during-collisions)
10. [Adding a new collision process: the playbook](#adding-a-new-collision-process-the-playbook)
11. [Pitfalls specific to collision code](#pitfalls-specific-to-collision-code)

## Where collision code lives

```
Source/Particles/Collision/
├── CollisionBase.{H,cpp}                  # Abstract base for all collision types
├── CollisionHandler.{H,cpp}               # Owner of all active collisions (in MultiParticleContainer)
├── ScatteringProcess.{H,cpp}              # Generic representation of a scattering process
├── BackgroundMCC/
│   ├── BackgroundMCCCollision.{H,cpp}     # The MCC driver class
│   ├── MCCProcess.{H,cpp}                 # Per-process descriptor (cross section, energy threshold, products)
│   ├── ImpactIonization.H                 # Impact ionization scattering
│   ├── BackgroundMCCCollision_K.H         # On-device kernels
│   └── ...
├── BackgroundStopping/
│   └── ...                                # Continuous slowing-down for fast particles in a background
├── BinaryCollision/
│   ├── BinaryCollision.{H,cpp}            # Templated driver for binary collisions
│   ├── ParticleCreationFunc.H             # Helpers for products of binary reactions
│   ├── ShuffleFisherYates.H               # Random shuffle for collision pairing
│   ├── Coulomb/                           # Perez et al. 2012 Coulomb collisions
│   ├── DSMC/                              # Direct Simulation Monte Carlo for neutrals
│   └── NuclearFusion/                     # D+D, D+T, p+B, ...
└── ScatteringProcess.{H,cpp}              # Shared scattering process abstraction
```

The exact file names shift across versions — verify in the current source tree before relying on a specific name.

## The collision-handler pattern

All active collisions in a simulation are owned by a `CollisionHandler` instance, which is itself a member of `MultiParticleContainer`. The handler holds a `std::vector<std::unique_ptr<CollisionBase>>` — one entry per `<collision_name>` declared in the input file.

`CollisionBase` is the polymorphic base. Each concrete collision type derives from it and implements:

```cpp
class CollisionBase {
public:
    virtual void doCollisions (amrex::Real cur_time, amrex::Real dt,
                                MultiParticleContainer* mypc) = 0;
    // ... other virtuals for diagnostics, restart, ...
};
```

`MultiParticleContainer::doCollisions()` is called once per step from `WarpX::OneStep_nosub`, and it dispatches to every active collision's `doCollisions()` method. Each collision then runs its own per-species loop and applies the appropriate scattering, particle removal, and particle creation.

When you add a new collision type, you create a class derived from `CollisionBase` (or `BinaryCollision<T>` if it's binary), implement `doCollisions()`, and register it in the handler's factory in `CollisionHandler::CollisionHandler(...)` so that the appropriate type is instantiated when a user declares it in the input file.

## Background MCC: against a thermal neutral background

`BackgroundMCCCollision` (in `Source/Particles/Collision/BackgroundMCC/`) handles collisions between simulated particles and an unsimulated thermal neutral background gas. This is the most important machinery for low-temperature-plasma work.

### Key design choices

- **Background is not simulated.** No neutral macroparticles. The background is described by a density profile and temperature, both supplied via input parameters. This is appropriate when the ionization fraction is small (<< 1) and the neutral gas is approximately equilibrated.
- **Null collision strategy** (Birdsall 1991) for efficiency. The maximum total collision cross section ν_max over a relevant energy range is precomputed. At each time step, only a fraction (1 - exp(-ν_max · dt)) of macroparticles is selected for further consideration; for those, a specific process is sampled. This reduces the per-particle work from "test every process for every particle" to "test only the pre-selected particles."
- **Relativistic energy calculation.** Collision energy is computed in the COM frame using `ParticleUtils::getCollisionEnergy()` (see [scattering primitives](#scattering-primitives-shared)).
- **Processes available**: elastic scattering (isotropic), back scattering, charge exchange, excitation (with energy loss), impact ionization (creates a new electron and ion).

### Configuration via input file

The user specifies:
- The simulated species: `<collision>.species` (one or more).
- The background species name: `<collision>.background_density`, `<collision>.background_temperature`, `<collision>.background_mass`.
- Processes: `<collision>.scattering_processes = elastic charge_exchange excitation1 ionization`, then per-process parameters (cross-section file paths, product species names, excitation energies, ionization energies).

### How `doCollisions` works

Per step, per species in the collision:

1. Compute `n_max = max(n_background)` over the relevant region; compute `ν_max = n_max * σ_total_max * v_max`.
2. Compute `P_max = 1 - exp(-ν_max * dt)`.
3. For each particle, draw a uniform `[0, 1]` random number. If < `P_max`, the particle is *considered* for a collision.
4. For considered particles, sample a background neutral velocity from a Maxwellian at the local temperature.
5. Boost to the COM frame, compute collision energy.
6. Evaluate σ for each enabled process at this energy. Sample which process occurs (or none, via the null-collision augmentation).
7. Apply the scattering: update particle velocity (elastic/back/excitation/charge exchange) or create products (ionization).

### Kernels

The hot work happens inside `ParallelFor` over particles. Cross-section interpolation is via `ScatteringProcess::getCrossSection(energy)` — a `[[nodiscard]] AMREX_GPU_HOST_DEVICE` function that does linear interpolation on a tabulated (energy, σ) array stored in device-accessible memory.

When adding a new MCC process (say, dissociative attachment, or a second excitation channel):

1. If the kinematics match an existing process (e.g., another elastic-like scattering), you can sometimes just add another `MCCProcess` instance to the species' process list — no new C++ code needed, just a different cross-section file and product specification.
2. If the kinematics are different, add a new scattering routine in a new header (e.g., `DissociativeAttachment.H`) following the pattern of `ImpactIonization.H`, and wire it up in `BackgroundMCCCollision::doCollisions`'s dispatch on process type.
3. Update `MCCProcessType` enum.
4. Update parser code in `MCCProcess::ReadParameters` (or equivalent) to accept the new process name as input.

## DSMC: binary collisions between simulated particles

DSMC handles binary collisions between two species that are *both simulated* as macroparticles. The standard use case: when neutral gas is dense enough that neutral-neutral collisions matter and you need to track neutrals as particles (e.g., near a wall, in a low-Kn region).

The algorithm (in `Source/Particles/Collision/BinaryCollision/DSMC/`):

1. **Sort particles by cell.** Use per-cell sorting helpers in `Source/Particles/Sorting/`.
2. **Shuffle within each cell.** `ShuffleFisherYates` provides a GPU-portable random shuffle.
3. **Pair particles within each cell.** If species A has more particles than species B in a cell, B's particles each pair with one A particle, and the surplus A particles are split (each pairs with one B partner counted multiple times). The pairing handles particle-weight imbalance via the standard Nanbu/Takizuka-Abe weighted-pair logic.
4. **For each pair, evaluate collision probability** based on the relative velocity, cross section, and macroparticle weights. Decide whether the pair collides.
5. **If colliding, scatter using the appropriate process** (elastic, charge exchange, excitation, ...). The scattering primitives are the same as MCC.

DSMC reuses MCC's `ScatteringProcess` infrastructure for cross sections and process descriptors. The differences are entirely in step 1–4 (pairing) and the bookkeeping for weight-mismatched pairs.

### Use cases

- **Atmospheric-pressure discharges** with explicit fast-neutral tracking: after charge exchange creates a fast neutral, you may want it in the DSMC system to thermalize via neutral-neutral collisions over a few hundred mean free paths.
- **Penning ionization** between simulated metastables.
- **Recombination** if you want to model it as a particle-particle process rather than a fluid sink.

### When to use MCC vs DSMC

| Use MCC when... | Use DSMC when... |
|---|---|
| The neutral gas is dilute enough that you can treat it as an unsimulated thermal background | Neutral-neutral collisions are frequent and matter (Kn < ~1) |
| Ionization fraction is small (no neutral depletion) | You need to resolve neutral depletion or non-Maxwellian neutral distributions |
| You don't need spatially resolved neutral dynamics | Fast neutrals from charge exchange need to be tracked |
| Wall interactions don't shape the bulk neutral distribution | You're near a wall or in a region with strong gradients |

The two are not mutually exclusive — a simulation can use MCC for electron-neutral processes (collisions are dominated by electrons sampling the background) and DSMC for neutral-neutral relaxation simultaneously.

## Coulomb collisions

Implemented in `Source/Particles/Collision/BinaryCollision/Coulomb/` following the Perez et al. (Phys. Plasmas 19, 083104, 2012) algorithm — relativistic, energy-and-momentum conserving binary Coulomb scattering. Input file declaration:

```
collisions.collision_names = coul1
coul1.species = electron ion
coul1.CoulombLog = 10.0   # or compute from local plasma parameters by omitting
coul1.ndt = 1             # apply every N time steps
```

The implementation reuses the binary-collision pairing infrastructure (shuffle, pair, scatter per pair) shared with DSMC. The scattering kernel itself is different — it computes a screened Rutherford deflection angle from the impact parameter distribution rather than using a discrete process selection.

When working on Coulomb collisions, key things to know:
- The Coulomb log can be specified or computed self-consistently from the local plasma parameters.
- Energy and momentum conservation are enforced per pair (recent PRs have improved this; see git log for details).
- Same-species collisions are supported (electron-electron, ion-ion).
- The collision frequency scales like 1/v³ at high energies, which means hot tails of distributions undergo few collisions per step — usually fine, but be alert for under-resolved relaxation.

## Nuclear fusion

Implemented in `Source/Particles/Collision/BinaryCollision/NuclearFusion/`. Supports D-D, D-T, D-He3, p-B11, and (as of recent versions) generalized two-product reactions including photon production. Uses the binary-pairing infrastructure plus reaction-specific cross sections (often from ENDF tables).

Products are created via the `ParticleCreationFunc` machinery — a templated functor that takes the reactant pair and produces new particles in the appropriate product species containers. The mechanism is GPU-portable (the creation happens inside the collision kernel).

When adding a new fusion reaction, the pattern is:
1. Implement the cross-section function (energy-dependent, possibly with multiple energy regimes).
2. Implement the reaction kinematics (product velocities given reactant velocities and Q-value).
3. Wire it into the fusion type enum and dispatch.
4. Add product-species mapping.

See `MultiParticleContainer::mapSpeciesProduct` for how the species name in the input file becomes a pointer to the correct destination container.

## Scattering primitives (shared)

These are in `Source/Particles/ParticleUtils.H` (header-only or with thin .cpp). They are device-callable and used by both MCC and DSMC:

### `ParticleUtils::getCollisionEnergy(m1, m2, u1, u2)`
Returns the relativistic collision energy in the COM frame:

```
E_coll = sqrt((γ m c² + M c²)² - (m u)²) - (m c² + M c²)
       = 2 M m u² / ((M + m + sqrt(M² + m² + 2 γ m M)) (γ + 1))
```

where `u = γv` is the normalized momentum and `m`, `M` are rest masses. In the classical limit (γ → 1), this reduces to the familiar `(1/2) (Mm/(M+m)) u²`. **Always use this function for collision-energy calculations**; the classical limit is wrong even for moderately fast ions.

### `ParticleUtils::doLorentzTransform(...)`
Boosts a 4-velocity between two frames via the general (non-axis-aligned) Lorentz transform. Used to go to/from the COM frame for scattering. The COM-frame boost velocity is computed from the relative velocity of the two collision partners.

### `ParticleUtils::RandomizeVelocity(...)`
Isotropic scattering in the COM frame. Given the magnitude of the COM velocity, generates a random orientation uniformly over the sphere and returns the new velocity vector.

### Scattering processes

Each of the four standard scattering types is implemented as a small kernel that takes a particle's pre-collision velocity (in the COM frame, in some processes) and a few process-specific parameters:

- **Elastic**: COM velocity magnitude preserved; direction randomized isotropically via `RandomizeVelocity`.
- **Back scattering**: COM velocity reversed (`-u_COM`).
- **Excitation**: COM velocity magnitude reduced according to the excitation energy threshold, then direction randomized isotropically.
- **Charge exchange**: For an ion-neutral collision producing a fast neutral and a slow ion, the labframe velocities are simply swapped — the ion becomes a neutral with the original neutral's velocity, and vice versa. No COM-frame rotation needed.
- **Impact ionization** (creates products): the projectile loses kinetic energy equal to the ionization threshold, then standard elastic scattering. A new electron is created with the appropriate energy partitioning (some fraction of the remaining energy after threshold), and a new ion is created with the neutral's pre-collision velocity.

When adding a new scattering kinematics, follow the pattern: keep the COM-frame transformation logic in the driver (so all processes share it) and write only the process-specific energy/direction update.

## Cross-section data: warpx-data and the LXCat format

Cross-section data is **not** stored in the main WarpX repository. It lives in `BLAST-WarpX/warpx-data/`, a separate repo that contains large auxiliary files (cross sections, reference fields, etc.). For MCC/DSMC work, the relevant subdirectory is `MCC_cross_sections/`, organized per atomic species:

```
warpx-data/MCC_cross_sections/
├── argon/
│   ├── electron_scattering.dat
│   ├── ionization.dat
│   ├── excitation_1.dat
│   └── ...
├── krypton/
├── xenon/
├── helium/
└── ...
```

File format: two-column ASCII, energy (eV) and cross section (m²), one row per data point, sorted by increasing energy. Linear interpolation is used between points. The data comes from various sources (Hayashi, Biagi, the LXCat database) — check the file headers for provenance.

When the user specifies `<process>.cross_section_file = path/to/data.dat` in the input deck, WarpX reads and tabulates the data at startup and stores it in device memory for fast access during the collision kernel.

For new species not in `warpx-data`, the appropriate place to contribute is via a PR to `warpx-data`, not by hardcoding cross sections in the WarpX source.

## Particle creation during collisions

When a collision produces new particles (ionization → new electron + ion; fusion → new products), creation happens *inside* the GPU kernel. The pattern uses a small set of utility functors:

- **`SmartCopy`** — copies an existing particle's attributes to a new particle, with optional per-attribute transforms.
- **`SmartCreate`** — creates a new particle from scratch with specified attributes.
- **`FilterCopyTransform`** — composable: filter which particles get copied, copy them, then transform the result.

These are defined in `Source/Particles/ParticleCreation/`. For collision use, the typical pattern is:

```cpp
// Inside the collision kernel, after deciding a pair undergoes ionization:
// 1. The original projectile loses energy and is scattered (modify in place)
// 2. A new electron is created (via the product species' SmartCreate machinery)
// 3. A new ion is created (via the product species' SmartCreate machinery)
```

The mechanics of finding the destination tiles for the new particles, atomically reserving slots, and writing the new particle data are handled by helpers in `BinaryCollision/ParticleCreationFunc.H`. You don't need to write the bookkeeping; you just specify *what* the new particles look like.

Particle creation is the trickiest part of writing a new collision process. **Use existing examples — ionization in `BackgroundMCC/ImpactIonization.H` and fusion product creation in `NuclearFusion/`.**

### Weight handling

Macroparticles in WarpX carry a weight `w` representing the number of real particles. When a collision occurs between particles of different weights:

- **MCC**: the simulated particle has weight `w`; the background is treated as a continuous medium, so the collision happens at the simulated particle's weight.
- **DSMC and other binary collisions**: with weight-mismatched pairs, the lower-weight particle dominates the collision rate. The standard Nanbu/Takizuka-Abe correction adjusts the collision probability per pair: if `w_A < w_B`, the pair has effective rate proportional to `min(w_A, w_B) / max(w_A, w_B)`. The implementation handles this in `BinaryCollision::doCoulombCollisionsWithinTile` (or similar) — check the current code for the exact factor.

When you create a product in an ionization or fusion event, its weight is typically set to the weight of the lower-weight parent particle, to conserve total physical particles. Charge conservation requires that the ion product weight equals the new electron weight.

## Adding a new collision process: the playbook

Putting it all together, here's the rough order of operations to add a new physical collision process:

1. **Decide which category the process is.** Single-species (slowing down)? Particle-background (MCC)? Pair (DSMC/Coulomb/fusion)? This determines which directory and which base class.

2. **Find the closest existing analog.** For low-temp plasma, e.g.:
   - Dissociative attachment → resembles ionization (creates a new particle, loses one)
   - Penning ionization → resembles fusion (two simulated particles → products)
   - Electron-electron Coulomb → already exists; same-species version is the template
   - Vibrational excitation → resembles existing excitation

3. **Cross-section data**: confirm it exists in `warpx-data`, or contribute it there. Add the path as an input file parameter for the user.

4. **Implement the kernel**: a small `[[nodiscard]] AMREX_GPU_HOST_DEVICE` function (or set of functions) that takes pre-collision state and process parameters, and returns post-collision state plus any product specifications.

5. **Wire into the driver** (`BackgroundMCCCollision::doCollisions` or `BinaryCollision::doCollisions`): add a case to the dispatch on process type. For new process types, extend the relevant `enum class`.

6. **Update parser** (`MCCProcess::ReadParameters`, `BinaryCollision::ReadParameters`, or similar) to accept the new process name and any new parameters.

7. **Document**: update `Docs/source/theory/multiphysics/collisions.rst` and `Docs/source/usage/parameters.rst`.

8. **Test**: write at least one unit-style test that exercises the new process. The test should verify a quantitative property — e.g., for an excitation process, energy loss per collision matches the threshold; for an ionization process, ionization rate matches the analytic rate computed from cross section and density. Examples: `Tests/Tests/test_*_mcc_*` test the existing MCC routines.

9. **Profile**: check that the new process doesn't dominate the per-step cost in a typical simulation. If it does, see whether the work can be amortized over multiple steps (via `ndt`) or moved out of the hot path.

## Pitfalls specific to collision code

- **Null collision sampling failure.** If `ν_max * dt > 1`, the null-collision approximation breaks down (the probability of no collision becomes negative). Check this assertion at startup, especially for high-energy beams or when adding processes with large cross sections.

- **Sample energy out of range.** Tabulated cross sections have a finite energy range. Particles with energies outside this range will use the clamped endpoint values. For very high-energy or very low-energy particles, this can silently give wrong rates. Add range checks if precision matters.

- **Frame of velocity arrays.** `ux`/`uy`/`uz` in WarpX are normalized momenta (γv/c is *not* the convention — it's γv, with v not divided by c, i.e., units of m/s). Always check whether a kernel expects `u` (proper velocity) or `v` (3-velocity), and convert correctly. `ParticleUtils::getCollisionEnergy` takes the WarpX-convention `u`.

- **Particle weight propagation through reactions.** When a parent particle ionizes, the resulting electron and ion must have the parent's weight to conserve charge. If you forget to set weight on a newly created particle, you get a default (often 0 or 1), which produces phantom physics.

- **Random number reproducibility.** For tests, you may need deterministic random numbers. WarpX uses `amrex::ParallelForRNG` with per-thread state — see `amrex::ResetRandomSeed`. Beware: GPU runs and CPU runs will not produce bit-identical results, so tests should check statistical properties, not exact trajectories.

- **Cross-check against known limits.** For any new process, derive an analytic rate in some limit (e.g., low-density, monoenergetic beam) and verify the simulation reproduces it before trusting more complex setups. Many collision bugs are quantitatively wrong by factors of order unity and go undetected without this check.

- **`ndt` (skip-steps) interaction with energy conservation.** If a collision module runs every N steps, the effective rate is correct on average but individual events transfer N times more energy than a per-step event would. For energetically violent events (impact ionization, fusion), this can produce spurious sub-grid heating. Use `ndt = 1` for these unless you have a specific reason.

- **MR and collisions.** Particles on different MR levels generally don't collide with each other. If your simulation has refinement near a collision-relevant region, verify the collision module behaves correctly across level boundaries. The existing handlers are tested for level-0 simulations; multi-level testing has been more sparse historically.
