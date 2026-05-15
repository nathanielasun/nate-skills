# Collisions and reactions: MCC, DSMC, Coulomb, ionization

This is the C++ side of WarpX's collision and reaction infrastructure. For PICMI configuration of these from Python, see the `warpx-python` skill. Read this file when extending collision physics, adding scattering processes, or debugging collision code.

## Decision tree: which collision type for what

| Physical setting | Collision type | Subdirectory |
|---|---|---|
| Charged species in low-density thermal neutral gas, no neutral depletion | **MCC** (Background Monte Carlo) | `Source/Particles/Collision/BackgroundMCC/` |
| Binary collisions between two simulated species (incl. self-collisions) | **DSMC** | `Source/Particles/Collision/BinaryCollision/DSMC/` |
| Charged-charged scattering (Coulomb collisions within or across species) | **Coulomb** (Perez et al. 2012) | `Source/Particles/Collision/BinaryCollision/Coulomb/` |
| Nuclear fusion reactions (D-T, D-D, etc.) | **NuclearFusion** | `Source/Particles/Collision/BinaryCollision/NuclearFusion/` |
| Strong-field ionization (laser) | **FieldIonization** (ADK only) | `Source/Particles/Collision/FieldIonization/` (or similar — see source) |

These can coexist. Typical atmospheric-pressure microwave discharge: MCC for electron-neutral + Coulomb for charged-charged + sometimes DSMC for neutral dynamics if neutrals are simulated.

## Reuse scattering primitives

Don't reinvent kinematics. Use the existing utilities in `ParticleUtils::`:

```cpp
#include "Particles/ParticleUtils.H"

// Energy in COM frame (relativistic):
amrex::ParticleReal E_coll = ParticleUtils::getCollisionEnergy(...);

// Transform to/from COM frame:
ParticleUtils::doLorentzTransform(...);

// Isotropic scattering in COM:
ParticleUtils::RandomizeVelocity(ux, uy, uz, v_mag, engine);
```

Cross-section file loading is handled by the existing `MCCProcess` / `ScatteringProcess` classes — don't write your own file readers.

## Adding a new MCC process

Background MCC architecture: each `BackgroundMCCCollision` (a concrete `CollisionBase` subclass) owns a list of `MCCProcess` objects, one per scattering channel (elastic, ionization, etc.).

To add a new process kind (the existing kinds are: elastic, back-scattering, charge_exchange, excitation, ionization):

1. **Read existing precedent first.** Open `Source/Particles/Collision/BackgroundMCC/MCCProcess.H` and look at how existing process types are dispatched.
2. **Add an enum value** to whatever process-kind enum exists in `MCCProcess.H` (e.g., `MCCProcessType::Elastic`, `MCCProcessType::Ionization`).
3. **Implement the scattering kinematics** in a new function or method. Reuse `ParticleUtils` primitives.
4. **Wire it into the dispatcher** so that the input-file specification with the new process name routes to your code.
5. **Add a cross-section file** to `BLAST-WarpX/warpx-data` (separate PR).
6. **Write a test** — an MCC test in `Examples/Tests/`.

### Cross-section file format

Cross-section data is loaded at runtime from paths specified in the input file. The format is two-column ASCII:

```
# header lines start with #
# Reference: Hayashi M., 2003, ...
# Energy [eV]   Cross-section [m^2]
1.000000e-02    7.500000e-21
1.500000e-02    7.700000e-21
2.000000e-02    8.100000e-21
...
```

Linear interpolation between points. Energies outside the tabulated range use clamped endpoint values.

Files live in the separate `BLAST-WarpX/warpx-data` repo under `MCC_cross_sections/<species>/`:

```
MCC_cross_sections/argon/electron_scattering.dat
MCC_cross_sections/argon/ionization.dat
MCC_cross_sections/argon/excitation_1.dat
MCC_cross_sections/argon/ion_scattering.dat
...
```

Available species (check current repo): argon, helium, krypton, xenon, oxygen, nitrogen.

## Binary collision architecture

`BinaryCollision<ScatteringProcess, FunctorType>` is a template class. Each concrete binary collision (DSMC, Coulomb, fusion) plugs in:

- **`ScatteringProcess`**: defines the cross section and kinematics per pair (e.g., `PairWiseCoulombScatteringFunc`, `DSMCFunc`, `FusionFunc`).
- **`FunctorType`**: handles product creation (only for reactive collisions, e.g., fusion).

For non-reactive binary collisions (Coulomb, elastic DSMC), `FunctorType` is `NoCollisionFunctor`.

For collisions that create products (DSMC ionization, fusion):
- `FunctorType` defines how to compute number and weights of products
- Products are added via `SmartCopy`/`SmartCreate` from `Source/Particles/ParticleCreation/`

### Pair construction (cell-level pairing)

Binary collisions need to pair particles within a cell. The pairing logic lives in `Source/Particles/Collision/BinaryCollision/ShuffleFisherYates.H` (or similar — check current source). The standard algorithm:

1. Sort particles into cells (`Source/Particles/Sorting/`).
2. Within each cell, shuffle the indices.
3. Pair consecutive shuffled indices: (0,1), (2,3), …
4. Apply the scattering process to each pair.

Weight handling: when paired particles have different weights, the collision probability is scaled, and only some pairs undergo the collision (Nanbu-Yonemura or Takizuka-Abe variants). This is in `BinaryCollisionUtils.H`.

## Adding a new DSMC process

If adding a new DSMC channel (e.g., a new charge-exchange reaction):

1. Look at the existing `DSMCFunc` in `Source/Particles/Collision/BinaryCollision/DSMC/DSMCFunc.{H,cpp}`.
2. Extend the process-type dispatcher.
3. For reactive DSMC (creating products), define the product creation in the corresponding `FunctorType`.
4. Cross-section file in `warpx-data/MCC_cross_sections/<species>/` (same format as MCC).

The DSMC ionization product creation is non-trivial — products inherit the appropriate parent particle's properties; the new electron gets the projectile's leftover momentum; the new ion takes the neutral target's pre-collision velocity.

## Coulomb collisions

`Source/Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombScatteringFunc.{H,cpp}` implements the Perez 2012 algorithm.

Key implementation details:
- **Coulomb log**: can be a fixed value or computed from local plasma parameters (density, temperature). The local computation is in `ComputeCoulombLogarithm`.
- **Relativistic**: the algorithm handles relativistic particles correctly via 4-momentum.
- **Weight imbalance**: handled by the standard Nanbu-Yonemura weight correction.

To modify Coulomb scattering: edit `PairWiseCoulombScatteringFunc::operator()`. The function receives a particle pair and applies the velocity rotation in their COM frame.

## Nuclear fusion

`Source/Particles/Collision/BinaryCollision/NuclearFusion/` implements several reactions:
- D-T → He-4 + n
- D-D → He-3 + n
- D-D → T + p
- p-B → 3 He-4

Cross sections come from compiled tabulated data (in the source, not in `warpx-data`). The `NuclearFusionFunc` dispatches based on the reaction type.

To add a new fusion reaction: add the reaction type to the dispatcher, supply cross-section tabulation, implement product creation. Look at how D-T is handled as the template.

## How collisions plug into the time loop

In `WarpX::OneStep_nosub`:

1. Particle push happens (`mypc->Evolve(...)`)
2. **Then** collisions run: `mypc->doCollisions(...)`
3. **Then** field solve, deposition, etc.

The `mypc->doCollisions` iterates over all `CollisionBase` instances and calls each one's `doCoulombCollisions` / `doMCC` / `doBinaryCollisions` / etc.

Each collision class has an `ndt` parameter — apply every N steps (default 1). Useful for expensive collisions (Coulomb in dense plasmas).

**As of late 2025+**: there is a `collisions_split_position_push` option that places binary collisions in the *middle* of the position push for better energy conservation in electrostatic PIC. If touching the push/collision interface, check current `OneStep_nosub` for ordering.

## Particle creation from collisions

New particles from ionization, charge exchange, or fusion are created via `Source/Particles/ParticleCreation/`:

- **`SmartCopy`**: copy a particle from one species to another, with transformations.
- **`SmartCreate`**: create a new particle from a source set with derived properties.
- **`FilterCopyTransform`**: filter source particles, copy to a new species, transform attributes.

Example: MCC ionization in `Source/Particles/Collision/BackgroundMCC/`:
- Source particle (electron) is preserved, with reduced energy.
- A new ion is created from the (unsimulated) background neutral, at the source particle's position.
- A new electron is created with the ionization-fragment momentum.
- All new particles inherit the parent particle's weight.

The bookkeeping for cross-species creation (which species receives the ion product) lives in `MultiParticleContainer`, set up from the input parameters during `ReadParameters`.

## Common pitfalls

- **Energy conservation**: large `ndt` values batch collision events; energy errors per event scale as `ndt`. Use `ndt=1` for energy-sensitive simulations.
- **Cross-section file paths**: relative paths are interpreted relative to the WarpX launch directory, not the input file. Use absolute paths or document the expected layout clearly.
- **Maximum background density**: when `background_density` is an expression in MCC, `max_background_density` must be set so the null-collision-method probability precomputation has an upper bound.
- **Weight conservation**: when a collision creates products with different weights than the parent, total charge/momentum must still be conserved. Standard `SmartCopy`/`SmartCreate` handle this; custom product-creation code must too.
- **Product species not registered**: if your process refers to a product species name that isn't in `particles.species_names`, you'll either crash at init or silently lose products.
- **Cell vs MFP**: at high pressure, mean free path may be much smaller than cell size; binary collisions assume colliding pairs are physically near each other (within a cell). Refine the grid if needed.
- **Self-species pairing**: when pairing within one species with non-uniform weights, the standard pairing algorithm has small statistical biases. Resample to even weights before relying on self-collision results.
- **`Redistribute()` after creation**: new particles may end up in wrong boxes; call `Redistribute()` after creation events to fix.

## Testing collision changes

Every collision PR needs a test. Patterns:

1. **Benchmark against literature**: e.g., the Turner et al. (2013) capacitive discharge benchmarks for MCC validation. See `Examples/Physics_applications/capacitive_discharge/`.
2. **Energy conservation**: total energy before vs. after, for elastic collisions.
3. **Detailed balance**: under equilibrium conditions, rate of A→B should equal rate of B→A.
4. **Conservation laws**: total charge, total momentum within numerical error.

Add the test under `Examples/Tests/` with a small problem (1D or 2D, few-cell, short runtime) and a CTest that compares against a reference.

## Quick links

- Cross sections: `https://github.com/BLAST-WarpX/warpx-data`
- LXCat (canonical low-temp plasma cross-section database): `https://lxcat.net`
- Perez et al. 2012 (Coulomb algorithm): Phys. Plasmas 19, 083104
- Turner et al. 2013 (capacitive discharge benchmark): Phys. Plasmas 20, 013507
- Nanbu 1997 / Takizuka-Abe 1977 (weight handling): cited in `BinaryCollisionUtils.H`
