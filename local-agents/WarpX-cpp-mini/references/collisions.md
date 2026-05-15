# Collisions: MCC, DSMC, Coulomb, ionization

For low-temperature plasma C++ work. To configure collisions from PICMI, see `warpx-python-mini`.

## Decision table

| Physical setting | Class | Subdirectory |
|---|---|---|
| Charged species in unsimulated thermal neutral gas | `BackgroundMCCCollision` | `Source/Particles/Collision/BackgroundMCC/` |
| Binary collisions between simulated species | DSMC (`BinaryCollision`-based) | `Source/Particles/Collision/BinaryCollision/DSMC/` |
| Charged-charged scattering | Coulomb (`BinaryCollision`-based) | `Source/Particles/Collision/BinaryCollision/Coulomb/` |
| Nuclear fusion | `BinaryCollision`-based | `Source/Particles/Collision/BinaryCollision/NuclearFusion/` |
| Strong-field ionization (laser) | `FieldIonization` (ADK only) | `Source/Particles/Collision/FieldIonization/` |

These coexist. Typical low-temp plasma: MCC alone at low ionization. Add Coulomb at high ionization. Add DSMC if neutrals are simulated.

## Always reuse scattering primitives

```cpp
#include "Particles/ParticleUtils.H"

// COM-frame energy (relativistic):
amrex::ParticleReal E_coll = ParticleUtils::getCollisionEnergy(...);

// Lorentz transform to/from COM:
ParticleUtils::doLorentzTransform(...);

// Isotropic scattering in COM:
ParticleUtils::RandomizeVelocity(ux, uy, uz, v_mag, engine);
```

Don't write your own kinematics — use these.

## MCC architecture

`BackgroundMCCCollision` owns a list of `MCCProcess` objects, one per channel (elastic, ionization, etc.). Each `MCCProcess` reads a cross-section file and applies the kinematics.

To add a new process kind:

1. **Read existing precedent in `Source/Particles/Collision/BackgroundMCC/MCCProcess.{H,cpp}`.**
2. Add an enum value to the process-kind enum.
3. Implement kinematics using `ParticleUtils` primitives.
4. Wire the input-file-name → process-kind dispatcher.
5. Add cross-section data to `BLAST-WarpX/warpx-data` (separate PR).
6. Add a test in `Examples/Tests/`.

## Cross-section file format

Two-column ASCII, in `BLAST-WarpX/warpx-data` repo under `MCC_cross_sections/<species>/`:

```
# Source: ...
# Energy [eV]   Cross-section [m^2]
1.000000e-02    7.500000e-21
1.500000e-02    7.700000e-21
...
```

Linear interpolation between points. Endpoints clamped outside the range.

Available species (check current repo): argon, helium, krypton, xenon, oxygen, nitrogen.

## Binary collision architecture

`BinaryCollision<ScatteringProcess, FunctorType>` template:

- **`ScatteringProcess`**: per-pair cross section + kinematics. Examples: `PairWiseCoulombScatteringFunc`, `DSMCFunc`, `NuclearFusionFunc`.
- **`FunctorType`**: product creation for reactive collisions. `NoCollisionFunctor` for non-reactive.

Pair construction (cell-level):
1. Sort particles into cells (`Source/Particles/Sorting/`)
2. Shuffle within each cell
3. Pair consecutive shuffled indices: (0,1), (2,3), ...
4. Apply scattering to each pair

Weight handling for unequal-weight pairs: Nanbu-Yonemura or Takizuka-Abe (in `BinaryCollisionUtils.H`).

## Adding a new DSMC process

1. Read `Source/Particles/Collision/BinaryCollision/DSMC/DSMCFunc.{H,cpp}`.
2. Extend the process-type dispatcher.
3. For reactive DSMC: define product creation in the `FunctorType`.
4. Cross-section file in `warpx-data/MCC_cross_sections/<species>/`.

## Coulomb collisions

`Source/Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombScatteringFunc.{H,cpp}`. Perez 2012 algorithm. Edit `operator()` to modify the velocity rotation.

The Coulomb log can be fixed or computed locally (`ComputeCoulombLogarithm`). Relativistic via 4-momentum. Weight imbalance handled by standard correction.

## Nuclear fusion

`Source/Particles/Collision/BinaryCollision/NuclearFusion/`. Supports D-T, D-D (→ He-3 + n), D-D (→ T + p), p-B. Dispatcher in `NuclearFusionFunc`.

To add a new reaction: extend dispatcher, supply cross-section, implement product creation. Use D-T as the template.

## Where collisions plug in

`WarpX::OneStep_nosub` sequence:
1. Particle push (`mypc->Evolve(...)`)
2. **Collisions** (`mypc->doCollisions(...)`)
3. Field solve, deposition, etc.

`mypc->doCollisions` iterates over all `CollisionBase` instances and calls each one's per-type method.

Each collision has `ndt` (apply every N steps). Use `ndt=1` for energy-sensitive simulations.

**Late 2025+**: option `collisions_split_position_push` places binary collisions in the middle of the position push for energy conservation in ES PIC. Check current `OneStep_nosub` before touching the push/collision interface.

## Particle creation from collisions

```cpp
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartCreate.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
```

- **`SmartCopy`**: copy particle from one species to another with transformations.
- **`SmartCreate`**: create new particle from a source set.
- **`FilterCopyTransform`**: filter + copy + transform.

After collision-driven particle creation, call `pc.Redistribute()`.

Cross-species bookkeeping (which species receives ionization products) is set up in `MultiParticleContainer::ReadParameters`.

## Common pitfalls

- **Energy conservation**: large `ndt` batches collision events; use `ndt=1` if energy matters.
- **File path interpretation**: relative paths are vs launch dir, not script. Use absolute.
- **Product species not registered**: name in `particles.species_names` must exist or you'll crash or lose products.
- **Cell vs MFP**: at atmospheric pressure with mm cells, DSMC pairing logic breaks (mfp ~μm). Refine grid.
- **Self-species pairing with unequal weights**: resample first to even out weights.
- **Forgotten `Redistribute()`** after creating new particles in callbacks → wrong-box particles.

## Test pattern

Every collision PR needs a test:

1. **Benchmark literature** (e.g., Turner et al. 2013 for MCC validation: `Examples/Physics_applications/capacitive_discharge/`)
2. **Energy conservation** for elastic collisions
3. **Detailed balance** at equilibrium
4. **Conservation laws** (charge, momentum) within numerical error

Add to `Examples/Tests/` with small problem (1D or 2D, few cells, short runtime) and a CTest comparing to a reference.

## Quick links

- Cross sections: `https://github.com/BLAST-WarpX/warpx-data`
- LXCat (cross-section database): `https://lxcat.net`
- Perez 2012 (Coulomb): Phys. Plasmas 19, 083104
- Turner 2013 (MCC benchmark): Phys. Plasmas 20, 013507
