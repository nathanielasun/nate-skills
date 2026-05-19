# Particles and species — source layer

The `Particles` class (SoA storage), `Species`, particle pushers, particle boundary conditions, injection, and the per-patch `dynamics` orchestration. The C++-side counterpart to namelist `Species` blocks.

## Table of contents

1. `Particles` — SoA layout
2. `Species` — the per-species container
3. Per-patch `dynamics` — the call chain
4. Pushers — `PusherBoris`, `PusherVay`, etc.
5. Interpolator — fields to particle
6. Projector — particle to current grid
7. Particle boundary conditions
8. `ParticleInjector`
9. Test particles
10. Particle merging
11. Adding a new pusher
12. Common pitfalls

---

## 1. `Particles` — SoA layout

`src/Particles/Particles.h` defines the SoA particle container. The key fields:

```cpp
class Particles {
public:
    std::vector<double> position[3];          // [0]=x, [1]=y, [2]=z (only used dims)
    std::vector<double> momentum[3];          // [0]=px, [1]=py, [2]=pz
    std::vector<double> weight;               // statistical weight per macro-particle
    std::vector<short>  charge;               // current charge state (integer, in units of e)
    std::vector<double> chi;                  // quantum parameter (radiation reaction)
    std::vector<uint64_t> id;                 // particle ID (for tracking)

    // Mass and other species-level constants are NOT per-particle —
    // they live on Species.

    // Accessors
    size_t   size()                  const { return weight.size(); }
    double*  getPtrPosition(int d)         { return position[d].data(); }
    double*  getPtrMomentum(int d)         { return momentum[d].data(); }
    short*   getPtrCharge()                { return charge.data(); }
    double*  getPtrWeight()                { return weight.data(); }
    double*  getPtrChi()                   { return chi.data(); }

    // Modifiers
    void resize(size_t n);                    // resize all arrays
    void reserve(size_t n);                   // pre-allocate
    void clear();
    void eraseParticle(size_t i);             // tombstone + swap with last; O(1)
};
```

**Key invariants:**

- All vectors share the same size — `weight.size() == position[0].size() == charge.size()` etc.
- Particles are indexed `0..size()-1`. The order is **not stable** — collisions, sorting, and migration reorder particles arbitrarily.
- `chi` is allocated only when radiation reaction is active; otherwise zero-sized.
- `id` is allocated only when particle tracking is active.

**Why SoA:**

- Vectorization: physics loops `for (i = 0; i < N; ++i) px[i] += ...` are SIMD-friendly.
- Cache: a typical Boris push touches `position[0..2]` and `momentum[0..2]` — 6 contiguous streams.
- GPU: SoA maps cleanly to one-thread-per-particle kernels.

**Adding a new per-particle attribute:**

1. Add a `std::vector<T>` field to `Particles`.
2. Update `resize()`, `clear()`, `reserve()`, `eraseParticle()` to handle the new vector.
3. Add a getter `getPtr*()`.
4. Update HDF5 dump/restart paths in `Diagnostic/` and `Checkpoints/`.
5. Document the new attribute in the doc.

This is intrusive — most physics modules don't need new attributes; they use existing ones with new physics meanings or new derived classes.

## 2. `Species` — the per-species container

`src/Particles/Species.h` defines `Species` — one instance per species per patch.

```cpp
class Species {
public:
    std::string name;                          // "electron", "Ar1", etc.
    double mass;                               // in m_e
    double charge;                             // initial charge state (can be incremented by ionization)
    short  atomic_number;                      // for ionization

    Particles* particles;                      // the SoA container for this species in this patch

    Pusher* Push;                              // species-specific pusher
    Interpolator* Interp;                      // field → particle
    Projector* Proj;                           // particle → current grid

    Ionization* Ionize;                        // field ionization (may be nullptr)
    // ...

    ParticleBoundaryConditions partBoundCond;  // BCs per direction/face

    // Methods
    void dynamics(double time, ...);           // per-step push + projection + boundaries
    void initializeParticles(...);              // namelist-driven initialization
    void sortParticles(...);                    // optional cell-based sorting
    // ...
};
```

**Critical**: `Species` is per-patch. A simulation with 4 patches and 3 species has 12 `Species` instances on each rank. Don't store global state on `Species`.

**Cross-references**:

- `Species*` instances are stored in `Patch::vecSpecies` (a `std::vector<Species*>`).
- Across patches, species share name/mass/atomic_number but have independent `particles` arrays.
- For diagnostics that aggregate across patches (e.g., `DiagScalar`), the aggregation walks `vecPatches[i]->vecSpecies[s]` for each rank.

## 3. Per-patch `dynamics` — the call chain

`Species::dynamics` (called from `Patch::dynamics`, called from `VectorPatch::dynamics`) is the main per-step physics for a species in a patch:

```cpp
void Species::dynamics(double time, ...) {
    // Resize current arrays for this species' contribution
    // ...

    // Field-ionization step (may add particles to ionization_electrons species)
    if (Ionize) {
        Ionize->apply(particles, EMfields);
    }

    // Field interpolation: read E, B at each particle's position
    Interp->fieldsAtParticle(EMfields, particles);

    // Pusher: update momentum and position
    (*Push)(*particles, ...);

    // Particle boundary conditions
    partBoundCond.apply(particles, ...);

    // Current deposition
    Proj->depositCurrents(JEM, *particles);
}
```

**Order matters.** Field ionization runs *before* the pusher so freed electrons see the field they would have at their birth time. Boundary conditions run *before* projection because particles outside the box don't deposit current.

Modifications to this order require care — many physics invariants depend on the sequence.

## 4. Pushers — `PusherBoris`, `PusherVay`, etc.

`src/Pusher/`. The base class:

```cpp
class Pusher {
public:
    virtual void operator()(Particles& particles, ...) = 0;
};
```

Concrete pushers:

| File | Pusher | Use |
|---|---|---|
| `PusherBoris.cpp` | Boris (default) | All-purpose relativistic |
| `PusherBorisNR.cpp` | Non-relativistic Boris | `a₀ ≪ 1`, faster |
| `PusherVay.cpp` | Vay | Better E×B at high `γ` |
| `PusherHigueraCary.cpp` | Higuera-Cary | Best accuracy at high `γ` |
| `PusherPhoton.cpp` | Photon (mass=0) | Ballistic at `c` |

**Boris pusher inner loop** (sketch from `PusherBoris.cpp`):

```cpp
void PusherBoris::operator()(Particles& particles, ...) {
    double dts2 = 0.5 * dt;
    const double charge_over_mass = charge / mass;

    #pragma omp simd
    for (size_t ip = ipart_min; ip < ipart_max; ++ip) {
        double Ex = local_E[0][ip], Ey = local_E[1][ip], Ez = local_E[2][ip];
        double Bx = local_B[0][ip], By = local_B[1][ip], Bz = local_B[2][ip];

        double umx = momentum[0][ip] + dts2 * charge_over_mass * Ex;
        double umy = momentum[1][ip] + dts2 * charge_over_mass * Ey;
        double umz = momentum[2][ip] + dts2 * charge_over_mass * Ez;

        double gm = std::sqrt(1.0 + umx*umx + umy*umy + umz*umz);
        double tx = dts2 * charge_over_mass * Bx / gm;
        double ty = dts2 * charge_over_mass * By / gm;
        double tz = dts2 * charge_over_mass * Bz / gm;
        // ... Boris rotation ...

        momentum[0][ip] = up[0] + dts2 * charge_over_mass * Ex;
        momentum[1][ip] = up[1] + dts2 * charge_over_mass * Ey;
        momentum[2][ip] = up[2] + dts2 * charge_over_mass * Ez;

        // Position update
        double inv_gamma = 1.0 / std::sqrt(1.0 + momentum[0][ip]*momentum[0][ip] + ...);
        position[0][ip] += dt * momentum[0][ip] * inv_gamma;
        position[1][ip] += dt * momentum[1][ip] * inv_gamma;
        position[2][ip] += dt * momentum[2][ip] * inv_gamma;
    }
}
```

The `#pragma omp simd` is essential — without it the inner loop is scalar. For GPU paths, a `#pragma omp target teams distribute parallel for` replaces or wraps the simd.

**`charge` here is the species charge, not the per-particle charge.** Wait — this is subtle. For species that ionize, charge varies per particle (the `charge[]` SoA field). The pusher must read `charge[ip]` per particle in that case. Look at `PusherBoris.cpp` for the exact handling — different versions of Smilei have shifted this.

## 5. Interpolator — fields to particle

`src/Interpolator/` — interpolates `E` and `B` from the staggered Yee grid to particle positions.

Base class:

```cpp
class Interpolator {
public:
    virtual void fieldsAtParticle(ElectroMagn* EMfields, Particles& particles, ...) = 0;
};
```

Concrete: `Interpolator2D4Order.cpp`, `Interpolator3D2Order.cpp`, etc. The order corresponds to `Main.interpolation_order` (2 default, 4 high-order).

**Output**: thread-local arrays `local_E[3][n_particles]`, `local_B[3][n_particles]` indexed by particle. These are then read by the pusher.

The interpolation uses staggered field locations per Yee convention — `Ex` at cell centers along x, `Ey` along y, `Bx` at face centers, etc. The Esirkepov-style staggering preserves charge conservation.

## 6. Projector — particle to current grid

`src/Projector/` — deposits particle currents onto the grid.

Smilei uses the **Esirkepov scheme** for current deposition — it preserves charge conservation `∂ρ/∂t + ∇·J = 0` exactly on the Yee grid. The trade-off: more arithmetic per particle than a naive direct projection.

Base class:

```cpp
class Projector {
public:
    virtual void depositCurrents(double* Jx, double* Jy, double* Jz, Particles& particles, ...) = 0;
};
```

Concrete: `Projector2D4Order.cpp`, `Projector3D2Order.cpp`, etc.

**Critical invariant**: `depositCurrents` writes to per-patch current arrays. After the per-patch loop, `VectorPatch::sumDensities` performs inter-patch communication to handle the overlap at patch boundaries (ghost cells exchanged via MPI).

## 7. Particle boundary conditions

`src/ParticleBC/` — actions when a particle crosses a patch boundary.

```cpp
class ParticleBoundaryConditions {
public:
    void apply(Particles& particles, ...);
};
```

Per-direction, per-face BCs are stored in this class (set from the namelist `Species(boundary_conditions=...)`).

**Implementation per BC type:**

- `"remove"` — `eraseParticle(ip)` (swap with last + decrement size)
- `"reflective"` — flip position and momentum sign on the perpendicular axis
- `"periodic"` — wrap position (`x -= grid_length[0]` or `x += grid_length[0]`); patch-internal periodic is handled by `SmileiMPI::exchangeParticles`
- `"thermalize"` — flip momentum, then redraw from a Maxwellian at `thermal_boundary_temperature`
- `"stop"` — leave at boundary, set momentum to zero

The `partBoundCond.apply` is called per-axis after the pusher. Each axis applies independently — a particle that exits along both `x` and `y` is reflected by `y` first then `x` (or vice versa, depending on order).

**Patch-boundary handling** (different from physical BC):

When a particle moves to a neighboring patch (still inside the simulation box), it's handed off via `SmileiMPI::exchangeParticles`. The particle is removed from the source patch and added to the destination patch — physical BCs do not apply at patch-internal boundaries.

## 8. `ParticleInjector`

`src/Particles/ParticleInjector.h` — continuous injection from a boundary.

Implemented as a separate per-patch object (not a member of `Species`). At each timestep, the injector:

1. Determines the flux to inject based on `time_envelope`.
2. Computes the number of macro-particles to add this step.
3. Initializes their positions, momenta, charges, weights according to the injector's namelist parameters.
4. Appends to the target species' `Particles`.

```cpp
class ParticleInjector {
    Species* species;                          // target species (pointer to existing Species)
    int box_side;                              // injection face
    Profile* density_profile;
    Profile* mean_velocity_profile;
    // ...

    void inject(double time, ...);
};
```

The injector reuses the target species' `Species` initialization machinery — it doesn't reimplement particle creation. This means a `ParticleInjector` and its target species share `mass`, `charge`, `pusher`, etc.

## 9. Test particles

`Species(is_test = True)` in the namelist sets `species->is_test = true` on the C++ side. Test particles:

- Are pushed by the fields (so they move correctly)
- Do **not** deposit current (the projector skips them)
- Do not participate in collisions

This requires explicit checks in `Projector::depositCurrents` (skip if `is_test`), and in `Collisions::apply` (skip if either species is test). The convention is centralized in `src/Tools/SmileiPP.h` macros.

When adding a new physics module, decide whether test particles should participate. Skipping them is usually the right default — test particles are meant to sample the dynamics without back-reaction.

## 10. Particle merging

`src/Merging/` — for high-density regimes where macro-particle count grows uncontrollably (ionization producing more particles than the grid can hold).

Merging combines two macro-particles into one with the merged weight, momentum-conserving and energy-conserving (within statistical tolerance). The algorithm:

1. Bin particles in momentum space.
2. Within a bin, pair particles randomly.
3. For each pair, replace with one particle carrying summed weight.

Triggered via the namelist `Species(merging_method=..., merge_every=...)`. Adds compute cost but bounds the particle count.

Implementation: `src/Merging/Merging.cpp`. Multiple algorithms available (`"vranic_cartesian"`, `"vranic_spherical"`).

For ionization-heavy excimer runs, merging can be essential to keep memory in bounds. Without it, a 10-stage ionization ladder can multiply particle count by `2¹⁰ = 1024`.

## 11. Adding a new pusher

If you want to add a pusher (e.g., a relativistic Drift-Kinetic pusher, a guiding-center pusher, a relativistic-corrected non-relativistic pusher):

1. **Read** `src/Pusher/PusherBoris.cpp` for the canonical structure.
2. **Create** `src/Pusher/PusherMyNew.h` and `.cpp` deriving from `Pusher`.
3. **Add** the case in `src/Pusher/PusherFactory.h`:
   ```cpp
   else if (params.pusher == "mynew") {
       return new PusherMyNew(params);
   }
   ```
4. **Add** the string to the validated list in `src/Params/Params.cpp::validatePusherName(...)`.
5. **Document** in `doc/Sphinx/Source/namelist.rst`.
6. **Test** — add a `validation/tst_*.py` namelist that exercises the new pusher, and a reference output.

For GPU support, follow the pattern in `PusherBoris.cpp` for the `_GPU_KERNEL_` brackets. CPU-first development with a GPU TODO is acceptable for an initial PR.

## 12. Common pitfalls

**Modifying particle order**: many code paths assume particles are ordered as inserted. `Particles::sortParticles` reorders by cell — call sites must tolerate this. If you add a new physics module that relies on particle order, document it and consider whether `sortParticles` invalidates your invariant.

**Forgetting to update both CPU and GPU paths**: if you add a new code path to `PusherBoris::operator()`, the GPU version (under `#ifdef _GPU`) must match. CI builds may or may not exercise GPU paths.

**Per-particle vs species-level data confusion**: per-particle charge (after ionization) is in `Particles::charge[ip]`. Species-level charge (initial state) is `Species::charge`. They diverge after the first ionization event. Pushers and projectors must read `Particles::charge[ip]` for ionizing species.

**Allocation in hot loops**: don't `new` inside `Species::dynamics`. Use member arrays sized at construction and `reserve()`'d to expected maximum.

**MPI/OpenMP mismatch in physics modules**: don't add MPI calls inside the per-patch loop. The OMP-parallel section can't synchronously call MPI safely. Aggregate data per-rank, then communicate after the parallel region.

**Test-particle bypasses**: when adding a module, check whether `is_test` particles should participate. The wrong default is silent — test particles affecting current density would be a physics bug.

**Photon pushers**: photons (`mass = 0`) use `PusherPhoton`, which is ballistic. Other pushers will produce NaN for mass-zero particles. Always dispatch on `mass == 0.` in the species factory.

**`SMILEI_ASSERT` performance**: assertions only fire in debug builds. Don't use them for production-critical checks; use explicit `if (condition) ERROR(...)` for those.
