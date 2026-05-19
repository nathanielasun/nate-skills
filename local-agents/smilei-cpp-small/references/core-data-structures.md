# Core data structures (small tier)

`Particles` SoA layout, `Species`, pushers, Maxwell solvers, fields, interpolation, current deposition. Merged from frontier's two largest source-layer references.

## `Particles` — SoA layout

`src/Particles/Particles.h`. Structure-of-arrays:

```cpp
class Particles {
public:
    std::vector<double> position[3];           // [0]=x, [1]=y, [2]=z
    std::vector<double> momentum[3];           // [0]=px, [1]=py, [2]=pz
    std::vector<double> weight;
    std::vector<short>  charge;                // current charge state (changes with ionization)
    std::vector<double> chi;                   // quantum parameter (radiation reaction)
    std::vector<uint64_t> id;                  // tracking

    // Mass is per-species, not per-particle.
    size_t size() const { return weight.size(); }
    double* getPtrPosition(int d) { return position[d].data(); }
    double* getPtrMomentum(int d) { return momentum[d].data(); }

    void resize(size_t n);
    void reserve(size_t n);
    void eraseParticle(size_t i);              // O(1) swap with last
};
```

**Invariants**: all vectors same size. Particles indexed `0..size()-1`. Order is NOT stable across calls (collisions/sorting/migration reorder). `chi` allocated only with radiation reaction; `id` allocated only with tracking.

**Adding a per-particle attribute**: add `std::vector<T>` field, update `resize()`/`clear()`/`reserve()`/`eraseParticle()`, add getter, update HDF5 dump/restart in `Diagnostic/` and `Checkpoints/`. Intrusive — usually unnecessary.

## `Species` — per-species container

`src/Particles/Species.h`. One instance per species **per patch**:

```cpp
class Species {
public:
    std::string name;                          // "electron", "Ar1", etc.
    double mass;                               // in m_e
    double charge;                             // initial charge state
    short  atomic_number;

    Particles* particles;                      // SoA for this species in this patch

    Pusher* Push;
    Interpolator* Interp;
    Projector* Proj;
    Ionization* Ionize;                        // may be nullptr

    ParticleBoundaryConditions partBoundCond;

    void dynamics(double time, ...);
    void initializeParticles(...);
};
```

`Species*` stored in `Patch::vecSpecies`. A simulation with 4 patches × 3 species has 12 `Species` instances per rank.

## Per-patch dynamics — call chain

`Species::dynamics()`:

```cpp
void Species::dynamics(double time, ...) {
    if (Ionize) Ionize->apply(particles, EMfields);    // field ionization
    Interp->fieldsAtParticle(EMfields, particles);     // E, B at particle positions
    (*Push)(*particles, ...);                          // momentum + position update
    partBoundCond.apply(particles, ...);               // particle BCs
    Proj->depositCurrents(JEM, *particles);            // current deposition
}
```

**Order matters**: ionization before pusher (freed electrons see field at birth). BCs before projection (out-of-box particles don't deposit).

## Pushers

`src/Pusher/`. Base:

```cpp
class Pusher {
    virtual void operator()(Particles& particles, ...) = 0;
};
```

| File | Pusher | Use |
|---|---|---|
| `PusherBoris.cpp` | Boris (default) | Relativistic all-purpose |
| `PusherBorisNR.cpp` | Non-relativistic | `a₀ << 1`, faster |
| `PusherVay.cpp` | Vay | High-γ E×B accuracy |
| `PusherHigueraCary.cpp` | Higuera-Cary | Highest-γ accuracy |
| `PusherPhoton.cpp` | Photon (mass=0) | Ballistic at c |

Inner Boris loop:

```cpp
void PusherBoris::operator()(Particles& particles, ...) {
    double dts2 = 0.5 * dt;
    double charge_over_mass = charge / mass;

    #pragma omp simd
    for (size_t ip = ipart_min; ip < ipart_max; ++ip) {
        // half-step E push, magnetic rotation, half-step E push, position update
        // Read momentum[d][ip], local_E[d][ip], local_B[d][ip]
        // Write momentum[d][ip], position[d][ip]
    }
}
```

`#pragma omp simd` essential. GPU paths replace simd with `target teams distribute parallel for`.

**Per-particle charge for ionizing species**: read `particles.charge[ip]`, not `Species::charge` (which is initial state only). Different Smilei versions handle this slightly differently — check `PusherBoris.cpp` for the current convention.

## Interpolator and Projector

**`Interpolator`** (`src/Interpolator/`): E, B from Yee grid to particle position. Output: thread-local arrays `local_E[3][n]`, `local_B[3][n]`. Concrete: `Interpolator2D4Order.cpp` etc. Order matches `Main.interpolation_order`.

**`Projector`** (`src/Projector/`): particle currents → grid (Esirkepov scheme — preserves `∂_t ρ + ∇·J = 0` exactly). ~50 flops/particle in 2D, ~100 in 3D. Hottest loop in many simulations. `#pragma omp simd` essential.

Per-patch deposition writes to per-patch current arrays. After per-patch loop, `VectorPatch::sumDensities` does inter-patch communication for boundary overlap.

GPU: atomic adds (particles in same cell collide). Main GPU bottleneck.

## Particle boundary conditions

`src/ParticleBC/`. Per-axis, per-face from `Species(boundary_conditions=...)`:

- `"remove"` — `eraseParticle(ip)` (swap with last + decrement)
- `"reflective"` — flip position and momentum on perpendicular axis
- `"periodic"` — wrap position; patch-internal periodic via `SmileiMPI::exchangeParticles`
- `"thermalize"` — flip + redraw momentum from Maxwellian at `thermal_boundary_temperature`
- `"stop"` — leave at boundary, zero momentum

Patch-internal boundaries handle particles moving to neighbor patches via `exchangeParticles` — physical BCs do not apply at patch-internal boundaries.

## `ParticleInjector`

`src/Particles/ParticleInjector.h`. Continuous injection from boundary:

```cpp
class ParticleInjector {
    Species* species;
    int box_side;
    Profile* density_profile;
    void inject(double time, ...);
};
```

At each timestep: determines flux from `time_envelope`, computes macro-particles to add, initializes positions/momenta/charges/weights, appends to target species' `Particles`. Reuses species' init machinery — doesn't re-implement.

## Test particles and merging

**Test particles** (`is_test = True`): pushed by fields, no current deposition, no collisions. Centralized checks in `Projector::depositCurrents` (skip if `is_test`) and `Collisions::apply`.

**Merging** (`src/Merging/`): for ionization-heavy runs where particle count grows uncontrollably. Combines two macro-particles into one (conserving momentum and energy within statistics). Triggered by namelist `Species(merging_method=..., merge_every=...)`. Algorithms: `"vranic_cartesian"`, `"vranic_spherical"`.

For 10-stage ionization ladders, merging can be essential — without it particle count multiplies by `2¹⁰ = 1024`.

## Adding a new pusher

1. Read `PusherBoris.cpp` for structure.
2. Create `src/Pusher/PusherMyNew.{h,cpp}` deriving from `Pusher`.
3. Add case in `PusherFactory.h`:
   ```cpp
   else if (params.pusher == "mynew") return new PusherMyNew(params);
   ```
4. Add string to `Params.cpp::validatePusherName(...)`.
5. Document in `doc/Sphinx/Source/namelist.rst`.
6. Test in `validation/tst_*.py` with reference output.

GPU support: follow `_GPU_KERNEL_` pattern in `PusherBoris.cpp`. CPU-first with GPU TODO acceptable.

## `ElectroMagn` — field container

`src/ElectroMagn/ElectroMagn.h`. One instance per patch:

```cpp
class ElectroMagn {
    Field* Ex, *Ey, *Ez;
    Field* Bx, *By, *Bz;
    Field* Jx, *Jy, *Jz;
    Field* rho;

    std::map<std::string, Field*> rho_s;            // per-species (for diagnostics)
    std::map<std::string, Field*> Jx_s, Jy_s, Jz_s;

    std::vector<Field*> PML_fields;                 // PML aux fields (lazy)

    virtual void solveMaxwell(...) = 0;
    virtual void centerMagneticFields() = 0;
    virtual void applyEMBoundaryConditions(...) = 0;
};
```

Concrete derivations: `ElectroMagn1D`, `ElectroMagn2D`, `ElectroMagn3D`, `ElectroMagnAM`.

AM cylindrical stores `Ex_m` (complex) per azimuthal mode `m = 0..M-1`. Indexing non-trivial — read `ElectroMagnAM.cpp` carefully.

## Yee staggered grid

Fields staggered by half a cell:

- `Ex` at integer-y, half-x
- `Ey` at half-y, integer-x
- `Bz` at half-x, half-y (cell centers in 2D)
- `rho` at integer-integer

`Ex(i,j)` array index `[i + nx*j]` corresponds to physical position `(i*dx + dx/2, j*dy)`. Offset implicit in code.

## Maxwell solvers

`Main.maxwell_solver` strings:

| String | File pattern | Ghost | Use |
|---|---|---|---|
| `"Yee"` | `MaxwellSolver_Yee_*.cpp` | 2 | Default — robust, 2nd order |
| `"Lehe"` | `MaxwellSolver_Lehe_*.cpp` | 3 | Relativistic LWFA — kills numerical Cherenkov |
| `"Cowan"` | `MaxwellSolver_Cowan_*.cpp` | 3 | Reduced grid anisotropy |
| `"Bouchard"` | `MaxwellSolver_Bouchard_*.cpp` | 4–5 | Higher-order, ~4× Yee cost |

Base class:

```cpp
class MaxwellSolver {
    virtual void operator()(ElectroMagn* EMfields) = 0;
};
```

Called once per timestep per patch. CFL: `dt² × (1/dx² + 1/dy² + 1/dz²) ≤ 1` for Yee. `timestep_over_CFL = 0.95` gives margin.

Standard Yee 2D update (sketch):

```cpp
// Bz from curl E
for (i,j) Bz(i,j) -= dt_dx * (Ey(i+1,j) - Ey(i,j))
                   - dt_dy * (Ex(i,j+1) - Ex(i,j));
// Ex, Ey from curl B and J
for (i,j) Ex(i,j) += dt_dy * (Bz(i,j) - Bz(i,j-1)) - dt * Jx(i,j);
// ... Ey similarly
```

## Ghost cells

Halo of cells outside each patch interior. Refreshed after each Maxwell step via `SmileiMPI::exchangeFields()` (non-blocking + waitall). Patches at simulation boundary fill ghosts via boundary conditions, not communication.

Per-patch ghost count is checked at startup — patch smaller than ghost extent → refuses to start.

## EM boundary conditions

`src/ElectroMagnBC/`. Per geometry:

| BC | Effect |
|---|---|
| `"silver-muller"` | First-order absorbing |
| `"PML"` | Perfectly Matched Layer (~−60 dB reflection) |
| `"periodic"` | Wraps (handled at patch-exchange level) |
| `"reflective"` | Perfect conductor |

Called from `Patch::applyEMBoundaryConditions()` after Maxwell solve.

**Laser injection** is a special BC on top of Silver-Müller or PML. Injected `B` at boundary = laser profile; absorbing BC handles reflected/outgoing waves.

## PML

```python
Main(
    EM_boundary_conditions = [["PML"], ["PML"]],
    number_of_pml_cells = [[12, 12], [12, 12]],
)
```

`ElectroMagnBC2D_PML.cpp` etc. PML allocates auxiliary fields and uses stretched coordinates that grow inside the layer. Damping profile typically `σ(d) = σ_max (d/L_PML)^n`.

Use over Silver-Müller for: long-time simulations (reflections accumulate), oblique-incidence lasers, dense plasmas with overdense reflections.

## Adding a new Maxwell solver

1. Create `MaxwellSolver_MySolver_2D.{h,cpp}` etc. per geometry.
2. Derive from `MaxwellSolver`.
3. Implement `operator()(ElectroMagn*)`.
4. Register in `MaxwellSolverFactory.h`.
5. Declare ghost-cell count in `PicParams::setMaxwellSolver(...)`.
6. Update CFL in `computeTimeStep()` if it differs.
7. Document and test.

## Common pitfalls

- **Order of particles not stable**: collisions/sort/migration reorder. Don't cache indices.
- **CPU/GPU paths diverge**: when modifying one, must update the other.
- **Per-particle vs species-level charge**: post-ionization use `particles.charge[ip]`, not `Species::charge`.
- **MPI inside `omp parallel`**: forbidden without `single`/`master` guards.
- **Allocations in hot loops**: don't `new` in `Species::dynamics`. Use member arrays, pre-reserve.
- **Wrong field location**: `Ex(i,j)` is at `(i*dx + dx/2, j*dy)`, not `(i*dx, j*dy)`.
- **Forgetting ghost-cell refresh**: stale ghosts produce silent errors growing with patch-boundary length.
- **Patch smaller than ghost extent**: refuses startup. Fix: larger patches or fewer patches.
- **GPU atomic missing in current deposition**: silent races, wrong physics.
- **`Field*` resize at runtime**: forbidden. Field arrays sized at patch construction.
- **Mixing solver order with interpolation order**: `interpolation_order = 2` with `"Bouchard"` (high-order) is inconsistent — physical accuracy limited by lower order.
- **Test particles bypass**: `is_test` skips current deposition AND collisions. New modules must respect.
- **Mass-zero with non-Photon pusher**: produces NaN. Dispatch on `mass == 0` in factory.
