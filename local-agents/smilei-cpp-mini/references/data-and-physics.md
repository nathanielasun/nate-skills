# Data structures and physics modules (mini tier)

`Particles` SoA, `Species`, pushers, Maxwell solvers, collisions, ionization — all source-layer essentials merged for the mini tier.

## `Particles` — SoA layout

`src/Particles/Particles.h`:

```cpp
class Particles {
    std::vector<double> position[3];   // x,y,z arrays
    std::vector<double> momentum[3];   // px,py,pz arrays
    std::vector<double> weight;
    std::vector<short>  charge;        // per-particle, changes with ionization
    std::vector<double> chi;           // radiation reaction
    std::vector<uint64_t> id;          // tracking

    size_t size() const { return weight.size(); }
    void resize(size_t n);
    void eraseParticle(size_t i);      // O(1) swap with last
};
```

All vectors same size. Order NOT stable across calls. `chi` allocated only with radiation reaction; `id` only with tracking.

## `Species` — per-patch container

`src/Particles/Species.h`. One instance per species **per patch**:

```cpp
class Species {
    std::string name;
    double mass, charge;               // species-level (initial)
    short atomic_number;
    Particles* particles;              // SoA in this patch
    Pusher* Push;
    Interpolator* Interp;
    Projector* Proj;
    Ionization* Ionize;                // may be nullptr
    ParticleBoundaryConditions partBoundCond;

    void dynamics(double time, ...);
};
```

Stored in `Patch::vecSpecies`. 4 patches × 3 species → 12 `Species` instances per rank.

## Per-patch dynamics call chain

```cpp
void Species::dynamics(double time, ...) {
    if (Ionize) Ionize->apply(particles, EMfields);    // field ionization
    Interp->fieldsAtParticle(EMfields, particles);     // E, B at particle
    (*Push)(*particles, ...);                          // momentum + position
    partBoundCond.apply(particles, ...);               // particle BCs
    Proj->depositCurrents(JEM, *particles);            // current deposition
}
```

Order matters: ionization before pusher (freed electrons see field at birth); BCs before projection (out-of-box doesn't deposit).

## Pushers

`src/Pusher/`. Base: `virtual void operator()(Particles&, ...) = 0;`

| File | Pusher | Use |
|---|---|---|
| `PusherBoris.cpp` | Boris | Default, relativistic |
| `PusherBorisNR.cpp` | Non-rel Boris | `a₀ << 1`, faster |
| `PusherVay.cpp` | Vay | High-γ E×B |
| `PusherHigueraCary.cpp` | Higuera-Cary | Highest-γ accuracy |
| `PusherPhoton.cpp` | Photon | mass=0, ballistic |

Inner loop has `#pragma omp simd` (CPU) or `#pragma omp target teams distribute parallel for` (GPU under `#ifdef _GPU`).

**Per-particle charge** (post-ionization): read `particles.charge[ip]`, not `Species::charge`.

## Interpolator / Projector

**Interpolator** (`src/Interpolator/`): E, B from Yee grid → particle position. Output thread-local arrays `local_E[3][n]`, `local_B[3][n]`. Order matches `Main.interpolation_order`.

**Projector** (`src/Projector/`): particles → current grid. Esirkepov scheme — preserves `∂_t ρ + ∇·J = 0` exactly. Hottest loop in many sims (~50 flops/particle 2D, ~100 3D). `#pragma omp simd` essential. GPU: atomic adds (race when particles in same cell).

## Particle BCs

`src/ParticleBC/`. Per-axis, per-face:
- `"remove"` — `eraseParticle(ip)`
- `"reflective"` — flip position + momentum sign perpendicular
- `"periodic"` — wrap; patch-internal handled by `SmileiMPI::exchangeParticles`
- `"thermalize"` — flip + redraw from Maxwellian
- `"stop"` — leave at boundary, zero momentum

## Adding a new pusher

1. Read `PusherBoris.cpp`.
2. Create `src/Pusher/PusherMyNew.{h,cpp}` deriving from `Pusher`.
3. Register in `PusherFactory.h`.
4. Validate string in `Params.cpp::validatePusherName(...)`.
5. Document, test in `validation/`.

## `ElectroMagn` — field container

`src/ElectroMagn/`. One per patch:

```cpp
class ElectroMagn {
    Field* Ex, *Ey, *Ez, *Bx, *By, *Bz;
    Field* Jx, *Jy, *Jz, *rho;
    std::map<std::string, Field*> rho_s, Jx_s, Jy_s, Jz_s;  // per-species
    std::vector<Field*> PML_fields;
    virtual void solveMaxwell(...) = 0;
};
```

Geometry variants: `ElectroMagn1D`, `2D`, `3D`, `AM`. AM stores `Ex_m` (complex) per azimuthal mode.

## Yee staggered grid

- `Ex` at integer-y, half-x
- `Ey` at half-y, integer-x
- `Bz` at half-x, half-y (cell centers in 2D)
- `rho` at integer-integer

`Ex(i,j)` physical position is `(i*dx + dx/2, j*dy)` — offset implicit in code. Bites when implementing physics needing absolute positions.

## Maxwell solvers

| Solver | File pattern | Ghost cells | CFL/cost | Use |
|---|---|---|---|---|
| Yee | `MaxwellSolver_Yee_*.cpp` | 2 | Standard | Default |
| Lehe | `MaxwellSolver_Lehe_*.cpp` | 3 | ~1.3× Yee | Relativistic LWFA (kills num Cherenkov) |
| Cowan | `MaxwellSolver_Cowan_*.cpp` | 3 | ~2× Yee | Reduced grid anisotropy |
| Bouchard | `MaxwellSolver_Bouchard_*.cpp` | 4–5 | ~4× Yee | Higher-order spatial |

CFL Yee: `dt² × (1/dx² + 1/dy² + 1/dz²) ≤ 1`. `timestep_over_CFL = 0.95` for margin.

## EM boundary conditions

`src/ElectroMagnBC/`. Per geometry:
- `"silver-muller"` — first-order absorbing
- `"PML"` — Perfectly Matched Layer (~−60 dB reflection)
- `"periodic"` — handled at patch-exchange
- `"reflective"` — perfect conductor

Laser injection = special BC on top of Silver-Müller or PML.

PML setup:
```python
Main(EM_boundary_conditions = [["PML"], ["PML"]],
     number_of_pml_cells = [[12, 12], [12, 12]])
```

12 cells typical. PML allocates aux fields, uses stretched coordinates.

## Adding a new Maxwell solver

1. Create `MaxwellSolver_MySolver_2D.{h,cpp}` per geometry.
2. Derive from `MaxwellSolver`.
3. Register in `MaxwellSolverFactory.h`.
4. Declare ghost count in `PicParams::setMaxwellSolver(...)`.
5. Update CFL in `computeTimeStep()` if differs.

## `Collisions` / `BinaryProcesses`

`src/Collisions/`:

```cpp
class Collisions {
    std::vector<unsigned int> species_group1, species_group2;
    double clog;                           // 0=auto, >0=fixed
    int every;
    bool intra_collisions, ionizing;
    int ionization_electrons_index;
    CollisionalNuclearReaction* nuclear_reaction;
    void apply(Patch* patch, int itime);
};

class BinaryProcesses {
    std::vector<Collisions*> coll;
    void apply(Patch* patch, int itime);
};
```

One `Collisions` per namelist block. `BinaryProcesses` dispatches. Called from `Patch::dynamics` inside OMP-parallel loop — each patch processes its particles independently.

## Pérez 2012 algorithm

`Collisions::apply()` per cell:

```cpp
for (each cell c in patch) {
    std::vector<size_t> idx1 = particles_in_cell(group1, c);
    std::vector<size_t> idx2 = particles_in_cell(group2, c);
    patch->rand->shuffle(idx1);
    patch->rand->shuffle(idx2);

    if (intra_collisions) n_pairs = idx1.size() / 2;          // pair idx1[2k] with idx1[2k+1]
    else                  n_pairs = min(idx1.size(), idx2.size());

    for (each pair p) {
        double angle = sample_perez_angle(...);
        update_momenta_post_collision(p, angle);
    }
}
```

Scattering angle in COM frame; momenta rotated preserving magnitudes → exact pair-by-pair conservation.

`Collisions.cpp::apply()` is ~500 lines (geometry variants, intra/inter, auto/fixed log, optional ionization, optional nuclear reactions).

## Weight handling

```cpp
double prob = base_prob * std::min(w1, w2) / std::max(w1, w2);
if (patch->rand->uniform() < prob) { apply_collision(); }
```

Preserves physical collision frequency under mixed weights.

## Coulomb log

Auto mode (`clog = 0`): requires `reference_angular_frequency_SI`. Computes from local T, n, masses, charges via NRL formulary + Lee-More for dense regimes + quantum-classical handling.

Fixed: returns `clog` directly.

## CollisionalIonization

When `ionizing = "electron"` (species name), per collision:

```cpp
if (ionizing && particles2_is_ion) {
    double E_inc = compute_incident_KE(...);
    double U_i = lookup_ionization_potential(atomic_number, current_charge_state);
    if (E_inc > U_i) {
        double sigma = lotz_drawin_cross_section(E_inc, U_i, ...);
        double prob = sigma * v_rel * dt;
        if (patch->rand->uniform() < prob) {
            particles2->charge[i2] += 1;
            int new_idx = electron_species->particles->size();
            electron_species->particles->resize(new_idx + 1);
            electron_species->particles->position[d][new_idx] = particles2->position[d][i2];
        }
    }
}
```

Lotz-Drawin cross-sections. Same `U_i` table as field ionization.

## Nuclear reactions

When `nuclear_reaction` is set:

```cpp
double E_com = compute_COM_energy(...);
double sigma = nuclear_reaction->crossSection(E_com);
double prob = sigma * v_rel * dt;
if (patch->rand->uniform() < prob) {
    particles1->eraseParticle(i1);
    particles2->eraseParticle(i2);
    nuclear_reaction->createProducts(patch, COM_pos, COM_p, ...);
}
```

ENDF tables. D-D, D-T, D-He³, p-B supported.

## `Ionization` class hierarchy

`src/Ionization/`:

| File | Model | String |
|---|---|---|
| `IonizationTunnel.cpp` | ADK | `"tunnel"` |
| `IonizationTunnelFullPPT.cpp` | Full PPT | `"tunnel_full_PPT"` |
| `IonizationTunnelBSI.cpp` | ADK + BSI | `"tunnel_BSI"` |
| `IonizationTunnelEnvelopeAveraged.cpp` | Cycle-averaged | `"tunnel_envelope_averaged"` |
| `IonizationFromRate.cpp` | User callable | `"from_rate"` |

```cpp
class Ionization {
    int atomic_number, max_charge_state;
    Species* ionization_electrons;
    virtual void apply(Particles*, ElectroMagn*, ...) = 0;
};
```

Called from `Species::dynamics()` before pusher.

## ADK implementation

```cpp
void IonizationTunnel::apply(Particles* particles, ElectroMagn* EMfields, ...) {
    for (size_t ip = 0; ip < particles->size(); ++ip) {
        short z = particles->charge[ip];
        if (z >= max_charge_state) continue;

        double E_mag = interp_field_magnitude(EMfields, particles->position, ip);
        double U_i = ionization_potential[atomic_number][z];
        double n_star = z / std::sqrt(2 * U_i / U_H);
        double rate = adk_formula(E_mag, U_i, n_star);

        double prob = 1.0 - std::exp(-rate * dt);
        if (patch->rand->uniform() < prob) {
            particles->charge[ip] = z + 1;
            create_electron(particles, ip, ionization_electrons, ...);
        }
    }
}
```

Same MC structure for PPT (replace `adk_formula` with `ppt_formula`), BSI (cap above threshold), envelope-averaged (cycle-averaged rate).

**For excimer regimes**, prefer `"tunnel_full_PPT"` — `γ_K ~ 1` makes ADK inaccurate.

## Ionization potential tables

`src/Ionization/IonizationTables.cpp`:

```cpp
double IonizationTables::ionization_potential[Z_MAX+1][Z_MAX+1];
//                                              Z         charge_state
```

`[18][2]` = 40.74 eV (Ar²⁺ → Ar³⁺), NOT cumulative. Shared between field and collisional ionization.

## `IonizationFromRate` — Python callable

```cpp
PyObject* result = PyObject_CallOneArg(rate_function, py_E);
double rate = PyFloat_AsDouble(result);
```

Slow — Python in hot loop. **GPU incompatible.** Used for excimer photoionization or tabulated rates not fit by ADK form.

## Adding a new ionization model

1. Create `IonizationMyNew.{h,cpp}` deriving from `Ionization`.
2. Implement `apply(Particles*, ElectroMagn*, ...)` with new rate formula.
3. Register in `IonizationFactory.h`.
4. Add string to `Params.cpp::validateIonizationModel(...)`.
5. Document, test in `validation/`.

## Field vs collisional ionization

**Independent additive channels.** Both increment `particles->charge[ip]` and create electrons. Order:

```cpp
Species::dynamics() {
    Ionize->apply(...);                    // field ionization
    Interp/Push/BCs/Proj ...
}
BinaryProcesses::apply(...);                // collisional ionization (in Collisions::apply)
```

Different points in timestep but same charge field.

## Common pitfalls

- **Order of particles** not stable across calls. Don't cache indices.
- **CPU/GPU paths diverge** — when modifying one, update the other.
- **Per-particle vs species charge**: post-ionization use `particles.charge[ip]`.
- **MPI inside `omp parallel`** without `single`/`master` guards.
- **Allocations in hot loops** — pre-reserve, don't `new` in `dynamics`.
- **Wrong field location**: `Ex(i,j)` at `(i*dx + dx/2, j*dy)`.
- **GPU atomic missing in current deposition** → silent race.
- **Patch smaller than ghost extent** → refuses startup.
- **Missing `reference_angular_frequency_SI` with auto Coulomb log** → NaN.
- **`ionizing = True`** (Boolean) vs **`ionizing = "electron"`** (species name string). Use name.
- **Missing `max_charge_state` check** in ionization → table out-of-bounds.
- **`from_rate` callable returning SI units** instead of normalized → off by `omega_r ≈ 10¹⁵`.
- **Calling Python from GPU kernel** — doesn't work. `IonizationFromRate` is CPU-only.
- **Race on freed-electron `resize(n+1)`** — patch-level OMP is the only safe layer.
- **Test particles** (`is_test`) skip current deposition AND collisions. New modules must respect.
- **Low PPC with collisions** (< 32) — noisy rates, oscillatory transfer.
- **Only e-i collisions** in dense plasma — electrons stay non-Maxwellian (missing e-e).
- **`Field*` resize at runtime** — forbidden.
- **Forgetting energy bookkeeping for new ionization model** — `Ubal` drifts in `DiagScalar`.
