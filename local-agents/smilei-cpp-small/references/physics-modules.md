# Physics modules — collisions and ionization (small tier)

Pérez 2012 binary collision implementation, collisional ionization, field ionization (ADK/PPT/BSI/envelope-averaged/from_rate). Merged from frontier's two physics-module references.

## `Collisions` and `BinaryProcesses`

`src/Collisions/`. Layout:

```cpp
class Collisions {
public:
    std::vector<unsigned int> species_group1, species_group2;
    double clog;                                // 0 = auto, >0 = fixed
    int every;
    bool intra_collisions;
    bool ionizing;
    int ionization_electrons_index;
    CollisionalNuclearReaction* nuclear_reaction;

    void apply(Patch* patch, int itime);
};

class BinaryProcesses {
    std::vector<Collisions*> coll;
    void apply(Patch* patch, int itime);
};
```

`BinaryProcesses` is a thin wrapper owning a list of `Collisions` instances (one per namelist block). Per-patch invocation from `Patch::dynamics` (OMP-parallel layer).

## Pérez 2012 algorithm in code

`Collisions::apply()` pseudocode:

```cpp
void Collisions::apply(Patch* patch, int itime) {
    if (itime % every != 0) return;

    for (each_cell c in patch) {
        std::vector<size_t> idx1 = particles_in_cell(species_group1, c);
        std::vector<size_t> idx2 = particles_in_cell(species_group2, c);

        patch->rand->shuffle(idx1);
        patch->rand->shuffle(idx2);

        // Pairing
        if (intra_collisions) n_pairs = idx1.size() / 2;     // pair idx1[2k] with idx1[2k+1]
        else                  n_pairs = min(idx1.size(), idx2.size());

        for (each pair p) {
            double scattering_angle = sample_perez_angle(...);
            update_momenta_post_collision(p, scattering_angle);
        }
    }
}
```

`Collisions.cpp::apply()` is ~500 lines because of geometry variants, intra/inter detection, auto/fixed Coulomb log, optional ionization, optional nuclear reactions. Find the section relevant to your work and treat others as orthogonal.

Energy and momentum are conserved pair-by-pair (scattering angle sampled in COM frame, momenta rotated preserving magnitudes).

## Intra vs inter detection

```cpp
intra_collisions = false;
for (int s1 : species_group1)
    for (int s2 : species_group2)
        if (s1 == s2) { intra_collisions = true; break; }
```

Pure code-path optimization — physics is identical.

## Weight handling

Macro-particles can have different weights. Pérez algorithm adjusts collision probability:

```cpp
double w1 = particles1->weight[i1];
double w2 = particles2->weight[i2];
double prob = base_prob * std::min(w1, w2) / std::max(w1, w2);
if (patch->rand->uniform() < prob) {
    // Apply collision
}
```

Preserves the underlying physical collision frequency.

## Coulomb logarithm

Auto mode (`clog = 0`): formula uses local plasma parameters. **Requires `reference_angular_frequency_SI`**:

```cpp
double computeCoulombLog(double T_norm, double n_norm, double mass1, ...) {
    double n_SI = n_norm * N_r_SI(omega_r_SI);
    double T_SI = T_norm * 511e3 * 1.602e-19;          // Joules
    double lambda_D = std::sqrt(eps0 * T_SI / (n_SI * e2));
    double b_min = compute_min_impact_parameter(...);
    return std::log(lambda_D / b_min);
}
```

NRL formulary baseline + Lee-More for dense regimes + relativistic/quantum-classical handling.

Fixed mode (`clog > 0`): returns `clog` directly.

## CollisionalIonization

When `ionizing` is set to an electron species name, extra step after standard collision:

```cpp
if (ionizing && particles2_is_ion) {
    double E_inc = compute_incident_kinetic_energy(...);
    double U_i = lookup_ionization_potential(atomic_number, current_charge_state);

    if (E_inc > U_i) {
        double sigma = lotz_drawin_cross_section(E_inc, U_i, ...);
        double prob_ionize = sigma * v_rel * dt;

        if (patch->rand->uniform() < prob_ionize) {
            particles2->charge[i2] += 1;
            // Create electron at ion position
            int new_idx = electron_species->particles->size();
            electron_species->particles->resize(new_idx + 1);
            electron_species->particles->position[d][new_idx] = particles2->position[d][i2];
            // Initialize momentum from energy budget
        }
    }
}
```

`ionization_electrons_index` resolved at construction. Lookup: `electron_species = patch->vecSpecies[ionization_electrons_index]`.

## Nuclear reactions

When `nuclear_reaction` is set, per-pair:

```cpp
if (nuclear_reaction) {
    double E_com = compute_COM_energy(...);
    double sigma = nuclear_reaction->crossSection(E_com);
    double prob = sigma * v_rel * dt;
    if (patch->rand->uniform() < prob) {
        particles1->eraseParticle(i1);
        particles2->eraseParticle(i2);
        nuclear_reaction->createProducts(patch, COM_pos, COM_p, ...);
    }
}
```

Cross-sections from ENDF tables. Supported: D-D, D-T, D-He³, p-B fusion, similar standard channels.

## Performance considerations for collisions

- O(N) per cell when pairing is random, O(cost × N) practically. Cost per cell linear in particle count.
- For dense plasma with PPC = 128, collisions can be 30–50% of step time.
- **Always enable `LoadBalancing`** for collisional runs.
- `every > 1` for slow channels (where physically justified).
- Avoid multi-species groups when not needed — `species1=["Ar1","Ar2","Ar3"], species2=["Ar1","Ar2","Ar3"]` is 9 pair combinations.
- GPU paths partial — CPU is current performance target.

## `Ionization` class hierarchy

`src/Ionization/`. Base:

```cpp
class Ionization {
    int atomic_number;
    int max_charge_state;
    Species* ionization_electrons;
    virtual void apply(Particles* particles, ElectroMagn* EMfields, ...) = 0;
};
```

| File | Model | `ionization_model` string |
|---|---|---|
| `IonizationTunnel.cpp` | ADK | `"tunnel"` |
| `IonizationTunnelFullPPT.cpp` | Full PPT | `"tunnel_full_PPT"` |
| `IonizationTunnelBSI.cpp` | ADK + BSI | `"tunnel_BSI"` |
| `IonizationTunnelEnvelopeAveraged.cpp` | DC-ADK cycle-averaged | `"tunnel_envelope_averaged"` |
| `IonizationFromRate.cpp` | User callable | `"from_rate"` |

One `Ionization*` per `Species` with `ionization_model != None`. Called from `Species::dynamics()` before pusher.

## `IonizationTunnel` (ADK)

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

ADK formula: Ammosov-Delone-Krainov 1986. Per-patch RNG (`patch->rand`). Energy `U_i` debited to `DiagScalar.Uexp` accumulator (bookkeeping, not field back-reaction per event).

## `IonizationTunnelFullPPT`

Same structure as ADK, replacing `adk_formula` with `ppt_formula(E, U_i, n_star, omega_r)`. PPT interpolates tunnel ↔ multiphoton without adiabatic approximation. For `γ_K → 0` reduces to ADK; for `γ_K → ∞` gives multiphoton power-law.

Cost ~2× ADK per particle. **Recommended for excimer regimes** with `γ_K ~ 1`.

## `IonizationTunnelBSI`

ADK rate capped at BSI threshold:

```cpp
double rate_adk = adk_formula(E_mag, U_i, n_star);
double E_bsi = U_i * U_i / (4 * z);                    // atomic units
if (E_mag > E_bsi) rate = geometric_bsi_rate(...);
else               rate = rate_adk;
```

Smooth interpolation in `IonizationTunnelBSI.cpp`. Use at intensities where field strength approaches `E_bsi` — for `H` (U_i = 13.6 eV), `I ≈ 2.6 × 10¹⁵ W/cm²` and up.

## `IonizationTunnelEnvelopeAveraged`

For envelope-laser simulations. Uses cycle-averaged rate `Γ_AC` instead of instantaneous. Linear vs circular polarization handled differently:
- Linear: `Γ_AC = (3/π) × sqrt(E/E_atomic) × Γ_DC` approximate; exact prefactor computed in `IonizationTunnelEnvelopeAveraged.cpp`.
- Circular: `Γ_AC = Γ_DC` (constant amplitude over cycle).

Freed-electron transverse momentum: initialized to zero by default (cycle-quiver momentum lost in averaging).

## `IonizationFromRate` — Python callable

```cpp
class IonizationFromRate : public Ionization {
    PyObject* rate_function;

    void apply(Particles* particles, ElectroMagn* EMfields, ...) override {
        for (each particle ip) {
            double E = compute_field_magnitude(...);
            PyObject* py_E = PyFloat_FromDouble(E);
            PyObject* result = PyObject_CallOneArg(rate_function, py_E);
            double rate = PyFloat_AsDouble(result);
            Py_DECREF(py_E); Py_DECREF(result);

            double prob = 1.0 - std::exp(-rate * dt);
            if (rand_uniform() < prob) increment_charge_and_create_electron(...);
        }
    }
};
```

**Slow**: millions of Python calls per second per rank. For dense plasma, 10× slowdown. **GPU incompatible** (Python doesn't run on GPU); factory falls back to CPU.

Path for excimer photoionization (rate from wavelength-specific cross-section, not ADK closed form).

## Ionization potential tables

`src/Ionization/IonizationTables.cpp`:

```cpp
double IonizationTables::ionization_potential[Z_MAX+1][Z_MAX+1];
//                                              Z         charge_state
```

`U_i` in eV. NIST atomic-spectra-database sourced. `[18][2]` = 40.74 eV (Ar²⁺ → Ar³⁺), not cumulative.

Same table used by `IonizationTunnel`, `IonizationTunnelFullPPT`, `IonizationTunnelBSI`, `CollisionalIonization`.

## Adding a new ionization model

1. Create `src/Ionization/IonizationMyNew.{h,cpp}` deriving from `Ionization`.
2. Implement `apply(Particles*, ElectroMagn*, ...)`:
   - Loop over particles, check `max_charge_state`.
   - Interpolate field magnitude.
   - Lookup `U_i`, compute rate with new formula.
   - MC step: `1 - exp(-rate × dt)`, draw uniform, ionize if true.
3. Register in `IonizationFactory.h`.
4. Add string to `Params.cpp::validateIonizationModel(...)`.
5. Document and test in `validation/`.

For GPU support, `#ifdef _GPU` path with same logic.

## Interaction with collisions

**Field and collisional ionization are independent additive channels.** Both increment `particles->charge[ip]` and create electrons; rates additive at low rates; no double-counting.

Order in per-patch loop:
```cpp
Species::dynamics() {
    Ionize->apply(particles, EMfields);    // field ionization
    Interp->fieldsAtParticle(...);
    Push->operator()(...);
    partBoundCond.apply(...);
    Proj->depositCurrents(...);
}
// later:
BinaryProcesses::apply(patch, itime);      // collisional ionization (inside Collisions::apply)
```

Both modify same `particles->charge[ip]` field but at different points.

## Common pitfalls

- **Missing `reference_angular_frequency_SI` with auto Coulomb log**: silent NaN or infinity. Hard rule 1.
- **`ionizing = True` (Boolean) vs `ionizing = "electron"` (species name string)**: current Smilei expects species name. Legacy True/False may not work.
- **Missing ion species in ladder**: must define `Ar1`, `Ar2` etc. explicitly. Smilei doesn't auto-create.
- **Wrong table lookup**: `ionization_potential[Z][z]` is energy from `z` → `z+1`, NOT cumulative. Off-by-one common.
- **Forgetting `max_charge_state` check**: rate formula evaluates with invalid `U_i` (out-of-table). Always guard.
- **Calling Python from GPU kernel**: doesn't work. `IonizationFromRate` is CPU-only.
- **Per-event field interpolation cost**: re-interpolating from scratch is comparable to pusher cost. Currently not shared with pusher's interpolation.
- **Forgetting energy bookkeeping**: `U_i` should be debited to `DiagScalar.Uexp`. Otherwise `Ubal` drifts.
- **Race on freed-electron species**: `resize(n+1)` from multiple threads in a patch would race. Don't parallelize within `Ionization::apply` — patch-level OMP is the only safe layer.
- **`from_rate` callable returning SI units instead of normalized**: user error, no auto-check. Symptom: rate off by `omega_r ≈ 10¹⁵`.
- **Is_test bypass**: test particles skip ionization in current code. New models should follow convention.
- **`coulomb_log_factor` left at non-1 from sensitivity study**: silent wrong physics in production.
- **Low PPC with collisions**: PPC = 8 with collisions → 35% per-cell variance, oscillatory transfer. Use ≥ 32.
- **Only e-i collisions**: in dense plasma, e-e drives electron thermalization. Without it, electrons stay non-Maxwellian.
- **Modifying `Collisions::apply` without GPU sync**: GPU paths partial; updates may diverge.
- **Storing state on `Collisions` shared across patches**: one `Collisions` per namelist block, shared. Use `Patch` for per-patch state.
