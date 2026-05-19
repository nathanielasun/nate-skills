# Ionization — source layer

Field-ionization implementation in Smilei: `Ionization` base class, ADK/PPT/BSI/envelope-averaged code paths, the Monte-Carlo machinery, the ionization-potential tables, and how to add a new model (including the `from_rate` Python-callable hook).

## Table of contents

1. `Ionization` class hierarchy
2. `IonizationTunnel` (ADK) — the canonical case
3. `IonizationTunnelFullPPT`
4. `IonizationTunnelBSI`
5. `IonizationTunnelEnvelopeAveraged`
6. `IonizationFromRate` — Python-callable hook
7. The Monte-Carlo step
8. Ionization potential tables
9. Adding a new ionization model
10. Interaction with collisions and the freed-electron species
11. Common pitfalls

---

## 1. `Ionization` class hierarchy

`src/Ionization/` contains the field-ionization module. Base class:

```cpp
class Ionization {
public:
    int atomic_number;
    int max_charge_state;
    Species* ionization_electrons;             // species receiving freed electrons

    virtual void apply(Particles* particles, ElectroMagn* EMfields, ...) = 0;
};
```

Concrete implementations:

| File | Model | Namelist `ionization_model` |
|---|---|---|
| `IonizationTunnel.cpp` | ADK | `"tunnel"` |
| `IonizationTunnelFullPPT.cpp` | Full PPT | `"tunnel_full_PPT"` |
| `IonizationTunnelBSI.cpp` | ADK + BSI | `"tunnel_BSI"` |
| `IonizationTunnelEnvelopeAveraged.cpp` | DC-ADK cycle-averaged | `"tunnel_envelope_averaged"` |
| `IonizationFromRate.cpp` | User-supplied callable | `"from_rate"` |

An `Ionization*` instance lives on each `Species` that has `ionization_model != None`. Allocated in `Species::Species()` from the factory `IonizationFactory::create(ionization_model, params, species)`.

**Call site**: `Species::dynamics()` calls `Ionize->apply(particles, EMfields, ...)` at the start of the per-step physics, before the pusher. This is so freed electrons see the field they would have at their birth time.

## 2. `IonizationTunnel` (ADK) — the canonical case

`IonizationTunnel.cpp::apply()` (sketch):

```cpp
void IonizationTunnel::apply(Particles* particles, ElectroMagn* EMfields, ...) {
    for (size_t ip = 0; ip < particles->size(); ++ip) {
        short z = particles->charge[ip];
        if (z >= max_charge_state) continue;        // ladder cap

        // 1. Interpolate local field at particle position
        double Ex = interp_Ex(EMfields, particles->position[0][ip], ...);
        double Ey = interp_Ey(EMfields, ...);
        double Ez = interp_Ez(EMfields, ...);
        double E_mag = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);

        // 2. ADK rate for state z → z+1
        double U_i = ionization_potential[atomic_number][z];   // table lookup
        double n_star = z / std::sqrt(2 * U_i / U_H);          // effective quantum number
        double rate = adk_formula(E_mag, U_i, n_star);

        // 3. Monte-Carlo sample
        double prob_ionize = 1.0 - std::exp(-rate * dt);
        if (patch->rand->uniform() < prob_ionize) {
            // 4. Increment charge, create electron
            particles->charge[ip] = z + 1;
            create_electron(particles, ip, ionization_electrons, ...);
        }
    }
}
```

The `adk_formula` is the standard ADK expression for the instantaneous tunnel rate as a function of `E`, `U_i`, and effective quantum numbers. Look at `IonizationTunnel.cpp` for the exact form (multiple equivalent expressions exist in the literature; Smilei uses Ammosov-Delone-Krainov 1986).

**Random-number stream**: per-patch RNG (`patch->rand`) — reproducibility with `Main.random_seed`.

**Energy bookkeeping**: when an electron is freed, the energy `U_i` is taken from the EM field (conservation). The implementation does this by recording the absorbed energy in a diagnostic accumulator (`Uexp` in `DiagScalar`). The field itself is not modified per-event — the absorbed energy is a bookkeeping quantity, not a back-reaction term.

## 3. `IonizationTunnelFullPPT`

`IonizationTunnelFullPPT.cpp::apply()` follows the same structure as `IonizationTunnel`, but replaces `adk_formula` with the Perelomov-Popov-Terent'ev expression:

```cpp
double rate_ppt = ppt_formula(E_mag, U_i, n_star, omega_r);
```

PPT interpolates between the tunnel and multiphoton limits without the adiabatic approximation. For `γ_K → 0` it reduces to ADK; for `γ_K → ∞` it gives the multiphoton power-law `~ (E)^(2N)` with `N` being the photon-count threshold.

**Cost**: ~2× ADK per particle, due to the slightly more complex special-function evaluation. Still negligible compared to the pusher/projector.

**Recommended** for excimer regimes with `γ_K ~ 1`. Smilei users coming from UHI work tend to default to `"tunnel"`; when this skill detects a UV regime, it should suggest `"tunnel_full_PPT"`.

## 4. `IonizationTunnelBSI`

`IonizationTunnelBSI.cpp::apply()`: same as ADK, with a cap that kicks in when the field exceeds the Barrier-Suppression Ionization threshold:

```cpp
double rate_adk = adk_formula(E_mag, U_i, n_star);

// BSI threshold field
double E_bsi = (U_i)^2 / (4 * z);                  // in atomic units

if (E_mag > E_bsi) {
    // Above-barrier — rate is geometric, not exponential
    rate = geometric_bsi_rate(E_mag, U_i, ...);
} else {
    rate = rate_adk;
}
```

The smooth interpolation between ADK and BSI is handled by the formula in `IonizationTunnelBSI.cpp` — look there for the exact form. Krainov-Ristic and other variants exist; Smilei's choice is documented in `doc/Sphinx/Source/highlights/Ionization.rst`.

**When to use**: at intensities where the field strength approaches `E_bsi` for the ions of interest. For `H` (U_i = 13.6 eV), `E_bsi ≈ 1.4 × 10¹⁰ V/m`, corresponding to `I ≈ 2.6 × 10¹⁵ W/cm²`. For higher-Z ions, `E_bsi` is larger.

## 5. `IonizationTunnelEnvelopeAveraged`

For envelope-laser simulations. The instantaneous field amplitude doesn't exist (the carrier oscillation is averaged out); we instead have the envelope `|E|_envelope`. The cycle-averaged ADK rate `Γ_ADK,AC` is the relevant quantity.

```cpp
double rate_ac = envelope_averaged_adk(E_envelope_mag, U_i, n_star, polarization_ellipticity);
```

The formula differs between linear and circular polarization:

- **Linear**: `Γ_AC = (3 / π) × sqrt(E / E_atomic) × Γ_DC` (approximate)
- **Circular**: `Γ_AC = Γ_DC` (constant amplitude over the cycle)

Smilei computes the exact prefactor using the cycle integral. See `IonizationTunnelEnvelopeAveraged.cpp` for the implementation.

**Freed-electron momentum initialization**: with envelope models, the carrier-cycle quiver momentum is lost in averaging. The freed electron's transverse momentum is initialized to zero by default. This is approximate — for physics depending on the freed-electron initial drift, run a full-EM benchmark.

## 6. `IonizationFromRate` — Python-callable hook

For arbitrary user-specified rates (photoionization, tabulated cross-sections, custom multiphoton):

```python
def my_rate(E):
    return some_rate_in_normalized_units

Species(
    ionization_model = "from_rate",
    ionization_rate = my_rate,
    ionization_electrons = "electron",
    ...
)
```

The Python callable is preserved as a `PyObject*` on the `IonizationFromRate` instance:

```cpp
class IonizationFromRate : public Ionization {
    PyObject* rate_function;                     // Python callable

    void apply(Particles* particles, ElectroMagn* EMfields, ...) override {
        for (each particle ip) {
            double E = compute_field_magnitude(...);
            // Call the Python function
            PyObject* py_E = PyFloat_FromDouble(E);
            PyObject* result = PyObject_CallOneArg(rate_function, py_E);
            double rate = PyFloat_AsDouble(result);
            Py_DECREF(py_E);
            Py_DECREF(result);

            // Same MC step as ADK
            double prob = 1.0 - std::exp(-rate * dt);
            if (rand_uniform() < prob) increment_charge_and_create_electron(...);
        }
    }
};
```

**Performance**: calling Python from a hot loop is *slow* — millions of Python calls per second per rank. For dense plasmas, this can be a 10× slowdown.

**Optimization**: when the user provides `ionization_rate`, Smilei tries to detect cases where the callable is a closed-form expression (e.g., `lambda E: c * E**6`) and JIT-compile it. The detection is conservative — only simple closures work.

**GPU**: incompatible with `from_rate` because Python can't run on GPU. The factory falls back to CPU even on GPU builds.

For excimer photoionization, this is the natural path — the rate is given by a wavelength-specific cross-section curve, not by an ADK-style closed form.

## 7. The Monte-Carlo step

All ionization models share the same Monte-Carlo step:

```cpp
double prob = 1.0 - std::exp(-rate * dt);

double u = patch->rand->uniform();

if (u < prob) {
    // Ionize
    particles->charge[ip] += 1;
    create_freed_electron(...);
}
```

The `1.0 - exp(-rate × dt)` form is mathematically equivalent to per-time-step probability `rate × dt` for small `dt × rate`, but is robust for fast rates (where `rate × dt > 1` would overflow).

**For very fast rates** (`rate × dt > 0.1`), the per-step model underestimates ionization (you can't ionize twice in one step). The remedy is to take smaller `dt` — if your ionization timescale is comparable to `dt`, your overall timestep is too coarse.

**Reproducibility**: `Main.random_seed` controls all randomness, including ionization. Two simulations with the same seed produce identical ionization histories.

## 8. Ionization potential tables

`src/Ionization/IonizationTables.cpp` holds the ionization-potential lookup:

```cpp
double IonizationTables::ionization_potential[ATOMIC_NUMBER_MAX+1][ATOMIC_NUMBER_MAX+1];
//                                                Z              charge_state
```

For each element `Z = 1..103` and each charge state `0..Z-1`, the table holds `U_i` in eV. Values come from NIST atomic-spectra-database and are documented in `doc/Sphinx/Source/highlights/Ionization.rst`.

**Looking up `U_i(Z=18, charge_state=2)` returns the energy to ionize from Ar²⁺ to Ar³⁺**, i.e., 40.74 eV. Not the cumulative energy.

**Updating the table**: rarely necessary. If new spectroscopic measurements update a value, edit `IonizationTables.cpp`, rebuild, and re-run validation to ensure no regressions.

## 9. Adding a new ionization model

To add a model (e.g., a more accurate intermediate-`γ_K` formula, a tabulated experimental rate, an ML-derived rate):

1. **Create** `src/Ionization/IonizationMyNew.h` and `.cpp` deriving from `Ionization`.

2. **Implement** `apply(Particles*, ElectroMagn*, ...)`:
   ```cpp
   void IonizationMyNew::apply(Particles* particles, ElectroMagn* EMfields, ...) {
       for (size_t ip = 0; ip < particles->size(); ++ip) {
           short z = particles->charge[ip];
           if (z >= max_charge_state) continue;

           double E = interpolate_field_magnitude(EMfields, particles->position, ip);
           double U_i = IonizationTables::ionization_potential[atomic_number][z];

           double rate = my_new_rate_formula(E, U_i, ...);

           double prob = 1.0 - std::exp(-rate * dt);
           if (patch->rand->uniform() < prob) {
               particles->charge[ip] = z + 1;
               create_freed_electron(ionization_electrons, particles, ip, ...);
           }
       }
   }
   ```

3. **Register** in `src/Ionization/IonizationFactory.h`:
   ```cpp
   else if (model_name == "my_new_model") {
       return new IonizationMyNew(params, species);
   }
   ```

4. **Add** the string to the validated list in `src/Params/Params.cpp::validateIonizationModel(...)`.

5. **Document** in `doc/Sphinx/Source/namelist.rst` with the formula and validity range.

6. **Test** — `validation/tst_ionization_mynew.py` exercising a known limit. Common choices: hydrogen at 10¹⁵ W/cm² (analytical ADK limit), or a published experiment with measured charge-state distribution.

For GPU support, add the same code in an `#ifdef _GPU` path. For models that read tabulated data (not pure formulas), the table must be uploaded to the GPU at startup.

## 10. Interaction with collisions and the freed-electron species

**Field ionization and collisional ionization are independent additive channels.** They both increment the charge state and create electrons; their rates are simply additive at low rates and they do not double-count physically.

**Order in the per-patch loop:**

```cpp
Species::dynamics() {
    Ionize->apply(particles, EMfields, ...);     // field ionization (this module)
    Interp->fieldsAtParticle(...);
    Push->operator()(particles, ...);
    partBoundCond.apply(...);
    Proj->depositCurrents(...);
}

// Later in the timestep:
BinaryProcesses::apply(patch, itime) {
    Collisions::apply(...);                       // includes collisional ionization
}
```

Field ionization runs *before* the pusher; collisional ionization runs *after* the full per-species dynamics (in the binary-processes phase). They modify the same `particles->charge[ip]` field but at different points in the timestep.

**Freed-electron species**: identified by `ionization_electrons_index` (set at construction from the namelist string). Each ionization event:

1. Increments `particles->charge[ip]` (the parent ion).
2. Resizes `ionization_electrons->particles` to add one entry.
3. Copies the parent ion's position to the new electron.
4. Initializes the new electron's momentum (typically zero for tunnel; can have a small kick consistent with ADK theory).
5. Sets the new electron's weight = parent ion's weight (per macro-particle, one electron is freed).

**Critical**: `ionization_electrons->particles->resize(n+1)` is not atomic. If you parallelize over particles within a patch (rather than across patches), this becomes a race. The current implementation is per-patch single-threaded for ionization, relying on patch-level OpenMP.

## 11. Common pitfalls

**Wrong table lookup**: `ionization_potential[Z][z]` is the energy to go from state `z` to state `z+1`. Off-by-one bugs are common when implementing new models.

**Forgetting to check `max_charge_state`**: without the check, the rate formula evaluates with invalid `U_i` (out-of-table or negative). Always:
```cpp
if (z >= max_charge_state) continue;
```

**Calling Python from a GPU kernel**: doesn't work. The `IonizationFromRate` path is CPU-only — factory must refuse GPU builds.

**Per-event field interpolation cost**: interpolating `E` at each particle position has cost comparable to the pusher's interpolation step. For light ionization rates, this can dominate the ionization runtime. The current code re-interpolates rather than reusing the pusher's interpolation — there's a TODO to share it.

**Forgetting energy bookkeeping**: when an electron is freed, `U_i` energy should be debited from the EM field. The current implementation records it in `DiagScalar.Uexp` for accounting; if you add a new model that doesn't update this, `Ubal` will drift.

**Per-particle momentum initialization for freed electrons**: ADK theory predicts a small initial momentum from the tunneling process. Smilei implements this; if you add a new model, decide whether the freed electron should have momentum set to zero (most cases) or to a model-specific kick.

**Race on freed-electron species**: `ionization_electrons->particles->resize(n+1)` from multiple threads within a patch would race. Don't parallelize the per-particle loop in `Ionization::apply` — patch-level parallelism is the only safe layer.

**`from_rate` callable returning wrong units**: the callable should return rate in *normalized* units (`1/T_r`), not SI (`1/s`). User code often forgets the conversion. Smilei doesn't check; the only symptom is wrong ionization rate (typically off by `omega_r ≈ 10¹⁵`).

**Ionization disabled by `is_test = True`**: test particles should not back-react on the EM field, but should they ionize? The current Smilei convention skips ionization for test particles. If you add a new model, follow the same convention.
