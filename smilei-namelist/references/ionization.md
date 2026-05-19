# Ionization

This reference covers Smilei's field-ionization models (ADK tunnel, Barrier-Suppression, PPT, envelope-averaged, user-defined rate), the Keldysh-parameter validity boundaries, charge-state ladder setup, and the interaction with the laser envelope model. Collisional (impact) ionization is covered in `collisions-and-reactions.md`. For UV-regime considerations and the ADK validity question at excimer wavelengths see `excimer-and-LTP-recipes.md` §9.

## Table of contents

1. Available ionization models
2. The Keldysh parameter and validity of each model
3. ADK tunnel ionization — the default
4. Barrier-Suppression Ionization
5. PPT (Perelomov-Popov-Terent'ev) — tunnel-full
6. Envelope-averaged ionization
7. User-defined ionization rate (`from_rate`)
8. Charge-state ladder configuration
9. The Monte-Carlo scheme
10. Interaction with collisions
11. Diagnostics for ionization
12. Common mistakes

---

## 1. Available ionization models

`ionization_model` on a `Species` block selects the model:

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    mass = 39.948*1836,
    ionization_model = "tunnel",                    # ← model selection
    ionization_electrons = "electron",              # name of receiving electron species
    maximum_charge_state = 18,
    # ...
)
```

| `ionization_model` | Description | Use when |
|---|---|---|
| `None` (default) | No field ionization | Pre-ionized plasma, collisionless cold plasma |
| `"tunnel"` | ADK tunneling, instantaneous field | Standard case for `γ_K ≲ 0.5` |
| `"tunnel_BSI"` | Tunnel + Barrier Suppression | When field approaches the BSI threshold (above 1 a.u.) |
| `"tunnel_full_PPT"` | Full PPT formula (multiphoton + tunnel) | When `γ_K ~ 1` (excimer/UV intermediate regime) |
| `"tunnel_envelope_averaged"` | DC-ADK rate averaged over a laser cycle | Use ONLY with envelope laser models |
| `"from_rate"` | User-supplied callable rate `Γ(E)` | Custom physics, tabulated rates, photoionization |

The default `"tunnel"` is what most published UHI laser-plasma papers use. The PPT and BSI variants exist for regimes where ADK underestimates or overestimates the rate.

## 2. The Keldysh parameter and validity of each model

The Keldysh adiabaticity parameter:

```
γ_K = ω × sqrt(2 m_e U_i) / (e E)
```

where `U_i` is the ionization potential (the energy needed to free the next bound electron) and `E` is the peak laser electric field amplitude.

| Regime | `γ_K` | Physics | Smilei model |
|---|---|---|---|
| Tunnel | `γ_K ≪ 1` (often `< 0.5`) | Quasi-static field, electron tunnels through suppressed barrier | `"tunnel"` |
| Intermediate | `γ_K ~ 1` | Mixed tunnel/multiphoton | `"tunnel_full_PPT"` |
| Multiphoton | `γ_K ≫ 1` | Few-photon absorption | Not natively supported; use `from_rate` |
| BSI | Field exceeds barrier threshold | Free electron without tunneling | `"tunnel_BSI"` |

### Worked example: KrF on Ar at varying intensity

`λ = 248 nm → ω = 7.59 × 10¹⁵ rad/s`, `U_i(Ar) = 15.76 eV`.

| Intensity [W/cm²] | E [V/m] | `γ_K` | Recommended model |
|---|---|---|---|
| 10¹² | 2.7 × 10⁹ | 18.8 | `from_rate` (multiphoton) or skip ionization |
| 10¹³ | 8.7 × 10⁹ | 5.9 | `"tunnel_full_PPT"` |
| 10¹⁴ | 2.7 × 10¹⁰ | 1.9 | `"tunnel_full_PPT"` |
| 10¹⁵ | 8.7 × 10¹⁰ | 0.6 | `"tunnel"` (acceptable) or PPT |
| 10¹⁶ | 2.7 × 10¹¹ | 0.2 | `"tunnel"` |

This is the canonical reason UV regimes need PPT, not ADK. At Ti:Sa (800 nm) intensities, `γ_K` is typically already small at relevant intensities; at UV the same intensity gives `γ_K` 3-5× larger because of the `ω` factor.

## 3. ADK tunnel ionization — the default

The Ammosov-Delone-Krainov rate uses the instantaneous laser electric field amplitude `|E(t)|` (instantaneous, not envelope-averaged) and the bound-electron quantum numbers `(n*, l*, m)` to compute the ionization rate per unit time. The macro-particle's charge state increments stochastically based on `Γ_ADK Δt`.

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    ionization_model = "tunnel",
    ionization_electrons = "electron",
    maximum_charge_state = 8,           # cap the ladder — see §8
)
```

ADK assumptions:
1. **Adiabatic field**: laser frequency much less than orbital frequency of the bound electron (`γ_K ≪ 1`)
2. **Hydrogenic effective potential**: bound electron in a Coulomb-like potential modified by quantum defects
3. **Single-active-electron**: only one electron is ionized per event; inner-shell electrons are frozen
4. **No magnetic field effects**: ignored

For excimer regimes, assumption 1 is the one that fails first. Use ADK only when `γ_K < 0.5` *throughout* the pulse — at the peak intensity. If `γ_K ~ 1` at any time during the pulse, switch to PPT.

ADK requires `reference_angular_frequency_SI` (hard-rule 1) because ionization potentials are in eV.

## 4. Barrier-Suppression Ionization

When the laser field is strong enough to suppress the binding potential below the bound-state energy, the electron is essentially free. ADK doesn't capture this — it predicts an exponentially increasing rate that diverges. BSI caps the rate at the geometric (over-the-barrier) limit.

```python
Species(
    # ...
    ionization_model = "tunnel_BSI",
)
```

BSI is the recommended model when peak intensity exceeds the BSI threshold for the ion of interest. For typical excimer regimes (`I < 10¹⁶ W/cm²`), BSI is rarely needed — the ions are not driven into the over-the-barrier regime. For UHI Ti:Sa runs above `10²⁰ W/cm²`, BSI matters for low-Z ions.

## 5. PPT (Perelomov-Popov-Terent'ev) — `tunnel_full_PPT`

The PPT formula interpolates between the multiphoton and tunnel limits without the adiabatic approximation. In the tunnel limit it reduces to ADK; in the multiphoton limit it gives the appropriate power-law `(E)^(2N)` with `N` being the minimum number of photons needed to ionize.

```python
Species(
    # ...
    ionization_model = "tunnel_full_PPT",
)
```

For excimer regimes with `γ_K` between 0.5 and 5, PPT is the most accurate model. Cost is comparable to ADK (a few extra trig/log evaluations per ionization check). It does not solve the multiphoton-dominated case (`γ_K > 5`) — for that, use `from_rate` with experimentally-tabulated rates.

## 6. Envelope-averaged ionization

Envelope laser models propagate the slowly-varying envelope rather than the carrier oscillation. ADK uses the *instantaneous* field, which the envelope model doesn't have. The correct rate is the cycle-averaged rate `Γ_ADK,AC`:

```python
# In the namelist where you use an envelope laser:
Main(...)
LaserEnvelopeGaussianAM(...)

Species(
    name = "Ar0",
    atomic_number = 18,
    ionization_model = "tunnel_envelope_averaged",   # ← matched to envelope laser
    ionization_electrons = "electron",
)
```

The `Γ_ADK,AC` formula differs from `Γ_ADK,DC` for linear polarization but they coincide for circular (`ellipticity = 1`) because the field amplitude is constant over the cycle.

For linear polarization, the relation is:

```
Γ_ADK,AC = (3/π) × (|E|/E_a)^(1/2) × Γ_ADK,DC          (approximate)
```

where `E_a` is the atomic field. Smilei computes the exact prefactor; the user just sets the model name. See Chen et al. 2013 for derivation.

**Newly-ionized electron momentum**: with envelope models, the freed electron's quiver momentum at ionization is lost in the average — Smilei initializes the electron with zero transverse momentum by default. For accurate downstream dynamics this can be wrong; the Smilei docs (Massimo 2020a) discuss the appropriate momentum initialization.

If your physics depends on the freed-electron initial momentum (e.g., subsequent acceleration), benchmark envelope-ionization runs against a short full-EM run.

## 7. User-defined ionization rate (`from_rate`)

For physics outside the ADK/BSI/PPT family — photoionization at specific wavelengths, tabulated experimental rates, multiphoton with a specific cross-section — provide a callable:

```python
def my_rate(E):
    # E is the local electric field magnitude in normalized units
    # return ionization probability per unit time in normalized units (1/T_r)
    import math
    # Example: 6-photon ionization with cross-section sigma_6
    sigma_6 = 1e-50                                          # cm^12 s^-1, hypothetical
    # ... convert E to intensity, compute rate ...
    return rate_in_normalized_units

Species(
    name = "Ar0",
    ionization_model = "from_rate",
    ionization_rate = my_rate,
    ionization_electrons = "electron",
    maximum_charge_state = 1,                                # for this example, just one step
)
```

This is the right pathway for excimer photoionization (single-photon at ArF 6.4 eV / KrF 5.0 eV) and for any tabulated experimental rate. The callable receives the instantaneous local electric field and returns a probability per unit time. The Monte-Carlo sampling (§9) takes care of the actual ionization events.

`ionization_rate` can be one callable for the whole ladder, or a list of callables (one per charge state) for arbitrarily detailed rate models.

## 8. Charge-state ladder configuration

For ionization that goes beyond a single step, you have two options:

### Option A — single `Species` walks the ladder

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    mass = 39.948 * 1836,
    ionization_model = "tunnel",
    ionization_electrons = "electron",
    maximum_charge_state = 18,                # cap at fully stripped
)
```

The same macro-particle's `charge` field increments as the electron is freed. The mass is unchanged (the freed electron mass is negligible). The diagnostic species name remains `"Ar0"` even when individual particles are in the +5 state — visualizing requires `DiagParticleBinning` with `weight_charge` to see the actual charge distribution.

This is the standard pattern for ionization with the `"tunnel"`, `"tunnel_BSI"`, `"tunnel_full_PPT"`, or `"tunnel_envelope_averaged"` models.

### Option B — separate `Species` per charge state (collisional ionization)

For collisional ionization, separate species are needed because the collision module groups by species. See `collisions-and-reactions.md` §5 for the pattern. The two systems can coexist: field ionization within one species walking the ladder, plus collisional ionization advancing through separate per-charge-state species.

### `maximum_charge_state`

Caps the ladder. If unset, Smilei goes up to `atomic_number` (fully stripped). For practical runs, set `maximum_charge_state` to the highest state you expect to populate to save compute:

```python
maximum_charge_state = 8                   # don't bother computing states 9-18 of Ar
```

This is an optimization, not a physics control — the rate is computed for each step, but states above the cap are not allowed.

## 9. The Monte-Carlo scheme

At each timestep, for each macro-particle of the ionizing species in each cell, Smilei:

1. Computes the local electric field magnitude at the particle position (via interpolation from grid).
2. Evaluates `Γ(|E|, current_charge_state) × Δt` — the ionization probability for this step.
3. Draws a uniform random number `r ∈ [0,1]`.
4. If `r < Γ × Δt`, the particle's `charge` is incremented, a new electron is added to the receiving species at the same position, and energy is removed from the EM field to conserve total energy (the freed electron's "binding energy" is taken from the field).

For the deposited electron's momentum: by default it inherits the parent ion's momentum (which is small for cold ions). For tunnel models, the freed electron also receives the small longitudinal kick consistent with ADK theory.

Random-number streams are reproducible if `Main.random_seed` is set.

## 10. Interaction with collisions

Field ionization and collisional ionization are independent channels:

```python
# Field ionization on Ar0
Species(name="Ar0", atomic_number=18, ionization_model="tunnel", ...)

# Collisional ionization on Ar0
Collisions(species1=["electron"], species2=["Ar0"], ionizing="electron")
```

Both increment the charge state and produce electrons. Their rates are simply additive at low rates, and they do not double-count physically (one is laser-driven, one is collision-driven). For excimer-class plasmas at moderate intensity, **both should be enabled** because each dominates in different parts of the pulse-plasma interaction:

- **Early in the pulse**: laser is on, plasma is cold and weakly ionized → field ionization dominates the first few % of ionization
- **Later**: plasma is hot, laser is past peak → collisional ionization drives the rest

Disabling either gives wrong dynamics. The natural assumption ("the laser is strong, collisions are slow, ignore collisional ionization") often fails because field ionization saturates at modest charge states (low rate at large `U_i`) while electron impact ionization keeps the ladder climbing.

## 11. Diagnostics for ionization

To watch the ionization process:

```python
# Distribution over charge state via ParticleBinning
DiagParticleBinning(
    deposited_quantity = "weight_charge",   # number-weighted charge
    every = 200,
    species = ["Ar0"],
    axes = [["charge", -0.5, 18.5, 19]],    # 19 bins covering 0 to 18
)

# Or total ionized electron density
DiagFields(every = 200, fields = ["Rho_electron"])

# Local ionization state map
DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = 200,
    species = ["Ar0"],
    axes = [["x", 0., 60., 200], ["y", 0., 40., 100]],
)
```

For events-level diagnostics on newly created electrons (where, when, energy):

```python
DiagNewParticles(
    species = "electron",
    every = 100,
    attributes = ["x", "y", "px", "py", "pz", "w"],
)
```

This records every electron creation event during the simulation. File size grows with ionization rate; use sparingly or with sampling.

## 12. Common mistakes

### Mistake 1: ADK at UV intensities where `γ_K > 1`

**BAD**:
```python
# KrF at 10^13 W/cm² — γ_K ≈ 6, ADK is in the wrong regime
Species(name="Ar0", ionization_model="tunnel", ...)
```

**GOOD**:
```python
Species(name="Ar0", ionization_model="tunnel_full_PPT", ...)
```

ADK underestimates the rate by orders of magnitude in the multiphoton regime. Always compute `γ_K` at peak intensity first.

### Mistake 2: missing `ionization_electrons`

**BAD**:
```python
Species(name="Ar0", ionization_model="tunnel", ...)
# ionization_electrons not set → no electron species to receive freed electrons
```

The simulation parses but no electrons are produced — the freed charge goes into the void and energy conservation breaks. Always specify `ionization_electrons` matching the name of an existing electron `Species`.

### Mistake 3: mismatched envelope model

**BAD**:
```python
# Envelope laser
LaserEnvelopeGaussianAM(...)
# But species uses the instantaneous-field model
Species(name="Ar0", ionization_model="tunnel", ...)
```

The instantaneous-field ADK rate has no instantaneous field to evaluate (the envelope discards the carrier). Smilei may issue a warning or compute with an averaged field; either way the result is wrong. Match the models:

```python
LaserEnvelopeGaussianAM(...)
Species(name="Ar0", ionization_model="tunnel_envelope_averaged", ...)
```

### Mistake 4: forgetting to define electron species

**BAD**:
```python
Species(name="Ar0", ionization_electrons="electron", ...)
# No Species(name="electron", ...) block defined
```

Smilei raises a parse error here. Easy to spot but easy to do in iterative namelist edits.

### Mistake 5: ladder cap too low

```python
maximum_charge_state = 1
# but the laser drives ionization to Ar+8 in practice
```

Macro-particles stop ionizing at the cap, so the plasma doesn't reach the right charge state, which underestimates collisional rates, mean-free-paths, etc. Set the cap based on the expected peak charge state with margin, or set it to `atomic_number` (no cap) for safety.

### Mistake 6: missing `reference_angular_frequency_SI`

Hard-rule 1. The ionization potential lookup requires SI units. Without it, the ADK formula evaluates with garbage `U_i` and either produces no ionization or saturates immediately.

### Mistake 7: testing ionization with `momentum_initialization = "cold"`

If the ion species is initialized cold and the field is too weak to ionize it, no electrons appear and the user concludes "ionization doesn't work." Increase intensity, switch model (try `from_rate` with a known nonzero rate for sanity), or add a small temperature.

### Mistake 8: assuming photoionization is included

Smilei's `"tunnel"` / `"tunnel_full_PPT"` / `"tunnel_BSI"` models cover only field ionization — they do NOT include direct photoionization (absorption of one or a few photons with energy above the ionization potential). For excimer wavelengths where photoionization is a real channel (ArF at 6.4 eV/photon, vacuum UV at higher energies), you must implement it as a `from_rate` model using tabulated cross-sections.

### Mistake 9: stochastic noise from low PPC

Ionization is a Monte-Carlo process — at low PPC the ionized-electron density is noisy. For initial-condition-sensitive physics (instabilities seeded by inhomogeneity), low PPC injects noise that's hard to distinguish from real seeds. Raise PPC for ionization-heavy runs.

### Mistake 10: comparing ADK and PPT without realizing they may give different `maximum_charge_state` effective behavior

A run that with ADK saturates at Ar+5 may with PPT saturate at Ar+8 simply because PPT has a different rate dependence on field strength. When changing models, the entire ionization dynamics may differ. Document model choice in any production paper.
