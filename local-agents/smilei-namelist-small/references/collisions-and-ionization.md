# Collisions and ionization

Binary collisions (Pérez 2012), collisional ionization, field ionization (ADK / PPT / BSI / envelope-averaged / `from_rate`), and the Keldysh-parameter regime selection. Merged from frontier's two physics-module references.

## Table of contents

1. The Pérez 2012 collision algorithm (one paragraph)
2. `Collisions` block syntax
3. Coulomb logarithm
4. `CollisionalIonization`
5. Multiple `Collisions` blocks — typical pattern
6. PPC requirements for collisional runs
7. Field-ionization models — overview
8. The Keldysh parameter — which model when
9. ADK tunnel ionization
10. PPT — `tunnel_full_PPT`
11. BSI — `tunnel_BSI`
12. Envelope-averaged ionization
13. User-defined rate — `from_rate`
14. Charge-state ladder configuration
15. Diagnostics for collisions and ionization
16. Common mistakes

---

## 1. The Pérez 2012 collision algorithm (one paragraph)

Smilei implements binary collisions via Pérez et al. (Phys. Plasmas 2012). Within each cell at each timestep, macro-particles of the two specified species are paired and a single Monte-Carlo binary collision is computed per pair, drawing the scattering angle from the relativistic Coulomb cross-section. The pairing is randomized each timestep. Different macro-particle weights are handled via probability rescaling. The scheme is statistical — collision rates have variance scaling as `1/sqrt(PPC)`. Practical floor for collisional runs: `PPC ≥ 32`. Safe default for excimer/LTP: `PPC ≥ 64`.

## 2. `Collisions` block syntax

```python
Collisions(
    species1 = ["electron"],
    species2 = ["ion"],
    coulomb_log = 0.,                 # 0 = auto, >0 = fixed
    debug_every = 0,                  # diagnostic cadence
    ionizing = False,                 # or species name string for collisional ionization
    every = 1,                        # apply every N timesteps
)
```

Both arguments are lists, so you can group:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1", "Ar2", "Ar3"],     # electrons collide against all three states
)
```

**Intra-species (collision of a species with itself):**

```python
Collisions(species1 = ["electron"], species2 = ["electron"], coulomb_log = 0.)
```

Smilei detects this and uses single-list pairing (no double-counting).

## 3. Coulomb logarithm

```python
coulomb_log = 0.        # automatic (computed from local plasma parameters)
coulomb_log = 5.        # fixed
```

**Automatic mode** (`coulomb_log = 0.`): Smilei computes `ln Λ` each timestep from local temperature, density, and species masses/charges. **Requires `reference_angular_frequency_SI`** (hard rule 1).

**Fixed mode**: for benchmarking against analytical results, reproducing legacy results, or sensitivity studies. Typical values: 3–10 for laser-plasma, 10–20 for fusion.

`coulomb_log_factor` (default 1.0) multiplies the auto-computed value. Use only for sensitivity studies; reset to 1.0 in production.

## 4. `CollisionalIonization`

Electron-impact ionization via Lotz-Drawin cross-sections. Enabled by setting `ionizing` to the **name of the electron species** that receives freed electrons:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1"],
    coulomb_log = 0.,
    ionizing = "electron",          # name of electron species, not a Boolean
)
```

For an ionization ladder, define separate `Species` per charge state and chain the collision blocks:

```python
Species(name = "Ar0", atomic_number = 18, charge = 0., ...)
Species(name = "Ar1", atomic_number = 18, charge = 1., ...)
Species(name = "Ar2", atomic_number = 18, charge = 2., ...)

Collisions(species1=["electron"], species2=["Ar0"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar1"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar2"], ionizing="electron")
```

Collisional ionization and field ionization are **independent additive channels**. Both should be enabled for excimer-class moderate-intensity plasmas — field dominates early (laser strong, plasma cold), collisional dominates later (laser past peak, plasma hot).

## 5. Multiple `Collisions` blocks — typical pattern

```python
# Hydrogen plasma — three channels
Collisions(species1=["electron"], species2=["electron"])
Collisions(species1=["electron"], species2=["ion"])
Collisions(species1=["ion"],      species2=["ion"])

# Argon plasma with charge-state ladder
Collisions(species1=["electron"], species2=["electron"])
Collisions(species1=["electron"], species2=["Ar1","Ar2","Ar3"])
Collisions(species1=["Ar1","Ar2","Ar3"], species2=["Ar1","Ar2","Ar3"])
# Plus collisional ionization
Collisions(species1=["electron"], species2=["Ar1"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar2"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar3"], ionizing="electron")
```

**Omitting e-e collisions** is the most common mistake — in dense plasma, e-e drives electron thermalization faster than e-i. Without it, electrons stay artificially non-Maxwellian.

## 6. PPC requirements for collisional runs

| Channel | Minimum PPC | Recommended |
|---|---|---|
| Collisionless | — | 8–16 |
| e-i only | 16 | 32 |
| Full (e-e, e-i, i-i) | 32 | 64 |
| Collisional ionization | 32 | 64+ per state |
| Nuclear reactions | 128 | 256+ |

If you see noisy or oscillatory collision-driven energy transfer, double PPC and retry — that's the diagnostic.

**Always enable `LoadBalancing`** for collisional runs. Cells with dense plasma have proportionally more collision work; without load balancing, MPI ranks idle waiting for dense-plasma ranks.

## 7. Field-ionization models — overview

`ionization_model` on a `Species` block:

| Model | When to use |
|---|---|
| `None` (default) | Pre-ionized plasma, no ionization |
| `"tunnel"` | ADK, instantaneous field. For `γ_K ≲ 0.5`. |
| `"tunnel_BSI"` | ADK + Barrier Suppression. When field approaches BSI threshold. |
| `"tunnel_full_PPT"` | PPT (full multiphoton+tunnel). For `γ_K ~ 1` (excimer intermediate). |
| `"tunnel_envelope_averaged"` | Cycle-averaged ADK. Use ONLY with envelope laser. |
| `"from_rate"` | User-supplied callable rate. For photoionization, tabulated rates. |

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    mass = 39.948*1836,
    ionization_model = "tunnel",
    ionization_electrons = "electron",          # name of receiving electron species
    maximum_charge_state = 8,                   # cap on the ladder
)
```

## 8. The Keldysh parameter — which model when

```
γ_K = ω × sqrt(2 m_e U_i) / (e E)
```

`U_i` = ionization potential, `E` = peak field amplitude.

| Regime | `γ_K` | Model |
|---|---|---|
| Tunnel | `≪ 1` (often `< 0.5`) | `"tunnel"` |
| Intermediate | `~ 1` | `"tunnel_full_PPT"` |
| Multiphoton | `≫ 1` | `"from_rate"` (not natively supported) |
| BSI | Field exceeds barrier | `"tunnel_BSI"` |

**Worked example: KrF on Ar (`U_i = 15.76 eV`, `ω = 7.59 × 10¹⁵ rad/s`):**

| Intensity [W/cm²] | `γ_K` | Recommended model |
|---|---|---|
| 10¹² | 18.8 | `from_rate` (multiphoton) or no ionization |
| 10¹³ | 5.9 | `"tunnel_full_PPT"` |
| 10¹⁴ | 1.9 | `"tunnel_full_PPT"` |
| 10¹⁵ | 0.6 | `"tunnel"` (acceptable) or PPT |
| 10¹⁶ | 0.2 | `"tunnel"` |

This is the canonical reason UV regimes need PPT, not ADK — `γ_K` is 3–5× higher at the same intensity due to the `ω` factor.

## 9. ADK tunnel ionization

Default for `γ_K ≪ 1`. Uses instantaneous laser field amplitude. Macro-particle's `charge` field increments stochastically based on `Γ_ADK × Δt`.

**Assumptions that fail at UV:**
- Adiabatic field (`γ_K ≪ 1`) — fails at moderate UV intensities
- Hydrogenic effective potential
- Single-active-electron

For excimer regimes, ADK is correct only at high intensities (`I > 10¹⁵ W/cm²` typically). Below, use PPT.

## 10. PPT — `tunnel_full_PPT`

Full Perelomov-Popov-Terent'ev formula. Interpolates between multiphoton and tunnel limits without the adiabatic approximation. Reduces to ADK for `γ_K ≪ 1`; gives correct multiphoton power-law `(E)^(2N)` for `γ_K ≫ 1`.

**Recommended for excimer regimes** with `γ_K` between 0.5 and 5. Computational cost comparable to ADK.

```python
Species(..., ionization_model = "tunnel_full_PPT", ...)
```

## 11. BSI — `tunnel_BSI`

ADK extended with Barrier-Suppression cap. When the laser field exceeds the BSI threshold (the electron is essentially free regardless of tunnel time), ADK's exponential rate is capped at the geometric limit.

Rarely needed for excimer regimes (`I < 10¹⁶ W/cm²`). Mandatory for UHI Ti:Sa above 10²⁰ W/cm².

## 12. Envelope-averaged ionization

Envelope laser models propagate the slowly-varying envelope, not the carrier. ADK uses the *instantaneous* field, which envelope doesn't have. Use the cycle-averaged rate `Γ_ADK,AC`:

```python
LaserEnvelopeGaussianAM(...)         # envelope laser

Species(
    name = "Ar0",
    atomic_number = 18,
    ionization_model = "tunnel_envelope_averaged",    # matched to envelope laser
    ionization_electrons = "electron",
)
```

The DC and AC rates differ for linear polarization but coincide for circular. Smilei computes the exact prefactor; the user just sets the model name.

**Caveat**: with envelope models, the freed electron's quiver momentum at ionization is lost in averaging. Smilei initializes freed electrons with zero transverse momentum. For physics depending on freed-electron initial momentum, benchmark against a short full-EM run.

## 13. User-defined rate — `from_rate`

For physics outside ADK/PPT — photoionization at specific wavelengths, tabulated rates, multiphoton at specific cross-section:

```python
def my_rate(E):
    # E in normalized units; return probability per unit normalized time
    import math
    rate_SI = 1e10  # 1/s, hypothetical
    return rate_SI / omega_r          # convert to normalized

Species(
    name = "Ar0",
    ionization_model = "from_rate",
    ionization_rate = my_rate,
    ionization_electrons = "electron",
    maximum_charge_state = 1,
)
```

**Right pathway for excimer photoionization** (single-photon at ArF 6.4 eV / KrF 5.0 eV) and for tabulated experimental rates. Smilei has no native photoionization model; this is the workaround.

`ionization_rate` can be one callable (same rate for all states) or a list of callables (one per charge state).

## 14. Charge-state ladder configuration

### For field ionization: single `Species` walks the ladder

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    mass = 39.948 * 1836,
    ionization_model = "tunnel",
    ionization_electrons = "electron",
    maximum_charge_state = 18,        # cap at fully stripped
)
```

The same macro-particle's `charge` field increments as electrons are freed. Mass is unchanged. The diagnostic species remains `"Ar0"` even when particles are in higher states — use `DiagParticleBinning` with `axes=[["charge", -0.5, 18.5, 19]]` to see the actual distribution.

### For collisional ionization: separate `Species` per charge state

Necessary because `Collisions` groups by species name. The two systems can coexist: field ionization walking within one species, plus collisional ionization advancing through separate per-charge-state species.

### `maximum_charge_state`

Caps the ladder. Set to expected peak charge state + margin, or to `atomic_number` (no cap) for safety.

## 15. Diagnostics for collisions and ionization

```python
# Charge-state distribution
DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = 200, species = ["Ar0"],
    axes = [["charge", -0.5, 18.5, 19]],
)

# Ionized electron density
DiagFields(every = 200, fields = ["Rho_electron"])

# Spatial map of ionization state
DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = 200, species = ["Ar0"],
    axes = [["x", 0., Lx, 200], ["y", 0., Ly, 100]],
)

# Per-event creation tracking
DiagNewParticles(
    species = "electron",
    every = 100,
    attributes = ["x", "y", "px", "py", "pz", "w"],
)

# Collision-rate debug
Collisions(..., debug_every = 100)
```

Watch `DiagScalar` outputs `Ukin` and `Uelm` for energy conservation. `Utot = Ukin + Uelm` should be conserved (modulo absorbed laser energy `Uexp`) — if not, simulation is under-resolved or PPC too low.

## 16. Common mistakes

**Missing `reference_angular_frequency_SI` with `Collisions`**: hard rule 1. Auto Coulomb log computes nonsense without it.

**Only e-i collisions, no e-e or i-i**: e-e drives electron thermalization. Dense plasmas without e-e stay artificially non-Maxwellian.

**Low PPC with collisions**: PPC = 8 with collisions enabled gives ~35% per-cell variance in collision rate. Symptoms: spurious temperature spikes, oscillatory energy transfer. Fix: `PPC ≥ 32`.

**`ionizing = True` (Boolean) instead of species name**: current Smilei expects `ionizing = "electron"` (the species name). The Boolean syntax is legacy.

**Missing ions in the ionization ladder**: must define `Ar1`, `Ar2`, etc. explicitly. Smilei doesn't auto-create them.

**ADK at UV `γ_K > 1`**: KrF at 10¹³ W/cm² has `γ_K ≈ 6` — ADK underestimates by orders of magnitude. Use PPT.

**Missing `ionization_electrons`**: freed electrons have no destination. Energy conservation breaks.

**Mismatched envelope-laser + ionization model**: envelope laser + `"tunnel"` (instantaneous-field) is inconsistent. Use `"tunnel_envelope_averaged"`.

**Ladder cap too low**: `maximum_charge_state = 1` when laser drives to Ar+8. Set with margin or use `atomic_number`.

**Forgetting `LoadBalancing` on collisional run**: idle MPI ranks waiting for dense-plasma ranks. Enable `LoadBalancing(every=50)`.

**`coulomb_log_factor ≠ 1` left over from sensitivity study**: silently wrong physics in production. Always check.
