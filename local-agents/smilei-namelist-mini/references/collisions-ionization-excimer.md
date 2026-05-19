# Collisions, ionization, and excimer recipes

Binary collisions, field/collisional ionization, the Keldysh regime, and excimer-specific worked namelists. Three frontier-tier references merged for the mini tier.

## Collisions

### `Collisions` block syntax

```python
Collisions(
    species1 = ["electron"],
    species2 = ["ion"],
    coulomb_log = 0.,                 # 0 = auto (needs reference_angular_frequency_SI), >0 = fixed
    debug_every = 0,
    ionizing = False,                 # or species name string for collisional ionization
    every = 1,                        # apply every N timesteps
)
```

Both arguments are lists, so you can group:

```python
Collisions(species1=["electron"], species2=["Ar1","Ar2","Ar3"], coulomb_log=0.)
```

Intra-species (same species in both lists) auto-detects single-list pairing.

### Typical multi-block pattern

```python
# Pure hydrogen plasma — three channels
Collisions(species1=["electron"], species2=["electron"])
Collisions(species1=["electron"], species2=["ion"])
Collisions(species1=["ion"],      species2=["ion"])
```

**Always enable all three.** Omitting e-e is the most common mistake; in dense plasma, e-e drives electron thermalization faster than e-i. Without it, electrons stay artificially non-Maxwellian.

### Coulomb logarithm

- `coulomb_log = 0.` — auto-compute from local plasma. **Requires `reference_angular_frequency_SI`.**
- `coulomb_log = 5.` — fixed. For benchmarking or sensitivity studies.
- `coulomb_log_factor` — multiplier on auto value (default 1.0). Reset to 1.0 in production.

### Collisional ionization

Set `ionizing` to the name of the receiving electron species:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1"],
    coulomb_log = 0.,
    ionizing = "electron",        # name of electron species (NOT True/False)
)
```

For an ionization ladder, define separate `Species` per charge state and chain blocks:

```python
Species(name="Ar0", atomic_number=18, charge=0., ...)
Species(name="Ar1", atomic_number=18, charge=1., ...)
Species(name="Ar2", atomic_number=18, charge=2., ...)

Collisions(species1=["electron"], species2=["Ar0"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar1"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar2"], ionizing="electron")
```

Collisional ionization and field ionization are **independent additive channels**. For excimer-class moderate-intensity, enable both.

### PPC for collisional runs

| Channel | Minimum PPC |
|---|---|
| Collisionless | 8–16 |
| e-i only | 16 |
| Full (e-e, e-i, i-i) | 32 (64 recommended) |
| Collisional ionization | 32+ per state |

**Always enable `LoadBalancing(every=50)` for collisional runs** — dense regions stall MPI ranks otherwise.

## Field ionization

### Models

```python
Species(
    name = "Ar0",
    atomic_number = 18,
    charge = 0.,
    mass = 39.948*1836,
    ionization_model = "tunnel",                    # or "tunnel_full_PPT", "tunnel_BSI",
                                                    # "tunnel_envelope_averaged", "from_rate"
    ionization_electrons = "electron",              # name of receiving electron species
    maximum_charge_state = 8,                       # cap on the ladder
)
```

| Model | Use when |
|---|---|
| `None` | Pre-ionized plasma, no ionization |
| `"tunnel"` (ADK) | For `γ_K ≲ 0.5` (deep tunneling) |
| `"tunnel_BSI"` | Field exceeds Barrier Suppression threshold |
| `"tunnel_full_PPT"` | For `γ_K ~ 1` (excimer intermediate regime) |
| `"tunnel_envelope_averaged"` | ONLY with envelope laser model |
| `"from_rate"` | User-supplied rate callable (photoionization, tabulated) |

### Keldysh parameter

```
γ_K = ω × sqrt(2 m_e U_i) / (e E)
```

Higher `ω` (UV) and higher `U_i` raise `γ_K`. ADK is valid only for `γ_K ≪ 1`.

**KrF on Ar (U_i = 15.76 eV):**

| Intensity [W/cm²] | `γ_K` | Recommended |
|---|---|---|
| 10¹² | 18.8 | `from_rate` or skip ionization |
| 10¹³ | 5.9 | `"tunnel_full_PPT"` |
| 10¹⁴ | 1.9 | `"tunnel_full_PPT"` |
| 10¹⁵ | 0.6 | `"tunnel"` (acceptable) or PPT |
| 10¹⁶ | 0.2 | `"tunnel"` |

**This is the canonical reason UV regimes need PPT** — `γ_K` is 3–5× higher at the same intensity due to the ω factor.

### Envelope-averaged ionization

For envelope laser models, use the cycle-averaged rate:

```python
LaserEnvelopeGaussianAM(...)

Species(
    name = "Ar0",
    atomic_number = 18,
    ionization_model = "tunnel_envelope_averaged",   # matched to envelope laser
    ionization_electrons = "electron",
)
```

DC and AC rates differ for linear polarization, coincide for circular. Smilei handles the prefactor.

### `from_rate` for photoionization

Smilei has **no native photoionization model**. For ArF (6.4 eV/photon) or KrF (5.0 eV/photon) one-photon ionization of low-IP atoms:

```python
def my_rate(E):
    # E in normalized units; return probability per unit normalized time
    rate_SI = 1e10  # 1/s
    return rate_SI / omega_r

Species(
    name = "Ar0",
    ionization_model = "from_rate",
    ionization_rate = my_rate,
    ionization_electrons = "electron",
)
```

### Charge-state ladder

For field ionization, a single `Species` walks the ladder — the macro-particle's `charge` field increments as electrons are freed. `maximum_charge_state` caps the ladder. Mass is unchanged.

For collisional ionization, you need separate `Species` per charge state.

## Excimer-laser plasma recipes

### Why excimer is different from default Smilei use

| Property | Typical UHI | Typical excimer |
|---|---|---|
| Wavelength | 800–1064 nm | 157–351 nm |
| `n_c` | 10²¹ cm⁻³ | 10²²–4×10²² cm⁻³ |
| Pulse duration | 10–100 fs | 1–20 ns (ps for ICF) |
| Intensity | 10¹⁸–10²² W/cm² | 10⁸–10¹³ W/cm² |
| `a₀` | ≳ 1 | 10⁻⁵–10⁻² |
| Collisions | Often negligible | Dominant |
| Field ionization | Tunnel (γ_K ≪ 1) | Multiphoton-tunnel boundary (γ_K ~ 1) |

**Default UHI namelists copied to excimer regimes produce nonsense.**

### Resolution constraints — Debye dominates at typical excimer densities

```python
import math
dx_laser = 2*math.pi/16                          # 16 cells/λ
dx_skin  = 0.3 / math.sqrt(n_norm)               # 0.3 c/ω_pe at peak
dx_debye = 2.0 * math.sqrt(T_norm/n_norm)        # 2 × λ_D

dx = min(dx_laser, dx_skin, dx_debye)            # take the tightest
```

For typical excimer (10²¹–10²² cm⁻³, 30–100 eV), Debye is binding — much tighter than laser-resolution.

### Worked namelist — KrF on overdense Ar

```python
import math

# === inputs (SI) ===
lambda_SI      = 248e-9
intensity_SI   = 1.0e12
pulse_FWHM_SI  = 100e-15
waist_SI       = 8e-6
n_e_SI         = 2.0e22                          # cm⁻³ ≈ 5 × n_c(248 nm)
T_e_eV         = 30.0
box_x_SI, box_y_SI = 25e-6, 16e-6
sim_time_SI    = 400e-15

# === derived ===
c, e_, m_e, eps0 = 2.99792458e8, 1.602e-19, 9.11e-31, 8.85e-12
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r     = eps0*m_e*omega_r**2/e_**2

a0      = 0.85*math.sqrt(intensity_SI*1e-18)*(lambda_SI*1e6)
n_e     = (n_e_SI*1e6)/N_r
T_e     = T_e_eV/511e3
pulse   = pulse_FWHM_SI*omega_r
waist   = waist_SI/L_r
Lx, Ly  = box_x_SI/L_r, box_y_SI/L_r
T_tot   = sim_time_SI*omega_r

lambda_D = math.sqrt(T_e/n_e)
dx = min(2*math.pi/20, 2.0*lambda_D)

# === namelist ===
Main(
    geometry = "2Dcartesian", interpolation_order = 2, timestep_over_CFL = 0.95,
    cell_length = [dx, dx], grid_length = [Lx, Ly],
    number_of_patches = [32, 16], simulation_time = T_tot,
    EM_boundary_conditions = [["silver-muller"]]*2,
    reference_angular_frequency_SI = omega_r, print_every = 200,
)

LoadBalancing(every = 50)

LaserGaussian2D(
    box_side = "xmin", a0 = a0, omega = 1.0,
    focus = [Lx*0.5, Ly*0.5], waist = waist,
    time_envelope = tgaussian(fwhm = pulse, center = 3*pulse),
)

def n_profile(x, y):
    ramp_start = 0.4*Lx
    ramp_width = 1e-6/L_r
    if x < ramp_start: return 0.0
    elif x < ramp_start + ramp_width: return n_e*(x-ramp_start)/ramp_width
    else: return n_e

Species(name="electron", position_initialization="regular",
        momentum_initialization="maxwell-juettner", temperature=[T_e]*3,
        particles_per_cell=64, mass=1.0, charge=-1.0,
        number_density=n_profile,
        boundary_conditions=[["remove"], ["periodic"]])

Species(name="Ar1", position_initialization="regular",
        momentum_initialization="maxwell-juettner", temperature=[1e-5]*3,
        particles_per_cell=64, mass=39.948*1836, charge=1.0, atomic_number=18,
        number_density=n_profile,
        boundary_conditions=[["remove"], ["periodic"]])

Collisions(species1=["electron"], species2=["electron"], coulomb_log=0.)
Collisions(species1=["electron"], species2=["Ar1"],      coulomb_log=0.)
Collisions(species1=["Ar1"],      species2=["Ar1"],      coulomb_log=0.)

DiagScalar(every=50)
DiagFields(every=200, fields=["Ex","Ey","Bz","Rho_electron","Rho_Ar1"])
DiagProbe(every=100, origin=[0., Ly*0.5], corners=[[Lx, Ly*0.5]], number=[400],
          fields=["Ex","Ey","Bz","Rho"])
DiagParticleBinning(deposited_quantity="weight_ekin", every=200, species=["electron"],
                    axes=[["ekin", 1e-4, 1e-1, 100, "logscale"]])
```

### When Smilei is wrong tool

| Use case | Better choice |
|---|---|
| Industrial UV ablation (LASIK, polymer cutting) | Hydro codes (FLASH, HELIOS) |
| Long-pulse ICF ablation (full ns history) | Radiation-hydro (DRACO, FCI2) |
| Plasma discharge / glow discharge | LTP codes (HPEM) |
| Photolithography plasma | Specialized lithography codes |
| Pre-formed plasma + short-pulse excimer | Smilei is correct |
| LPI in DUV ICF (TPD, SBS, SRS) | Smilei is correct |

For long-pulse excimer ablation, the right approach is: hydro code for plume evolution + Smilei for kinetic sub-process in a representative window.

## Mistakes summary

- Missing `reference_angular_frequency_SI` with Collisions → silent OOM rate error (rule 1)
- Only e-i collisions → electrons stay non-Maxwellian
- PPC = 8 with collisions → noisy rates, spurious temperature spikes
- `ionizing = True` (Boolean) → use species name string instead
- Missing ions in ladder → must define `Ar1`, `Ar2`, etc. explicitly
- ADK at UV `γ_K > 1` → underestimates rate by orders of magnitude (use PPT)
- Missing `ionization_electrons` → freed electrons have no destination
- Envelope laser + `tunnel` (instantaneous-field model) → mismatched, wrong rate (use `tunnel_envelope_averaged`)
- `maximum_charge_state = 1` when laser drives higher → ladder caps too early
- No `LoadBalancing` on collisional run → idle MPI ranks
