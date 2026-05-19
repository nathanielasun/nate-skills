# Basics and units

This reference covers Smilei's normalization scheme, the reference angular frequency, how to convert between normalized and SI units in both directions, and the common unit-related mistakes that produce subtle physics errors.

## Table of contents

1. Why Smilei normalizes the way it does
2. The reference quantities
3. `reference_angular_frequency_SI` — when it's optional, when it's mandatory
4. Critical density, plasma frequency, skin depth
5. Converting SI inputs to normalized values inside the namelist
6. Converting normalized outputs to SI in post-processing
7. Spatial profiles and the normalization trap
8. Temperature, momentum, and energy units
9. Common unit mistakes (BAD/GOOD pairs)
10. Quick reference card

---

## 1. Why Smilei normalizes the way it does

PIC codes that solve Maxwell-Vlasov are scale-invariant in a useful way: if you express all quantities in units derived from the speed of light `c`, the elementary charge `e`, the electron mass `m_e`, and one chosen frequency `ω_r`, the equations themselves become dimensionless. Smilei exploits this. The user picks `ω_r` *implicitly* by their choice of laser wavelength or plasma frequency, and the code never needs to know its SI value — *unless* a physics module breaks the scale invariance.

This has three practical consequences:

1. **Most simulations don't need to set `reference_angular_frequency_SI` at all.** A pure laser-plasma run can stay fully normalized and the user scales results afterward.
2. **Anything that breaks scale invariance requires a fixed `ω_r` in SI.** Collisions, ionization, radiation reaction, and certain QED processes all involve physical constants (Coulomb logarithm, atomic ionization potentials, Compton wavelength) that don't normalize away. Smilei needs `ω_r` in rad/s to compute these.
3. **The user's mental model and the code's units can diverge silently.** A user thinking "I'm setting density to 0.1 cm⁻³" who writes `number_density = 0.1` has actually set it to `0.1 × n_c` where `n_c` depends on their (possibly unset) `ω_r`. This is the single most common source of confusion.

## 2. The reference quantities

Once `ω_r` is chosen (explicitly or implicitly), the rest follow:

| Quantity | Symbol | Definition | SI expression |
|---|---|---|---|
| Angular frequency | `ω_r` | chosen | — |
| Time | `T_r` | `1/ω_r` | seconds |
| Length | `L_r` | `c/ω_r` | meters |
| Velocity | — | `c` | m/s |
| Mass | — | `m_e` | kg |
| Charge | — | `e` | C |
| Electric field | `E_r` | `m_e c ω_r / e` | V/m |
| Magnetic field | `B_r` | `m_e ω_r / e` | T |
| Number density | `N_r` | `ε_0 m_e ω_r² / e²` | m⁻³ (this is `n_c` for the chosen `ω_r`) |
| Current density | `J_r` | `c e N_r` | A/m² |
| Energy | — | `m_e c²` | J (= 511 keV) |
| Temperature | — | `m_e c² / k_B` | K (relativistic), or `m_e c²` in energy units |
| Pressure | — | `m_e c² N_r` | Pa |

The key identity: **`N_r = n_c(ω_r)`** — the reference number density is exactly the critical density at the reference frequency. So `number_density = 0.5` always means "half the critical density at the reference frequency you picked."

## 3. `reference_angular_frequency_SI` — when it's optional, when it's mandatory

### When it's optional

Pure Maxwell-Vlasov simulations: laser propagation in a cold or warm plasma, wakefield acceleration without ionization, RT/Weibel instabilities, basic plasma waves, magnetic reconnection, anything where the only physics is the Vlasov-Maxwell system with collisionless particles.

In these cases the user is free to scale results a posteriori by choosing any `ω_r` that makes sense for their problem. Leaving it unset is fine.

### When it's mandatory

Any of these blocks/options in the namelist forces `reference_angular_frequency_SI` to be set:

- `Collisions(...)` — uses Coulomb logarithm, which depends on absolute densities
- `CollisionalIonization(...)` — uses ionization potentials in eV
- `Species(..., ionization_model="tunnel")` — uses atomic field in V/m
- `Species(..., ionization_model="tunnel_envelope_averaged")` — same
- `Species(..., ionization_model="tunnel_full_PPT")` — same
- `Species(..., ionization_model="from_rate")` — typically uses SI-tabulated rates
- `RadiationReaction(...)` — uses Compton wavelength
- `DiagRadiationSpectrum(...)` — outputs photon energies in keV/MeV
- `MultiphotonBreitWheeler` — uses fundamental constants

### What it actually does

`reference_angular_frequency_SI` only affects the modules listed above. It does **not** convert simulation outputs to SI — they remain in normalized units regardless. Conversion happens at post-processing time via happi (`units=['um','fs','GV/m']` and similar).

### How to compute it for common lasers

```python
import math
c = 2.99792458e8                     # speed of light

# Excimer lasers
omega_F2   = 2*math.pi*c/157e-9      # ≈ 1.200e16  (F₂ 157 nm)
omega_ArF  = 2*math.pi*c/193e-9      # ≈ 9.760e15  (ArF 193 nm)
omega_KrF  = 2*math.pi*c/248e-9      # ≈ 7.593e15  (KrF 248 nm)
omega_XeCl = 2*math.pi*c/308e-9      # ≈ 6.114e15  (XeCl 308 nm)
omega_XeF  = 2*math.pi*c/351e-9      # ≈ 5.368e15  (XeF 351 nm)

# Common solid-state and gas lasers
omega_Nd_3w = 2*math.pi*c/351e-9     # 3ω Nd (matches XeF — useful coincidence)
omega_Ti_Sa = 2*math.pi*c/800e-9     # ≈ 2.355e15  (Ti:Sapphire)
omega_Nd    = 2*math.pi*c/1064e-9    # ≈ 1.770e15  (Nd:YAG fundamental)
omega_CO2   = 2*math.pi*c/10.6e-6    # ≈ 1.777e14  (CO₂)

Main(
    # ...
    reference_angular_frequency_SI = omega_KrF,
)
```

## 4. Critical density, plasma frequency, skin depth

The critical density at the reference frequency is **always** 1 in normalized units, but its SI value scales as `1/λ²`:

| λ | n_c (cm⁻³) | n_c (m⁻³) |
|---|---|---|
| 10.6 μm (CO₂) | 9.9 × 10¹⁸ | 9.9 × 10²⁴ |
| 1.064 μm (Nd) | 9.7 × 10²⁰ | 9.7 × 10²⁶ |
| 800 nm (Ti:Sa) | 1.7 × 10²¹ | 1.7 × 10²⁷ |
| 351 nm (XeF / 3ω Nd) | 8.9 × 10²¹ | 8.9 × 10²⁷ |
| 308 nm (XeCl) | 1.1 × 10²² | 1.1 × 10²⁸ |
| 248 nm (KrF) | 1.8 × 10²² | 1.8 × 10²⁸ |
| 193 nm (ArF) | 2.9 × 10²² | 2.9 × 10²⁸ |
| 157 nm (F₂) | 4.5 × 10²² | 4.5 × 10²⁸ |

In normalized units:
- Electron plasma frequency: `ω_pe / ω_r = sqrt(n_e / n_c) = sqrt(number_density)` for non-relativistic, fully-ionized plasma
- Electron skin depth: `c/ω_pe = 1 / sqrt(number_density)` in units of `L_r`
- Debye length: `λ_D / L_r = sqrt(T_e / number_density)` where `T_e` is in units of `m_e c²`

**The resolution rule of thumb** for collisionless laser-plasma: `cell_length ≤ λ/16` AND `cell_length ≤ 0.5 × λ_D` AND `cell_length ≤ 0.3 × c/ω_pe`. For collisional regimes, the Debye constraint typically dominates and can be much tighter than the laser-resolution constraint.

## 5. Converting SI inputs to normalized values inside the namelist

The cleanest pattern: keep the user's physical parameters at the top of the namelist in SI, then derive the normalized values immediately. This makes the choices visible and reviewable.

```python
import math

# ============================================================
# PHYSICAL PARAMETERS (SI) — edit these
# ============================================================
lambda_SI      = 248e-9          # laser wavelength, m
intensity_SI   = 1.0e15          # peak intensity, W/cm² → convert below
pulse_FWHM_SI  = 50e-15          # pulse FWHM duration, s
waist_SI       = 5e-6            # 1/e² waist radius, m
n_e_SI         = 5e21            # electron density, cm⁻³
T_e_SI         = 100.0           # electron temperature, eV
T_i_SI         = 10.0            # ion temperature, eV
box_x_SI       = 30e-6           # x-extent, m
box_y_SI       = 20e-6           # y-extent, m
sim_time_SI    = 200e-15         # total simulation time, s

# ============================================================
# DERIVED REFERENCE QUANTITIES — do not edit
# ============================================================
c              = 2.99792458e8
e              = 1.602176634e-19
m_e            = 9.1093837015e-31
eps0           = 8.8541878128e-12

omega_r        = 2*math.pi*c/lambda_SI
L_r            = c/omega_r                                # ≈ λ/(2π)
T_r            = 1/omega_r
N_r_m3         = eps0*m_e*omega_r**2/e**2                 # critical density in m⁻³
N_r_cm3        = N_r_m3 * 1e-6
E_r            = m_e*c*omega_r/e                          # V/m

# a₀ from intensity (linear polarization, in W/cm²):
a0             = 0.85 * math.sqrt(intensity_SI*1e-18) * (lambda_SI*1e6)

# ============================================================
# NORMALIZED VALUES
# ============================================================
n_e_norm       = (n_e_SI*1e6) / N_r_m3                    # cm⁻³ → m⁻³ → /N_r
T_e_norm       = T_e_SI / 511e3                           # eV → m_e c² units
T_i_norm       = T_i_SI / 511e3
pulse_norm     = pulse_FWHM_SI * omega_r
waist_norm     = waist_SI / L_r
box_x_norm     = box_x_SI / L_r
box_y_norm     = box_y_SI / L_r
sim_t_norm     = sim_time_SI * omega_r

# ============================================================
# Main block uses the derived normalized values
# ============================================================
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [2*math.pi/16, 2*math.pi/16],
    grid_length        = [box_x_norm, box_y_norm],
    number_of_patches  = [16, 8],
    simulation_time    = sim_t_norm,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_r,
    print_every = 200,
)
```

This pattern scales: a parameter scan becomes a wrapper script that edits the top block, and the namelist body is unchanged.

## 6. Converting normalized outputs to SI in post-processing

happi handles this through the `units` argument. The simulation HDF5 files are always in normalized units — happi reads `omega_r` from the namelist (which is why hard-rule 1 matters here too) and converts on plot:

```python
import happi
S = happi.Open("/path/to/simulation")

# Default: normalized units
S.Probe(0, "Ex").slide()

# Convert: x in μm, t in fs, E in GV/m
S.Probe(0, "Ex", units=["um", "fs", "GV/m"]).slide()

# Mixed units fine; happi uses Pint internally
S.Field(0, "Rho_electron", units=["um", "fs", "cm**-3"]).slide()

# To get the conversion factor explicitly:
omega_r_SI = S.namelist.Main.reference_angular_frequency_SI
L_r_SI     = 2.99792458e8 / omega_r_SI                    # meters
N_r_SI     = 8.854e-12 * 9.11e-31 * omega_r_SI**2 / 1.602e-19**2  # m⁻³
```

If `reference_angular_frequency_SI` was not set in the namelist, happi will refuse to convert and will raise. The fix is to set it (in the namelist, then re-run) or pass `reference_angular_frequency_SI` to `happi.Open()` directly:

```python
S = happi.Open("...", reference_angular_frequency_SI = 2*3.14159*3e8/248e-9)
```

This is a post-hoc fix and only affects post-processing — it does NOT retroactively correct collisions/ionization in the run itself.

## 7. Spatial profiles and the normalization trap

Density and temperature profiles are arbitrary Python callables. They receive coordinates in normalized units and must return values in normalized units. This is where most "my plasma is in the wrong place" bugs originate.

```python
# Coordinates passed to the callable are in c/omega_r.
# A 5 μm waist at 248 nm: L_r = 248e-9/(2π) m ≈ 39.5 nm
# → 5 μm in normalized units ≈ 127

def density_profile(x, y):
    # x, y are in units of L_r
    if 20.0 < x < 40.0:                          # ~785 nm to ~1570 nm — probably not what you want
        return 0.1                               # 0.1 * n_c
    return 0.0
```

The right pattern: convert SI distances to normalized once at the top, then use them:

```python
ramp_start_SI = 1e-6    # 1 μm
ramp_end_SI   = 5e-6    # 5 μm
ramp_start    = ramp_start_SI / L_r
ramp_end      = ramp_end_SI   / L_r

def density_profile(x, y):
    if ramp_start < x < ramp_end:
        return 0.1
    return 0.0
```

Smilei also provides built-in profile helpers that handle the boundary logic correctly: `trapezoidal`, `gaussian`, `polygonal`, `cosine`. Prefer these over hand-written callables for standard shapes.

## 8. Temperature, momentum, and energy units

- **Temperature** in `Species(temperature=...)` is in units of `m_e c² = 511 keV`. So `temperature = [1e-3, 1e-3, 1e-3]` is ~511 eV. For low-temperature plasma users this is counterintuitive — note it and convert explicitly.
- **Momentum** is in units of `m_e c`. So `mean_velocity = [0.1, 0., 0.]` is `0.1 c` ≈ 30,000 km/s — a non-relativistic but fast velocity. For a thermal sub-eV electron, this is enormous.
- **Energy** in diagnostic outputs is in `m_e c²` for relativistic species. ParticleBinning with `deposited_quantity = "weight_ekin"` returns kinetic energy in these units.
- **Field amplitudes** in laser blocks: `a0` is dimensionless and is the standard relativistic-laser amplitude.
- **Time profile arguments** like `ttrapezoidal(start=, plateau=, ...)` take time in normalized units (`T_r`). A 50 fs pulse at 248 nm: `50e-15 * omega_KrF ≈ 380`.

## 9. Common unit mistakes (BAD/GOOD pairs)

### Mistake 1: density in cm⁻³

**BAD**:
```python
Species(..., number_density = 5e21)   # 5×10²¹ * n_c — way too dense
```
**GOOD**:
```python
n_e_SI = 5e21                          # cm⁻³
n_e_norm = (n_e_SI*1e6) / N_r_m3       # convert via reference N_r
Species(..., number_density = n_e_norm)
```

### Mistake 2: temperature in eV

**BAD**:
```python
Species(..., temperature = [100., 100., 100.])   # 100 × 511 keV = 51 MeV electrons
```
**GOOD**:
```python
T_e_eV = 100.0
T_e_norm = T_e_eV / 511e3
Species(..., temperature = [T_e_norm]*3)
```

### Mistake 3: position in μm

**BAD**:
```python
LaserGaussian2D(..., focus = [15., 10.])   # 15 c/omega_r ≈ 590 nm at 248 nm — way too close
```
**GOOD**:
```python
focus_x_SI = 15e-6
focus_y_SI = 10e-6
LaserGaussian2D(..., focus = [focus_x_SI/L_r, focus_y_SI/L_r])
```

### Mistake 4: forgetting `reference_angular_frequency_SI` with Collisions

This is hard-rule 1. Collisions silently run with `ω_r = 1 rad/s` (no, it's actually undefined behavior in older versions; in v5+ it raises, but the error is sometimes only at first collision step, well into the run).

### Mistake 5: setting `cell_length` from a habit

**BAD** (mixing wavelengths):
```python
# Was working for 800 nm; now running at 248 nm with same numerical setup
cell_length = [0.125, 0.125]   # 0.125 c/omega_r ≈ 5 nm at 248 nm — OK for laser
# But forgot to check Debye: dx > 0.5*λ_D triggers numerical heating
```
**GOOD**: compute both constraints and take the minimum:
```python
dx_laser = 2*math.pi/16                                # 16 cells per λ
dx_debye = 0.4 * math.sqrt(T_e_norm/n_e_norm)          # 0.4 of Debye, safe margin
cell_length = [min(dx_laser, dx_debye)]*2
```

### Mistake 6: simulation_time in fs

**BAD**:
```python
simulation_time = 100   # 100 / omega_r ≈ 13 fs at 248 nm — probably too short
```
**GOOD**:
```python
sim_time_fs = 200
simulation_time = sim_time_fs * 1e-15 * omega_r
```

## 10. Quick reference card

```
omega_r [rad/s]         = 2π c / λ_SI
T_r     [s]             = 1/omega_r
L_r     [m]             = c/omega_r  ≈ λ/(2π)
N_r     [m⁻³]           = ε₀ m_e ω_r² / e²    (= n_c)
E_r     [V/m]           = m_e c ω_r / e

In the namelist:
- distance:    x_norm = x_SI / L_r
- time:        t_norm = t_SI * omega_r
- density:     n_norm = n_SI[/m³] / N_r
- temperature: T_norm = T_eV / 511e3
- momentum:    p_norm = v_SI / c     (non-relativistic limit)
- a0 from I:   a0 = 0.85 * sqrt(I[W/cm²] * 1e-18) * λ[μm]

Mandatory `reference_angular_frequency_SI` modules:
- Collisions
- CollisionalIonization
- Ionization (any tunnel/from_rate model)
- RadiationReaction / DiagRadiationSpectrum
- MultiphotonBreitWheeler

In happi:
- units=["um","fs","GV/m","cm**-3", ...]  — Pint-style
- happi.Open(..., reference_angular_frequency_SI=...)  — override if namelist missing it
```
