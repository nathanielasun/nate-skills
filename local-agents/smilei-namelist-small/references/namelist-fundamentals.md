# Namelist fundamentals

Covers Smilei normalization and units, laser configuration, and species/profile setup. The three foundational namelist topics merged for the small tier.

## Table of contents

1. Normalization scheme — the reference quantities
2. `reference_angular_frequency_SI` — when mandatory
3. Critical density and resolution constraints
4. Converting SI → normalized inside the namelist
5. Spatial profiles and the normalization trap
6. Lasers — two interfaces
7. `LaserGaussian2D` / `LaserGaussianAM` / `LaserGaussian3D`
8. Time envelopes and polarization
9. Envelope-model lasers
10. `Species` block at a glance
11. Position and momentum initialization
12. Density profiles
13. Boundary conditions and pushers
14. Common mistakes

---

## 1. Normalization scheme — the reference quantities

Smilei uses Vlasov-Maxwell-natural units. The user picks `ω_r` implicitly via the laser wavelength. The derived units:

| Quantity | Symbol | Definition |
|---|---|---|
| Time | `T_r` | `1/ω_r` |
| Length | `L_r` | `c/ω_r` (≈ `λ/(2π)`) |
| Velocity | — | `c` |
| Mass | — | `m_e` |
| Charge | — | `e` |
| Electric field | `E_r` | `m_e c ω_r / e` |
| Magnetic field | `B_r` | `m_e ω_r / e` |
| Number density | `N_r` | `ε₀ m_e ω_r² / e²` (= `n_c(ω_r)`) |
| Energy | — | `m_e c²` (= 511 keV) |

Key identity: **`N_r = n_c`** — the reference density is the critical density at the chosen frequency. `number_density = 0.5` always means `0.5 × n_c`.

## 2. `reference_angular_frequency_SI` — when mandatory

**Optional** for pure Maxwell-Vlasov (laser-plasma without collisions/ionization/QED).

**Mandatory** when any of these appear:
- `Collisions`, `CollisionalIonization`
- `Species(ionization_model="tunnel" | "tunnel_envelope_averaged" | "tunnel_full_PPT" | "tunnel_BSI" | "from_rate")`
- `RadiationReaction`, `DiagRadiationSpectrum`
- `MultiphotonBreitWheeler`

It only affects modules that break scale invariance. Simulation output stays in normalized units regardless; conversion to SI happens at post-processing.

```python
import math
c = 2.99792458e8

omega_F2   = 2*math.pi*c/157e-9      # ≈ 1.200e16
omega_ArF  = 2*math.pi*c/193e-9      # ≈ 9.760e15
omega_KrF  = 2*math.pi*c/248e-9      # ≈ 7.593e15
omega_XeCl = 2*math.pi*c/308e-9      # ≈ 6.114e15
omega_XeF  = 2*math.pi*c/351e-9      # ≈ 5.368e15
omega_TiSa = 2*math.pi*c/800e-9      # ≈ 2.355e15
omega_Nd   = 2*math.pi*c/1064e-9     # ≈ 1.770e15
```

## 3. Critical density and resolution constraints

`n_c` in SI scales as `1/λ²`. For excimer wavelengths, `n_c` is 30× higher than at 1 μm, so resolution must be tighter at the same physical density.

**Three constraints — take the minimum:**

```python
dx_laser = 2*math.pi/16                              # 16 cells per laser λ
dx_skin  = 0.3 * (1/math.sqrt(n_e_norm))             # 0.3 × c/ω_pe at peak n_e
dx_debye = 2.0 * math.sqrt(T_e_norm/n_e_norm)        # 2 × λ_D

cell_length = [min(dx_laser, dx_skin, dx_debye)]*2
```

For excimer regimes at moderate density (10²¹–10²² cm⁻³, 30–100 eV), the **Debye constraint typically wins** — laser-resolution is satisfied automatically.

## 4. Converting SI → normalized inside the namelist

The clean pattern: SI inputs at the top, derived normalized values below. Make choices visible.

```python
import math

# === PHYSICAL INPUTS (SI) ===
lambda_SI      = 248e-9
intensity_SI   = 1.0e13          # W/cm²
pulse_FWHM_SI  = 100e-15
waist_SI       = 5e-6
n_e_SI         = 5e21            # cm⁻³
T_e_eV         = 50.0
box_x_SI       = 30e-6
box_y_SI       = 20e-6
sim_time_SI    = 400e-15

# === DERIVED ===
c, e, m_e, eps0 = 2.99792458e8, 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r_SI  = eps0*m_e*omega_r**2/e**2

a0      = 0.85 * math.sqrt(intensity_SI*1e-18) * (lambda_SI*1e6)
n_e     = (n_e_SI*1e6) / N_r_SI
T_e     = T_e_eV / 511e3
pulse   = pulse_FWHM_SI * omega_r
waist   = waist_SI / L_r
Lx      = box_x_SI / L_r
Ly      = box_y_SI / L_r
T_tot   = sim_time_SI * omega_r

# === Main block uses the derived values ===
Main(
    geometry = "2Dcartesian",
    cell_length = [2*math.pi/16]*2,
    grid_length = [Lx, Ly],
    timestep_over_CFL = 0.95,
    number_of_patches = [16, 8],
    simulation_time = T_tot,
    EM_boundary_conditions = [["silver-muller"]]*2,
    reference_angular_frequency_SI = omega_r,
)
```

**Quick conversions:**
- distance: `x_norm = x_SI / L_r`
- time: `t_norm = t_SI * omega_r`
- density: `n_norm = n_SI[/m³] / N_r`
- temperature: `T_norm = T_eV / 511e3`
- momentum (non-rel): `p_norm = v_SI / c`

## 5. Spatial profiles and the normalization trap

Profile callables receive coordinates in `c/ω_r` and return density in `n_c`. Pre-compute SI distances:

```python
ramp_start = 1e-6 / L_r            # 1 μm
ramp_end   = 5e-6 / L_r            # 5 μm

def n_profile(x, y):
    if ramp_start < x < ramp_end:
        return 0.1                  # 0.1 × n_c
    return 0.0
```

Smilei provides built-in helpers: `trapezoidal`, `gaussian`, `polygonal`, `cosine`. Prefer these over hand-written callables for standard shapes:

```python
number_density = trapezoidal(
    0.1,                            # plateau
    xvacuum=5., xplateau=50., xslope1=2., xslope2=2.,
)

number_density = gaussian(
    0.1, xfwhm=10., xcenter=25., xorder=1,
)
```

## 6. Lasers — two interfaces

**Generic `Laser`**: direct control of `By(y,t)` and `Bz(y,t)` at the boundary. Use for arbitrary chirp, non-Gaussian transverse modes, tabulated profiles.

**Predefined helpers** (`LaserGaussian2D` etc.): the common cases pre-wired. Use these unless you have a reason not to.

Both inject from a boundary specified by `box_side` — only `"silver-muller"` or `"PML"` faces support laser injection.

## 7. `LaserGaussian2D` / `LaserGaussianAM` / `LaserGaussian3D`

```python
LaserGaussian2D(
    box_side = "xmin",
    a0 = 1.0,
    omega = 1.0,                                # frequency in ω_r units
    focus = [30., 20.],                         # [x_focus, y_focus] in c/ω_r
    waist = 5.0,                                # 1/e² waist radius
    incidence_angle = 0.,                       # radians
    polarization_phi = 0.,                      # 0 → E along z (out of 2D plane)
    ellipticity = 0.,                           # 0 linear, 1 circular
    time_envelope = tgaussian(fwhm = 30., center = 50.),
)
```

`LaserGaussianAM` is the cylindrical (r, θ) analogue: focus is 1-element list, transverse focus is the axis. Set `Main.number_of_AM = 2` for linearly polarized fundamental.

`LaserGaussian3D` takes a 3-element focus and a 2-element `incidence_angle = [angle_y, angle_z]`.

## 8. Time envelopes and polarization

```python
time_envelope = tgaussian(fwhm = 30., center = 50.)
time_envelope = ttrapezoidal(start = 0., plateau = 100., slope1 = 10., slope2 = 10.)
time_envelope = tconstant(start = 0.)
time_envelope = tpolygonal(points = [0., 20., 80., 100.], values = [0., 1., 1., 0.])
time_envelope = tsin2plateau(start = 0., fwhm = 30., plateau = 50.)
```

Custom callable also accepted: `lambda t: math.exp(-((t-50.)/20.)**4)`.

**Always set `center` explicitly** for `tgaussian` — default 0 means pulse peaks at t=0 (already past peak when simulation starts). Use `center ≥ 3 × fwhm`.

**Polarization:**
- `polarization_phi = 0` → E along z (out of plane in 2D)
- `polarization_phi = π/2` → E along y (in plane)
- `ellipticity = 0` linear, `1` circular

For circular at fixed peak intensity, `a₀_circ = a₀_lin / sqrt(2)`.

**`a₀` ↔ intensity:**
- Linear: `a₀ = 0.85 × sqrt(I_18 × λ_μm²)`
- Circular: `a₀ = 0.60 × sqrt(I_18 × λ_μm²)`

| Wavelength | Intensity | `a₀` |
|---|---|---|
| 248 nm | 10¹² W/cm² | 2.1 × 10⁻³ |
| 248 nm | 10¹⁴ W/cm² | 2.1 × 10⁻² |
| 193 nm | 10¹⁴ W/cm² | 1.6 × 10⁻² |
| 800 nm | 10¹⁸ W/cm² | 0.68 |
| 800 nm | 10²⁰ W/cm² | 6.8 |

## 9. Envelope-model lasers

Averages out the laser oscillation; allows much coarser longitudinal cells. Use only when:
- Pulse much longer than laser period
- Slowly-varying-envelope approximation holds
- Plasma response is ponderomotive-dominated (caution at excimer collisional intensities)

```python
LaserEnvelopeGaussianAM(
    a0 = 1.0, omega = 1.0,
    focus = [100.], waist = 20.,
    time_envelope = tgaussian(fwhm=200., center=300.),
    polarization_phi = 0.,
    envelope_solver = "explicit",
    Envelope_boundary_conditions = [["PML"]],
)
```

**Must accompany an envelope laser:**
- `Main.cell_length` along propagation can be ~5–10 cells per envelope width (not per laser period)
- For ionization, use `tunnel_envelope_averaged` (see `collisions-and-ionization.md`)
- `Envelope_boundary_conditions` must be set explicitly (default reflective is almost always wrong)

## 10. `Species` block at a glance

```python
Species(
    name = "electron",
    mass = 1.0,                                 # in m_e
    charge = -1.0,                              # in e
    atomic_number = 0,                          # required for ionization

    # Spatial
    number_density = 0.1,
    position_initialization = "regular",
    particles_per_cell = 16,

    # Momentum
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-3]*3,                     # [Tx, Ty, Tz] in m_e c²
    mean_velocity = [0.]*3,

    # Boundary
    boundary_conditions = [["remove"]]*2,

    # Optional / advanced
    pusher = "boris",
    ionization_model = None,
    ionization_electrons = None,
    is_test = False,
)
```

**Common masses:**

| Species | `mass` |
|---|---|
| electron | 1.0 |
| H⁺ / proton | 1836.15 |
| D⁺ | 3671.5 |
| He²⁺ | 7294.3 |
| Ar | 73357 (= 39.948 × 1836) |
| Au | 361688 (= 196.97 × 1836) |
| photon | 0.0 |

**`atomic_number`** is required when `ionization_model` is set or the species participates in `CollisionalIonization`. Smilei uses it to look up ionization potentials from internal atomic tables. Same `atomic_number` for all states in an ionization ladder (e.g., 18 for Ar0, Ar1, ..., Ar18).

## 11. Position and momentum initialization

```python
position_initialization = "regular"    # default; uniform sub-grid in cell (low noise)
position_initialization = "random"     # uniform random
position_initialization = "centered"   # all at cell center (rarely useful)
position_initialization = "ion"        # take from another species (e.g., "ion" → quasi-neutral)
```

```python
momentum_initialization = "maxwell-juettner"   # relativistic Maxwellian, always correct
momentum_initialization = "cold"               # all particles at mean_velocity
momentum_initialization = "rectangular"        # uniform within ±sqrt(3T)
```

**Particles per cell:**

| Use case | PPC |
|---|---|
| Quick test/debug | 4–8 |
| Standard collisionless UHI | 8–16 |
| Standard collisional | 32–64 |
| Excimer collisional + ionization | 64–128 |
| High-fidelity LPI | 64–256 |
| Test-particle subset | 1–4 (with `is_test = True`) |

Collisional runs need higher PPC because the collision module pairs macro-particles within a cell — too few gives noisy and biased rates.

## 12. Density profiles

`number_density` accepts a scalar, a callable, or a built-in helper:

```python
number_density = 0.1                                          # constant
number_density = trapezoidal(0.1, xvacuum=5., xplateau=50., xslope1=2., xslope2=2.)
number_density = gaussian(0.1, xfwhm=10., xcenter=25., xorder=1)
number_density = polygonal(xpoints=[0., 10., 30., 40.], xvalues=[0., 0.1, 0.1, 0.])
number_density = cosine(0.1, xamplitude=0.05, xvacuum=5., xlength=30., xnumber=1)

# Or custom callable
def n_profile(x, y):
    if x < 5.: return 0.0
    if x < 25.: return 0.1
    return 0.0
```

Returning 0.0 → empty cell (no macro-particles, no compute cost). Returning negative → undefined behavior; test profiles standalone first.

**Temperature** in `Species(temperature=...)` is in `m_e c²` (= 511 keV). `temperature = [1e-3]*3` is ~511 eV. Convert eV → normalized: `T_norm = T_eV / 511e3`.

**Mean velocity** in `c` for sub-relativistic, or normalized momentum (γβ) for relativistic.

## 13. Boundary conditions and pushers

```python
boundary_conditions = [
    ["xmin_bc", "xmax_bc"],
    ["ymin_bc", "ymax_bc"],     # third inner list for 3D
]
```

| BC | Behavior |
|---|---|
| `"remove"` | Particle deleted on contact (open boundary) |
| `"reflective"` | Specular reflection |
| `"periodic"` | Wraps to opposite face (both faces of axis must match) |
| `"thermalize"` | Reflected with new thermal momentum; requires `thermal_boundary_temperature` |
| `"stop"` | Particle stays at boundary |

Shortcut: `boundary_conditions = [["remove"]]*2` (2D) or `*3` (3D).

**Pushers:**

| Pusher | Use when |
|---|---|
| `"borisnr"` | Non-relativistic Boris. Faster, valid for `a₀ ≪ 1`. Recommended for excimer regimes. |
| `"boris"` | Default. Standard relativistic Boris. Fine for moderate `a₀`. |
| `"vay"` | High-`γ`. Better E×B drift. |
| `"higueracary"` | Very high-`γ`. Best accuracy. |

For excimer/collisional regimes, `"borisnr"` is appropriate but `"boris"` is safe everywhere.

## 14. Common mistakes

**Density in SI**: `number_density = 5e21` is interpreted as `5×10²¹ × n_c` (way too dense). Convert: `n_norm = n_SI[/m³] / N_r`.

**Temperature in eV**: `temperature = [100.]*3` is 100 × 511 keV = 51 MeV electrons. Convert: `T_norm = T_eV / 511e3`.

**Position in μm**: `focus = [15., 10.]` at 248 nm is ~590 nm (15 × `c/ω_r`). Convert: `focus = [x_SI/L_r, y_SI/L_r]`.

**Missing `reference_angular_frequency_SI` with `Collisions`**: hard rule 1. Silent physics error.

**Pulse center earlier than rise time**: `tgaussian(fwhm=30., center=5.)` already peaks before t=0. Set `center ≥ 3 × fwhm`.

**Periodic on only one face of an axis**: both faces of an axis must be `"periodic"` together. Raises parse error.

**Re-using a species name**: must be unique. For multiple electron populations, use `"bulk_electron"`, `"beam_electron"`, etc.

**Laser launched from a non-supporting face**: `EM_boundary_conditions` must be `"silver-muller"` or `"PML"` on the laser's `box_side`. Periodic → laser silently ignored.

**Patch count not power of 2**: `number_of_patches = [12, 8]` → opaque error at startup. Pick `[2^a, 2^b]`.

**Mismatched ionization-laser pair**: envelope laser + `ionization_model="tunnel"` → wrong rate. Use `"tunnel_envelope_averaged"` (see `collisions-and-ionization.md`).
