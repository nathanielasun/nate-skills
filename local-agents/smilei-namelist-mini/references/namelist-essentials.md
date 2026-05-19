# Namelist essentials

Units, `Main` block, lasers, species, profiles, boundary conditions, and pushers — the foundational namelist topics for the mini tier.

## Units and `Main` block

Smilei uses natural units. `ω_r` is the reference frequency (typically `2π c / λ_laser`).

```python
Main(
    geometry = "1Dcartesian" | "2Dcartesian" | "3Dcartesian" | "AMcylindrical",
    interpolation_order = 2,                     # 2 standard; 4 for high-order accuracy
    timestep_over_CFL = 0.95,                    # ≤ 1 for Yee (default)
    cell_length = [dx, dy],                      # in c/ω_r
    grid_length = [Lx, Ly],                      # in c/ω_r
    number_of_patches = [16, 8],                 # MUST be power of 2 per direction
    simulation_time = T_tot,                     # in 1/ω_r
    EM_boundary_conditions = [["silver-muller"]]*2,
    reference_angular_frequency_SI = omega_r,    # required for Collisions / Ionization
    maxwell_solver = "Yee",                      # default; alternatives: "Lehe", "Cowan", "Bouchard"
    print_every = 200,
    random_seed = 0,
)
```

For AM cylindrical: add `number_of_AM = 2` (or more for higher azimuthal harmonics).

**Patch count rule**: power of 2 per direction, total ≥ 4 × MPI ranks.

## SI → normalized conversions

```python
import math
c, e, m_e, eps0 = 2.99792458e8, 1.602e-19, 9.11e-31, 8.85e-12

omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r     = eps0*m_e*omega_r**2/e**2

x_norm  = x_SI / L_r
t_norm  = t_SI * omega_r
n_norm  = n_SI_per_m3 / N_r              # for cm⁻³, multiply by 1e6 first
T_norm  = T_eV / 511e3
a0      = 0.85 * math.sqrt(I_per_cm2 * 1e-18) * (lambda_SI*1e6)
```

## Lasers

### Two interfaces

**Generic `Laser`**: direct control of `By(y,t)` and `Bz(y,t)`. Use for non-Gaussian profiles, chirp, tabulated.

**Helpers** (`LaserGaussian2D`, etc.): pre-wired Gaussian. Use these unless you have a reason.

Only `silver-muller` or `PML` boundaries support laser injection.

### `LaserGaussian2D` (most common)

```python
LaserGaussian2D(
    box_side = "xmin",                       # "xmin"/"xmax"/"ymin"/"ymax"
    a0 = 1.0,                                # dimensionless amplitude
    omega = 1.0,                             # frequency in ω_r units
    focus = [30., 20.],                      # [x_focus, y_focus] in c/ω_r
    waist = 5.0,                             # 1/e² waist radius
    incidence_angle = 0.,                    # radians
    polarization_phi = 0.,                   # 0 → E along z (out of 2D plane)
    ellipticity = 0.,                        # 0 linear, 1 circular
    time_envelope = tgaussian(fwhm = 30., center = 50.),
)
```

`LaserGaussianAM` (cylindrical): `focus = [x_focus]` (1 element); requires `Main.number_of_AM ≥ 2`.

`LaserGaussian3D`: `focus = [x,y,z]`, `incidence_angle = [angle_y, angle_z]`.

### Time envelopes

```python
time_envelope = tgaussian(fwhm = 30., center = 50.)
time_envelope = ttrapezoidal(start=0., plateau=100., slope1=10., slope2=10.)
time_envelope = tconstant(start = 0.)
time_envelope = tpolygonal(points=[0., 20., 80., 100.], values=[0., 1., 1., 0.])
time_envelope = tsin2plateau(start=0., fwhm=30., plateau=50.)
```

**Always set `center` explicitly for `tgaussian`** — default 0 means pulse peaks at t=0 (past peak when sim starts). Use `center ≥ 3 × fwhm`.

### Polarization

- `polarization_phi = 0` → E along z (out of 2D plane)
- `polarization_phi = π/2` → E along y (in plane)
- `ellipticity = 0` linear, `1` circular
- For circular at fixed peak intensity, `a₀_circ = a₀_lin / sqrt(2)`

### `a₀` ↔ intensity (linear polarization)

```
a₀ = 0.85 × sqrt(I_18 × λ_μm²)
```

where `I_18 = I/(10¹⁸ W/cm²)` and `λ_μm` is wavelength in μm.

| Wavelength | Intensity | `a₀` |
|---|---|---|
| 248 nm | 10¹² W/cm² | 2.1e-3 |
| 248 nm | 10¹⁴ W/cm² | 2.1e-2 |
| 193 nm | 10¹⁴ W/cm² | 1.6e-2 |
| 800 nm | 10¹⁸ W/cm² | 0.68 |
| 800 nm | 10²⁰ W/cm² | 6.8 |

### Envelope-model lasers

For slowly-varying pulses where the carrier oscillation can be averaged:

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

Use only when: pulse much longer than laser period, SVEA holds, ponderomotive-dominated, not laser-cycle-resolved physics. For ionization with envelope laser, use `tunnel_envelope_averaged`.

## Species

```python
Species(
    name = "electron",
    mass = 1.0,                                  # in m_e
    charge = -1.0,                               # in e
    atomic_number = 0,                           # required for ionization

    # Spatial
    number_density = 0.1,
    position_initialization = "regular",         # or "random", "centered", "ion" (copy from another species)
    particles_per_cell = 16,

    # Momentum
    momentum_initialization = "maxwell-juettner", # or "cold", "rectangular"
    temperature = [1e-3]*3,                      # [Tx, Ty, Tz] in m_e c²
    mean_velocity = [0.]*3,                      # in c (sub-rel) or γβ (rel)

    # Boundary
    boundary_conditions = [["remove"]]*2,

    # Pusher
    pusher = "boris",                            # or "borisnr" (non-rel, faster for a₀ << 1)

    # Optional
    is_test = False,
    ionization_model = None,
    ionization_electrons = None,
)
```

### Common masses

| Species | `mass` |
|---|---|
| electron | 1.0 |
| H⁺ | 1836.15 |
| D⁺ | 3671.5 |
| He²⁺ | 7294.3 |
| Ar | 73357 |
| Au | 361688 |
| photon | 0.0 |

### Particles per cell (PPC)

| Use case | PPC |
|---|---|
| Quick test / debug | 4–8 |
| Standard collisionless UHI | 8–16 |
| Standard collisional | 32–64 |
| Excimer collisional + ionization | 64–128 |

Collisional runs need higher PPC because the collision module pairs particles within a cell. Too few → noisy and biased rates.

## Density profiles

```python
number_density = 0.1                                          # constant
number_density = trapezoidal(0.1, xvacuum=5., xplateau=50., xslope1=2., xslope2=2.)
number_density = gaussian(0.1, xfwhm=10., xcenter=25., xorder=1)
number_density = polygonal(xpoints=[0.,10.,30.,40.], xvalues=[0.,0.1,0.1,0.])
number_density = cosine(0.1, xamplitude=0.05, xvacuum=5., xlength=30., xnumber=1)

# Or custom callable — receives (x), (x,y), or (x,y,z) in c/ω_r, returns n_norm
def n_profile(x, y):
    if x < 5.: return 0.0
    if x < 25.: return 0.1
    return 0.0
```

Return 0.0 → empty cell (no compute cost). Return negative → undefined behavior; test profiles standalone first.

## Boundary conditions

```python
boundary_conditions = [
    ["xmin_bc", "xmax_bc"],
    ["ymin_bc", "ymax_bc"],                # third for 3D
]
```

| BC | Behavior |
|---|---|
| `"remove"` | Particle deleted (open boundary) |
| `"reflective"` | Specular reflection |
| `"periodic"` | Wraps to opposite face (both faces of axis must match) |
| `"thermalize"` | Reflected with new thermal momentum (set `thermal_boundary_temperature`) |
| `"stop"` | Particle stays at boundary |

Shortcut: `[["remove"]]*2` for 2D or `*3` for 3D.

For periodic-y with absorbing-x: `[["remove"], ["periodic"]]`.

## Pushers

| Pusher | Use when |
|---|---|
| `"borisnr"` | Non-relativistic Boris. Fastest for `a₀ ≪ 1`. **Recommended for excimer regimes.** |
| `"boris"` | Standard Boris. Default. Safe for any `a₀`. |
| `"vay"` | High-γ. Better E×B drift. |
| `"higueracary"` | Very high-γ. Best accuracy. |

## Common mistakes summary

- Density in cm⁻³: convert to `n_norm = n_SI[/m³] / N_r`
- Temperature in eV: convert to `T_norm = T_eV / 511e3`
- Position in μm: convert to `x_norm = x_SI / L_r`
- Missing `reference_angular_frequency_SI` with Collisions/Ionization (rule 1)
- `tgaussian` without `center` → pulse past peak at t=0
- Non-power-of-2 patches → opaque startup error
- Laser launched from periodic boundary → silently ignored
- Re-using species names → conflict
- Periodic on only one face → parse error (both faces of axis must match)
