# Species and profiles

This reference covers the `Species` block: identity (mass/charge/atomic_number), position and momentum initialization, density and temperature profiles, particle pushers, boundary conditions, and particle injection. Ionization is covered briefly here (just enough to know which fields belong on a `Species` block) — full treatment is in `ionization.md`. Collisions are covered in `collisions-and-reactions.md`.

## Table of contents

1. The `Species` block at a glance
2. Identity: name, mass, charge, atomic_number
3. Position initialization
4. Momentum initialization
5. `particles_per_cell` and the convergence question
6. Density profiles
7. Temperature
8. Mean velocity
9. Boundary conditions
10. Particle pushers
11. `ParticleInjector` — continuous injection
12. Relativistic field initialization
13. Ionization-related fields (cross-reference)
14. Common mistakes

---

## 1. The `Species` block at a glance

```python
Species(
    name = "electron",                            # used to reference this species elsewhere
    mass = 1.0,                                   # in m_e
    charge = -1.0,                                # in e
    atomic_number = 0,                            # only for ionization

    # Spatial layout
    number_density = 0.1,                         # constant, callable, or profile helper
    position_initialization = "regular",          # or "random", "centered"
    particles_per_cell = 16,

    # Momentum layout
    momentum_initialization = "maxwell-juettner", # or "cold", "rectangular"
    temperature = [1e-3]*3,                       # [Tx, Ty, Tz] in m_e c²
    mean_velocity = [0.]*3,                       # [vx, vy, vz] in c

    # Boundaries
    boundary_conditions = [["remove"]]*2,         # one entry per axis

    # Numerics
    pusher = "boris",                             # default

    # Optional / advanced
    ionization_model = None,
    ionization_electrons = None,
    is_test = False,                              # test particles don't feedback to fields
    track_every = 0,
)
```

Most fields have sensible defaults. The mandatory ones in a basic run are `name`, `mass`, `charge`, `number_density`, `position_initialization`, `momentum_initialization`, and `boundary_conditions`.

## 2. Identity: name, mass, charge, atomic_number

- **`name`** is the user-facing label. It appears in diagnostic filenames (`Rho_<name>`, `Jx_<name>`), is the string passed to `Collisions(species1=[...])`, and identifies the species in `ParticleBinning(species=[...])`. Names are case-sensitive and must be unique. Conventional names: `"electron"`, `"ion"`, `"Ar1"`, `"H+"`, `"proton"`, `"positron"`.

- **`mass`** is in units of electron mass. Some common values:

| Species | `mass` |
|---|---|
| electron | 1.0 |
| positron | 1.0 |
| proton, H⁺ | 1836.15 |
| deuteron, D⁺ | 3671.5 |
| He²⁺ | 7294.3 |
| Ar | 39.948 × 1836.15 ≈ 73357 |
| Au | 196.97 × 1836.15 ≈ 361688 |
| photon | 0.0 |

For arbitrary ion isotopes, `mass = atomic_mass * 1836.15` where atomic_mass is in atomic mass units (1 u ≈ 1836 m_e is a useful shortcut). Photons (mass = 0) are special — they cannot be pushed by the Lorentz force and exist only as targets for QED processes.

- **`charge`** is in units of `e` and is the *initial* charge state. For species that will be ionized, this is the starting state and Smilei tracks the running charge per macro-particle. For ions explicitly set up with multiple charge states, define separate `Species` blocks for `Ar1`, `Ar2`, `Ar3`, etc.

- **`atomic_number`** is the nuclear charge Z. Required when `ionization_model` is set or when this species participates in `CollisionalIonization`. Smilei uses it to look up ionization potentials from internal atomic data tables. For ions in a ladder, set `atomic_number` to the actual element (e.g., 18 for argon at every stage `Ar1`, `Ar2`, ..., `Ar18`).

## 3. Position initialization

Macro-particle positions inside each cell. Three options:

```python
position_initialization = "regular"        # default; uniform grid inside each cell
position_initialization = "random"         # uniform random; reseed with Main.random_seed
position_initialization = "centered"       # all macro-particles at cell center (rarely useful)
```

`regular` is the default and the right choice for most physics. It minimizes Poisson noise in initial density. `random` is appropriate when the user wants stochasticity, or when `regular` would produce artificial patterns (e.g., diffraction off the grid).

You can also seed positions from another species:

```python
position_initialization = "ion"            # take positions from species named "ion"
```

This guarantees charge neutrality at `t = 0` — useful for quasi-neutral cold plasmas. Note: the source species must be defined earlier in the namelist.

### `regular_number` — controlling the regular grid

By default, `particles_per_cell` particles are arranged as a regular sub-grid (e.g., 4×4 for `particles_per_cell = 16` in 2D). To override the arrangement:

```python
particles_per_cell = 16
regular_number = [4, 4]                    # explicit sub-grid; product must equal particles_per_cell
```

For elongated species distributions, an anisotropic `regular_number` can reduce noise in one direction.

## 4. Momentum initialization

Sets the initial momentum distribution.

```python
momentum_initialization = "maxwell-juettner"   # relativistic Maxwellian
momentum_initialization = "cold"               # all particles have p = mean_velocity
momentum_initialization = "rectangular"        # uniform within ±sqrt(3T) per axis
```

`maxwell-juettner` is the right choice for thermal plasmas at any temperature. For sub-relativistic plasmas it reduces to the non-relativistic Maxwell-Boltzmann distribution; for relativistic temperatures (`T > m_e c²`) it correctly handles the relativistic Maxwellian.

`cold` is appropriate for beams (`mean_velocity ≠ 0`, `temperature = 0`), test particles, and pre-formed plasmas where the temperature is set by subsequent dynamics rather than initial state.

`rectangular` is a simpler distribution with the same variance as Maxwell-Boltzmann; rarely the right physics choice but occasionally useful for debugging.

For loaded test particles (e.g., a witness beam in LWFA), use `cold` and set `mean_velocity` to the design momentum.

## 5. `particles_per_cell` and the convergence question

The number of macro-particles per cell controls statistical noise. The trade-off is convergence vs. compute cost.

| Use case | Recommended PPC |
|---|---|
| Quick test, debugging | 4–8 |
| Standard collisionless UHI | 8–16 |
| Standard collisional run | 32–64 |
| Excimer-class collisional + ionization | 64–128 |
| High-fidelity LPI study | 64–256 |
| Test-particle subset | 1–4 (with `is_test = True`) |

The cost scales linearly with PPC. Memory scales with PPC × number of cells. For 3D runs PPC is usually the first thing to reduce when scaling down to small clusters.

**Collisional runs need more PPC** than collisionless because the collision algorithm pairs macro-particles within a cell; with too few, collision rates have high variance and the binary scheme can produce nonphysical momentum bursts.

**Adaptive PPC** can be set per-species using a callable for `particles_per_cell` (function of position) — useful when some regions need higher resolution than others. Not commonly used; standard practice is uniform PPC per species.

## 6. Density profiles

`number_density` accepts:

- A scalar (constant density everywhere): `number_density = 0.1`
- A Python callable `(x), (x,y), or (x,y,z) → float`: arbitrary profile
- One of Smilei's built-in profile constructors

### Built-in profile helpers

```python
# Trapezoidal — flat top with linear ramps
number_density = trapezoidal(
    0.1,                                # plateau value
    xvacuum = 5.0,                      # length of vacuum region at the start
    xplateau = 50.0,                    # length of flat plateau
    xslope1 = 2.0,                      # length of ramp up
    xslope2 = 2.0,                      # length of ramp down
    yvacuum = 0., yplateau = 1e20, yslope1 = 0., yslope2 = 0.,    # transverse (effectively constant if plateau is huge)
)

# Gaussian
number_density = gaussian(
    0.1,                                # peak value
    xfwhm = 10., xcenter = 25., xorder = 1,        # Gaussian along x
    yfwhm = 10., ycenter = 20., yorder = 1,
)
# xorder > 1 → super-Gaussian

# Polygonal — piecewise-linear
number_density = polygonal(
    xpoints = [0., 10., 30., 40.],
    xvalues = [0., 0.1, 0.1, 0.],
)

# Cosine — half-cosine profile
number_density = cosine(
    0.1,                                # peak value
    xamplitude = 0.05, xvacuum = 5., xlength = 30., xnumber = 1,
)

# Constant — equivalent to scalar, but works with the helper-driven workflow
number_density = constant(0.1)
```

For arbitrary profiles, use a callable:

```python
import math

def n_profile(x, y):
    # Step-up region with exponential decay outside
    if x < 5.: return 0.0
    if x < 15.: return 0.1
    return 0.1 * math.exp(-(x - 15.)/5.)
```

The callable receives coordinates in normalized units (`c/omega_r`) and returns density in normalized units (fraction of `N_r`).

### Density of `0.0`

Returning `0.0` from the profile leaves no macro-particles in that cell. This is the correct way to make a vacuum region — Smilei skips empty cells in particle loops, so vacuum is free.

### Negative density

Smilei does not raise an error on negative density returned by a callable, but the behavior is undefined. Always sanity-check profiles return non-negative values. Test the profile standalone before plugging in:

```python
import numpy as np
import matplotlib.pyplot as plt
xs = np.linspace(0, 50, 500)
plt.plot(xs, [n_profile(x, 25.) for x in xs])
plt.show()
```

## 7. Temperature

```python
temperature = [Tx, Ty, Tz]              # per-axis temperatures in m_e c²
```

Three entries, one per Cartesian axis. Setting `[T, T, T]` gives an isotropic Maxwellian. Anisotropic temperatures (`Tx ≠ Ty`) are physical (Weibel instability seeds, magnetized plasmas) and supported.

Temperatures can also be callables — useful for steep gradients or plasma with spatially varying preheat:

```python
def T_e_profile(x, y):
    return 1e-3 if x > 10. else 1e-5    # hot region at x > 10
temperature = [T_e_profile]*3
```

For `momentum_initialization = "cold"`, `temperature` is ignored (the particles are mono-energetic at `mean_velocity`).

## 8. Mean velocity

```python
mean_velocity = [vx, vy, vz]            # in c (sub-relativistic) or normalized momentum (relativistic)
```

For sub-relativistic species, `mean_velocity` is in units of c — drift velocity directly. For relativistic species (e.g., a witness beam in LWFA), the field actually represents normalized momentum `p/(m c)`, which equals `γβ`. Smilei interprets this consistently as the drift four-velocity component.

Can also be a callable, useful for sheared flows or beam profiles:

```python
def shear_vy(x, y):
    return 0.001 * (x - 20.)            # linear shear
mean_velocity = [None, shear_vy, None]
```

## 9. Boundary conditions

```python
boundary_conditions = [
    ["xmin_bc", "xmax_bc"],
    ["ymin_bc", "ymax_bc"],
    # ["zmin_bc", "zmax_bc"],            # for 3D
]
```

One inner list per axis, two strings per axis (for the two faces). Available options per face:

| BC | Behavior |
|---|---|
| `"remove"` | Particle is deleted on contact (open boundary) |
| `"reflective"` | Specular reflection |
| `"periodic"` | Wraps to opposite face; both faces of the axis must be `"periodic"` |
| `"thermalize"` | Particle is reflected with thermal momentum redrawn from a thermal distribution |
| `"stop"` | Particle stays at the boundary (rarely useful) |

For `"thermalize"`, two extra fields are needed:

```python
Species(
    # ...
    boundary_conditions = [["thermalize", "remove"], ["periodic", "periodic"]],
    thermal_boundary_temperature = [1e-4]*3,
    thermal_boundary_velocity = [0., 0., 0.],
)
```

`"thermalize"` is useful for simulating walls held at fixed temperature (e.g., chamber walls in a gas-discharge simulation).

### Shortcut: same BC on both faces

```python
boundary_conditions = [["remove"]]*2    # 2D: both axes set to ["remove","remove"]
boundary_conditions = [["remove"]]*3    # 3D
```

## 10. Particle pushers

```python
pusher = "boris"          # default; symplectic, second-order, fast
pusher = "vay"            # relativistic, slightly more accurate at high γ
pusher = "higueracary"    # high-accuracy relativistic
pusher = "borisnr"        # non-relativistic Boris; faster, valid for v << c
```

Pusher selection rules:

| Regime | Recommendation |
|---|---|
| `a₀ < 0.1`, no relativistic dynamics expected | `"borisnr"` — faster, no relativistic overhead |
| `0.1 ≤ a₀ ≲ 5` | `"boris"` — default, fine for moderate relativity |
| `a₀ > 5` or pair plasma at high `γ` | `"vay"` or `"higueracary"` — preserve E×B drift at high `γ` |
| Excimer regimes typically | `"borisnr"` is appropriate; `"boris"` works but is slightly slower |
| Photons (mass = 0) | No pusher applies — photons advance via Maxwell only |

### Special pusher: photon pusher

When `mass = 0`, the species is a photon and Smilei advances it ballistically at `c` along its initial momentum vector. Photons are used as targets for radiation reaction (Monte-Carlo photon emission) and multiphoton Breit-Wheeler pair generation. They don't participate in ordinary species dynamics.

## 11. `ParticleInjector` — continuous injection

Standard `Species` blocks initialize particles only at `t = 0`. For continuous injection from a boundary (e.g., a particle beam, a steady gas inflow, a plasma source), use `ParticleInjector`:

```python
ParticleInjector(
    name = "injected_electrons",
    species = "electron",                # links to a Species block defining the type
    box_side = "xmin",
    number_density = 0.05,
    mean_velocity = [0.1, 0., 0.],       # drifting in +x
    temperature = [1e-4]*3,
    time_envelope = tconstant(),         # injection rate vs time
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = 16,
)
```

The injector adds particles at each timestep along the specified boundary. Multiple injectors can target the same species. Useful for:
- Gas inflow in discharge simulations
- Beam injection in beam-plasma codes
- Steady-state plasma source (where boundary thermalization isn't enough)

## 12. Relativistic field initialization

For dense relativistic beams (e.g., a witness beam in LWFA, a self-injected bunch with `γ > 100`), the static fields of the beam must be initialized self-consistently. Smilei provides a relativistic Poisson solver triggered by a flag on the `Species`:

```python
Species(
    name = "witness_beam",
    relativistic_field_initialization = True,
    # ...
)
```

When this is set, Smilei runs a relativistic Poisson solve at `t = 0` to set up the electromagnetic fields co-moving with the beam. Without this, a high-density relativistic species at `t = 0` will radiate spuriously as the fields catch up.

Use only for species with `mean_velocity` close to `c` (`γ > ~5`) and high density. For thermal plasmas at rest, this is unnecessary and the default Poisson initialization (or no initialization at all) is correct.

## 13. Ionization-related fields (cross-reference)

These fields belong on a `Species` block but are documented in `ionization.md`:

```python
Species(
    # ...
    ionization_model = "tunnel",                  # or "tunnel_envelope_averaged", "tunnel_full_PPT", "from_rate"
    ionization_electrons = "electron",            # name of the Species receiving freed electrons
    maximum_charge_state = 18,                    # cap on the ladder
    ionization_rate = None,                       # callable for "from_rate" model
)
```

The neutral species must have `atomic_number > 0`. Each ionization step produces one electron (added to `ionization_electrons` species) and increments the ion's charge state. The same `Species` macro-particle is reused as it climbs the ladder. See `ionization.md` for which model to choose at which intensity, the Keldysh γ threshold, and envelope-averaged considerations.

## 14. Common mistakes

### Mistake 1: density in physical units

**BAD**:
```python
number_density = 5e21                       # interpreted as 5×10²¹ × n_c
```
**GOOD**: see `basics-and-units.md` §5 for the conversion.

### Mistake 2: forgetting to set `boundary_conditions` per axis

**BAD**:
```python
Species(..., boundary_conditions = "remove")        # type error
Species(..., boundary_conditions = [["remove"]])    # only one axis specified in 2D
```
**GOOD**: one inner list per axis, two strings per inner list.

### Mistake 3: mismatched `atomic_number` and `mass` for ions

```python
Species(name = "Ar1", mass = 39.948*1836, charge = 1.0, atomic_number = 18, ...)  # consistent
Species(name = "Ar1", mass = 4.0*1836,    charge = 1.0, atomic_number = 18, ...)  # implausible — He-mass Ar
```
Smilei doesn't validate physical consistency. `atomic_number` is only used for ionization; `mass` controls dynamics. Both must be set correctly.

### Mistake 4: zero `particles_per_cell` in a non-empty cell

**BAD**:
```python
Species(..., particles_per_cell = 0)
```
The simulation runs but the species has no representation. No error.

**GOOD**: minimum is 1; recommended is 8+ for any species that affects fields.

### Mistake 5: `"periodic"` on only one face of an axis

**BAD**:
```python
boundary_conditions = [["periodic", "remove"], ["remove"]*2]
```
Periodicity requires both faces of an axis to be `"periodic"`. Smilei raises a parse error here.

### Mistake 6: forgetting `relativistic_field_initialization` for a beam

A dense `γ = 1000` beam set up without `relativistic_field_initialization = True` will produce a burst of radiation in the first few timesteps as the fields catch up. Often diagnosed as "I see a spurious EM pulse at `t = 0`."

### Mistake 7: re-using a species name

```python
Species(name = "electron", ...)
Species(name = "electron", ...)   # conflict — Smilei takes the second one or errors
```
Species names must be unique. For multiple electron populations, use distinct names: `"bulk_electron"`, `"beam_electron"`, etc.

### Mistake 8: thermalize without temperature

```python
boundary_conditions = [["thermalize", "thermalize"], ...]
# but no thermal_boundary_temperature set
```
The default is zero temperature, which makes thermalize equivalent to a cold reflective wall. Always set `thermal_boundary_temperature` explicitly when using `"thermalize"`.

### Mistake 9: `mean_velocity` as a list of three numbers when one is a callable

```python
mean_velocity = [some_func, 0., 0.]            # mixed types — fragile
```
Either all three are scalars, or all three are callables (use `lambda x, y: 0.` for the zero entries).
