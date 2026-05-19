# Laser configuration

This reference covers every way to inject a laser in a Smilei namelist: the generic `Laser` block, the predefined `LaserGaussian*` helpers, time and space profiles, polarization, envelope-model variants, oblique incidence, and multi-laser setups. Boundary-condition coupling and the relationship between `a₀`, intensity, and wavelength are covered here; for the meaning of `omega = 1` in normalized units see `basics-and-units.md`.

## Table of contents

1. The two ways to specify a laser
2. The generic `Laser` block
3. `LaserGaussian2D`, `LaserGaussianAM`, `LaserGaussian3D`
4. Time envelopes
5. Space profiles
6. Polarization
7. Envelope-model lasers
8. Oblique incidence and offset lasers
9. Multi-laser setups
10. Boundary conditions for lasers
11. `a₀` ↔ intensity conversion
12. Common mistakes

---

## 1. The two ways to specify a laser

Smilei represents a laser as an oscillating boundary condition on the magnetic field at one face of the simulation box. There are two interfaces:

**Generic `Laser`** — direct control of the `By` and `Bz` profiles (in 2D/3D) as functions of space and time. Maximum flexibility, more code to write. Use when no helper covers your case (arbitrary chirp, custom spatial mode, non-Gaussian transverse profile).

**Predefined helpers** — `LaserGaussian2D`, `LaserGaussianAM`, `LaserGaussian3D`, `LaserEnvelopeGaussianAM`, etc. The common cases pre-wired. Use these unless you have a reason not to.

In both cases, the laser is launched from the boundary specified by `box_side`. Only `"silver-muller"` and `"PML"` boundary conditions support laser injection — `"periodic"` and `"reflective"` cannot launch lasers.

## 2. The generic `Laser` block

```python
Laser(
    box_side = "xmin",                          # one of "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"
    space_time_profile = [By_profile, Bz_profile],
)
```

The two profiles are callables of `(y, t)` for `box_side="xmin"` or `"xmax"`, of `(x, t)` for `"ymin"` or `"ymax"`, etc. They return the field value in normalized units. For cylindrical geometry with azimuthal-mode decomposition, use `space_time_profile_AM` instead, which takes a list of (real, imaginary) component pairs per mode.

```python
import math

def By(y, t):
    # transverse Gaussian, longitudinal Gaussian envelope, p-polarized
    y0     = 20.0                               # center, in c/omega_r
    waist  = 5.0
    t_peak = 30.0
    fwhm_t = 20.0
    sigma_y = waist
    sigma_t = fwhm_t/(2*math.sqrt(2*math.log(2)))
    env = math.exp(-(y - y0)**2/(2*sigma_y**2)) * math.exp(-(t - t_peak)**2/(2*sigma_t**2))
    return 0.01 * env * math.sin(t)             # a0 = 0.01

def Bz(y, t):
    return 0.0                                  # linear polarization in y

Laser(
    box_side = "xmin",
    space_time_profile = [By, Bz],
)
```

The generic block is the right choice for:
- Chirped pulses with custom frequency-time relation
- Spatial modes that aren't Hermite-Gaussian (Bessel beams, super-Gaussian)
- Pulses where the user has tabulated `B(y, t)` from a separate calculation
- Two-color or mixed-frequency configurations

For anything else, prefer a predefined helper.

## 3. `LaserGaussian2D`, `LaserGaussianAM`, `LaserGaussian3D`

These are the workhorses. Almost all production namelists use one of them.

### `LaserGaussian2D`

```python
LaserGaussian2D(
    box_side = "xmin",
    a0 = 1.0,                                   # dimensionless laser amplitude
    omega = 1.0,                                # frequency in units of omega_r (= 1 by default)
    focus = [30., 20.],                         # focus position [x_focus, y_focus] in normalized units
    waist = 5.0,                                # 1/e² waist radius at focus
    incidence_angle = 0.,                       # radians, w.r.t. box_side normal
    polarization_phi = 0.,                      # polarization angle (linear if ellipticity=0)
    ellipticity = 0.,                           # 0 = linear, 1 = circular, in between = elliptical
    time_envelope = tgaussian(fwhm = 30.),     # see §4
    phase_offset = 0.,                          # carrier-envelope phase
)
```

Key conventions:
- `focus[0]` is the longitudinal focal-plane position. The laser is launched from `box_side` (e.g., `x=0` for `xmin`) and propagates to focus at `focus[0]`. Pick `focus[0]` such that the Rayleigh range covers the target.
- `waist` is the radius where the intensity drops to `1/e²` of peak (= the e-folding radius of the field).
- `omega = 1.0` corresponds to the reference frequency (the wavelength encoded in `reference_angular_frequency_SI`). Set `omega ≠ 1` only for second-harmonic generation, mixed wavelengths, etc.
- The default `polarization_phi = 0` and `ellipticity = 0` gives a wave linearly polarized along `y` (i.e., `By = 0`, `Bz ≠ 0` for `xmin` launch).
- `incidence_angle` rotates the launch direction in the simulation plane.

### `LaserGaussianAM`

For `AMcylindrical` geometry — laser propagates along `x`, transverse extent decomposed into azimuthal modes.

```python
LaserGaussianAM(
    box_side = "xmin",
    a0 = 2.0,
    omega = 1.0,
    focus = [50.],                              # in AM, focus is a 1-element list (just x_focus)
    waist = 10.0,
    polarization_phi = 0.,
    ellipticity = 0.,
    time_envelope = tgaussian(fwhm = 25.),
)
```

In AM geometry the transverse coordinates (r, θ) replace (y, z), and the focus argument is one-dimensional — the transverse focus position is by construction the axis. The laser is decomposed into azimuthal modes; a linearly-polarized fundamental Gaussian needs modes 0 and 1 (set `number_of_AM = 2` in `Main`).

### `LaserGaussian3D`

```python
LaserGaussian3D(
    box_side = "xmin",
    a0 = 1.0,
    omega = 1.0,
    focus = [30., 20., 20.],                    # x_focus, y_focus, z_focus
    waist = 5.0,
    incidence_angle = [0., 0.],                 # [angle_y, angle_z]
    polarization_phi = 0.,
    ellipticity = 0.,
    time_envelope = tgaussian(fwhm = 30.),
)
```

Same conventions as 2D, with the third axis added. 3D is expensive; use AM cylindrical if azimuthal symmetry holds approximately.

## 4. Time envelopes

The `time_envelope` argument accepts any callable `t → amplitude`, but Smilei provides four shape constructors that cover most needs:

```python
# Gaussian — most common
time_envelope = tgaussian(fwhm = 30., center = 50.)
# fwhm in T_r units; center is the time of peak amplitude

# Trapezoidal — flat top, useful for "long pulse" approximations
time_envelope = ttrapezoidal(start = 0., plateau = 100., slope1 = 10., slope2 = 10.)
# linear ramp up, flat plateau, linear ramp down

# Constant — CW or DC test
time_envelope = tconstant(start = 0.)

# Polygonal — arbitrary piecewise-linear envelope
time_envelope = tpolygonal(
    points = [0., 20., 80., 100.],              # times
    values = [0., 1., 1., 0.],                  # amplitudes
)

# sin² ramp + plateau
time_envelope = tsin2plateau(start = 0., fwhm = 30., plateau = 50.)

# Or a custom callable
def custom_env(t):
    import math
    return math.exp(-((t-50.)/20.)**4)          # super-Gaussian
time_envelope = custom_env
```

For excimer ns pulses simulated as representative windows, `tgaussian` with FWHM matching the resolved window is typical. For ICF-style "shaped pulses," `tpolygonal` lets you reproduce the engineered envelope.

## 5. Space profiles

For the predefined `LaserGaussian*` helpers, the transverse profile is Gaussian by construction (set via `waist`). To use a non-Gaussian transverse profile, fall back to the generic `Laser` block (§2) and provide your own spatial dependence inside `space_time_profile`.

Note: Smilei also provides spatial profile helpers for use in `Species(number_density=...)` etc. — `gaussian`, `polygonal`, `trapezoidal`, `cosine`. These are documented in `species-and-profiles.md` because they're primarily used there, but they can also appear inside custom laser callables.

## 6. Polarization

In `LaserGaussian*` helpers, polarization is controlled by two parameters:

- `polarization_phi` — angle of the linear polarization in the transverse plane, in radians. For `xmin` launch, `phi = 0` aligns E along z (in 2D) or in the (y,z) plane (in 3D).
- `ellipticity` — 0 for linear, 1 for circular. Intermediate values give elliptical.

```python
# Pure linear, E along z
LaserGaussian2D(..., polarization_phi = 0.,     ellipticity = 0.)

# Pure linear, E along y
LaserGaussian2D(..., polarization_phi = math.pi/2, ellipticity = 0.)

# Circular (right-hand for ellipticity > 0)
LaserGaussian2D(..., polarization_phi = 0.,     ellipticity = 1.)
```

For circular polarization at fixed peak intensity, `a₀_circ = a₀_lin / sqrt(2)`. If specifying intensity, set `a₀` accordingly.

For the generic `Laser` block, polarization is controlled implicitly by the relative phases of `By_profile` and `Bz_profile`:

```python
# Linear E_z
space_time_profile = [lambda y, t: 0.0, lambda y, t: amplitude(y, t)]

# Circular
space_time_profile = [
    lambda y, t: amplitude(y, t) * math.cos(t),
    lambda y, t: amplitude(y, t) * math.sin(t),
]
```

## 7. Envelope-model lasers

The envelope solver averages out laser oscillations and lets you use a coarser longitudinal cell. Use when the pulse is many laser periods long and the slowly-varying-envelope approximation holds.

```python
LaserEnvelopeGaussianAM(
    a0 = 1.0,
    omega = 1.0,                                # carrier frequency
    focus = [100.],
    waist = 20.,
    time_envelope = tgaussian(fwhm = 200., center = 300.),
    polarization_phi = 0.,
    envelope_solver = "explicit",               # "explicit" (Benedetti) or "explicit_reduced"
    Envelope_boundary_conditions = [["PML"]],   # absorbing for the envelope at both ends
)
```

```python
LaserEnvelopeGaussian2D(                        # 2D cartesian envelope
    # same arguments as AM version, but focus = [x, y]
    focus = [100., 20.],
)
```

Settings that must accompany an envelope laser:
- `Main.maxwell_solver` is the standard solver for the residual field; the envelope has its own solver chosen via `envelope_solver`.
- `Main.cell_length` along the propagation direction can be much coarser (typically 5–10 cells per envelope width, not per laser period).
- `Main.timestep_over_CFL` may need adjustment based on envelope-solver dispersion.
- For ionization, set `ionization_model = "tunnel_envelope_averaged"` on the relevant species — the DC-ADK rate is wrong when integrated against the envelope (see `ionization.md`).

**Envelope is wrong when:**
- Plasma response varies on laser-period timescales (some resonance phenomena)
- Polarization rotation matters (envelope assumes carrier polarization is fixed)
- Pulse contains few laser cycles (the SVEA breaks down)
- The user wants instantaneous laser-field visualization (envelope `DiagFields` does not include the laser)

## 8. Oblique incidence and offset lasers

For non-normal incidence, two options:

**`incidence_angle` argument** of `LaserGaussian*` rotates the laser within the simulation plane. Simplest if the angle is small and the focus stays inside the box.

**`LaserOffset`** projects a laser focus that sits outside the simulation box back to the boundary. Useful for grazing-incidence or when the natural geometry has the focus far from the launch boundary.

```python
LaserOffset(
    box_side = "xmin",
    space_time_profile = [By_profile, Bz_profile],
    offset = 20.0,                              # distance from box to virtual focus, in normalized units
    extra_envelope = lambda y, t: 1.0,          # additional modulation
    keep_n_strongest_modes = 100,               # truncate angular spectrum
    angle = 0.0,                                # tilt
)
```

This computes the appropriate phase ramp on the boundary so the laser focuses at the offset position. Internally it projects the user's profile through vacuum diffraction.

## 9. Multi-laser setups

Multiple `Laser*` blocks in the same namelist are additive. Each block injects its own contribution; the boundary sums them.

```python
# Pump
LaserGaussian2D(
    box_side = "xmin",
    a0 = 2.0, omega = 1.0, focus = [50., 30.], waist = 8.,
    time_envelope = tgaussian(fwhm = 30., center = 50.),
)

# Probe — same boundary, different timing, different frequency (2ω)
LaserGaussian2D(
    box_side = "xmin",
    a0 = 0.1, omega = 2.0, focus = [50., 30.], waist = 4.,
    time_envelope = tgaussian(fwhm = 15., center = 80.),
)

# Counter-propagating pulse from xmax (e.g., for two-beam SRS or photon-photon studies)
LaserGaussian2D(
    box_side = "xmax",
    a0 = 1.0, omega = 1.0, focus = [50., 30.], waist = 8.,
    time_envelope = tgaussian(fwhm = 30., center = 80.),
)
```

For cross-beam energy transfer studies (e.g., Oudin et al. 2022), launch two or more lasers from the same side with controlled relative angle and frequency offset. Spatial smoothing can be added by perturbing the `space_time_profile` with stochastic phase plates.

## 10. Boundary conditions for lasers

Only `"silver-muller"` and `"PML"` support laser injection. The choice matters:

| BC | Pros | Cons |
|---|---|---|
| `silver-muller` | Cheap, simple, well-tested. Adequate for normal incidence and moderate angles. | Reflection grows with incidence angle. Imperfect absorption of out-going radiation. |
| `PML` | Near-perfect absorption at all angles. Stable for long simulations. | Expensive (extra cells in each PML direction); requires `number_of_pml_cells`. |

```python
# Silver-Müller everywhere (cheap, default)
Main(
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
)

# PML in x (laser direction), reflective in y (e.g., simulating a slab)
Main(
    EM_boundary_conditions = [["PML"], ["reflective"]],
    number_of_pml_cells = [[12, 12], [0, 0]],
)

# PML everywhere — robust but expensive
Main(
    EM_boundary_conditions = [["PML"], ["PML"]],
    number_of_pml_cells = [[12, 12], [12, 12]],
)
```

For excimer/dense-plasma runs where reflected light off the overdense surface needs to leave cleanly, use PML on the laser boundary. 8–12 cells of PML are typical.

## 11. `a₀` ↔ intensity conversion

The dimensionless laser strength `a₀` and the peak intensity `I` are related by:

```
a₀ = (e E) / (m_e c ω) = 0.85 × sqrt(I_18 × λ_μm²)        (linear polarization)
a₀ = (e E) / (m_e c ω) = 0.60 × sqrt(I_18 × λ_μm²)        (circular polarization)
```

where `I_18 = I / (10^18 W/cm²)` and `λ_μm` is wavelength in micrometers.

Worked examples (linear polarization):

| Wavelength | Intensity [W/cm²] | `a₀` |
|---|---|---|
| 248 nm (KrF) | 10¹² | 2.1 × 10⁻³ |
| 248 nm (KrF) | 10¹³ | 6.7 × 10⁻³ |
| 248 nm (KrF) | 10¹⁴ | 2.1 × 10⁻² |
| 193 nm (ArF) | 10¹⁴ | 1.6 × 10⁻² |
| 193 nm (ArF) | 10¹⁶ | 1.6 × 10⁻¹ |
| 800 nm (Ti:Sa) | 10¹⁸ | 0.68 |
| 800 nm (Ti:Sa) | 10²⁰ | 6.8 |
| 1064 nm (Nd) | 10¹⁵ | 2.7 × 10⁻² |

For excimer regimes `a₀ ≪ 1` is the norm. This is sub-relativistic and the ponderomotive force is small compared to collisional and ionization effects. Always set `a₀`, not intensity, in the namelist — convert in the SI-inputs block at the top of the namelist.

## 12. Common mistakes

### Mistake 1: launching from a boundary that doesn't support lasers

**BAD**:
```python
Main(EM_boundary_conditions = [["periodic"], ["silver-muller"]])
LaserGaussian2D(box_side = "xmin", ...)        # periodic on x — laser ignored
```
The laser block parses but injects nothing. No error is raised.

**GOOD**: ensure `box_side` matches a `silver-muller` or `PML` boundary.

### Mistake 2: focus outside the box

**BAD**:
```python
Main(grid_length = [50., 40.])
LaserGaussian2D(focus = [200., 20.], waist = 5.)
```
The laser is technically valid (it'll focus to a virtual point past the box), but the Rayleigh range may be much larger than the simulation, so the laser appears unconverged on entry. This is sometimes intentional (cylindrical-symmetry test cases) but usually a mistake.

**GOOD**: place `focus[0]` inside the box, typically at the target plane.

### Mistake 3: pulse center earlier than rise time

**BAD**:
```python
time_envelope = tgaussian(fwhm = 30., center = 5.)
```
The Gaussian extends back to `t < 0`, so the pulse is already significantly on at simulation start. Causes an abrupt onset and spurious radiation.

**GOOD**: set `center ≥ 3 × fwhm` so the leading edge has decayed to negligible values before `t = 0`.

### Mistake 4: setting `omega ≠ 1` without understanding what it does

`omega` in a `Laser*` block is the carrier frequency in units of `omega_r`. Setting `omega = 2` means the laser oscillates at twice the reference frequency — useful for second-harmonic studies or when the reference is set to a different anchor frequency. Setting it casually to "2.0" because you want a 2 μm laser instead of 1 μm will not produce a 2 μm laser unless `reference_angular_frequency_SI` is also changed.

### Mistake 5: pmu confusing focus + waist + Rayleigh range

The Rayleigh range `z_R = π w₀² / λ` (in SI) or `z_R = waist² / 2` (in `c/omega_r` for `omega = 1`). For a 5-unit waist, `z_R = 12.5` units. The laser converges from waist-times-sqrt(2) (the boundary of the Rayleigh region) over `z_R` distance. If your box is much shorter than `z_R`, the laser is essentially collimated; if much longer, the focus is sharp.

### Mistake 6: forgetting to match polarization in 2D

In 2D Cartesian, `polarization_phi = 0` gives E in the `z` direction (out of plane). `polarization_phi = π/2` gives E in `y` (in plane). The choice matters for what physics the 2D simulation can capture:
- In-plane E (`phi = π/2`): captures ponderomotive force in the simulation plane, x-y particle motion.
- Out-of-plane E (`phi = 0`): the "z-mode," cleanest for laser-plasma instabilities like SRS in slab geometry.

Choose deliberately.

### Mistake 7: setting `time_envelope = tgaussian(fwhm=...)` but no `center`

The default `center` is 0, meaning the Gaussian peaks at `t = 0`. The pulse is already past its peak before the simulation starts. Always set `center` explicitly.

### Mistake 8: envelope laser without setting `Envelope_boundary_conditions`

The envelope is a separate field with its own boundary handling. Omitting `Envelope_boundary_conditions` defaults to reflective, which is almost always wrong — the envelope reflects off the front face and interferes with itself. Use PML on the launch face at minimum.
