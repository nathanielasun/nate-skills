---
name: smilei-namelist-mini
description: Use for Smilei Python namelist authoring, editing, debugging, running, or post-processing with happi. Trigger on "Smilei", `smilei_test`, `happi.Open`, `Main(...)`, `Species`, `Laser`, `LaserGaussian2D`, `Collisions`, `Ionization`, `DiagScalar`, `DiagFields`, `DiagProbe`, `DiagParticleBinning`, or any PIC namelist file. Includes excimer-laser plasma content (ArF 193, KrF 248, XeCl 308, XeF 351, F2 157 nm) where collisions and ionization dominate over the ponderomotive regime Smilei is most often used in. Use this skill instead of guessing — Smilei has version-sensitive details and silent physics errors (missing reference frequency, non-power-of-2 patches, wrong ionization model) only surface deep into runs. For C++ source work, use smilei-cpp-mini instead.
---

# Smilei Namelist Skill (mini tier)

Smilei is a Particle-In-Cell (PIC) code for laser-plasma physics. This skill covers writing Python namelist files, running Smilei, and post-processing with `happi`. Tested against Smilei v5.1.

**Smilei's defaults assume UHI Ti:Sapphire (800 nm, fs, relativistic, collisionless).** For excimer-laser regimes (UV, ns, sub-relativistic, collisional, ionizing), these defaults are wrong. This skill carries both.

## Hard rules — never violate

1. **`reference_angular_frequency_SI` is mandatory** when `Collisions`, `CollisionalIonization`, `Ionization`, or `RadiationReaction` appears. Set it to `2*math.pi*c/lambda_SI`. Missing → silent order-of-magnitude collision-rate error.

2. **`number_of_patches` must be a power of 2 per direction.** `[16, 8]` valid; `[12, 8]` not. Aim for 4–16 patches per MPI rank.

3. **CFL ≤ 1 for Yee solver.** Use `timestep_over_CFL = 0.95` for margin.

4. **At UV wavelengths (< 400 nm), resolve at actual `n_c`, not at 1 μm.** `n_c ∝ 1/λ²` — at 193 nm `n_c ≈ 30×` higher than at 1 μm. Cell must resolve `c/ω_pe` at peak density.

5. **`number_density = 0.1` means `0.1 × n_c`, not `0.1 cm⁻³`.** Convert SI explicitly.

6. **Define electron `Species` before referencing it via `ionization_electrons`.**

7. **Enable `LoadBalancing(every=50)` for any collisional or laser-driven run.**

8. **Don't use `time_average > 1` in `DiagFields` for instantaneous quantities** (laser, polarization). It smears physics.

9. **Test density-profile callables standalone before plugging in.** Negative or NaN return → undefined behavior.

10. **Always run `smilei_test <namelist>.py` first.** Catches ~90% of errors in seconds.

## Normalization scheme

Smilei uses natural units. With `ω_r` chosen via the laser wavelength:
- Time: `T_r = 1/ω_r`
- Length: `L_r = c/ω_r` (≈ `λ/(2π)`)
- Density: `N_r = ε₀ m_e ω_r² / e²` = **the critical density** at that wavelength
- Field: `E_r = m_e c ω_r / e`
- Energy: `m_e c²` = 511 keV

**Conversions:**
- distance: `x_norm = x_SI / L_r`
- time: `t_norm = t_SI * ω_r`
- density: `n_norm = n_SI[/m³] / N_r`
- temperature: `T_norm = T_eV / 511e3`
- `a₀` from intensity (linear pol.): `a₀ = 0.85 × sqrt(I_W_per_cm² × 1e-18) × λ_μm`

## Excimer-laser quick reference

| λ | Laser | `n_c` [cm⁻³] | `reference_angular_frequency_SI` |
|---|---|---|---|
| 157 nm | F₂ | 4.45e22 | `2*pi*3e8/157e-9` ≈ 1.20e16 |
| 193 nm | ArF | 2.94e22 | `2*pi*3e8/193e-9` ≈ 9.76e15 |
| 248 nm | KrF | 1.78e22 | `2*pi*3e8/248e-9` ≈ 7.59e15 |
| 308 nm | XeCl | 1.15e22 | `2*pi*3e8/308e-9` ≈ 6.11e15 |
| 351 nm | XeF | 8.92e21 | `2*pi*3e8/351e-9` ≈ 5.37e15 |
| 800 nm | Ti:Sa | 1.74e21 | `2*pi*3e8/800e-9` ≈ 2.35e15 |
| 1064 nm | Nd:YAG | 9.69e20 | `2*pi*3e8/1064e-9` ≈ 1.77e15 |

**Excimer regime defaults (λ < 400 nm):**
- Enable `Collisions` for all charged species pairs
- Enable `CollisionalIonization` if neutrals present
- Enable `Ionization` (`tunnel_full_PPT` typically) if `I > 10¹³ W/cm²`
- For pulses > 1 ps, use envelope mode or simulate a representative window
- `cell_length` such that `dx < min(skin_depth, λ_Debye) / 10`
- `particles_per_cell ≥ 32`, ideally 64

## Operating modes — choose one

| Mode | Laser | Plasma | Required physics |
|---|---|---|---|
| **UHI / relativistic** | `a₀ ≳ 1`, fs, NIR | Underdense → overdense | None mandatory |
| **Collisional moderate** (excimer default) | `a₀ ≪ 1`, ns–ps, UV–visible | Often overdense, partially ionized | Collisions + CollisionalIonization + Field Ionization mandatory |
| **Laser-envelope** | Long, slowly varying | Underdense LWFA or long pulse | Envelope solver; `tunnel_envelope_averaged` |
| **Wakefield / LWFA** | `a₀ ≳ 1`, fs, NIR | Underdense, pre-ionized | Moving window |

## Routing — read references when needed

| For details on... | Read |
|---|---|
| Units, lasers, species, profiles, pushers, boundaries | `references/namelist-essentials.md` |
| Collisions, field/collisional ionization, Keldysh regime, excimer recipes | `references/collisions-ionization-excimer.md` |
| Diagnostic blocks and happi post-processing | `references/happi-and-diagnostics.md` |

## The complete minimum template (KrF on Ar — excimer default)

This is a working namelist. Adapt rather than retype.

```python
import math

# === physical scales (SI) ===
lambda_SI = 248e-9                                  # KrF wavelength
omega_SI  = 2*math.pi*3e8/lambda_SI                 # reference frequency

# === Main block ===
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [2*math.pi/16, 2*math.pi/16],   # 16 cells / λ
    grid_length        = [60.,  40.],                    # in units of c/ω_r
    number_of_patches  = [16, 8],
    simulation_time    = 200.,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_SI,           # ← rule 1
    print_every = 200,
)

LoadBalancing(every = 50)                                # ← rule 7

# === Laser ===
LaserGaussian2D(
    box_side  = "xmin",
    a0        = 0.01,                                    # I ~ 4.6e15 W/cm² at 248 nm
    omega     = 1.,
    focus     = [30., 20.],
    waist     = 5.,
    time_envelope = ttrapezoidal(start=0., plateau=80., slope1=10., slope2=10.),
)

# === Electron species ===
Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-4]*3,                              # ~50 eV
    particles_per_cell = 32,                             # ≥ 32 for collisional
    mass = 1.0, charge = -1.0,
    number_density = trapezoidal(0.5, xvacuum=20., xplateau=20., xslope1=2., xslope2=2.),
    boundary_conditions = [["remove"]]*2,
)

# === Ion species ===
Species(
    name = "ion",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-5]*3,
    particles_per_cell = 32,
    mass = 39.948*1836,                                  # Ar mass in m_e
    charge = 1.0,
    number_density = trapezoidal(0.5, xvacuum=20., xplateau=20., xslope1=2., xslope2=2.),
    boundary_conditions = [["remove"]]*2,
    atomic_number = 18,                                  # for ionization
)

# === Collisions (rule 1 already satisfied) ===
Collisions(species1=["electron"], species2=["ion"],     coulomb_log = 0.)
Collisions(species1=["electron"], species2=["electron"], coulomb_log = 0.)
Collisions(species1=["ion"],      species2=["ion"],     coulomb_log = 0.)

# === Diagnostics ===
DiagScalar(every = 50)
DiagFields(every = 200, fields = ["Ex","Ey","Bz","Rho_electron","Rho_ion"])
DiagProbe(every = 100, origin=[0.,20.], corners=[[60.,20.]], number=[300],
          fields=["Ex","Ey","Bz","Rho"])
DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = 200, species = ["electron"],
    axes = [["x", 0., 60., 200], ["px", -0.5, 0.5, 100]],
)
```

## SI → normalized in the namelist — top-of-file pattern

```python
import math

# === physical inputs ===
lambda_SI    = 248e-9
intensity_SI = 1.0e13                  # W/cm²
n_e_SI       = 5e21                    # cm⁻³
T_e_eV       = 50.0
waist_SI     = 5e-6
pulse_FWHM_SI = 100e-15

# === derived ===
c, e, m_e, eps0 = 2.99792458e8, 1.602e-19, 9.11e-31, 8.85e-12
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r     = eps0*m_e*omega_r**2/e**2

# === normalized values ===
a0    = 0.85 * math.sqrt(intensity_SI*1e-18) * (lambda_SI*1e6)
n_e   = (n_e_SI*1e6) / N_r
T_e   = T_e_eV / 511e3
waist = waist_SI / L_r
pulse = pulse_FWHM_SI * omega_r
```

## Workflow

**New namelist:** identify operating mode → set `reference_angular_frequency_SI` (rule 1) → start from the template above → adjust species/laser/diagnostics → check all 10 hard rules → run `smilei_test`.

**Debug a non-running namelist:** check hard rules 1–3 first (~80% of bugs). Then 5, 9. Then `smilei_test` output.

**Wrong physics:** check operating-mode match (collisions disabled in collisional regime? ADK at high `γ_K`?) → re-verify rule 1 and 5.

**Post-process:** use `happi` (see `happi-and-diagnostics.md`).

## Running Smilei

```bash
./smilei_test my_namelist.py             # test parse only
./smilei my_namelist.py                  # single rank

export OMP_NUM_THREADS=4
mpirun -np 4 ./smilei my_namelist.py     # 4 MPI × 4 OMP = 16 cores
```

On an N-core node, choose `MPI × OMP = N`. Patch count rule: `≥ 4 × MPI`, power of 2 per direction.

## Common mistakes

- **Density in cm⁻³**: `number_density = 5e21` is 5×10²¹ × n_c. Convert: `n_norm = n_SI[/m³] / N_r`.
- **Temperature in eV**: `temperature = [100.]*3` is 100 × 511 keV = 51 MeV. Convert: `T_norm = T_eV / 511e3`.
- **Position in μm**: `focus = [15., 10.]` at 248 nm is ~590 nm. Convert: `[x_SI/L_r, y_SI/L_r]`.
- **Missing `reference_angular_frequency_SI` with `Collisions`**: silent order-of-magnitude error. Rule 1.
- **`tgaussian` without `center`**: pulse already past peak at t=0. Set `center ≥ 3 × fwhm`.
- **Non-power-of-2 patches**: opaque error. Use `[2^a, 2^b]`.
- **Laser from periodic boundary**: silently ignored. Use `silver-muller` or `PML`.
- **ADK at UV `γ_K > 1`**: underestimates rate by orders of magnitude. Use `tunnel_full_PPT`.
- **PPC = 8 with collisions**: rates noisy and biased. Use ≥ 32.
- **Only e-i collisions**: e-e drives thermalization too. Enable all three (e-e, e-i, i-i).

## Closing notes

- When given intensity in W/cm², convert to `a₀` explicitly: `a₀ = 0.85 × sqrt(I_18 × λ_μm²)`. State the conversion.
- For collisional/ionization-heavy regimes, read `collisions-ionization-excimer.md` before answering.
- For unusual physics regimes, suggest the Smilei publications page (smileipic.github.io/Smilei/Overview/material.html) as the precedent source.
- For Smilei C++ work, defer to `smilei-cpp-mini`.
