---
name: smilei-namelist
description: Use for writing, editing, debugging, or running a Smilei Python namelist, post-processing with happi, or planning a Smilei PIC simulation. Trigger on "Smilei", `smilei_test`, `happi.Open`, `Main(geometry=...)`, `Species`, `Laser`, `LaserGaussian2D`, `LaserGaussianAM`, `LaserEnvelope`, `Collisions`, `Ionization`, `DiagScalar`, `DiagFields`, `DiagProbe`, `DiagParticleBinning`, `DiagTrackParticles`, or any mention of a PIC namelist. Includes specialized content for excimer-laser plasmas (ArF 193, KrF 248, XeCl 308, XeF 351, F2 157 nm) where collisional and ionization physics dominate over the ponderomotive regime Smilei is most often used in. Use this skill instead of relying on memory — Smilei's API evolves across releases and silent physics errors (missing reference frequency for collisions, non-power-of-2 patches, wrong ionization model for the Keldysh regime) only manifest deep into a run. For Smilei C++ source work, use smilei-cpp instead.
---

# Smilei Namelist Skill (frontier tier)

This skill covers authoring, debugging, and running Smilei Python namelist files, plus post-processing with the `happi` Python module. It does **not** cover modifying Smilei's C++ source (use `smilei-cpp` for that). Tested against Smilei v5.1 (current as of mid-2026); call out version-specific syntax when it appears.

Smilei was originally tuned for ultra-high-intensity (UHI) femtosecond laser-plasma interaction. Many users — especially those working on **excimer-laser plasmas, low-temperature plasmas, gas discharges, and inertial-fusion-energy contexts** — operate in regimes where Smilei's default assumptions (collisionless, ponderomotive-dominated, relativistic) are wrong. This skill carries excimer/LTP content alongside the standard UHI material so the model picks the right regime on the user's first message.

## Hard rules (do not violate)

These rules prevent the most common silent failures. They are repeated in the reference files but stated once here so the model carries them even if no reference is loaded.

1. **`reference_angular_frequency_SI` is mandatory whenever `Collisions`, `CollisionalIonization`, `Ionization`, `RadiationReaction`, or `RadiationSpectrum` appears in the namelist.** Compute it from the relevant laser wavelength: `reference_angular_frequency_SI = 2*math.pi*c/lambda_SI`. Missing this gives no error at parse time but produces collision rates that are off by orders of magnitude. For KrF 248 nm: `2*pi*3e8/248e-9 ≈ 7.59e15`. For ArF 193 nm: `2*pi*3e8/193e-9 ≈ 9.76e15`. For Ti:Sa 800 nm: `2*pi*3e8/800e-9 ≈ 2.35e15`.

2. **`number_of_patches` must be a power of 2 in each direction.** `[16, 8]` is valid; `[12, 8]` is not. Smilei will refuse to start if violated, but the error message is opaque. Pick patch counts so that each MPI rank owns at least one patch and ideally 4–16.

3. **The CFL ceiling depends on the Maxwell solver.** Default `"Yee"`: `timestep_over_CFL ≤ 1`. `"Cowan"`, `"Lehe"`, `"Bouchard"`: stricter. Set `timestep_over_CFL = 0.95` rather than 1.0 by default to leave margin.

4. **For UV wavelengths (< 400 nm), the resolution budget is set by the local skin depth at the working `n_c`, not at 1 μm reference.** `n_c ∝ 1/λ²`, so at 193 nm `n_c ≈ 2.9 × 10²² cm⁻³` — roughly 30× higher than at 1 μm. Cell size must resolve `c/ω_pe` at the *actual* electron density, not at the reference. Under-resolving causes numerical heating that swamps real physics.

5. **Density profiles are evaluated in normalized units.** `number_density = 0.1` means `0.1 × n_c` (the reference critical density), not `0.1 cm⁻³`. When users give SI numbers, convert explicitly inside the namelist using `reference_angular_frequency_SI` and document the conversion.

6. **`Species` ordering matters for `Ionization`.** When `ionization_model` is set on a neutral species, list the electron species that receives ionized electrons via `ionization_electrons = "elec"` BEFORE configuring `Collisions` blocks that involve those electrons. Smilei resolves species references at parse time.

7. **Patch arrangement and load balancing interact.** For laser-driven simulations with a moving feature, enable `LoadBalancing` — patches without it are statically distributed and you will see severe imbalance. Set `every = 20` to 100 typically.

8. **`DiagFields` with `time_average > 1` smooths in time and changes physical meaning.** Use it for envelope-like outputs; do NOT use it for instantaneous field measurements (e.g., laser polarization tracking).

9. **`number_density` callables must return a non-negative float for all `(x, y, z)`.** Returning `NaN` from a Python profile silently produces undefined behavior. Test profiles standalone before plugging into the namelist.

10. **Always run `smilei_test <namelist>.py` before launching a real simulation.** It executes the Python interpretation step without doing the PIC loop and catches ~90% of namelist errors in seconds.

## Four operating modes — pick one explicitly

Smilei is used across very different physical regimes and the right namelist structure differs. Identify the mode on first contact and proceed accordingly.

| Mode | Laser | Plasma | Required physics | Typical `omega_r` source |
|---|---|---|---|---|
| **UHI / relativistic** | a₀ ≳ 1, fs pulse, NIR (~800 nm) | Underdense → overdense | None mandatory; collisions/ionization optional | Laser λ |
| **Collisional moderate-intensity** | a₀ ≪ 1, ns–ps, UV–visible | Often overdense, often partially ionized | Collisions + CollisionalIonization + (Field) Ionization mandatory | Laser λ (excimer wavelength) |
| **Laser-envelope** | Long pulse, slowly varying | LWFA-relevant underdense, or long ns-class pulse | Envelope solver; ionization uses `tunnel_envelope_averaged` | Laser λ |
| **Wakefield / LWFA** | a₀ ≳ 1, fs, NIR | Underdense, pre-ionized | Moving window mandatory; envelope optional | Laser λ |

**Excimer-laser plasmas live in the second row by default.** ArF/KrF/XeCl/XeF lasers produce ns-class pulses at moderate intensities (10⁸ – 10¹³ W/cm²) at UV wavelengths. The plasma is dense (`n_c` at 193 nm is ~2.9 × 10²² cm⁻³), collision-dominated, and partially ionized — Smilei's collisionless default is wrong here. There are exceptions: short-pulse KrF amplifier output (sub-ps) for ICF (NRL Nike/Electra-class) can push into UHI territory, in which case the mode is closer to row 1 but with the higher `n_c` of row 2.

## Routing — when to read which reference

Read the relevant reference file(s) BEFORE answering, not after. Multiple references can apply to one question.

| If the user is asking about… | Read |
|---|---|
| Choosing units, normalizations, or `reference_angular_frequency_SI` | `references/basics-and-units.md` |
| Laser configuration: profile, polarization, focus, envelope, multi-laser | `references/laser-configuration.md` |
| Species definition, density profiles, particle initialization, pushers, boundary conditions | `references/species-and-profiles.md` |
| Collisions, collisional ionization, nuclear reactions, bound-electron screening | `references/collisions-and-reactions.md` |
| Field ionization (ADK, BSI, Tong-Lin), envelope-averaged ionization, charge-state ladder | `references/ionization.md` |
| All diagnostic blocks (Scalar/Fields/Probe/ParticleBinning/Track/Screen/etc.) | `references/diagnostics.md` |
| Opening simulations with happi, post-processing, unit conversion, plotting, VTK export | `references/happi-and-postprocessing.md` |
| **Excimer-specific recipes (KrF/ArF worked namelists), UV-regime gotchas, long-pulse strategies** | `references/excimer-and-LTP-recipes.md` |
| Running Smilei, MPI/OpenMP, GPU, machine files, debugging non-starting simulations | `references/running-and-troubleshoot.md` |

For excimer-laser plasma questions, **always read `excimer-and-LTP-recipes.md`**, then load the topic-specific reference (collisions / ionization / laser) as needed. The excimer file is the orientation; the others are the API.

## Workflow patterns

### Pattern A — New namelist from scratch
1. Identify the operating mode (the table above).
2. Read `basics-and-units.md` and confirm `reference_angular_frequency_SI` with the user if collisions/ionization will be used.
3. Read mode-specific references (e.g., `laser-configuration.md` + `species-and-profiles.md` + `excimer-and-LTP-recipes.md` for excimer).
4. Read `diagnostics.md` to choose minimum diagnostic set.
5. Draft the namelist using the minimal template at the bottom of this file as scaffolding.
6. Run `smilei_test` mentally: are all referenced species defined? Is CFL satisfied? Patches power-of-2? Reference frequency set?
7. Present to user with explicit notes on what was chosen and why.

### Pattern B — Modify existing namelist
1. Ask for or read the existing namelist in full before editing — Smilei namelists have implicit cross-references (species names, ionization targets, boundary types).
2. Identify which blocks need changes.
3. Make minimum changes consistent with the user's request.
4. Re-check hard rules 1–10 against the modified namelist.

### Pattern C — Debug a non-running or wrong-result simulation
1. Read `references/running-and-troubleshoot.md` for the symptom table.
2. Check the hard rules in order. ~80% of bugs are rule 1, 2, 3, or 9.
3. If physics seems wrong rather than crashing: check operating mode mismatch first (are collisions disabled when the regime demands them?).

### Pattern D — Post-process Smilei output
1. Read `references/happi-and-postprocessing.md`.
2. Confirm which diagnostic blocks were in the namelist (the user may not know).
3. Generate happi code with explicit unit conversion if SI output is wanted.

### Pattern E — Parameter scan
1. Identify the swept variable.
2. Build a wrapper Python script that generates namelists in a directory tree.
3. Use `preprocess()` callbacks inside the namelist for parameters derived from the sweep variable (e.g., timestep derived from λ).

## Excimer-laser plasmas quick reference

This section is intentionally embedded in `SKILL.md` so it loads on every excimer-related query, not deferred to a reference.

### Critical density at common excimer wavelengths

| Wavelength | Laser | `n_c` [cm⁻³] | `reference_angular_frequency_SI` |
|---|---|---|---|
| 157 nm | F₂ | 4.45 × 10²² | `2*math.pi*3e8/157e-9` ≈ 1.20e16 |
| 193 nm | ArF | 2.94 × 10²² | `2*math.pi*3e8/193e-9` ≈ 9.76e15 |
| 248 nm | KrF | 1.78 × 10²² | `2*math.pi*3e8/248e-9` ≈ 7.59e15 |
| 308 nm | XeCl | 1.15 × 10²² | `2*math.pi*3e8/308e-9` ≈ 6.11e15 |
| 351 nm | XeF | 8.92 × 10²¹ | `2*math.pi*3e8/351e-9` ≈ 5.37e15 |
| 1064 nm | Nd:YAG (compare) | 9.69 × 10²⁰ | `2*math.pi*3e8/1064e-9` ≈ 1.77e15 |

### Regime decision (excimer defaults)

```
If wavelength < 400 nm:
    enable Collisions for all charged species pairs
    enable CollisionalIonization where neutrals are present
    If laser intensity > 1e13 W/cm²: also enable Ionization (tunnel model)
    If pulse duration > 1 ps: prefer envelope model, or simulate a representative window
    cell_length such that dx < min(skin_depth, λ_Debye) / 10
```

### Common excimer-namelist mistakes

**BAD** — forgetting reference frequency (collisions will compute with `omega_r = 1`, giving nonsense rates):
```python
Main(
    geometry = "2Dcartesian",
    # ... no reference_angular_frequency_SI ...
)
Collisions(species1 = ["electron"], species2 = ["ion"], ...)
```

**GOOD**:
```python
import math
lambda_KrF = 248e-9
Main(
    geometry = "2Dcartesian",
    reference_angular_frequency_SI = 2*math.pi*3e8/lambda_KrF,
    # ...
)
Collisions(species1 = ["electron"], species2 = ["ion"], ...)
```

**BAD** — picking `cell_length` from a 1 μm habit at 248 nm:
```python
# At 800 nm this resolves the laser; at 248 nm it's 3× under-resolved
cell_length = [0.125, 0.125]  # in units of c/omega_r — fine for 800 nm
```

**GOOD** — recompute against the actual wavelength:
```python
# omega_r is tied to λ_KrF; cell ≤ λ/16 still resolves laser AND skin depth at n_e ~ n_c
cell_length = [2*math.pi/16, 2*math.pi/16]  # 16 cells per laser wavelength
```

## Minimum complete template (KrF on overdense Ar plasma)

This is a working starting point for the excimer regime. Adapt rather than retype.

```python
import math

# ---- physical constants and chosen scales ----
lambda_SI = 248e-9                      # KrF
omega_SI  = 2*math.pi*3e8/lambda_SI
n_c_SI    = 1.78e22                     # cm^-3, critical density at 248 nm

# ---- Main block ----
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [2*math.pi/16, 2*math.pi/16],   # 16 cells / lambda
    grid_length        = [60.,  40.],                    # in units of c/omega_r
    number_of_patches  = [16, 8],
    simulation_time    = 200.,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_SI,           # rule 1
    print_every = 200,
)

LoadBalancing(every = 50, cell_load = 1.0, frequency_load = 0.1)

# ---- laser ----
LaserGaussian2D(
    box_side  = "xmin",
    a0        = 0.01,                   # I ~ 4.6e15 W/cm² at 248 nm
    omega     = 1.,                     # in units of omega_r
    focus     = [30., 20.],
    waist     = 5.,
    time_envelope = ttrapezoidal(start=0., plateau=80., slope1=10., slope2=10.),
)

# ---- species ----
Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-4]*3,             # ~50 eV; in units of m_e c^2
    particles_per_cell = 16,
    mass = 1.0, charge = -1.0,
    number_density = trapezoidal(0.5, xvacuum=20., xplateau=20., xslope1=2., xslope2=2.),
    boundary_conditions = [["remove"]]*2,
)

Species(
    name = "ion",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-5]*3,
    particles_per_cell = 16,
    mass = 39.948*1836,                 # Ar
    charge = 1.0,
    number_density = trapezoidal(0.5, xvacuum=20., xplateau=20., xslope1=2., xslope2=2.),
    boundary_conditions = [["remove"]]*2,
    atomic_number = 18,
)

# ---- collisions (rule 1 satisfied above) ----
Collisions(species1 = ["electron"], species2 = ["ion"],  coulomb_log = 0.)
Collisions(species1 = ["electron"], species2 = ["electron"], coulomb_log = 0.)
Collisions(species1 = ["ion"],      species2 = ["ion"],      coulomb_log = 0.)

# ---- diagnostics ----
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

## Closing notes for the frontier model

- When the user provides intensity in W/cm², convert to `a0` using `a0 = 0.85 * sqrt(I_18 * lambda_um²)` where `I_18` is intensity in 10¹⁸ W/cm² and `lambda_um` is wavelength in μm. State the conversion explicitly.
- When in doubt about a Smilei API, prefer to verify against `references/` or the official docs over guessing. Smilei syntax has evolved across versions and confident wrong syntax is worse than asking.
- The user may have a specific Smilei version installed. Ask if you see version-sensitive features (GPU, AM cylindrical, envelope cylindrical).
- The publications page at https://smileipic.github.io/Smilei/Overview/material.html is the best source for finding precedent papers in a given regime — suggest it for users whose problem is unusual.
