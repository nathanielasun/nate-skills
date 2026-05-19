---
name: smilei-namelist-small
description: Use for writing, editing, debugging, or running a Smilei Python namelist, post-processing with happi, or planning a Smilei PIC simulation. Trigger on "Smilei", `smilei_test`, `happi.Open`, `Main(geometry=...)`, `Species`, `Laser`, `LaserGaussian2D`, `LaserEnvelope`, `Collisions`, `Ionization`, `DiagScalar`, `DiagFields`, `DiagProbe`, `DiagParticleBinning`, or any PIC namelist file. Includes excimer-laser plasma content (ArF 193, KrF 248, XeCl 308, XeF 351, F2 157 nm) where collisional and ionization physics dominate. Use this skill instead of guessing — Smilei's API has version-sensitive details and silent physics errors (missing reference frequency, non-power-of-2 patches, wrong ionization model for the Keldysh regime) only surface deep into runs. For Smilei C++ source work, use smilei-cpp-small instead.
---

# Smilei Namelist Skill (small tier)

This skill covers Smilei Python namelist authoring, debugging, running, and post-processing with happi. It does not cover Smilei's C++ source (use `smilei-cpp-small`). Tested against Smilei v5.1.

Smilei was designed for ultra-high-intensity femtosecond laser-plasma interaction. Excimer-laser and low-temperature plasma users sit in regimes where Smilei's collisionless / ponderomotive defaults are wrong. This skill carries excimer/LTP content alongside standard UHI material.

## Hard rules (do not violate)

1. **`reference_angular_frequency_SI` is mandatory when `Collisions`, `CollisionalIonization`, `Ionization`, or `RadiationReaction` appears.** Compute from laser wavelength: `2*math.pi*c/lambda_SI`. Missing this produces silent order-of-magnitude errors in collision rates.

2. **`number_of_patches` must be a power of 2 in each direction.** `[16, 8]` valid; `[12, 8]` not. Pick patch counts so each MPI rank owns 4–16 patches.

3. **CFL ceiling depends on the Maxwell solver.** Default Yee: `timestep_over_CFL ≤ 1`. Set 0.95 by default for margin.

4. **For UV wavelengths (< 400 nm), resolve at the actual `n_c`, not at 1 μm.** `n_c ∝ 1/λ²` — at 193 nm `n_c ≈ 2.9 × 10²² cm⁻³`, ~30× higher than at 1 μm. Cell size must resolve `c/ω_pe` at the local electron density.

5. **Density profiles are in normalized units.** `number_density = 0.1` means `0.1 × n_c`, not `0.1 cm⁻³`. Convert SI inputs explicitly using `reference_angular_frequency_SI`.

6. **`Species` ordering matters for `Ionization`.** Define the electron species before the neutral that references it via `ionization_electrons`.

7. **Enable `LoadBalancing` for any laser-driven or collisional run.** Without it, dense regions stall MPI ranks. `LoadBalancing(every=50)` is a fine default.

8. **`DiagFields` with `time_average > 1` smears physics.** Use only for envelope-like outputs; never for laser snapshots or polarization.

9. **`number_density` callables must return non-negative floats for all coordinates.** Test profiles standalone (`matplotlib`) before plugging in.

10. **Always run `smilei_test <namelist>.py` before launching.** It catches ~90% of namelist errors in seconds.

## Four operating modes — pick one explicitly

| Mode | Laser | Plasma | Required physics |
|---|---|---|---|
| **UHI / relativistic** | `a₀ ≳ 1`, fs, NIR ~800 nm | Underdense → overdense | None mandatory |
| **Collisional moderate-intensity** | `a₀ ≪ 1`, ns–ps, UV–visible | Often overdense, partially ionized | Collisions + CollisionalIonization + Field Ionization mandatory |
| **Laser-envelope** | Long pulse, slowly varying | Underdense LWFA or ns-class long pulse | Envelope solver; `tunnel_envelope_averaged` |
| **Wakefield / LWFA** | `a₀ ≳ 1`, fs, NIR | Underdense, pre-ionized | Moving window; envelope optional |

**Excimer-laser plasmas default to row 2.** UV wavelengths, dense `n_c`, ns pulses, collision- and ionization-dominated. Smilei's collisionless defaults are wrong. Exception: sub-ps KrF amplifier output for ICF can approach row 1 physics but with row-2 `n_c`.

## Routing — read references before answering

| If the question is about… | Read |
|---|---|
| Units, `reference_angular_frequency_SI`, density/temperature conversion, laser config, species definition, profiles, pushers, boundary conditions | `references/namelist-fundamentals.md` |
| Binary collisions, collisional ionization, field ionization (ADK/PPT/BSI/envelope-averaged), charge-state ladder | `references/collisions-and-ionization.md` |
| Diagnostic blocks (Scalar/Fields/Probe/ParticleBinning/Track/Screen) and post-processing with happi | `references/diagnostics-and-happi.md` |
| **Excimer-specific recipes (worked KrF/ArF namelists), UV gotchas, long-pulse strategies** | `references/excimer-recipes.md` |
| Running Smilei (MPI/OpenMP/GPU), checkpoints, troubleshooting symptom tables | `references/running-and-troubleshoot.md` |

For excimer-laser questions, **always read `excimer-recipes.md`**, then the topic-specific reference.

## Workflow patterns

**New namelist from scratch:** Identify operating mode → read `namelist-fundamentals.md` → confirm `reference_angular_frequency_SI` if collisions/ionization will run → read mode-specific references → draft using the template below → verify hard rules 1–10 → run `smilei_test`.

**Debug a non-running namelist:** Check hard rules in order. ~80% of bugs are rule 1, 2, 3, or 9. Then consult `running-and-troubleshoot.md` symptom table.

**Wrong-physics results:** Check operating mode mismatch first (collisions disabled when regime demands them, ADK at high Keldysh γ). Then re-verify rule 1 and rule 5.

**Post-process output:** Read `diagnostics-and-happi.md`. Confirm which diagnostic blocks the namelist had. Generate happi code with explicit unit conversion if SI output is wanted.

## Excimer-laser plasmas quick reference

| λ | Laser | `n_c` [cm⁻³] | `reference_angular_frequency_SI` |
|---|---|---|---|
| 157 nm | F₂ | 4.45 × 10²² | `2*pi*3e8/157e-9` ≈ 1.20e16 |
| 193 nm | ArF | 2.94 × 10²² | `2*pi*3e8/193e-9` ≈ 9.76e15 |
| 248 nm | KrF | 1.78 × 10²² | `2*pi*3e8/248e-9` ≈ 7.59e15 |
| 308 nm | XeCl | 1.15 × 10²² | `2*pi*3e8/308e-9` ≈ 6.11e15 |
| 351 nm | XeF | 8.92 × 10²¹ | `2*pi*3e8/351e-9` ≈ 5.37e15 |
| 1064 nm | Nd:YAG (compare) | 9.69 × 10²⁰ | `2*pi*3e8/1064e-9` ≈ 1.77e15 |

**Excimer regime defaults (λ < 400 nm):**
- Enable `Collisions` for all charged species pairs
- Enable `CollisionalIonization` where neutrals are present
- Enable `Ionization` (tunnel or PPT) if `I > 10¹³ W/cm²`
- For pulses > 1 ps, prefer envelope model or simulate a representative window
- `cell_length` such that `dx < min(skin_depth, λ_Debye) / 10`

`a₀` from intensity (linear polarization): `a₀ = 0.85 × sqrt(I_18 × λ_μm²)` with `I_18 = I/10¹⁸ W/cm²`.

## Minimum complete template (KrF on overdense Ar — excimer default)

```python
import math

# === physical scales ===
lambda_SI = 248e-9                                  # KrF
omega_SI  = 2*math.pi*3e8/lambda_SI

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [2*math.pi/16, 2*math.pi/16],   # 16 cells / λ
    grid_length        = [60.,  40.],                    # in c/ω_r
    number_of_patches  = [16, 8],
    simulation_time    = 200.,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_SI,           # rule 1
    print_every = 200,
)

LoadBalancing(every = 50)

LaserGaussian2D(
    box_side  = "xmin",
    a0        = 0.01,                                    # I ~ 4.6e15 W/cm² at 248 nm
    omega     = 1.,
    focus     = [30., 20.],
    waist     = 5.,
    time_envelope = ttrapezoidal(start=0., plateau=80., slope1=10., slope2=10.),
)

Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [1e-4]*3,                              # ~50 eV
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
    mass = 39.948*1836,                                  # Ar
    charge = 1.0,
    number_density = trapezoidal(0.5, xvacuum=20., xplateau=20., xslope1=2., xslope2=2.),
    boundary_conditions = [["remove"]]*2,
    atomic_number = 18,
)

Collisions(species1 = ["electron"], species2 = ["ion"],     coulomb_log = 0.)
Collisions(species1 = ["electron"], species2 = ["electron"], coulomb_log = 0.)
Collisions(species1 = ["ion"],      species2 = ["ion"],      coulomb_log = 0.)

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

## Closing reminders

- Intensity to `a₀`: `a₀ = 0.85 × sqrt(I_18 × λ_μm²)`. State the conversion explicitly when given W/cm².
- When uncertain about a Smilei API, prefer verifying against `references/` over guessing.
- Ask about Smilei version for version-sensitive features (GPU, AM cylindrical, envelope cylindrical).
- The Smilei publications page (smileipic.github.io/Smilei/Overview/material.html) is the best precedent source for unusual regimes.
