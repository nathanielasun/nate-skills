# Excimer-laser plasmas and recipes

Orientation for Smilei use outside its typical UHI regime — excimer lasers (ArF, KrF, XeCl, XeF, F₂), DUV/EUV, gas discharges, pre-formed plasmas at moderate intensities. Conventional Smilei wisdom assumes 800 nm Ti:Sapphire physics; almost every default is wrong for excimer wavelengths.

## Table of contents

1. Why excimer plasmas are different
2. Regime taxonomy
3. UV-specific hard checks
4. Resolution budget calculator
5. Worked namelist A — KrF on overdense Ar (collisional)
6. Worked namelist B — ArF on pre-ionized H (short-pulse ICF window)
7. Long-pulse strategy
8. Physics modules needing care at UV
9. When Smilei is the wrong tool

---

## 1. Why excimer plasmas are different

Smilei was built for UHI femtosecond Ti:Sapphire physics. Excimer regimes are very different:

| Property | Typical UHI | Typical excimer |
|---|---|---|
| Wavelength | 800 nm – 1.06 μm | 157 – 351 nm |
| `n_c` | 10²¹ cm⁻³ | 10²² – 4×10²² cm⁻³ |
| Pulse duration | 10 – 100 fs | 1 ns – 20 ns (ps for ICF) |
| Intensity | 10¹⁸ – 10²² W/cm² | 10⁸ – 10¹³ W/cm² (10¹⁴–10¹⁶ for ICF short-pulse) |
| `a₀` | ≳ 1 (relativistic) | 10⁻⁵ – 10⁻² (sub-relativistic) |
| Dominant absorption | Vacuum heating, J×B, ponderomotive | Inverse bremsstrahlung, collisional, photoionization |
| Collisions | Often negligible | Dominant |
| Field ionization | Tunnel (γ_K ≪ 1) | Multiphoton/tunnel boundary (γ_K ~ 1) |
| Plasma state | Fully ionized fast | Partially ionized, evolving |

**Default UHI namelists copied to excimer regimes produce nonsense** — under-resolved skin depth, missing collisional absorption, wrong ionization channel.

## 2. Regime taxonomy

| Regime | Intensity | Smilei applicability |
|---|---|---|
| I — Industrial / medical (low) | 10⁸–10¹⁰ W/cm² | Usually wrong tool. Use hydro codes. |
| II — Mid-intensity laboratory | 10¹¹–10¹³ W/cm² | Works, requires Collisions + Ionization, expensive |
| III — ICF / DUV high-energy | 10¹⁴–10¹⁶ W/cm² | Reasonable for short-pulse or windowed |
| IV — Short-pulse excimer | sub-ps, can exceed 10¹⁷ W/cm² | UHI-like, use standard templates with correct UV `n_c` |

**Excimer regime II is the default Smilei use case for this skill** — moderate-intensity collision-dominated plasmas.

## 3. UV-specific hard checks

- [ ] `reference_angular_frequency_SI = 2π c / λ_excimer` (rule 1)
- [ ] `cell_length` resolves `λ_laser/16` AND `c/ω_pe` at peak density AND `2× λ_D` minimum
- [ ] `Collisions` block for all charged species pairs
- [ ] `CollisionalIonization` if neutrals present
- [ ] `Ionization` (field) if `I > 10¹³ W/cm²` (Keldysh threshold at UV)
- [ ] `particles_per_cell ≥ 32` for collisional runs
- [ ] `LoadBalancing` enabled
- [ ] Boundary conditions match physical setup (gas-jet → reflective ions; solid → absorbing)
- [ ] `simulation_time` realistic given pulse and `n_c`

## 4. Resolution budget calculator

```python
import math

lambda_SI    = 248e-9
n_e_peak_SI  = 5e21              # cm⁻³
T_e_eV       = 50.0

c, e_, m_e, eps0 = 2.99792458e8, 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r     = eps0*m_e*omega_r**2/e_**2

# Three constraints (normalized)
n_norm = (n_e_peak_SI*1e6)/N_r
T_norm = T_e_eV/511e3

dx_laser = 2*math.pi/16
dx_skin  = 0.3 * 1/math.sqrt(n_norm)
dx_debye = 2.0 * math.sqrt(T_norm/n_norm)

dx = min(dx_laser, dx_skin, dx_debye)
print(f"laser={dx_laser:.4f}  skin={dx_skin:.4f}  debye={dx_debye:.4f}")
print(f"use dx = {dx:.4f} c/ω_r = {dx*L_r*1e9:.2f} nm")
```

For typical excimer cases (10²¹–10²² cm⁻³, 30–100 eV), **Debye dominates** — `dx_debye` is the binding constraint.

## 5. Worked namelist A — KrF on overdense Ar (collisional)

The canonical "moderate intensity, dense plasma, collision-dominated" case.

```python
"""
KrF 248 nm on overdense argon.
I ≈ 10¹² W/cm², n_e ≈ 5×n_c, partially ionized.
"""
import math

# === Physical inputs (SI) ===
lambda_SI      = 248e-9
intensity_SI   = 1.0e12
pulse_FWHM_SI  = 100e-15
waist_SI       = 8e-6
n_e_SI         = 2.0e22                       # cm⁻³ ≈ 5 × n_c(248 nm)
T_e_eV         = 30.0
T_i_eV         = 5.0
box_x_SI       = 25e-6
box_y_SI       = 16e-6
sim_time_SI    = 400e-15

# === Derived ===
c, e_, m_e, eps0 = 2.99792458e8, 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r     = eps0*m_e*omega_r**2/e_**2

a0      = 0.85 * math.sqrt(intensity_SI*1e-18) * (lambda_SI*1e6)
n_e     = (n_e_SI*1e6)/N_r
T_e     = T_e_eV/511e3
T_i     = T_i_eV/511e3
pulse   = pulse_FWHM_SI*omega_r
waist   = waist_SI/L_r
Lx, Ly  = box_x_SI/L_r, box_y_SI/L_r
T_tot   = sim_time_SI*omega_r

lambda_D = math.sqrt(T_e/n_e)
dx = min(2*math.pi/20, 2.0*lambda_D)

# === Namelist ===
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL = 0.95,
    cell_length = [dx, dx],
    grid_length = [Lx, Ly],
    number_of_patches = [32, 16],
    simulation_time = T_tot,
    EM_boundary_conditions = [["silver-muller"]]*2,
    reference_angular_frequency_SI = omega_r,
    print_every = 200,
)

LoadBalancing(every = 50, cell_load = 1.0, frequency_load = 0.1)

LaserGaussian2D(
    box_side = "xmin", a0 = a0, omega = 1.0,
    focus = [Lx*0.5, Ly*0.5], waist = waist,
    time_envelope = tgaussian(fwhm = pulse, center = 3*pulse),
)

ramp_start = 0.4*Lx
ramp_width = 1e-6/L_r

def n_profile(x, y):
    if x < ramp_start:               return 0.0
    elif x < ramp_start + ramp_width: return n_e*(x-ramp_start)/ramp_width
    else:                             return n_e

Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [T_e]*3,
    particles_per_cell = 64,
    mass = 1.0, charge = -1.0,
    number_density = n_profile,
    boundary_conditions = [["remove"], ["periodic"]],
)

Species(
    name = "Ar1",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    temperature = [T_i]*3,
    particles_per_cell = 64,
    mass = 39.948*1836, charge = 1.0, atomic_number = 18,
    number_density = n_profile,
    boundary_conditions = [["remove"], ["periodic"]],
)

# All charged pairs
Collisions(species1=["electron"], species2=["electron"], coulomb_log=0.)
Collisions(species1=["electron"], species2=["Ar1"],      coulomb_log=0.)
Collisions(species1=["Ar1"],      species2=["Ar1"],      coulomb_log=0.)

# Diagnostics
DiagScalar(every = 50)
DiagFields(every = 200, fields = ["Ex","Ey","Bz","Rho_electron","Rho_Ar1"])
DiagProbe(every = 100, origin=[0., Ly*0.5], corners=[[Lx, Ly*0.5]], number=[400],
          fields=["Ex","Ey","Bz","Rho"])
DiagParticleBinning(
    deposited_quantity = "weight_ekin", every = 200, species = ["electron"],
    axes = [["ekin", 1e-4, 1e-1, 100, "logscale"]],
)
```

**Key choices**: `PPC = 64` (collisional needs high PPC). `dx` from Debye, not laser. Pulse modeled as 100 fs window — full ns pulse is impractical.

## 6. Worked namelist B — ArF on pre-ionized H (short-pulse ICF window)

Regime III/IV: representative window from a Nike-class ArF experiment.

```python
"""
ArF 193 nm short pulse on pre-ionized H plasma.
I ≈ 5×10¹⁴ W/cm², a₀ ≈ 0.012, ICF-relevant.
"""
import math

lambda_SI    = 193e-9
intensity_SI = 5e14
pulse_FWHM_SI = 500e-15
waist_SI     = 30e-6
n_e_SI       = 1e21                          # cm⁻³ — pre-ablated corona
T_e_eV       = 800.0
box_x_SI     = 80e-6
box_y_SI     = 60e-6
sim_time_SI  = 2e-12

c, e_, m_e, eps0 = 2.99792458e8, 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12
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

# High T_e → Debye is large → laser-resolution dominates
dx = 2*math.pi/20

Main(
    geometry = "2Dcartesian", interpolation_order = 2, timestep_over_CFL = 0.95,
    cell_length = [dx, dx], grid_length = [Lx, Ly],
    number_of_patches = [32, 16], simulation_time = T_tot,
    EM_boundary_conditions = [["silver-muller"]]*2,
    reference_angular_frequency_SI = omega_r,
    maxwell_solver = "Yee", print_every = 500,
)

LoadBalancing(every = 100)

LaserGaussian2D(
    box_side = "xmin", a0 = a0, omega = 1.0,
    focus = [Lx*0.5, Ly*0.5], waist = waist,
    time_envelope = tgaussian(fwhm = pulse, center = 3*pulse),
)

def n_profile(x, y):
    scale = 5e-6/L_r
    if x < 0.3*Lx: return 0.0
    return n_e * math.exp(-(x - 0.3*Lx)/scale)

Species(name="electron", position_initialization="regular",
        momentum_initialization="maxwell-juettner", temperature=[T_e]*3,
        particles_per_cell=32, mass=1.0, charge=-1.0,
        number_density=n_profile,
        boundary_conditions=[["remove"], ["periodic"]])

Species(name="ion", position_initialization="regular",
        momentum_initialization="maxwell-juettner", temperature=[T_e*0.1]*3,
        particles_per_cell=32, mass=1836.0, charge=1.0,
        number_density=n_profile,
        boundary_conditions=[["remove"], ["periodic"]])

Collisions(species1=["electron"], species2=["ion"], coulomb_log=0.)
Collisions(species1=["electron"], species2=["electron"], coulomb_log=0.)

DiagScalar(every=100)
DiagFields(every=500, fields=["Ex","Ey","Bz","Rho_electron"])
DiagProbe(every=200, origin=[0., Ly*0.5], corners=[[Lx, Ly*0.5]],
          number=[500], fields=["Ex","Ey","Bz"])
```

Note: 500 fs window. A full Nike pulse is 4 ns — 8000× longer. Window-based simulation is essential.

## 7. Long-pulse strategy

A 20 ns excimer pulse is impractical to simulate end-to-end. Two options:

**Quasi-steady-state representative window:** Use hydrodynamic code to evolve the plasma to steady-state ablation profile → import as initial condition → simulate ~10 ps with PIC for kinetic physics.

**Chained segments with checkpoint/restart:**

```python
Checkpoints(
    restart_dir = None,                    # set to previous segment's dir on restart
    dump_step = 5000,
    exit_after_dump = True,                # exit after writing — controller resubmits
    keep_n_dumps = 2,
)
```

A wrapper bash/Python script handles the chain. Each segment adjusts the laser intensity to match the pulse envelope at that time.

## 8. Physics modules needing care at UV

**Field ionization**: Keldysh `γ_K` is borderline at UV. For KrF at `I = 10¹³ W/cm²`, `γ_K ≈ 1.2`. ADK strictly valid for `γ_K ≪ 1`. **Use `tunnel_full_PPT` for excimer regimes** with `γ_K` between 0.5 and 5.

**Photoionization**: Smilei has no native photoionization. ArF at 6.4 eV/photon can directly ionize low-IP atoms. Workaround: `from_rate` model with tabulated cross-section.

**Collisional ionization**: Lotz cross-sections, validated for moderate-Z. Lower accuracy for high-Z (W, Au, lanthanides).

**Radiation reaction**: Sub-relativistic excimer regimes have negligible radiation reaction. Do not enable; adds overhead.

**Maxwell solver**: `"Yee"` is fine at UV. Avoid `"Lehe"` (optimized for relativistic LWFA, may distort short-wavelength response).

**Particle pusher**: `"borisnr"` (non-relativistic Boris) is fastest and slightly more accurate for `a₀ ≪ 1`. Consider it for pure Regime I/II excimer runs.

## 9. When Smilei is the wrong tool

| Use case | Better choice |
|---|---|
| Industrial UV ablation (LASIK, polymer cutting) | Hydro codes (FLASH, HELIOS) with EOS + radiation transport |
| Photolithography plasma effects | Specialized lithography models |
| Long-pulse ICF ablation (full ns history) | Radiation-hydrodynamics (FLASH, MULTI, DRACO, FCI2) |
| Plasma discharge / glow discharge | Low-temperature plasma codes (HPEM, LXcat fluid models) |
| EUV source plasma (Sn, Xe) | Hybrid codes; PIC sub-modules |
| Pre-formed plasma + short-pulse excimer | Smilei is appropriate |
| LPI in DUV ICF (TPD, SBS, SRS) | Smilei is appropriate |
| Cross-beam energy transfer with smoothed beams | Smilei (Oudin et al. 2022) |
| Photoionization-dominated plasma formation | Smilei with `from_rate` workaround, or codes with native photoionization (PIConGPU) |

If asked "how do I model a 20 ns KrF pulse ablating Si," the most useful answer is **"Smilei is not the right primary tool — use a radiation-hydro code, and use Smilei only for the kinetic sub-process of interest within a representative window."** Be willing to say this.

If the user does want PIC physics in an excimer regime (laser-plasma instabilities, hot-electron transport, kinetic effects in pre-formed plasma), Smilei is well-suited — configure carefully along the lines of §5 and §6.
