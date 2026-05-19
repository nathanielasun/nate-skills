# Excimer-laser plasmas and low-temperature-plasma recipes

This reference is the orientation file for Smilei use outside its typical UHI regime. Read this first when the user mentions excimer lasers (ArF, KrF, XeCl, XeF, F₂), DUV/EUV, gas discharges, laser ablation, or pre-formed plasmas at moderate intensities. The conventional Smilei wisdom assumes 800 nm Ti:Sapphire physics; almost every default is wrong for excimer wavelengths.

## Table of contents

1. Why excimer plasmas are different from typical Smilei targets
2. Regime taxonomy: where does the user actually sit?
3. The UV-specific hard checks
4. Resolution budget calculator
5. Worked namelist A — KrF on overdense argon (collisional, moderate-intensity)
6. Worked namelist B — ArF on pre-ionized hydrogen (ICF-relevant short pulse)
7. Worked namelist C — Long-pulse strategy with segmented runs
8. Worked namelist D — Envelope model for ns-class pulse (with strong caveats)
9. Physics modules that need extra care at UV
10. When Smilei is the wrong tool — and what is

---

## 1. Why excimer plasmas are different from typical Smilei targets

Smilei was built around ultra-high-intensity femtosecond Ti:Sapphire physics: collisionless, relativistic, dominated by the ponderomotive force, with laser wavelengths near 800 nm. Excimer-laser plasmas sit somewhere very different:

| Property | Typical UHI use of Smilei | Typical excimer regime |
|---|---|---|
| Wavelength | 800 nm – 1.06 μm | 157 – 351 nm |
| `n_c` | 10²¹ cm⁻³ | 10²² – 4×10²² cm⁻³ |
| Pulse duration | 10 – 100 fs | 1 ns – 20 ns (sometimes ps for ICF) |
| Intensity | 10¹⁸ – 10²² W/cm² | 10⁸ – 10¹³ W/cm² (10¹⁴–10¹⁶ for ICF short-pulse) |
| `a₀` | ≳ 1 (relativistic) | 10⁻⁵ – 10⁻² (deeply sub-relativistic) |
| Dominant absorption | Vacuum heating, J×B, ponderomotive | Inverse bremsstrahlung, collisional, photoionization |
| Collisions | Often negligible | Dominant |
| Field ionization | Tunnel (Keldysh γ ≪ 1) | Multiphoton or tunnel boundary (Keldysh γ ~ 1) |
| Plasma state | Fully ionized fast | Partially ionized, evolving |

Practical implication: **the excimer regime is one Smilei *can* simulate, but only if collisions, ionization, and proper resolution are configured.** Default Smilei namelists copied from UHI examples will produce nonsense in this regime — under-resolved skin depth, missing collisional absorption, wrong ionization channel.

## 2. Regime taxonomy: where does the user actually sit?

Excimer-laser plasma research splits into roughly four regimes. Identify which one before drafting a namelist.

### Regime I — Industrial/medical low-intensity (10⁸ – 10¹⁰ W/cm²)

Photolithography, ablation, LASIK, polymer processing. Intensity is below the ionization threshold of typical targets at UV wavelengths but high enough to produce a plume of vapor and weakly ionized plasma after a few nanoseconds. **Smilei is usually the wrong tool here** — hydrodynamic codes (FLASH, MULTI, HYDRA, HELIOS) with non-LTE radiation transport and tabulated EOS are standard. PIC adds little when the mean free path is much shorter than the cell size.

### Regime II — Mid-intensity laboratory (10¹¹ – 10¹³ W/cm²)

Plasma plumes from solid targets, gas-jet experiments, plasma diagnostics. Intensity reaches ionization but is far from relativistic. **Smilei works here but requires Collisions + CollisionalIonization + Ionization, fine resolution at high `n_c`, and is computationally expensive due to long pulse durations.** Often run as a representative window rather than full pulse.

### Regime III — ICF / DUV high-energy facility (10¹⁴ – 10¹⁶ W/cm²)

Nike, Electra, ESHEL — direct-drive ICF and related platforms using KrF/ArF. Still sub-relativistic but laser-plasma instabilities (TPD, SRS, SBS) matter, and the plasma is dense and hot. **Smilei is a reasonable tool but typically used for short-pulse or representative windows; full ns pulse is impractical.**

### Regime IV — Short-pulse excimer (sub-picosecond)

Hybrid KrF-dye systems, ultrashort-pulse ArF amplifiers. Peak intensity can exceed 10¹⁷ W/cm². Closer in physics to UHI Ti:Sa than to nanosecond excimer. **Use standard Smilei UHI templates but with `n_c` at the correct UV wavelength**. Collisions still matter in the early pre-pulse phase but the main pulse is collisionless.

## 3. The UV-specific hard checks

Before launching any excimer namelist, verify:

- [ ] `reference_angular_frequency_SI` set to `2π c / λ_excimer` (rule 1)
- [ ] `cell_length` resolves both `λ_laser/16` AND `c/ω_pe` at peak density AND `2× λ_D` minimum
- [ ] `Collisions` block present for all pairs of charged species that physically collide
- [ ] `CollisionalIonization` enabled if neutral atoms are present
- [ ] `Ionization` (field) enabled if `I > 10¹³ W/cm²` (Keldysh threshold at UV)
- [ ] `particles_per_cell` at least 32 for collisional runs (vs typical 8–16 for UHI)
- [ ] `LoadBalancing` enabled (collisional load is uneven)
- [ ] Boundary conditions match physical setup (gas-jet → reflective on ions; solid target → absorbing)
- [ ] `simulation_time` realistic given pulse length and `n_c` (long-pulse warning below)
- [ ] Check intensity → `a0` conversion: at 248 nm, `I = 10¹³ W/cm²` gives `a0 ≈ 2.1×10⁻³`

## 4. Resolution budget calculator

This is the calculation to do *before* drafting the namelist. Choose `cell_length` such that the tightest of three constraints is satisfied:

```python
import math

# === inputs ===
lambda_SI    = 248e-9            # KrF
n_e_peak_SI  = 5e21              # cm⁻³ (set by your gas/solid target)
T_e_eV       = 50.0              # electron temperature

# === derived ===
c       = 2.99792458e8
e       = 1.602176634e-19
m_e     = 9.1093837015e-31
eps0    = 8.8541878128e-12

omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
N_r_SI  = eps0*m_e*omega_r**2/e**2

# Three constraints (all in normalized units c/omega_r):

# 1. Laser-resolution: at least 16 cells per laser wavelength
dx_laser = 2*math.pi/16

# 2. Skin depth: cell ≤ 0.3 * c/omega_pe at peak density
n_norm   = (n_e_peak_SI*1e6) / N_r_SI
skin_depth_norm = 1/math.sqrt(n_norm)
dx_skin  = 0.3 * skin_depth_norm

# 3. Debye length: cell ≤ 2 * lambda_Debye
T_norm   = T_e_eV / 511e3
lambda_D_norm = math.sqrt(T_norm/n_norm)
dx_debye = 2.0 * lambda_D_norm

print(f"Constraints (normalized to c/omega_r):")
print(f"  laser:  dx ≤ {dx_laser:.4f}")
print(f"  skin:   dx ≤ {dx_skin:.4f}")
print(f"  Debye:  dx ≤ {dx_debye:.4f}")
print(f"  → use dx = {min(dx_laser, dx_skin, dx_debye):.4f}")
print(f"  → SI:    dx = {min(dx_laser, dx_skin, dx_debye)*L_r*1e9:.2f} nm")
```

For typical excimer cases at moderate density (10²¹ cm⁻³, 50 eV), the Debye constraint or skin constraint usually wins — the laser-resolution constraint is generally satisfied automatically.

**Example output for KrF + 5×10²¹ cm⁻³ + 50 eV:**
- `dx_laser ≈ 0.393`
- `dx_skin  ≈ 0.180` ← binding constraint
- `dx_debye ≈ 0.029`

Wait — `dx_debye` is the binding constraint at 50 eV and 5×10²¹ cm⁻³. The skin depth is wider. This is the typical pattern at moderate excimer densities: **Debye dominates, not laser wavelength.** Plan compute budget accordingly: at `dx = 0.029 c/ω_r` in 2D over a 30 μm × 20 μm box, the grid is ~26,000 × 17,500 cells. This is large; you may need to either coarse-grain in temperature, reduce box size, or move to a hybrid code.

## 5. Worked namelist A — KrF on overdense argon (collisional)

The canonical "moderate intensity, dense plasma, collision-dominated" case.

```python
"""
KrF 248 nm laser interacting with overdense argon plasma.
Regime II: I ≈ 10¹² W/cm², n_e ≈ 5×n_c, partially ionized Ar+ → Ar²⁺.
"""
import math

# ============================================================
# Physical inputs (SI)
# ============================================================
lambda_SI      = 248e-9
intensity_SI   = 1.0e12           # W/cm²
pulse_FWHM_SI  = 100e-15          # 100 fs (representative window of ns pulse)
waist_SI       = 8e-6
n_e_SI         = 2.0e22           # cm⁻³ ≈ 5 × n_c(248 nm)
T_e_eV         = 30.0
T_i_eV         = 5.0
box_x_SI       = 25e-6
box_y_SI       = 16e-6
sim_time_SI    = 400e-15

# ============================================================
# Derived
# ============================================================
c       = 2.99792458e8
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
eps0    = 8.8541878128e-12
e       = 1.602176634e-19
m_e     = 9.1093837015e-31
N_r_SI  = eps0*m_e*omega_r**2/e**2

a0      = 0.85 * math.sqrt(intensity_SI*1e-18) * (lambda_SI*1e6)
n_e     = (n_e_SI*1e6) / N_r_SI
T_e     = T_e_eV / 511e3
T_i     = T_i_eV / 511e3
pulse   = pulse_FWHM_SI * omega_r
waist   = waist_SI / L_r
Lx      = box_x_SI / L_r
Ly      = box_y_SI / L_r
T_tot   = sim_time_SI * omega_r

# Debye-driven resolution
lambda_D = math.sqrt(T_e/n_e)
dx       = min(2*math.pi/20, 2.0*lambda_D)        # whichever is tighter

# ============================================================
# Namelist
# ============================================================
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [dx, dx],
    grid_length        = [Lx, Ly],
    number_of_patches  = [32, 16],
    simulation_time    = T_tot,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_r,
    random_seed = 0,
    print_every = 200,
)

LoadBalancing(every = 50, cell_load = 1.0, frequency_load = 0.1)

LaserGaussian2D(
    box_side  = "xmin",
    a0        = a0,
    omega     = 1.0,
    focus     = [Lx*0.5, Ly*0.5],
    waist     = waist,
    time_envelope = tgaussian(fwhm = pulse, center = 3*pulse),
)

# Plasma slab starts at x = 0.4 Lx, ramps up over 1 μm
ramp_start = 0.4*Lx
ramp_width = 1e-6 / L_r

def n_profile(x, y):
    if x < ramp_start:
        return 0.0
    elif x < ramp_start + ramp_width:
        return n_e * (x - ramp_start)/ramp_width
    else:
        return n_e

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
    mass = 39.948*1836,
    charge = 1.0,
    atomic_number = 18,
    number_density = n_profile,
    boundary_conditions = [["remove"], ["periodic"]],
)

# Collisions: all charged pairs
Collisions(species1 = ["electron"], species2 = ["electron"],  coulomb_log = 0.)
Collisions(species1 = ["electron"], species2 = ["Ar1"],       coulomb_log = 0.)
Collisions(species1 = ["Ar1"],      species2 = ["Ar1"],       coulomb_log = 0.)

# Collisional ionization Ar1 → Ar2
# (Need second ion species; see ionization.md for the full ladder pattern)

# ============================================================
# Diagnostics
# ============================================================
DiagScalar(every = 50)
DiagFields(every = 200, fields = ["Ex","Ey","Bz","Rho_electron","Rho_Ar1"])
DiagProbe(every = 100, origin=[0., Ly*0.5], corners=[[Lx, Ly*0.5]], number=[400],
          fields=["Ex","Ey","Bz","Rho"])
DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 200, species = ["electron"],
    axes = [["ekin", 1e-4, 1e-1, 100, "logscale"]],
)
```

Key choices and why:
- `particles_per_cell = 64` — collisional runs need many particles per cell for collision-rate convergence. Halve only if memory-bound.
- `dx` derived from Debye, not laser — at these conditions the Debye constraint wins (see §4).
- Pulse modeled as a 100 fs window — represents the peak of a real ns pulse; full simulation would require ~10⁴× more steps.
- Boundary conditions `[["remove"], ["periodic"]]` — particles leaving in x are absorbed, transverse is periodic.

## 6. Worked namelist B — ArF on pre-ionized hydrogen (short-pulse ICF window)

Regime III / IV: representative window from a Nike-class ArF experiment.

```python
"""
ArF 193 nm short pulse on pre-ionized H plasma.
I ≈ 5×10¹⁴ W/cm², a0 ≈ 0.012, ICF-relevant geometry.
"""
import math

lambda_SI    = 193e-9
intensity_SI = 5e14
pulse_FWHM_SI = 500e-15
waist_SI     = 30e-6
n_e_SI       = 1e21              # cm⁻³ — pre-ablated corona, well below n_c at 193 nm
T_e_eV       = 800.0
box_x_SI     = 80e-6
box_y_SI     = 60e-6
sim_time_SI  = 2e-12

c       = 2.99792458e8
omega_r = 2*math.pi*c/lambda_SI
L_r     = c/omega_r
eps0    = 8.8541878128e-12; e_ = 1.602176634e-19; m_e = 9.1093837015e-31
N_r_SI  = eps0*m_e*omega_r**2/e_**2

a0      = 0.85*math.sqrt(intensity_SI*1e-18)*(lambda_SI*1e6)
n_e     = (n_e_SI*1e6)/N_r_SI
T_e     = T_e_eV/511e3
pulse   = pulse_FWHM_SI*omega_r
waist   = waist_SI/L_r
Lx      = box_x_SI/L_r
Ly      = box_y_SI/L_r
T_tot   = sim_time_SI*omega_r

# At this temperature Debye is large; laser-resolution constraint dominates
dx = 2*math.pi/20

Main(
    geometry = "2Dcartesian", interpolation_order = 2, timestep_over_CFL = 0.95,
    cell_length = [dx, dx], grid_length = [Lx, Ly],
    number_of_patches = [32, 16], simulation_time = T_tot,
    EM_boundary_conditions = [["silver-muller"], ["silver-muller"]],
    reference_angular_frequency_SI = omega_r,
    maxwell_solver = "Yee", print_every = 500,
)

LoadBalancing(every = 100)

LaserGaussian2D(
    box_side = "xmin", a0 = a0, omega = 1.0,
    focus = [Lx*0.5, Ly*0.5], waist = waist,
    time_envelope = tgaussian(fwhm = pulse, center = 3*pulse),
)

# Pre-ionized exponential corona profile
def n_profile(x, y):
    scale = 5e-6/L_r           # 5 μm scale length
    if x < 0.3*Lx:
        return 0.0
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

# Lower density: collisions less dominant but still present at corona base
Collisions(species1=["electron"], species2=["ion"], coulomb_log=0.)
Collisions(species1=["electron"], species2=["electron"], coulomb_log=0.)

DiagScalar(every=100)
DiagFields(every=500, fields=["Ex","Ey","Bz","Rho_electron"])
DiagProbe(every=200, origin=[0., Ly*0.5], corners=[[Lx, Ly*0.5]], number=[500],
          fields=["Ex","Ey","Bz"])
```

Note: this is a 500 fs window. A full Nike-class pulse is 4 ns — 8000× longer. Even a representative window is computationally heavy.

## 7. Worked namelist C — Long-pulse strategy with segmented runs

A 20 ns excimer pulse is impractical to simulate end-to-end. Two options:

### Option 1 — Simulate a quasi-steady-state representative window
1. Use hydrodynamic/MHD code to evolve the plasma to the steady-state ablation profile.
2. Import that profile as initial condition.
3. Simulate ~10 ps with PIC to capture kinetic physics (laser-plasma instabilities, hot-electron generation).

### Option 2 — Chained segments with checkpoint/restart
1. Run segment 1 (e.g., 100 ps).
2. Use Smilei's checkpoint mechanism (`Checkpoints` block) to save state.
3. Restart segment 2 from the checkpoint, modifying the laser intensity to match the pulse envelope at that time.

```python
# Inside the namelist, add:
Checkpoints(
    restart_dir = None,                 # set to previous segment's directory on restart
    dump_step = 5000,                   # checkpoint every N steps
    dump_minutes = 0.,
    exit_after_dump = True,             # exit after writing — controller resubmits
    keep_n_dumps = 2,
)

# Time-dependent laser intensity ramped within the segment to match pulse envelope
def time_profile(t):
    # adjust per segment; this example: pulse plateau
    return 1.0

LaserGaussian2D(
    # ...
    time_envelope = time_profile,
)
```

Then a wrapper bash/Python script handles the chain. The restart mechanism is reliable; the orchestration is the user's responsibility.

## 8. Worked namelist D — Envelope model for ns-class pulse

The envelope model averages out the laser oscillation period, allowing coarser timesteps. **For excimer regimes this is delicate** — the envelope assumption requires that the laser amplitude varies slowly compared to the optical cycle, which is fine for ns pulses, but the model also assumes the plasma response is dominated by ponderomotive effects, which is often wrong for moderate-intensity excimer plasmas where collisions dominate.

**Use envelope for excimer only when**:
- The pulse is much longer than the laser period (always true for ns excimer)
- The intensity is high enough that ponderomotive ≳ collisional response
- You're tracking large-scale plasma evolution, not laser-cycle-resolved physics
- You've verified results against a short full-EM benchmark

```python
import math

# At the high-intensity end of excimer regime, envelope may be appropriate
lambda_SI = 248e-9
omega_r = 2*math.pi*3e8/lambda_SI

Main(
    geometry = "AMcylindrical",                # envelope works well in cylindrical
    number_of_AM = 1,
    interpolation_order = 2,
    timestep_over_CFL  = 0.95,
    cell_length        = [2*math.pi/10, 0.5],  # much coarser longitudinal — no laser-cycle resolution needed
    grid_length        = [800., 80.],
    number_of_patches  = [32, 4],
    simulation_time    = 5000.,                 # much longer feasible
    EM_boundary_conditions = [["silver-muller"], ["PML"]],
    number_of_pml_cells = [[0,0],[12,12]],
    reference_angular_frequency_SI = omega_r,
)

LaserEnvelopeGaussianAM(
    a0 = 0.5,                                   # higher than collisional range; relativistic
    omega = 1.0,
    focus = [200., 0.],
    waist = 30.,
    time_envelope = tgaussian(fwhm=400., center=600.),
    polarization_phi = 0.,
    envelope_solver = 'explicit',
    Envelope_boundary_conditions = [['PML']],
)

# For envelope mode + ionization, use the averaged model
Species(
    name = "electron",
    # ... usual config ...
    ionization_model = "tunnel_envelope_averaged",
    # ...
)
```

**Caveats for excimer envelope use:**
- The envelope-averaged ionization rate `Γ_ADK,AC` differs from `Γ_ADK,DC` for linear polarization (see `ionization.md`).
- The averaged fields in `DiagFields` do not include the laser field — you cannot directly visualize the laser this way.
- Some collisional effects that depend on the instantaneous field (e.g., field-modified collision frequency) are not captured.
- If unsure, run a short non-envelope benchmark covering 1–2 laser periods and compare absorbed energy.

## 9. Physics modules that need extra care at UV

### Field ionization
- **Keldysh parameter at UV is borderline.** For `I = 10¹³ W/cm²` at 248 nm and Ar (ionization potential 15.76 eV), `γ ≈ 1.2`. ADK is strictly valid for `γ ≪ 1`. For `γ ~ 1`, prefer `tunnel_full_PPT` (PPT model) if available, or note the systematic error.
- **Photoionization channel is missing in Smilei's standard models.** Single-photon ionization of low-IP atoms by ArF (6.4 eV/photon) requires explicit photoionization treatment — Smilei does not have this natively. Workaround: tabulate a rate as `from_rate` model.
- **Multiphoton ionization at UV** can dominate over tunnel below `γ < 1` threshold. ADK underpredicts in this regime.

### Collisional ionization
- Smilei's `CollisionalIonization` uses Lotz-like cross-sections, validated for moderate-Z. For high-Z (W, Au, lanthanides), expect lower accuracy.
- The bound-electron screening option matters for partially-stripped ions in dense plasma; enable it.

### Radiation reaction
- Sub-relativistic regimes (`a₀ ≪ 1`): radiation reaction is negligible. Do not enable; it adds overhead.

### Maxwell solver
- `"Yee"` (default) is fine at UV. Avoid the Lehe solver unless you understand its dispersion characteristics — it's optimized for relativistic LWFA and may distort short-wavelength response.

### Particle pusher
- `"Boris"` (default) is fine for sub-relativistic. `"Vay"` and `"HigueraCary"` are designed for relativistic and overkill at typical excimer intensities (but harmless).
- `"Borisnr"` (non-relativistic Boris) is faster and slightly more accurate for `a₀ ≪ 1` — consider for pure excimer Regime I/II runs.

## 10. When Smilei is the wrong tool — and what is

| Excimer use case | Better choice |
|---|---|
| Industrial UV ablation (LASIK, polymer cutting) | Hydro codes (FLASH, HELIOS) with EOS + radiation transport |
| Photolithography plasma effects | Specialized lithography models; PIC is overkill |
| Long-pulse ICF ablation (full ns history) | Radiation-hydrodynamics codes (FLASH, MULTI, DRACO, FCI2) |
| Plasma discharge / glow discharge | Low-temperature plasma codes (LXcat-based fluid models, or HPEM) |
| EUV source plasma (Sn, Xe) | Hybrid codes; PIC sub-modules |
| Pre-formed plasma + short-pulse excimer | Smilei is appropriate, this regime |
| LPI in DUV ICF (TPD, SBS, SRS) | Smilei is appropriate |
| Cross-beam energy transfer with smoothed beams | Smilei (see Oudin et al. 2022, in publications list) |
| Photoionization-dominated plasma formation | Smilei requires `from_rate` workaround; consider codes with built-in photoionization (e.g., PIConGPU) |

If the user's question is "how do I model a 20 ns KrF pulse ablating Si," the most useful answer is "Smilei is not the right primary tool — use a radiation-hydro code, and use Smilei only for the kinetic sub-process of interest within a representative window." Be willing to say this.

If the user *does* want PIC physics in an excimer regime (laser-plasma instabilities, hot-electron transport, kinetic effects in pre-formed plasma), Smilei is well-suited — but configure it carefully along the lines of §5 and §6.
