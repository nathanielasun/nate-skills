# Collisions and reactions in PICMI

This reference covers the WarpX collision and reaction Python interface — MCC, DSMC, Coulomb, and field ionization. Read this when configuring any process that scatters, ionizes, or reacts particles in a PICMI script. For low-temperature plasma work this is the central file.

## Table of contents

1. [Overview: choosing the right collision type](#overview-choosing-the-right-collision-type)
2. [Background MCC: against a thermal neutral background](#background-mcc-against-a-thermal-neutral-background)
3. [The scattering_processes dictionary](#the-scattering_processes-dictionary)
4. [Cross-section data: warpx-data and file format](#cross-section-data-warpx-data-and-file-format)
5. [DSMC: binary collisions between simulated species](#dsmc-binary-collisions-between-simulated-species)
6. [Coulomb collisions](#coulomb-collisions)
7. [Field ionization (ADK)](#field-ionization-adk)
8. [Putting collisions into a Simulation](#putting-collisions-into-a-simulation)
9. [Worked example: Turner et al. capacitive discharge](#worked-example-turner-et-al-capacitive-discharge)
10. [Worked example: argon DSMC in a microwave field](#worked-example-argon-dsmc-in-a-microwave-field)
11. [Performance and accuracy considerations](#performance-and-accuracy-considerations)
12. [Pitfalls](#pitfalls)

## Overview: choosing the right collision type

| You have... | You want... | Use |
|---|---|---|
| A simulated plasma in a low-density neutral gas | Electron-neutral and ion-neutral collisions; neutrals not depleted | `MCCCollisions` |
| Simulated neutral macroparticles + simulated charged species | Charge-exchange, elastic scattering between simulated species | `DSMCCollisions` |
| Charged-charged collisions within or across species | Relativistic Coulomb scattering | `CoulombCollisions` |
| Strong-field ionization (intense laser, beam wake) | ADK ionization model | `FieldIonization` |
| Photon-particle, multi-photon, or quantum processes | QED (Breit-Wheeler, Schwinger, photon emission) | QED interactions (separate machinery, not covered here) |

These are not exclusive — a typical low-temperature plasma simulation often combines MCC (electron-neutral) with DSMC (neutral-neutral or fast-neutral relaxation) and possibly Coulomb collisions (at higher ionization fraction). They all work simultaneously on the same simulation.

## Background MCC: against a thermal neutral background

`picmi.MCCCollisions` is the most-used collision class for low-temperature plasma. It models collisions between a simulated charged species and an *unsimulated* thermal neutral background (described only by density and temperature).

### Basic usage

```python
mcc_electrons = picmi.MCCCollisions(
    name='mcc_electrons',
    species=electrons,                     # the simulated species
    background_density=9.64e20,             # m^-3 (e.g., He at ~30 mTorr)
    background_temperature=300.0,           # K
    background_mass=4 * picmi.constants.m_p,  # He atomic mass (m_p ~ He/4)
    scattering_processes={
        'elastic': {'cross_section': 'path/to/electron_scattering.dat'},
        'excitation1': {
            'cross_section': 'path/to/excitation_1.dat',
            'energy': 19.82,                # eV, excitation threshold
        },
        'excitation2': {
            'cross_section': 'path/to/excitation_2.dat',
            'energy': 20.61,
        },
        'ionization': {
            'cross_section': 'path/to/ionization.dat',
            'energy': 24.55,                # eV, ionization threshold
            'species': 'he_ions',            # name of the product ion species
        },
    },
    ndt=1,                                  # apply every step
)
```

A separate `MCCCollisions` instance is typically declared for each simulated species that undergoes collisions:

```python
mcc_ions = picmi.MCCCollisions(
    name='mcc_ions',
    species=ions,
    background_density=9.64e20,
    background_temperature=300.0,
    background_mass=4 * picmi.constants.m_p,
    scattering_processes={
        'elastic': {'cross_section': 'path/to/ion_scattering.dat'},
        'back': {'cross_section': 'path/to/ion_back.dat'},
        'charge_exchange': {'cross_section': 'path/to/charge_exchange.dat'},
    },
    ndt=1,
)
```

### Constructor arguments in detail

```python
picmi.MCCCollisions(
    name,                          # required: unique string identifier
    species,                       # required: a picmi.Species instance (one)
    background_density,            # required: m^-3, or string expression
    background_temperature,        # required: K, or string expression
    scattering_processes,          # required: dict described below
    background_mass=None,          # kg; if None, inferred from process
    max_background_density=None,   # required if background_density is an expression
    ndt=1,                          # apply every N steps
)
```

- **`background_density`** can be a float (constant) or a string expression of `(x, y, z, t)`. For a spatially varying gas profile (e.g., a pressure gradient or a localized cooled region), use the expression form. When using an expression, you *must* also specify `max_background_density` for the null-collision-method precomputation.

  ```python
  background_density='1e20 * exp(-(z/0.01)**2)',     # Gaussian profile in z
  max_background_density=1e20,                       # the max value of the above
  ```

- **`background_temperature`** can similarly be a constant or expression. For non-isothermal gas (heating near electrodes, for example).

- **`background_mass`** in kg. For atomic species, use the atomic mass × `picmi.constants.m_p` as a rough approximation (or use scipy's atomic weights for accuracy). Helium ≈ `4*m_p`, argon ≈ `40*m_p`, krypton ≈ `84*m_p`.

- **`ndt`** is the skip-step interval. `ndt=1` (default) runs collisions every step. Higher values reduce cost but space out events; energy-violent processes (ionization, excitation) should use `ndt=1` unless cost demands otherwise.

## The scattering_processes dictionary

The `scattering_processes` argument is the heart of MCC configuration. Its structure:

```python
scattering_processes = {
    '<process_name>': {
        'cross_section': '<path/to/file>',
        # process-specific extra keys:
        'energy': <threshold_eV>,         # for excitation, ionization
        'species': '<product_species>',   # for ionization (the ion product)
    },
    ...
}
```

### Supported process names

- **`'elastic'`** — isotropic elastic scattering in COM frame. No threshold.
- **`'back'`** — back-scattering (180° in COM frame). No threshold.
- **`'charge_exchange'`** — ion + neutral → fast neutral + slow ion (velocities swap). No threshold.
- **`'excitation_N'`** (N integer) — excitation with energy loss; particle loses `energy` eV to internal state. Multiple excitation channels can be declared with `excitation_1`, `excitation_2`, etc.
- **`'ionization'`** — projectile loses `energy` eV, creates a new ion (in `species`) and a new electron. The product ion takes the neutral's pre-collision velocity; the new electron shares the projectile's remaining energy.

The process names are exact strings — `'excitation_1'`, not `'excitation 1'` or `'Excitation_1'`. Case-sensitive.

### Process-specific keys

| Key | Required for | Description |
|---|---|---|
| `cross_section` | all | Path to cross-section file (LXCat-style ASCII) |
| `energy` | excitation, ionization | Threshold energy in eV |
| `species` | ionization | Name of the product ion species (must be declared as a Species) |

Charge-exchange products are *not* declared via the dict — the products are implicit (the same simulated species with swapped velocity, plus an implicit fast neutral that is not represented).

### Example: a complete argon scattering set

```python
ar_cross_section_dir = 'path/to/warpx-data/MCC_cross_sections/argon/'
scattering_processes_e = {
    'elastic': {
        'cross_section': f'{ar_cross_section_dir}/electron_scattering.dat',
    },
    'excitation_1': {
        'cross_section': f'{ar_cross_section_dir}/excitation_1.dat',
        'energy': 11.5,
    },
    'excitation_2': {
        'cross_section': f'{ar_cross_section_dir}/excitation_2.dat',
        'energy': 13.0,
    },
    'ionization': {
        'cross_section': f'{ar_cross_section_dir}/ionization.dat',
        'energy': 15.76,
        'species': 'ar_ions',
    },
}

scattering_processes_i = {
    'elastic': {
        'cross_section': f'{ar_cross_section_dir}/ion_elastic.dat',
    },
    'back': {
        'cross_section': f'{ar_cross_section_dir}/ion_back.dat',
    },
}
```

Note: charge-exchange between Ar+ and Ar at low energies is typically dominant; the `'back'` process is sometimes used in place of explicit charge_exchange in benchmark sets.

## Cross-section data: warpx-data and file format

Cross-section files are *not* shipped with the main WarpX repository. They live in [`BLAST-WarpX/warpx-data`](https://github.com/BLAST-WarpX/warpx-data), a separate repo. For low-temperature plasma work, the relevant directory is `MCC_cross_sections/<species>/`.

To use, clone alongside your input:

```bash
git clone https://github.com/BLAST-WarpX/warpx-data.git
```

Available subdirectories (as of late 2025; check the current repo for additions):
- `MCC_cross_sections/argon/`
- `MCC_cross_sections/helium/`
- `MCC_cross_sections/xenon/`
- `MCC_cross_sections/krypton/`
- `MCC_cross_sections/oxygen/`
- `MCC_cross_sections/nitrogen/`

Each contains files like `electron_scattering.dat`, `ionization.dat`, `excitation_1.dat`, etc. Headers in each file document the source (Hayashi, Biagi, BSR, LXCat) and the process modeled.

### File format

Two-column ASCII text:

```
# Comment lines start with #
# Source: ...
# energy [eV]   cross section [m^2]
1.000000e-02    7.500000e-21
1.500000e-02    7.700000e-21
...
```

Linear interpolation is used between data points. Energies outside the tabulated range use the clamped endpoint values (so a 1 keV electron will use the highest-energy data point's cross section, which may be inaccurate). For high-energy applications, ensure the tabulated range covers the relevant energy span.

### Contributing new cross sections

If you need a species or process not in `warpx-data`, the correct fix is a PR to that repo, not hardcoding cross sections in the WarpX Python script. The convention is:
1. Find the data from a reputable source (LXCat database is the canonical reference).
2. Format as two-column ASCII with a header documenting the source.
3. Add to the appropriate species directory.
4. Open a PR.

## DSMC: binary collisions between simulated species

`picmi.DSMCCollisions` handles binary collisions between two simulated species (or a species with itself). Used when neutrals are simulated as macroparticles — typically at higher gas densities or when neutral dynamics matter (charge exchange producing fast neutrals, neutral-neutral relaxation at high collision frequency).

```python
dsmc_e_n = picmi.DSMCCollisions(
    name='dsmc_electron_neutral',
    species=[electrons, neutrals],          # exactly 2 species (or same species twice for self-coll)
    scattering_processes={
        'elastic': {
            'cross_section': 'path/to/electron_neutral_elastic.dat',
        },
        'ionization': {
            'cross_section': 'path/to/electron_impact_ionization.dat',
            'energy': 15.76,
            'species': ['ar_ions', 'electrons'],   # products (ion + electron)
        },
    },
    ndt=1,
)
```

Key differences from MCC:
- `species` is a list of two (or repeated species for self-collisions).
- Both species must be declared `picmi.Species` objects (DSMC pairs *simulated* particles).
- For `ionization` and similar reactions producing different products, `species` in the process dict is a *list* of product species names.

### When to choose DSMC over MCC

- **Density**: if neutral-neutral collision frequency is high enough that neutral kinetics matter, DSMC needs to track them.
- **Neutral depletion**: ionization in MCC doesn't remove background neutrals (the background is unsimulated). If ionization fraction exceeds ~1% and you care about local neutral density, switch to DSMC.
- **Wall interactions**: thermalized neutrals from wall recombination, neutral pumping, etc. — needs explicit neutral macroparticles, hence DSMC.
- **Fast-neutral physics**: charge exchange produces fast neutrals that may be relevant for sputtering, surface bombardment, or downstream transport. DSMC tracks them; MCC does not.

For typical atmospheric-pressure microwave plasma at low ionization fraction (the user's domain), MCC + frozen-neutral DSMC (neutrals declared but not pushed) is the standard pattern. See the example in `species-and-injection.md`.

### Self-species DSMC

For neutral-neutral relaxation:

```python
dsmc_nn = picmi.DSMCCollisions(
    name='neutral_neutral',
    species=[neutrals, neutrals],
    scattering_processes={
        'elastic': {'cross_section': 'path/to/ar_ar_vhs.dat'},
    },
    ndt=1,
)
```

Pairs within the same species. Used to thermalize neutrals after charge exchange or wall interactions inject non-Maxwellian neutrals.

## Coulomb collisions

`picmi.CoulombCollisions` handles relativistic binary Coulomb scattering (Perez et al. 2012 algorithm) between charged species.

```python
coul = picmi.CoulombCollisions(
    name='coulomb_ei',
    species=[electrons, ions],
    CoulombLog=10.0,             # or None to compute from local plasma parameters
    ndt=1,
)

# Self-species (electron-electron):
coul_ee = picmi.CoulombCollisions(
    name='coulomb_ee',
    species=[electrons, electrons],
    CoulombLog=None,             # auto-compute is usually fine for self-coll
    ndt=1,
)
```

For most low-temp plasma at low ionization fraction (< 0.1%), Coulomb collisions are negligible — electron-neutral collisions dominate. But for higher ionization, or for relativistic beam in plasma, or for thermal equilibration studies, Coulomb collisions matter.

### `CoulombLog` choice

- **Fixed value (e.g., 10)**: appropriate for simulations where the Coulomb log doesn't vary much; gives reproducible behavior; matches some published benchmarks.
- **`None` (auto-compute)**: WarpX computes the Coulomb log from local density and temperature each step. More physically accurate but slower and less reproducible.

In doubt, use a fixed value matching the bulk plasma conditions. For magnetic confinement or astrophysical plasmas with widely-varying density, auto-compute.

## Field ionization (ADK)

Strong-field ionization via the Ammosov-Delone-Krainov (ADK) model. The only ionization model currently in WarpX:

```python
field_ion = picmi.FieldIonization(
    model='ADK',                    # only ADK is currently implemented
    ionized_species=neutrals,       # the species being ionized
    product_species=ions,           # where new ions go (electrons are added to a default species)
)
sim.add_interaction(field_ion)     # or sim.interactions = [field_ion, ...]
```

Used for laser ionization of gas — the laser's electric field strips electrons from atoms. Cross-coupled with the EM solver: the laser's E field is read each step, the ionization probability computed, and new electrons/ions are created.

For low-temp plasma work in atmospheric pressure microwave discharges, field ionization is usually *not* the dominant ionization channel (impact ionization via MCC is). Use field ionization for laser-plasma, attosecond physics, or strong-field studies.

### Charge-state stripping

For multi-electron atoms, ionization can produce multiply-ionized states:

```python
ar_neutral = picmi.Species(particle_type='Ar', charge_state=0, name='ar_neutral', ...)
ar_1 = picmi.Species(particle_type='Ar', charge_state=1, name='ar_1', ...)
ar_2 = picmi.Species(particle_type='Ar', charge_state=2, name='ar_2', ...)

ion_step_1 = picmi.FieldIonization(model='ADK', ionized_species=ar_neutral, product_species=ar_1)
ion_step_2 = picmi.FieldIonization(model='ADK', ionized_species=ar_1, product_species=ar_2)
```

Each ionization step is a separate interaction. The ADK rates depend on the ionization energy of each charge state, which WarpX includes for common atoms.

## Putting collisions into a Simulation

The current convention is to attach collisions as the `collisions` attribute or via `warpx_collisions` keyword:

```python
sim = picmi.Simulation(...)

# Option 1: attribute assignment
sim.collisions = [mcc_electrons, mcc_ions, coul_ee]

# Option 2: keyword argument at construction
sim = picmi.Simulation(..., warpx_collisions=[mcc_electrons, mcc_ions])

# Option 3 (newer convention): add_interaction
sim.add_interaction(mcc_electrons)
sim.add_interaction(mcc_ions)
sim.add_interaction(coul_ee)
```

The exact attribute name and method have shifted across versions. Check the current `pywarpx.picmi.Simulation` docstring if you encounter issues. `sim.collisions = [...]` and `add_interaction()` have both been supported.

Field ionization, since it's an interaction not a collision, generally goes through `sim.add_interaction()`.

## Worked example: Turner et al. capacitive discharge

The canonical low-temperature-plasma PICMI benchmark. From `Examples/Physics_applications/capacitive_discharge/inputs_base_1d_picmi.py`. Abbreviated:

```python
from pywarpx import picmi

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# Case 1 parameters from Turner et al. (Phys. Plasmas 2013)
V_voltage = 450      # V
P_gas = 30e-3 * 133.322   # 30 mTorr in Pa
T_gas = 300          # K
n_gas = P_gas / (1.38e-23 * T_gas)
f_rf = 13.56e6       # Hz
length = 6.7e-2      # m, gap

T_e = 30000          # K, initial electron temperature
T_i = 300            # K
n0 = 2.56e14         # initial plasma density, m^-3

# 1D grid with applied potential
grid = picmi.Cartesian1DGrid(
    number_of_cells=[128],
    lower_bound=[0],
    upper_bound=[length],
    lower_boundary_conditions=['dirichlet'],
    upper_boundary_conditions=['dirichlet'],
    lower_boundary_conditions_particles=['absorbing'],
    upper_boundary_conditions_particles=['absorbing'],
    warpx_potential_hi_x=f'{V_voltage} * cos(2*pi*{f_rf}*t)',
)

solver = picmi.ElectrostaticSolver(grid=grid, method='Multigrid',
                                    required_precision=1e-6, maximum_iterations=1000)

# Species: electrons + He+ ions
v_th_e = (T_e * 1.38e-23 / m_e) ** 0.5
v_th_i = (T_i * 1.38e-23 / (4*1.673e-27)) ** 0.5

elec_dist = picmi.UniformDistribution(density=n0, rms_velocity=[v_th_e]*3,
                                       lower_bound=[0], upper_bound=[length])
ion_dist = picmi.UniformDistribution(density=n0, rms_velocity=[v_th_i]*3,
                                      lower_bound=[0], upper_bound=[length])

electrons = picmi.Species(particle_type='electron', name='electrons',
                           initial_distribution=elec_dist,
                           warpx_save_particles_at_xlo=True,
                           warpx_save_particles_at_xhi=True)
ions = picmi.Species(particle_type='He', charge_state=1, name='he_ions',
                     initial_distribution=ion_dist,
                     warpx_save_particles_at_xlo=True,
                     warpx_save_particles_at_xhi=True)

layout = picmi.GriddedLayout(n_macroparticle_per_cell=[200])

# Collisions
cs_dir = 'warpx-data/MCC_cross_sections/helium/'
mcc_e = picmi.MCCCollisions(
    name='mcc_e',
    species=electrons,
    background_density=n_gas,
    background_temperature=T_gas,
    background_mass=4 * 1.673e-27,
    scattering_processes={
        'elastic': {'cross_section': f'{cs_dir}electron_scattering.dat'},
        'excitation_1': {'cross_section': f'{cs_dir}excitation_1.dat', 'energy': 19.82},
        'ionization': {'cross_section': f'{cs_dir}ionization.dat',
                       'energy': 24.55, 'species': 'he_ions'},
    },
)
mcc_i = picmi.MCCCollisions(
    name='mcc_i',
    species=ions,
    background_density=n_gas,
    background_temperature=T_gas,
    background_mass=4 * 1.673e-27,
    scattering_processes={
        'elastic': {'cross_section': f'{cs_dir}ion_scattering.dat'},
        'back': {'cross_section': f'{cs_dir}ion_back.dat'},
    },
)

# Diagnostics
field_diag = picmi.FieldDiagnostic(name='diag1', grid=grid, period=1000,
                                    data_list=['rho', 'E', 'phi'], write_dir='./diags',
                                    warpx_format='openpmd')
part_diag = picmi.ParticleDiagnostic(name='diag1', period=1000,
                                      species=[electrons, ions],
                                      data_list=['position', 'momentum', 'weighting'],
                                      write_dir='./diags', warpx_format='openpmd')

# Simulation
sim = picmi.Simulation(solver=solver, time_step_size=1/(400*f_rf),
                        max_steps=50000, verbose=1)
sim.add_species(electrons, layout=layout)
sim.add_species(ions, layout=layout)
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)
sim.collisions = [mcc_e, mcc_i]

sim.step()
```

This script reproduces Turner et al. (2013) Case 1: a 13.56 MHz capacitive helium discharge between parallel plates. The runtime is hours on a single core; minutes on many cores.

## Worked example: argon DSMC in a microwave field

A sketch of a DSMC setup (full script would be longer):

```python
# Both species simulated as particles
electrons = picmi.Species(particle_type='electron', name='electrons',
                           initial_distribution=picmi.UniformDistribution(
                               density=1e16, ...))
ar_neutrals = picmi.Species(particle_type='Ar', charge_state=0, name='neutrals',
                             initial_distribution=picmi.UniformDistribution(
                                 density=2.45e25, ...),    # ~1 atm at 300 K
                             warpx_do_not_push=False,        # let them push
                             warpx_do_not_deposit=True,       # no charge contribution
                             warpx_do_not_gather=True)        # no E force
ar_ions = picmi.Species(particle_type='Ar', charge_state=1, name='ar_ions', ...)

# Electron-neutral DSMC
dsmc_en = picmi.DSMCCollisions(
    name='dsmc_en',
    species=[electrons, ar_neutrals],
    scattering_processes={
        'elastic': {'cross_section': 'warpx-data/MCC_cross_sections/argon/electron_scattering.dat'},
        'ionization': {'cross_section': 'warpx-data/MCC_cross_sections/argon/ionization.dat',
                       'energy': 15.76,
                       'species': ['ar_ions', 'electrons']},
    },
    ndt=1,
)

# Ion-neutral DSMC (charge exchange, elastic)
dsmc_in = picmi.DSMCCollisions(
    name='dsmc_in',
    species=[ar_ions, ar_neutrals],
    scattering_processes={
        'charge_exchange': {'cross_section': '...'},
        'elastic': {'cross_section': '...'},
    },
    ndt=1,
)

# Neutral-neutral DSMC (thermal relaxation)
dsmc_nn = picmi.DSMCCollisions(
    name='dsmc_nn',
    species=[ar_neutrals, ar_neutrals],
    scattering_processes={
        'elastic': {'cross_section': 'warpx-data/MCC_cross_sections/argon/ar_ar_vhs.dat'},
    },
    ndt=10,   # neutral-neutral is expensive; can skip steps
)

sim.collisions = [dsmc_en, dsmc_in, dsmc_nn]
```

The mix of MCC-style cross-section files and DSMC-style pairing handles atmospheric-pressure plasmas where neutral dynamics need explicit treatment.

## Performance and accuracy considerations

### Cost scaling

- MCC cost is approximately proportional to the number of simulated charged particles (one check per particle per step). For typical low-temp plasma, this is a small fraction of total cost.
- DSMC cost scales with particles per cell squared (pairing). With 500+ particles per cell, DSMC can dominate runtime. Use sorting (`warpx_sort_intervals`) and consider resampling to control particle counts.
- Coulomb cost scales similarly to DSMC.
- Field ionization cost is per-particle, similar to MCC.

### Accuracy: cell-size and time-step constraints

For collisions to be accurate:
- The cell size should be comparable to or smaller than the mean free path between collisions. Otherwise, intra-cell pairing is statistically misleading.
- The time step should be small enough that the probability of multiple collisions per step is small. Typically `ν * dt < 0.1` where `ν` is the collision frequency.

For a typical low-temperature argon discharge at 30 mTorr, 1 eV electrons: λ_mfp ≈ 1 cm, ν ≈ 1e8 Hz, so cells of 1 mm and dt of 1e-10 s comfortably satisfy both.

For 1 atm argon: λ_mfp ≈ 100 μm, ν ≈ 1e10 Hz — cells must be much finer (~ 30 μm), and dt smaller still.

### `ndt` tradeoffs

`ndt=N` runs the collision every N steps but with N times the per-particle probability, so the time-averaged rate is correct. But:
- Energy violently transferred per collision is N times larger than it would be per step. This can produce spurious heating events in cells.
- Resolution of transient phenomena (waves, sheaths) is degraded.

For low-temp plasma steady-state, `ndt=1` is the safe default. For Coulomb collisions in hot plasmas where the rate is slow, `ndt=10` or even higher can be acceptable.

## Pitfalls

- **Cross-section file paths**: relative paths are interpreted relative to where you launch WarpX from, not relative to the script. Use absolute paths or `os.path.abspath()` to avoid the "file not found" failure mode at startup.

- **Background mass not specified**: if you omit `background_mass` for MCC, WarpX uses a default. For ions colliding with their parent neutral (charge exchange) the default may be wrong. Always specify explicitly.

- **Ionization product species not registered**: if `scattering_processes['ionization']['species']` references a species name that wasn't added via `sim.add_species()`, the simulation will fail or silently lose the products. Always declare the product species before adding the collision.

- **Product weights and charge conservation**: when ionization creates new electrons and ions, both products should carry the same weight as the parent neutral (or the projectile, depending on convention). If you see net charge accumulating in the simulation that shouldn't, check the weights of ionization products.

- **Resampling + ionization interaction**: high ionization rates with no resampling lead to memory exhaustion in O(1000) steps. Enable resampling on the electron and ion species before running anything serious.

- **MCC null-collision validity**: if `max_background_density × σ_max × v_max × dt > 1`, the null-collision approximation fails. WarpX should warn about this; if it doesn't, manually check that the maximum collision probability is < 1.

- **Cell size vs mean free path**: if you run at atmospheric pressure with mm-scale cells, collisions are happening many times per cell-traversal, which is *fine* for MCC (each macroparticle is sampled independently). But for DSMC the pairing logic assumes that pairs are physically nearby, which becomes lossy at mm-cells with μm-mfp. Refine the grid.

- **Coulomb log auto-compute can NaN**: if local density is zero (a vacuum cell), the auto-compute can produce NaN. Either set a fixed Coulomb log or ensure all cells have non-zero density.

- **DSMC self-species with very different macroparticle weights**: the Nanbu/Takizuka-Abe weight correction handles weight mismatches between species pairs, but for self-species pairing, large weight variances within the species can produce statistical bias. Resample to even out weights before relying on DSMC.

- **MCC + applied periodic potential at high frequency**: if the RF period is short compared to the collision frequency, MCC may sample collision events at favorable phases and bias the result. Either use `ndt=1` and resolve the RF period with many time steps, or verify your benchmark against a known result before publishing.
