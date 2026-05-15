# Collisions: MCC, DSMC, Coulomb, ionization

The low-temp plasma core. Copy patterns and adapt.

## Decision table

| Setting | Use |
|---|---|
| Charged species in unsimulated thermal neutral gas | `picmi.MCCCollisions` |
| Binary collisions between two simulated species | `picmi.DSMCCollisions` |
| Charged-charged scattering | `picmi.CoulombCollisions` |
| Strong-field ionization (laser) | `picmi.FieldIonization` |

These coexist. Typical low-temp plasma at low ionization: MCC alone. Add Coulomb at high ionization. Add DSMC if neutrals are simulated.

## MCC: against thermal background

```python
mcc_e = picmi.MCCCollisions(
    name='mcc_electrons',
    species=electrons,                          # ONE species (not list)
    background_density=9.64e20,                  # m^-3 OR string expression
    background_temperature=300.0,                # K OR string expression
    background_mass=4 * picmi.constants.m_p,    # kg (ALWAYS set explicitly)
    scattering_processes={
        'elastic': {
            'cross_section': 'path/to/electron_scattering.dat',
        },
        'excitation_1': {                        # underscore_N for multiple channels
            'cross_section': 'path/to/excitation_1.dat',
            'energy': 19.82,                     # eV threshold
        },
        'ionization': {
            'cross_section': 'path/to/ionization.dat',
            'energy': 24.55,
            'species': 'he_ions',                # product ion species NAME (string)
        },
    },
    ndt=1,                                       # apply every N steps
)
```

Declare a separate `MCCCollisions` for each charged species.

## scattering_processes dict structure

```python
scattering_processes = {
    '<process_name>': {
        'cross_section': '<file_path>',
        # extras (process-specific):
        'energy': <threshold_eV>,                # for excitation, ionization
        'species': '<product_species>',          # for ionization
    },
}
```

## Process names (exact strings)

| Name | Description | Extra keys |
|---|---|---|
| `'elastic'` | Isotropic elastic in COM | — |
| `'back'` | Back-scattering (180° in COM) | — |
| `'charge_exchange'` | ion+neutral → fast neutral + slow ion | — |
| `'excitation_N'` (N=1,2,...) | Excitation with energy loss | `'energy'` |
| `'ionization'` | Creates new ion + new electron | `'energy'`, `'species'` |

Case-sensitive. `'excitation_1'`, not `'Excitation_1'`.

## Spatially-varying background

```python
mcc_e = picmi.MCCCollisions(
    ...,
    background_density='1e20 * exp(-(z/0.01)**2)',
    max_background_density=1e20,                # REQUIRED when density is expression
    background_temperature=300.0,
    ...
)
```

## DSMC: simulated species pairs

```python
dsmc_en = picmi.DSMCCollisions(
    name='dsmc_electron_neutral',
    species=[electrons, neutrals],               # list of TWO (or same species twice)
    scattering_processes={
        'elastic': {
            'cross_section': 'path/to/electron_neutral_elastic.dat',
        },
        'ionization': {
            'cross_section': 'path/to/electron_impact_ionization.dat',
            'energy': 15.76,
            'species': ['ar_ions', 'electrons'], # LIST for DSMC (ion + electron)
        },
    },
    ndt=1,
)
```

Differences from MCC:
- `species` is a **list of two** (or repeated for self-coll)
- Both species must be added to the simulation as `picmi.Species`
- Ionization `species` is a **list** of product species names

## Coulomb collisions

```python
coul_ei = picmi.CoulombCollisions(
    name='coulomb_ei',
    species=[electrons, ions],
    CoulombLog=10.0,                # or None to auto-compute (slower but more accurate)
    ndt=1,
)

# Self-species (electron-electron):
coul_ee = picmi.CoulombCollisions(
    name='coulomb_ee',
    species=[electrons, electrons],
    CoulombLog=None,
    ndt=1,
)
```

For typical low-temp plasma at < 0.1% ionization: skip Coulomb. Electron-neutral dominates.

## Field ionization (ADK only)

```python
field_ion = picmi.FieldIonization(
    model='ADK',                    # only ADK is implemented
    ionized_species=neutrals,
    product_species=ions,
)
sim.add_interaction(field_ion)
```

For laser ionization. Not typically needed for microwave low-temp plasma (impact ionization via MCC dominates).

### Multi-electron stripping

```python
ar_neutral = picmi.Species(particle_type='Ar', charge_state=0, name='ar_neutral', ...)
ar_1 = picmi.Species(particle_type='Ar', charge_state=1, name='ar_1', ...)
ar_2 = picmi.Species(particle_type='Ar', charge_state=2, name='ar_2', ...)

ion_step_1 = picmi.FieldIonization(model='ADK', ionized_species=ar_neutral,
                                    product_species=ar_1)
ion_step_2 = picmi.FieldIonization(model='ADK', ionized_species=ar_1,
                                    product_species=ar_2)
```

Each step a separate `add_interaction()`.

## Attaching to simulation

```python
# Most robust pattern:
sim.collisions = [mcc_e, mcc_i, coul_ee]

# Or via add_interaction (newer):
sim.add_interaction(mcc_e)
sim.add_interaction(mcc_i)

# FieldIonization always uses add_interaction:
sim.add_interaction(field_ion)
```

## Cross-section data

NOT shipped with WarpX. Lives in `BLAST-WarpX/warpx-data`:

```bash
git clone https://github.com/BLAST-WarpX/warpx-data.git
```

Subdirectories: `MCC_cross_sections/argon/`, `MCC_cross_sections/helium/`, `MCC_cross_sections/krypton/`, `MCC_cross_sections/xenon/`, `MCC_cross_sections/oxygen/`, `MCC_cross_sections/nitrogen/`.

### File format

Two-column ASCII:
```
# Source: ...
# Energy [eV]   Cross-section [m^2]
1.000000e-02    7.500000e-21
...
```

### Path conventions

```python
# BAD — relative path breaks if launching from elsewhere
'cross_section': 'electron_scattering.dat'

# GOOD — absolute path
import os
CS_DIR = os.path.abspath('./warpx-data/MCC_cross_sections/argon')
'cross_section': f'{CS_DIR}/electron_scattering.dat'
```

Paths are interpreted vs launch directory, NOT script directory. Use absolute.

## Worked example: Turner et al. capacitive discharge

The canonical low-temp plasma PICMI benchmark (Phys. Plasmas 20, 013507, 2013). Case 1: 13.56 MHz helium discharge.

```python
from pywarpx import picmi

c = picmi.constants.c
m_e = picmi.constants.m_e

# Turner Case 1 parameters
V_voltage = 450               # V
P_gas = 30e-3 * 133.322       # 30 mTorr in Pa
T_gas = 300                   # K
n_gas = P_gas / (1.38e-23 * T_gas)
f_rf = 13.56e6                # Hz
length = 6.7e-2               # m

T_e = 30000                   # K initial
T_i = 300
n0 = 2.56e14                  # m^-3
v_th_e = (T_e * 1.38e-23 / m_e) ** 0.5
v_th_i = (T_i * 1.38e-23 / (4 * 1.673e-27)) ** 0.5

# Grid + solver
grid = picmi.Cartesian1DGrid(
    number_of_cells=[128],
    lower_bound=[0], upper_bound=[length],
    lower_boundary_conditions=['dirichlet'],
    upper_boundary_conditions=['dirichlet'],
    lower_boundary_conditions_particles=['absorbing'],
    upper_boundary_conditions_particles=['absorbing'],
    warpx_potential_hi_x=f'{V_voltage} * cos(2*pi*{f_rf}*t)',
)
solver = picmi.ElectrostaticSolver(grid=grid, method='Multigrid',
                                    required_precision=1e-6)

# Species
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
cs = 'warpx-data/MCC_cross_sections/helium/'
mcc_e = picmi.MCCCollisions(
    name='mcc_e', species=electrons,
    background_density=n_gas, background_temperature=T_gas,
    background_mass=4 * 1.673e-27,
    scattering_processes={
        'elastic': {'cross_section': f'{cs}electron_scattering.dat'},
        'excitation_1': {'cross_section': f'{cs}excitation_1.dat', 'energy': 19.82},
        'ionization': {'cross_section': f'{cs}ionization.dat',
                       'energy': 24.55, 'species': 'he_ions'},
    },
)
mcc_i = picmi.MCCCollisions(
    name='mcc_i', species=ions,
    background_density=n_gas, background_temperature=T_gas,
    background_mass=4 * 1.673e-27,
    scattering_processes={
        'elastic': {'cross_section': f'{cs}ion_scattering.dat'},
        'back': {'cross_section': f'{cs}ion_back.dat'},
    },
)

# Simulation
sim = picmi.Simulation(solver=solver, time_step_size=1/(400*f_rf),
                        max_steps=50000, verbose=1)
sim.add_species(electrons, layout=layout)
sim.add_species(ions, layout=layout)
sim.collisions = [mcc_e, mcc_i]

sim.step()
```

Adapt for argon/krypton/xenon by swapping species + cross sections.

## Common pitfalls

| Pitfall | Fix |
|---|---|
| Cross-section path wrong | Use absolute paths |
| Missing `background_mass` | Always set explicitly (kg) |
| Ionization product species not registered | Add product species before declaring collision |
| No resampling with ionization | Enable on charged species (memory exhausts in ~1000 steps otherwise) |
| `max_background_density` missing for expression density | Required — set the max value |
| `ndt > 1` with ionization | Use `ndt=1` for energy-sensitive simulations |
