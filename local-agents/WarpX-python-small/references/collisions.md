# Collisions and reactions

The central reference for low-temperature plasma work. MCC, DSMC, Coulomb, and field ionization.

## Decision tree

| Setting | Use |
|---|---|
| Charged species in unsimulated thermal neutral gas | `picmi.MCCCollisions` |
| Two simulated species (incl. simulated neutrals) | `picmi.DSMCCollisions` |
| Charged-charged scattering | `picmi.CoulombCollisions` |
| Strong-field ionization (laser) | `picmi.FieldIonization` |

These can coexist. Typical low-temp plasma at low ionization: MCC alone. At atmospheric pressure: MCC + frozen-neutral DSMC. High ionization: + Coulomb.

## MCC: against thermal background

`picmi.MCCCollisions` is the most-used class for low-temp plasma. Models collisions between a simulated charged species and an *unsimulated* thermal background described by density + temperature.

```python
mcc_e = picmi.MCCCollisions(
    name='mcc_electrons',
    species=electrons,                          # ONE species (not list)
    background_density=9.64e20,                  # m^-3, or expression of (x,y,z,t)
    background_temperature=300.0,                # K, or expression
    background_mass=4 * picmi.constants.m_p,    # kg (always set explicitly)
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
            'species': 'he_ions',                # product ion species name
        },
    },
    ndt=1,                                       # apply every N steps
)
```

Declare a separate `MCCCollisions` for each charged species (e.g., one for electrons, one for ions).

### scattering_processes dict structure

```python
scattering_processes = {
    '<process_name>': {
        'cross_section': '<file_path>',
        # process-specific keys:
        'energy': <threshold_eV>,                # for excitation, ionization
        'species': '<product_species>',          # for ionization
    },
}
```

### Supported process names

| Name | Description | Extra keys |
|---|---|---|
| `'elastic'` | Isotropic elastic in COM | — |
| `'back'` | Back-scattering (180° in COM) | — |
| `'charge_exchange'` | ion+neutral → fast neutral + slow ion | — |
| `'excitation_N'` (N=1,2,...) | Excitation with energy loss | `'energy'` |
| `'ionization'` | Creates new ion + new electron | `'energy'`, `'species'` |

Names are exact case-sensitive strings. `'excitation_1'`, not `'excitation 1'` or `'Excitation_1'`.

### Spatially-varying background

```python
mcc_e = picmi.MCCCollisions(
    ...,
    background_density='1e20 * exp(-(z/0.01)**2)',
    max_background_density=1e20,                # REQUIRED when density is expression
    background_temperature=300.0,
    ...
)
```

`max_background_density` must be set so the null-collision precomputation has an upper bound.

## DSMC: simulated species pairs

For binary collisions between two simulated species (or self-collisions):

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
            'species': ['ar_ions', 'electrons'], # LIST for DSMC (product ion + electron)
        },
    },
    ndt=1,
)
```

Key differences from MCC:
- `species` is a list of **two** (or repeated for self-coll).
- Both species must be `picmi.Species` instances added to the simulation.
- Ionization product `species` is a **list** (DSMC creates both products explicitly).

### When to use DSMC over MCC

- Neutral kinetics matter (high pressure, neutral-neutral relaxation).
- Ionization fraction > 1% and you care about local neutral density.
- Wall recombination produces non-Maxwellian neutrals.
- Charge-exchange fast neutrals matter for sputtering/transport.

For low ionization at low pressure: MCC is sufficient and much cheaper.

### Frozen-neutral pattern (atmospheric-pressure plasma)

Common for low-temp plasma where neutrals are present in DSMC but otherwise stationary:

```python
neutrals = picmi.Species(
    particle_type='Ar', charge_state=0, name='neutrals',
    initial_distribution=picmi.UniformDistribution(density=2.5e25, ...),  # ~1 atm @ 300K
    warpx_do_not_push=True,        # frozen
    warpx_do_not_gather=True,       # no E force
    warpx_do_not_deposit=True,      # no charge contribution
)
```

The neutrals exist as collision targets for DSMC pairing but otherwise don't participate.

## Coulomb collisions

Perez et al. 2012 algorithm, between charged species:

```python
coul_ei = picmi.CoulombCollisions(
    name='coulomb_ei',
    species=[electrons, ions],
    CoulombLog=10.0,            # or None to compute from local plasma parameters
    ndt=1,
)

coul_ee = picmi.CoulombCollisions(
    name='coulomb_ee',
    species=[electrons, electrons],     # self-coll: same species twice
    CoulombLog=None,                     # auto-compute is fine for self-coll
    ndt=1,
)
```

For typical low-temp plasma at < 0.1% ionization: Coulomb collisions are negligible vs electron-neutral. Skip them.

For higher ionization, magnetic confinement, or hot plasma: include them.

Fixed `CoulombLog`: more reproducible, matches benchmarks. `None`: auto-computed each step, more accurate, can NaN in vacuum cells.

## Field ionization (ADK only)

```python
field_ion = picmi.FieldIonization(
    model='ADK',                    # only ADK is currently implemented
    ionized_species=neutrals,        # what's being ionized
    product_species=ions,            # where new ions go
)
sim.add_interaction(field_ion)
```

Used for laser-induced ionization. For low-temp microwave plasma where impact ionization dominates, FieldIonization is usually NOT needed.

### Charge-state stripping (multi-electron)

```python
ar_neutral = picmi.Species(particle_type='Ar', charge_state=0, name='ar_neutral', ...)
ar_1 = picmi.Species(particle_type='Ar', charge_state=1, name='ar_1', ...)
ar_2 = picmi.Species(particle_type='Ar', charge_state=2, name='ar_2', ...)

ion_step_1 = picmi.FieldIonization(model='ADK', ionized_species=ar_neutral, product_species=ar_1)
ion_step_2 = picmi.FieldIonization(model='ADK', ionized_species=ar_1, product_species=ar_2)
```

Each step a separate interaction.

## Attaching collisions to the simulation

Three equivalent ways (current versions support all; the exact preferred one shifts):

```python
# Option 1: attribute assignment
sim.collisions = [mcc_e, mcc_i, coul_ee]

# Option 2: keyword at construction
sim = picmi.Simulation(..., warpx_collisions=[mcc_e, mcc_i])

# Option 3 (newer): add_interaction
sim.add_interaction(mcc_e)
sim.add_interaction(mcc_i)
sim.add_interaction(coul_ee)
```

`FieldIonization` always uses `add_interaction`.

## Cross-section data

NOT shipped with WarpX. Lives in the separate `BLAST-WarpX/warpx-data` repo:

```bash
git clone https://github.com/BLAST-WarpX/warpx-data.git
```

Subdirectories (check repo for current list): `MCC_cross_sections/argon/`, `MCC_cross_sections/helium/`, `MCC_cross_sections/xenon/`, `MCC_cross_sections/krypton/`, `MCC_cross_sections/oxygen/`, `MCC_cross_sections/nitrogen/`.

### File format

Two-column ASCII:

```
# Comments start with #
# Source: ...
# Energy [eV]   Cross-section [m^2]
1.000000e-02    7.500000e-21
1.500000e-02    7.700000e-21
...
```

Linear interpolation between points. Energies outside the range use clamped endpoints.

### Path conventions

```python
# BAD — relative path breaks if you launch from elsewhere
'cross_section': 'electron_scattering.dat'

# BETTER — relative to a documented base
'cross_section': 'warpx-data/MCC_cross_sections/argon/electron_scattering.dat'

# BEST for production scripts — absolute
import os
CS_DIR = os.path.abspath('./warpx-data/MCC_cross_sections/argon')
'cross_section': f'{CS_DIR}/electron_scattering.dat'
```

Cross-section file paths are interpreted relative to where you launch WarpX, not relative to the script. Absolute paths avoid surprises.

## Worked example: Turner et al. capacitive discharge

The canonical low-temp plasma PICMI benchmark. Case 1 from Turner et al. (Phys. Plasmas 20, 013507, 2013) — 13.56 MHz helium discharge.

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
T_e = 30000                   # K (initial)
T_i = 300                     # K
n0 = 2.56e14                  # m^-3 (initial plasma density)

v_th_e = (T_e * 1.38e-23 / m_e) ** 0.5
v_th_i = (T_i * 1.38e-23 / (4 * 1.673e-27)) ** 0.5     # He, mass ~ 4*m_p

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

sim = picmi.Simulation(solver=solver, time_step_size=1/(400*f_rf),
                        max_steps=50000, verbose=1)
sim.add_species(electrons, layout=layout)
sim.add_species(ions, layout=layout)
sim.collisions = [mcc_e, mcc_i]

sim.step()
```

This script reproduces Turner Case 1. Adapt for argon/krypton/xenon by swapping species + cross sections; adapt for different driving frequencies/voltages by changing `f_rf` and `V_voltage`.

## Common pitfalls

- **Cross-section file path wrong** → "file not found" at startup. Use absolute paths.
- **Missing `background_mass`** → wrong default used; collision kinematics off. Always set explicitly.
- **Ionization `species` not registered** → product ions lost or crash. Add the product species before declaring the collision.
- **No resampling with ionization** → particle count grows unboundedly, memory exhausted in ~1000 steps. Enable resampling on charged species.
- **`max_background_density` missing** when `background_density` is an expression → null-collision probability unbounded → init error.
- **MCC + small cells at high pressure** → DSMC pairing logic breaks down (assumes pairs are physically near). Refine grid.
- **DSMC self-species with wildly different weights** → statistical bias. Resample first.

## Performance notes

- **MCC**: cost scales with charged particle count. Usually small fraction of total.
- **DSMC**: cost scales with particles-per-cell-squared (pairing). Can dominate at 500+ ppc.
- **Coulomb**: similar to DSMC.
- **Cell vs MFP**: cell size should be ≲ mean free path. For 30 mTorr Ar 1eV: λ ~ 1cm. For 1 atm Ar: λ ~ 100μm — cells must be much finer.
- **ndt > 1**: reduces cost but spaces out events; energy-violent processes (ionization) should use `ndt=1`.
