# PICMI basics: script recipes

Copy-paste templates. Replace placeholders; don't redesign structure.

## Canonical minimal script

```python
from pywarpx import picmi

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# 1. Grid (selects dimensionality)
grid = picmi.Cartesian2DGrid(
    nx=64, ny=64,
    xmin=-1e-3, xmax=1e-3,
    ymin=-1e-3, ymax=1e-3,
    lower_boundary_conditions=['periodic', 'periodic'],
    upper_boundary_conditions=['periodic', 'periodic'],
    lower_boundary_conditions_particles=['periodic', 'periodic'],
    upper_boundary_conditions_particles=['periodic', 'periodic'],
    warpx_max_grid_size=32,
)

# 2. Solver
solver = picmi.ElectrostaticSolver(grid=grid, method='Multigrid',
                                    required_precision=1e-7)

# 3. Species + distribution + layout
dist = picmi.UniformDistribution(density=1e18, rms_velocity=[c*0.01]*3)
electrons = picmi.Species(particle_type='electron', name='electrons',
                          initial_distribution=dist)
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2], grid=grid)

# 4. Diagnostic
field_diag = picmi.FieldDiagnostic(
    name='diag1', grid=grid, period=100,
    data_list=['rho', 'E', 'phi'],
    write_dir='./diags', warpx_format='openpmd',
)

# 5. Simulation
sim = picmi.Simulation(solver=solver, max_steps=200, verbose=1,
                        warpx_random_seed=42)
sim.add_species(electrons, layout=layout)
sim.add_diagnostic(field_diag)

# 6. Run
sim.step()
```

## Simulation class

```python
sim = picmi.Simulation(
    solver=my_solver,                       # required
    time_step_size=1e-12,                   # required if CFL not set on solver
    max_steps=10000,
    verbose=1,
    particle_shape='linear',                 # 'NGP', 'linear', 'quadratic', 'cubic'
    # Common warpx_* options:
    warpx_current_deposition_algo='direct',  # 'direct', 'esirkepov', 'vay'
    warpx_particle_pusher_algo='boris',      # 'boris', 'vay', 'higuera'
    warpx_random_seed=42,
)

sim.add_species(species, layout=layout)
sim.add_diagnostic(diagnostic)
sim.add_interaction(interaction)             # for collisions/ionization
sim.add_applied_field(applied_field)

sim.step(nsteps=100)
# or
sim.write_input_file('inputs')
```

## Grid classes

### Cartesian 1D / 2D / 3D

```python
# 3D
grid = picmi.Cartesian3DGrid(
    number_of_cells=[64, 64, 256],
    lower_bound=[-1e-3, -1e-3, 0.0],
    upper_bound=[1e-3, 1e-3, 1e-2],
    lower_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    upper_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    lower_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    upper_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    warpx_max_grid_size=64,
    warpx_potential_lo_z='V_RF * sin(2*pi*f_RF*t)',
    warpx_potential_hi_z='0',
)

# 2D: same args with nx/ny, xmin/xmax/ymin/ymax (no z)
# 1D: same with nx, xmin, xmax
```

### Cylindrical (RZ)

```python
grid = picmi.CylindricalGrid(
    nr=128, nz=512,
    n_azimuthal_modes=1,            # 1 = mode 0 only; >1 includes asymmetric modes
    rmin=0.0, rmax=1e-2,
    zmin=0.0, zmax=2e-2,
    bc_rmax='dirichlet',             # rmin at axis; no BC needed
    bc_zmin='dirichlet', bc_zmax='dirichlet',
)
```

### Boundary conditions

Field BCs: `'periodic'`, `'dirichlet'`, `'neumann'`, `'open'`, `'absorbing_silver_mueller'`, `'pec'`, `'pmc'`.

Particle BCs: `'periodic'`, `'absorbing'`, `'reflect'`, `'thermal'`, `'open'`.

Common RF discharge: field `'dirichlet'` + particle `'absorbing'` + `warpx_potential_*` for time-varying voltage.

**For periodic axes**: BOTH field AND particle BCs must be periodic.

## Field solvers

### Electrostatic (low-temp plasma workhorse)

```python
solver = picmi.ElectrostaticSolver(
    grid=grid,
    method='Multigrid',              # 'Multigrid' (MLMG) or 'FFT' (periodic only)
    required_precision=1e-7,
    maximum_iterations=200,
)
```

### Electromagnetic

```python
solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',                    # 'Yee', 'CKC', 'Lehe', 'PSATD', 'ECT'
    cfl=0.99,
)
```

Use `'Yee'` for general EM. `'PSATD'` for laser-plasma (needs `WarpX_FFT=ON`). `'ECT'` for embedded boundaries with internal conductors.

## Species class

```python
electrons = picmi.Species(
    particle_type='electron',          # openPMD species type
    name='electrons',                   # used in diagnostics/collisions
    initial_distribution=my_distribution,
    method='Boris',                     # 'Boris', 'Vay', 'Higuera-Cary', 'free-streaming'
    # WarpX-specific:
    warpx_do_not_push=False,
    warpx_do_not_deposit=False,
    warpx_do_not_gather=False,
    warpx_do_resampling=True,           # essential when ionization is on
    warpx_resampling_trigger_intervals='1000',
    warpx_resampling_trigger_max_avg_ppc=200,
    warpx_save_particles_at_xlo=True,    # record particles crossing xlo boundary
    warpx_save_particles_at_xhi=True,
)
```

### Common `particle_type` values

`'electron'`, `'positron'`, `'proton'`, `'photon'`, `'H'`, `'He'`, `'Ar'`, `'Kr'`, `'Xe'`, `'N'`, `'O'`, `'D'`, `'T'`, `'alpha'`.

### CRITICAL: ions REQUIRE `charge_state`

```python
# BAD — neutral (charge=0); ions won't feel E fields
bad = picmi.Species(particle_type='Ar', name='ions', ...)

# GOOD — Ar+ singly-ionized
good = picmi.Species(particle_type='Ar', charge_state=1, name='ions', ...)
```

This is the #1 Python-side bug. Always specify `charge_state` for ions.

### Frozen neutrals (DSMC pattern)

```python
neutrals = picmi.Species(
    particle_type='Ar', charge_state=0, name='neutrals',
    initial_distribution=picmi.UniformDistribution(density=2.5e25, ...),  # ~1 atm @ 300K
    warpx_do_not_push=True,        # frozen
    warpx_do_not_gather=True,       # no E force
    warpx_do_not_deposit=True,      # no charge contribution
)
```

## Initial distributions

### UniformDistribution

```python
dist = picmi.UniformDistribution(
    density=1e18,                              # m^-3
    lower_bound=[-1e-3, -1e-3, 0.0],
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th, v_th, v_th],          # thermal spread (m/s)
    directed_velocity=[0, 0, c*0.1],           # drift (m/s)
)
```

Thermal speed: `v_th = sqrt(T * q_e / m)` for T in eV.

### AnalyticDistribution

```python
dist = picmi.AnalyticDistribution(
    density_expression='n0 * exp(-(x**2 + y**2) / (2*r0**2))',
    n0=1e18,                                  # passed as kwargs
    r0=5e-4,
    lower_bound=[-1e-3, -1e-3, 0],            # ALWAYS SET — restricts init region
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th, v_th, v_th],
)
```

**Always set bounds** — default is the entire box, wasting memory for localized distributions.

### GaussianBunchDistribution (relativistic beams)

```python
beam = picmi.GaussianBunchDistribution(
    n_physical_particles=1e10,
    rms_bunch_size=[1e-6, 1e-6, 1e-5],
    rms_velocity=[c*1e-4, c*1e-4, c*1e-3],
    centroid_position=[0, 0, 0],
    centroid_velocity=[0, 0, c*0.9999],       # γv, NOT v
)
```

## Layouts

```python
# Deterministic placement:
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2, 4])    # 16 total in 3D

# Random:
layout = picmi.PseudoRandomLayout(n_macroparticles_per_cell=16, seed=42)
```

Particles per cell:
- Cold plasma: 4-16
- Hot plasma w/ thermal noise: 16-64
- Collisional (collision stats matter): 32-128
- Turner et al. benchmark: 200-500

## Diagnostics

### FieldDiagnostic

```python
field_diag = picmi.FieldDiagnostic(
    name='diag1', grid=grid, period=1000,
    data_list=['rho', 'rho_electrons', 'rho_ar_ions', 'E', 'phi'],
    write_dir='./diags',
    warpx_format='openpmd',
    warpx_openpmd_backend='h5',                # 'h5' or 'bp' (ADIOS2)
)
sim.add_diagnostic(field_diag)
```

`data_list` options:
- `'E'`, `'B'`, `'J'`: all three components
- `'rho'`: total charge density
- `'rho_<species>'`: per-species (e.g., `'rho_electrons'`)
- `'phi'`: scalar potential (ES only)

### ParticleDiagnostic

```python
part_diag = picmi.ParticleDiagnostic(
    name='diag1',                              # same name as Field = consolidated openPMD series
    period=1000,
    species=[electrons, ar_ions],               # None = all
    data_list=['position', 'momentum', 'weighting'],
    write_dir='./diags',
    warpx_format='openpmd',
    warpx_random_fraction=None,                 # e.g., 0.01 for 1% subsample
    warpx_plot_filter_function=None,            # e.g., '(uz > 1e7) and (z > 0)'
)
sim.add_diagnostic(part_diag)
```

Set both Field and Particle `name='diag1'` → consolidated openPMD series.

### Checkpoint (for restart)

```python
checkpoint = picmi.Checkpoint(period=10000, write_dir='./diags', name='chk')
sim.add_diagnostic(checkpoint)

# Restart later:
sim = picmi.Simulation(..., warpx_amr_restart='./diags/chk001000')
```

Restarts preserve fields, particles, RNG. They do NOT preserve Python state (callbacks, ML models).

## Constants

```python
picmi.constants.c        # speed of light
picmi.constants.ep0      # vacuum permittivity
picmi.constants.mu0      # vacuum permeability
picmi.constants.q_e      # elementary charge (absolute)
picmi.constants.m_e      # electron mass
picmi.constants.m_p      # proton mass
```

All SI.

## Applied external fields

```python
# Once at t=0:
sim.add_applied_field(picmi.AnalyticInitialField(
    Ez_expression='1e6 * exp(-((z - 0.005)/1e-3)**2)',
))

# Static, every step:
sim.add_applied_field(picmi.ConstantAppliedField(Bz=0.5))

# Time-varying:
sim.add_applied_field(picmi.AnalyticAppliedField(
    Bz_expression='0.5 * (1 + tanh((t - 1e-8) / 1e-9))',
    lower_bound=[-1, -1, 0],
    upper_bound=[1, 1, 0.01],
))
```

String expressions: `x`, `y`, `z`, `t` available; constant `pi`; functions `sin`, `cos`, `exp`, `sqrt`, `log`; ternary as `if(cond, t, f)`. Custom params via kwargs.

## Common pitfalls

| Pitfall | Fix |
|---|---|
| Ions without `charge_state` | Always set `charge_state` for ions |
| `AnalyticDistribution` without bounds | Always set `lower_bound`/`upper_bound` |
| Mismatched periodic BCs | Field + particle BC must both be periodic on the same axis |
| Forgetting `warpx_*` prefix | Kwargs silently ignored as unknown |
| `centroid_velocity` confusion | It's γv (momentum/mass), not raw v |
| Wrong dimensionality install | Install `pywarpx[<dim>]` matching your `Grid` class |

## Verify everything against the live docs

`https://warpx.readthedocs.io/en/latest/usage/python.html` is the ground truth for class names, kwargs, and `warpx_*` options. Always check there before relying on a remembered API.
