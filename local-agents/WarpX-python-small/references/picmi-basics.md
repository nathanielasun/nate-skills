# PICMI basics: Simulation, Grid, Solver, Species

Compact reference for setting up a PICMI script. For collisions see `collisions.md`. For diagnostics see `output-and-analysis.md`.

## A canonical minimal script

The skeleton. Adapt this for new scripts:

```python
from pywarpx import picmi

# 1. Constants
c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# 2. Grid (selects dimensionality)
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

# 3. Solver
solver = picmi.ElectrostaticSolver(grid=grid, method='Multigrid',
                                    required_precision=1e-7)

# 4. Species + distribution + layout
distribution = picmi.UniformDistribution(
    density=1e18,
    rms_velocity=[c*0.01, c*0.01, c*0.01],
)
electrons = picmi.Species(
    particle_type='electron', name='electrons',
    initial_distribution=distribution,
)
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2], grid=grid)

# 5. Simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=200,
    verbose=1,
    warpx_random_seed=42,
)
sim.add_species(electrons, layout=layout)

# 6. Run
sim.step()
```

## Simulation class

```python
sim = picmi.Simulation(
    solver=my_solver,                       # required
    time_step_size=1e-12,                   # required if CFL not set on solver
    max_steps=10000,                         # one of: max_steps, max_time, or step() loop
    max_time=None,
    verbose=1,
    particle_shape='linear',                # 'NGP', 'linear', 'quadratic', 'cubic'
    gamma_boost=None,
    # WarpX-specific (selection — many more exist):
    warpx_current_deposition_algo='direct',  # 'direct', 'esirkepov', 'vay'
    warpx_particle_pusher_algo='boris',      # 'boris', 'vay', 'higuera'
    warpx_grid_type='staggered',             # 'collocated', 'staggered', 'hybrid'
    warpx_random_seed=42,
    warpx_amrex_use_gpu_aware_mpi=True,
    warpx_load_balance_intervals='100',
)

sim.add_species(species, layout=layout)
sim.add_laser(laser, injection_method)
sim.add_diagnostic(diagnostic)
sim.add_interaction(interaction)             # for collisions/ionization
sim.add_applied_field(applied_field)

sim.step(nsteps=100)                          # advance
# or
sim.write_input_file('inputs')                # generate C++ input file
```

## Grid classes

Pick exactly one — the choice determines dimensionality.

### Cartesian3DGrid

```python
grid = picmi.Cartesian3DGrid(
    number_of_cells=[64, 64, 256],
    lower_bound=[-1e-3, -1e-3, 0.0],
    upper_bound=[1e-3, 1e-3, 1e-2],
    lower_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    upper_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    lower_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    upper_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    warpx_max_grid_size=64,
    warpx_blocking_factor=8,
    # Time-varying potential on Dirichlet boundary:
    warpx_potential_lo_z='V_RF * sin(2*pi*f_RF*t)',
    warpx_potential_hi_z='0',
)
```

Other classes: `Cartesian2DGrid` (use `nx`, `ny`, `xmin`/`xmax`, `ymin`/`ymax`), `Cartesian1DGrid` (`nx`, `xmin`/`xmax`), `CylindricalGrid` (`nr`, `nz`, `n_azimuthal_modes`, `rmin`/`rmax`, `zmin`/`zmax`).

### Boundary conditions

Field BCs: `'periodic'`, `'dirichlet'`, `'neumann'`, `'open'`, `'absorbing_silver_mueller'`, `'pec'`, `'pmc'`.

Particle BCs: `'periodic'`, `'absorbing'`, `'reflect'`, `'thermal'`, `'open'`.

For typical RF discharge: field `'dirichlet'` + particle `'absorbing'`, with `warpx_potential_lo_x` / `warpx_potential_hi_x` for time-varying voltage.

For periodic axes, BOTH field AND particle BCs must be periodic on that axis.

### Domain decomposition options

- `warpx_max_grid_size=32` (CPU) or `128` (GPU) — max tile size
- `warpx_blocking_factor=8` — tile sizes are multiples of this
- `warpx_numprocs=[2, 2, 4]` — force exact MPI decomposition (must multiply to total ranks)

## Field solvers

### Electromagnetic

```python
em_solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',          # 'Yee', 'CKC', 'Lehe', 'PSATD', 'PSTD', 'GPSTD', 'DS', 'ECT'
    cfl=0.99,
)
```

Use `'Yee'` for general EM, `'PSATD'` for laser-plasma (needs `WarpX_FFT=ON`), `'ECT'` for embedded boundaries with internal conductors.

### Electrostatic (workhorse for low-temp plasma)

```python
es_solver = picmi.ElectrostaticSolver(
    grid=grid,
    method='Multigrid',     # 'Multigrid' (MLMG) or 'FFT' (periodic only)
    required_precision=1e-7,
    maximum_iterations=200,
)
```

`'Multigrid'` handles arbitrary BCs and embedded boundaries. `'FFT'` is faster but periodic-only.

### Hybrid (Ohm's law)

```python
hybrid_solver = picmi.HybridPICSolver(
    grid=grid,
    Te=10.0,                # eV
    n0_reference_density=1e18,
    plasma_resistivity='0.0',
    substeps=10,
)
```

For ion-kinetic problems where electrons are a massless fluid.

## Embedded boundaries

```python
eb = picmi.EmbeddedBoundary(
    implicit_function='-(((x - 0)**2 + (y - 0)**2) - r_inner**2)',
    r_inner=2e-3,
    potential='V_bias',
    V_bias='100.0',
)
# Pass to simulation:
sim = picmi.Simulation(..., warpx_embedded_boundary=eb)
```

Positive `implicit_function` values are inside the simulation domain; negative are inside the solid.

For STL geometries:
```python
eb = picmi.EmbeddedBoundary(stl_file='/path/to/geometry.stl', potential='100.0')
```

**Not supported**: RZ/RCYLINDER/RSPHERE geometries, PSATD solver. Works with ES Multigrid, EM Yee/CKC, and ECT.

## Species class

```python
electrons = picmi.Species(
    particle_type='electron',                # use openPMD species type strings
    name='electrons',                         # used in diagnostics and references
    initial_distribution=my_distribution,
    method='Boris',                           # 'Boris', 'Vay', 'Higuera-Cary', 'free-streaming'
    particle_shape='linear',
    # WarpX-specific:
    warpx_do_not_push=False,
    warpx_do_not_deposit=False,
    warpx_do_not_gather=False,
    warpx_do_resampling=True,                 # essential for ionization runaway
    warpx_resampling_trigger_intervals='1000',
    warpx_resampling_trigger_max_avg_ppc=200,
    warpx_save_particles_at_xlo=True,
    warpx_save_particles_at_xhi=True,
)
```

### Common `particle_type` values

`'electron'`, `'positron'`, `'proton'`, `'photon'`, `'H'`, `'He'`, `'Ar'`, `'Kr'`, `'Xe'`, `'N'`, `'O'`, `'D'`, `'T'`, `'alpha'`.

### CRITICAL: ions need `charge_state`

```python
# BAD — neutral species (charge=0), doesn't respond to E fields
bad_ions = picmi.Species(particle_type='Ar', name='ions', ...)

# GOOD — singly-ionized Ar+
good_ions = picmi.Species(particle_type='Ar', charge_state=1, name='ions', ...)
```

This is the single most common Python-side bug. The simulation runs without error but ions don't feel electric fields.

### Explicit charge/mass (when not using particle_type)

```python
heavy_e = picmi.Species(
    name='heavy_electron',
    charge=-picmi.constants.q_e,
    mass=10 * picmi.constants.m_e,
    initial_distribution=...,
)
```

## Initial distributions

### UniformDistribution

```python
dist = picmi.UniformDistribution(
    density=1e18,                              # m^-3
    lower_bound=[-1e-3, -1e-3, 0.0],          # spatial bounds (default: full domain)
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th, v_th, v_th],          # thermal spread (m/s)
    directed_velocity=[0, 0, c*0.1],           # drift (m/s)
)
```

For thermal at temperature T (kelvin): `v_th = sqrt(k_B * T / m)`. For T in eV: `v_th = sqrt(T * q_e / m)`.

### AnalyticDistribution

```python
dist = picmi.AnalyticDistribution(
    density_expression='n0 * exp(-(x**2 + y**2) / (2*r0**2))',
    n0=1e18,                                  # passed as kwargs, substituted in expression
    r0=5e-4,
    momentum_expressions=[
        '0',
        '0',
        'uz_drift * (1 - (x**2 + y**2)/r0**2)',
    ],
    uz_drift=1e7,
    lower_bound=[-1e-3, -1e-3, 0],            # ALWAYS SET — restricts initialization region
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th, v_th, v_th],
)
```

**ALWAYS set `lower_bound`/`upper_bound`.** Default is the entire box; for a small high-density region in a large box, this wastes memory and time as the initializer creates particles everywhere before filtering.

### GaussianBunchDistribution (for relativistic beams)

```python
beam = picmi.GaussianBunchDistribution(
    n_physical_particles=1e10,
    rms_bunch_size=[1e-6, 1e-6, 1e-5],         # meters
    rms_velocity=[c*1e-4, c*1e-4, c*1e-3],
    centroid_position=[0, 0, 0],
    centroid_velocity=[0, 0, picmi.constants.c*0.9999],   # γ*v, NOT v
)
```

`centroid_velocity` is `γv` (normalized momentum), not raw `v`. For ultra-relativistic beams the third component is `c * γ` approximately.

### ParticleListDistribution

```python
import numpy as np
N = 1000
dist = picmi.ParticleListDistribution(
    x=np.random.uniform(-1e-3, 1e-3, N),
    y=np.random.uniform(-1e-3, 1e-3, N),
    z=np.zeros(N),
    ux=np.zeros(N), uy=np.zeros(N),
    uz=np.full(N, c * 0.1),
    weight=np.full(N, 1e10 / N),
)
```

For prescribed positions. Use for tests, debugging, importing from external codes.

## Particle layouts

```python
# Deterministic, equal placement per cell:
layout = picmi.GriddedLayout(
    n_macroparticle_per_cell=[2, 2, 4],       # per-dimension (16 total in 3D)
)

# Random placement:
layout = picmi.PseudoRandomLayout(
    n_macroparticles_per_cell=16,
    seed=42,
)
```

Particles per cell rules of thumb:
- Cold/quiescent plasma: 4-16
- Hot plasma with thermal noise: 16-64
- Collisional plasma (collision statistics matter): 32-128
- Turner et al. capacitive discharge benchmarks: 200-500

## Applied external fields

```python
# Applied once at t=0:
init_field = picmi.AnalyticInitialField(
    Ez_expression='1e6 * exp(-((z - 0.005)/1e-3)**2)',
)

# Static, every step:
b_field = picmi.ConstantAppliedField(Bz=0.5)

# Time-varying, every step:
applied = picmi.AnalyticAppliedField(
    Bz_expression='0.5 * (1 + tanh((t - 1e-8) / 1e-9))',
    lower_bound=[-1, -1, 0],
    upper_bound=[1, 1, 0.01],
)

sim.add_applied_field(init_field)
sim.add_applied_field(b_field)
sim.add_applied_field(applied)
```

String expressions use `x`, `y`, `z`, `t` as variables, support `sin`, `cos`, `exp`, `sqrt`, etc. Custom params passed as kwargs.

## Constants

```python
picmi.constants.c        # speed of light
picmi.constants.ep0      # vacuum permittivity
picmi.constants.mu0      # vacuum permeability
picmi.constants.q_e      # elementary charge (absolute value)
picmi.constants.m_e      # electron mass
picmi.constants.m_p      # proton mass
```

All SI. Use these instead of hardcoded values.

## Common pitfalls

- **Ions without `charge_state`** → neutral; bug runs silently. Always set `charge_state` for ions.
- **`AnalyticDistribution` without bounds** → particles created everywhere in the box, wasting memory.
- **Mismatched field/particle BCs on periodic axes** → either both periodic or neither.
- **Forgetting `warpx_*` prefix** for WarpX-specific args → silently ignored as unknown kwargs.
- **Mixing dim across the script** → import `picmi.Cartesian3DGrid` but install `pywarpx[2D]` → runtime error.
- **`centroid_velocity` confusion** → it's γv, not v. For relativistic beams this matters; for non-rel it's ~the same.
