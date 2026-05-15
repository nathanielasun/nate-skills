# PICMI fundamentals: Simulation, Grid, and Solver

This reference covers the skeleton of any WarpX PICMI script — the Simulation object, the Grid class (geometry), the Solver class, embedded boundaries, applied fields, and the lifecycle of a PICMI script. Read this when starting a new script, choosing a geometry, or configuring boundary conditions.

## Table of contents

1. [The Simulation class](#the-simulation-class)
2. [Grid classes (geometry)](#grid-classes-geometry)
3. [Boundary conditions](#boundary-conditions)
4. [Field solvers](#field-solvers)
5. [Embedded boundaries](#embedded-boundaries)
6. [Applied fields](#applied-fields)
7. [The two operating modes in detail](#the-two-operating-modes-in-detail)
8. [Constants](#constants)
9. [A canonical minimal script](#a-canonical-minimal-script)

## The Simulation class

`pywarpx.picmi.Simulation` is the central object. Every other PICMI class is passed *to* it or registered through one of its `add_*` methods.

```python
from pywarpx import picmi

sim = picmi.Simulation(
    solver=my_solver,                          # required
    time_step_size=1e-12,                      # required if CFL not set on solver
    max_steps=10000,                            # one of max_steps, max_time, or use step()
    verbose=1,
    particle_shape='linear',                    # 'NGP', 'linear', 'quadratic', 'cubic'
    gamma_boost=None,                           # for boosted-frame simulations
    # Plus many warpx_* options:
    warpx_current_deposition_algo='direct',     # 'direct', 'esirkepov', 'vay'
    warpx_charge_deposition_algo='standard',
    warpx_field_gathering_algo='energy-conserving',  # or 'momentum-conserving'
    warpx_particle_pusher_algo='boris',         # 'boris', 'vay', 'higuera'
    warpx_use_filter=None,
    warpx_grid_type='staggered',                # 'collocated', 'staggered', 'hybrid'
    warpx_random_seed=42,                       # int or 'random'
    warpx_serialize_initial_conditions=False,
    warpx_amrex_use_gpu_aware_mpi=True,
    warpx_load_balance_intervals='100',         # interval syntax string
    warpx_load_balance_costs_update='heuristic',
    warpx_sort_intervals=None,                  # default -1 on CPU, 4 on GPU
    warpx_break_signals=[],                     # ['SIGINT', 'SIGTERM'] etc.
    warpx_checkpoint_signals=[],
)
```

The full list of `warpx_*` arguments is in the [Simulation class reference](https://warpx.readthedocs.io/en/latest/usage/python.html#pywarpx.picmi.Simulation). Most defaults are sensible; touch them only if you have a specific reason.

### `add_*` methods

After construction:

```python
sim.add_species(species, layout, initialize_self_field=False)
sim.add_laser(laser, injection_method)
sim.add_diagnostic(diagnostic)
sim.add_interaction(interaction)            # for collisions, ionization
sim.add_applied_field(applied_field)
```

These accept the PICMI class instances you've created. The `layout` argument on `add_species` is critical — it specifies how particles are placed (gridded vs pseudo-random) and the count per cell. See `species-and-injection.md`.

For collisions, the WarpX-specific pattern is to either pass them via `warpx_collisions` on the Simulation object or to call `sim.collisions = [...]` directly. Recent versions standardize on `add_interaction`.

### Running

```python
# Interactive
sim.step(nsteps=100)       # advance 100 steps
sim.step()                  # advance to max_steps or max_time

# Or, preprocessor mode
sim.write_input_file('inputs')
```

After `step()` finishes, the simulation state is still accessible via `sim.extension` until the Python process ends.

## Grid classes (geometry)

Each PICMI script targets one geometry; the choice of `Grid` class determines the dimensionality.

### Cartesian3DGrid

```python
grid = picmi.Cartesian3DGrid(
    number_of_cells=[64, 64, 256],
    lower_bound=[-1e-3, -1e-3, 0.0],
    upper_bound=[1e-3, 1e-3, 1e-2],
    lower_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    upper_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
    # Particle BCs (default to field BCs if not set):
    lower_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    upper_boundary_conditions_particles=['periodic', 'periodic', 'absorbing'],
    warpx_max_grid_size=64,                # tile size for parallel decomposition
    warpx_blocking_factor=8,
    # Time-varying potential on Dirichlet boundaries:
    warpx_potential_lo_z='V_RF * sin(2*pi*f_RF*t)',
    warpx_potential_hi_z='0',
)
```

Equivalent classes for other dimensionalities:
- `picmi.Cartesian2DGrid` — `(x, y)` typically, with `z` translationally invariant. Use `nx`, `ny`, `xmin`, `xmax`, etc.
- `picmi.Cartesian1DGrid` — `(x,)`. Mainly for benchmarks and 1D PIC studies (capacitive discharge).
- `picmi.CylindricalGrid` — axisymmetric. Use `nr`, `nz`, `n_azimuthal_modes`, `rmin`, `rmax`, `zmin`, `zmax`.

### Cylindrical-specific notes

```python
grid = picmi.CylindricalGrid(
    nr=128, nz=512,
    n_azimuthal_modes=1,              # 0 = pure axisymmetric, >0 includes asymmetric modes
    rmin=0.0, rmax=1e-2,
    zmin=0.0, zmax=2e-2,
    bc_rmax='dirichlet',              # rmin is at axis; no BC needed
    bc_zmin='dirichlet', bc_zmax='dirichlet',
)
```

`n_azimuthal_modes=1` means mode 0 only (single azimuthal mode, pure axisymmetric). `n_azimuthal_modes=2` includes mode 1 (first asymmetric harmonic), etc. Higher modes capture asymmetry but multiply cost roughly linearly.

For low-temp plasma in RZ geometry, mode 0 is usually sufficient unless studying breakdown patterns or asymmetric driving.

### Mesh refinement via `refined_regions`

```python
grid = picmi.Cartesian2DGrid(
    nx=128, ny=256,
    xmin=0, xmax=1e-2, ymin=0, ymax=1e-2,
    refined_regions=[
        [1, [0, 0.4e-2], [0.5e-2, 0.6e-2], [2, 2]],  # level 1, region, refinement factor
    ],
    ...
)
```

The list format is `[level, lower_corner, upper_corner, refinement_factor]`. Default factor is `[2, 2, ...]`. Only one level beyond level 0 (i.e., level 1) is supported in some configurations (notably the hybrid PIC solver).

MR adds significant complexity. For low-temp plasma at near-uniform density, AMR is rarely needed and usually not worth the complication. For laser-plasma or beam-plasma with sharply localized regions, AMR is essential.

### `warpx_max_grid_size` and `warpx_blocking_factor`

These control the AMReX domain decomposition. The simulation domain is chopped into tiles of maximum size `max_grid_size`, with blocking factor `blocking_factor` enforcing that tile sizes are multiples of this number.

- For CPU runs: `max_grid_size=32` to `64` is typical. Smaller tiles enable better load balancing but increase MPI communication overhead.
- For GPU runs: `max_grid_size=128` or larger to expose enough parallel work per tile.
- `blocking_factor=8` is typical; powers of 2 work well.

If you see "particle deposition shape doesn't fit in tile guards" errors, your `max_grid_size` is too small for the particle shape. Reduce particle shape order or increase grid size.

### `warpx_numprocs` — explicit domain decomposition

```python
grid = picmi.Cartesian3DGrid(...,
    warpx_numprocs=[2, 2, 4],   # exactly 2*2*4 = 16 MPI ranks, decomposed this way
)
```

When set, this forces exact domain decomposition. Useful for reproducibility or to enforce a specific layout (e.g., long axis along z gets more ranks). When unset, AMReX chooses automatically.

## Boundary conditions

PICMI distinguishes between **field BCs** and **particle BCs**. Each can be different.

### Field BCs

Set via `lower_boundary_conditions` / `upper_boundary_conditions` vectors or `bc_xmin`, `bc_xmax`, etc. individually. Options:

- **`'periodic'`** — wraps around. Both lo and hi must be periodic on the same axis.
- **`'dirichlet'`** — fixed value (default 0, can be set via `warpx_potential_*` parameters for ES solvers).
- **`'neumann'`** — zero-gradient. ES solvers only.
- **`'open'`** — outgoing waves (for EM solvers, similar to PML).
- **`'absorbing_silver_mueller'`** — local first-order absorbing BC; works best at normal incidence.
- **`'pec'`** — perfect electric conductor (tangential E = 0).
- **`'pmc'`** — perfect magnetic conductor (tangential B = 0).

For typical low-temp plasma in 1D Cartesian: `bc_xmin='dirichlet', bc_xmax='dirichlet'` with time-varying `warpx_potential_lo_x` / `warpx_potential_hi_x` for the applied RF voltage.

### Particle BCs

Set via `lower_boundary_conditions_particles` / `upper_boundary_conditions_particles`. Options:

- **`'periodic'`** — wraps around (matches field BC).
- **`'absorbing'`** — particles leaving the boundary are deleted.
- **`'reflect'`** — particles bounce off elastically.
- **`'thermal'`** — particles re-injected from a thermal distribution at temperature `warpx_boundary_u_th[species]`.
- **`'open'`** — same as absorbing for particles.

Common combination for capacitive discharge: field BC `'dirichlet'`, particle BC `'absorbing'`. Particles hitting the electrode are removed.

For periodic axes, both field and particle BCs must be periodic together.

### `warpx_save_particles_at_*` on Species

To record particles as they cross a boundary (for downstream analysis of fluxes, energy spectra, etc.):

```python
electrons = picmi.Species(...,
    warpx_save_particles_at_xlo=True,
    warpx_save_particles_at_xhi=True,
)
```

This produces "boundary scraper" output in the diagnostic directory. See `diagnostics-and-output.md` for accessing this data.

## Field solvers

### Electromagnetic solver

```python
em_solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',           # 'Yee', 'CKC', 'Lehe', 'PSATD', 'PSTD', 'GPSTD', 'DS', 'ECT'
    cfl=0.99,
    divE_cleaning=False,
    divB_cleaning=False,
)
```

Method selection:
- **`Yee`** — standard FDTD on staggered grid. Default. Works in all geometries.
- **`CKC`** — Cole-Karkkainen-Cowan extended stencil. Better dispersion for relativistic plasmas. FDTD-like.
- **`Lehe`** — modified dispersion CKC variant. For laser-plasma in boosted frame.
- **`PSATD`** — pseudo-spectral. Requires `WarpX_FFT=ON` at build. Eliminates Numerical Cherenkov Instability, supports higher-order time integration.
- **`ECT`** — Enlarged Cell Technique. For internal conductors (embedded boundaries).

For low-temp plasma, EM solvers are rarely needed — electrostatic is usually appropriate. Use EM when modeling wave propagation, microwave coupling, or relativistic effects.

### Electrostatic solver

```python
es_solver = picmi.ElectrostaticSolver(
    grid=grid,
    method='Multigrid',         # 'Multigrid' (MLMG) or 'FFT' (IGF)
    required_precision=1e-7,
    maximum_iterations=200,
    warpx_relativistic=False,
    warpx_self_fields_verbosity=2,
)
```

- **`Multigrid`** — geometric multigrid via AMReX's MLMG. Handles arbitrary BCs and embedded boundaries. **Workhorse for low-temp plasma.**
- **`FFT`** — periodic-only (or with specific symmetric BCs). Fast, but limited applicability.

For an ES discharge simulation, this is your default solver. Pair with `warpx_potential_lo_x` etc. on the grid for time-varying applied voltage.

### Hybrid (Ohm's law) solver

```python
hybrid_solver = picmi.HybridPICSolver(
    grid=grid,
    Te=10.0,                              # eV, electron temperature (massless fluid)
    n0_reference_density=1e18,            # reference density for normalization
    plasma_resistivity='0.0',             # string expression of (rho, J)
    substeps=10,                          # field solve substeps per particle step
    ...
)
```

For ion-kinetic problems where electrons are treated as a massless fluid. See the WarpX `Examples/Tests/ohm_solver_*` examples.

## Embedded boundaries

WarpX supports embedded boundaries (EB) for representing arbitrary solid surfaces inside the simulation domain. Two forms:

### Implicit function (analytic)

```python
eb = picmi.EmbeddedBoundary(
    implicit_function='-(((x - 0)**2 + (y - 0)**2) - r_inner**2)',
    r_inner=2e-3,                  # parsed as user parameter
    potential='V_bias',            # time-varying potential on the surface
    V_bias='100.0',                # in volts
)
```

Positive values of `implicit_function` are *inside* the simulation domain (fluid region); negative are inside the solid. The example above carves out a cylindrical conductor of radius `r_inner` at the origin (in 2D Cartesian; would be axial in 3D).

Pass as `warpx_embedded_boundary=eb` on the Simulation.

### STL file

```python
eb = picmi.EmbeddedBoundary(
    stl_file='/path/to/geometry.stl',
    stl_scale=1.0,
    stl_center=[0, 0, 0],
    stl_reverse_normal=False,
    potential='100.0',
)
```

For complex geometries (CAD-derived). The STL is read at startup and converted to an implicit function internally.

### Restrictions

- Not currently supported in RZ/RCYLINDER/RSPHERE geometries (check current docs for status).
- ES Multigrid solver and EM Yee/CKC support EB. PSATD does not.
- The ECT EM solver is specifically designed for EB with internal conductors.
- Particles hitting EB are absorbed by default. Use `warpx_save_particles_at_eb=True` on a species to record them.

## Applied fields

For external fields not produced by the simulated charges (e.g., a static B field, a coil-produced field, a magnetic mirror lattice):

### `picmi.AnalyticInitialField`
Applied **once at t=0**:

```python
init_field = picmi.AnalyticInitialField(
    Ez_expression='1e6 * exp(-((z - 0.005)/1e-3)**2)',
)
sim.add_applied_field(init_field)
```

### `picmi.ConstantAppliedField`
Static, applied at every step:

```python
b_field = picmi.ConstantAppliedField(Bz=0.5)  # 0.5 Tesla uniform Bz
sim.add_applied_field(b_field)
```

### `picmi.AnalyticAppliedField`
Time-varying, applied at every step:

```python
applied = picmi.AnalyticAppliedField(
    Bz_expression='0.5 * (1 + tanh((t - 1e-8) / 1e-9))',
    lower_bound=[-1, -1, 0],         # region restriction
    upper_bound=[1, 1, 0.01],
)
sim.add_applied_field(applied)
```

### `picmi.PlasmaLens`
Periodic lattice of focusing elements (for accelerator-relevant simulations):

```python
lens = picmi.PlasmaLens(
    period=0.1,
    starts=[0.02],
    lengths=[0.005],
    strengths_E=[1e9],
    strengths_B=[10],
)
```

### `picmi.Mirror`
Zero-out fields in a slab (e.g., behind a laser to simulate a perfect reflector):

```python
mirror = picmi.Mirror(
    z_front_location=0.005,
    depth=1e-6,
    number_of_cells=4,
)
```

## The two operating modes in detail

### Interactive mode

```python
# script.py
from pywarpx import picmi
sim = picmi.Simulation(...)
# ... configure ...
sim.step()
```

Run with `python3 script.py` or `mpiexec -n N python3 script.py`. The simulation runs *inside* the Python process. Native callbacks and `Simulation.extension` are available. The process holds the entire WarpX state in memory.

After `step()` returns, you can:
- Access `sim.extension.warpx` for the C++ WarpX object
- Compute derived quantities from the final state
- Run additional analysis or save data

### Preprocessor mode

```python
# script.py
from pywarpx import picmi
sim = picmi.Simulation(...)
# ... configure ...
sim.write_input_file('inputs')
```

Run with `python3 script.py` (no MPI needed for this — it's pure Python). Produces a text file `inputs` in the AMReX input-file format. The simulation then runs as:

```bash
mpiexec -n N warpx.3d.MPI.OMP.DP.OPMD.QED inputs
```

(or whatever your built binary is called). No Python at runtime. Callbacks not available.

The generated input file is a useful debugging aid even when running interactively: read it to verify that your PICMI parameters mapped to C++ parameters as you intended. Many "why is my simulation doing X?" bugs become obvious when you see the input file.

## Constants

PICMI provides SI-unit constants for convenience:

```python
from pywarpx import picmi

picmi.constants.c         # 2.998e8, speed of light
picmi.constants.ep0       # 8.854e-12, vacuum permittivity
picmi.constants.mu0       # 1.257e-6, vacuum permeability
picmi.constants.q_e       # 1.602e-19, elementary charge (absolute value)
picmi.constants.m_e       # 9.109e-31, electron mass
picmi.constants.m_p       # 1.673e-27, proton mass
```

Use these instead of hardcoding values to keep scripts readable and self-documenting. For specialized constants (Boltzmann, atomic masses for other species), import from `scipy.constants` or define explicitly.

## A canonical minimal script

A complete-but-simple PICMI script. This is the template to start from:

```python
#!/usr/bin/env python3
"""
Minimal 2D PICMI example.
Free-streaming electrons in a periodic box, no fields, for sanity-checking
the install and the diagnostic output.
"""
from pywarpx import picmi

# 1. Constants and parameters (use SI throughout)
c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# 2. Grid
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
solver = picmi.ElectromagneticSolver(grid=grid, method='Yee', cfl=0.99)

# 4. Species
distribution = picmi.UniformDistribution(
    density=1e18,
    rms_velocity=[c*0.01, c*0.01, c*0.01],
)
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=distribution,
)
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2], grid=grid)

# 5. Diagnostics
field_diag = picmi.FieldDiagnostic(
    name='diag1',
    grid=grid,
    period=50,
    data_list=['rho', 'E', 'B', 'J'],
    write_dir='./diags',
    warpx_format='openpmd',
)
part_diag = picmi.ParticleDiagnostic(
    name='diag1',
    period=50,
    species=[electrons],
    data_list=['position', 'momentum', 'weighting'],
    write_dir='./diags',
    warpx_format='openpmd',
)

# 6. Simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=200,
    verbose=1,
    warpx_random_seed=42,
)
sim.add_species(electrons, layout=layout)
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)

# 7. Run
sim.step()
```

This template has all the structural elements but minimal physics. Extend it by:
- Replacing `ElectromagneticSolver` with `ElectrostaticSolver` for electrostatic problems.
- Replacing `UniformDistribution` with `AnalyticDistribution`, `GaussianBunchDistribution`, or `ParticleListDistribution` (see `species-and-injection.md`).
- Adding collision interactions via `sim.collisions = [...]` (see `collisions-and-reactions.md`).
- Adding callbacks via `pywarpx.callbacks.installafterstep` (see `runtime-extension.md`).
- Changing the diagnostic `period`, `data_list`, and `warpx_format` (see `diagnostics-and-output.md`).
