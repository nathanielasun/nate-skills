# Species, distributions, and injection

This reference covers particle species setup: declaring species, initial distributions, layouts, beam injection, runtime particle injection, and WarpX-specific species behaviors (reflection models, save-at-boundary, resampling). Read this when configuring particle species in a PICMI script.

## Table of contents

1. [The Species class](#the-species-class)
2. [Particle types and the openPMD species extension](#particle-types-and-the-openpmd-species-extension)
3. [Initial distributions](#initial-distributions)
4. [Layouts: how particles fill a cell](#layouts-how-particles-fill-a-cell)
5. [WarpX-specific species options](#warpx-specific-species-options)
6. [Boundary handling: reflection, absorption, thermal](#boundary-handling-reflection-absorption-thermal)
7. [Save-at-boundary diagnostics](#save-at-boundary-diagnostics)
8. [Particle resampling](#particle-resampling)
9. [Beam injection (continuous injection)](#beam-injection-continuous-injection)
10. [Adding particles at runtime](#adding-particles-at-runtime)
11. [Putting it together: example setups](#putting-it-together-example-setups)

## The Species class

```python
from pywarpx import picmi

electrons = picmi.Species(
    particle_type='electron',              # use OpenPMD species type strings
    name='electrons',                       # name used in inputs/diagnostics
    initial_distribution=my_distribution,
    method='Boris',                         # 'Boris', 'Vay', 'Higuera-Cary', 'free-streaming', 'LLRK4'
    particle_shape='linear',                # overrides Simulation default if set
    density_scale=1.0,                      # scale factor on density
    # Optional: explicit charge/mass override (rarely needed if particle_type is set)
    charge=None,
    mass=None,
    charge_state=None,                      # for atoms/ions
    # WarpX-specific:
    warpx_do_not_push=False,
    warpx_do_not_deposit=False,
    warpx_do_not_gather=False,
    warpx_do_resampling=False,
    warpx_self_fields_required_precision=1e-11,
    warpx_save_previous_position=False,
)
```

Register with the simulation:

```python
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2, 4], grid=grid)
sim.add_species(electrons, layout=layout, initialize_self_field=False)
```

The `layout` determines how many macroparticles per cell and how they're placed (gridded vs random). The `initial_distribution` determines the *physical* density and velocity statistics. Together they define the initial state.

### `name` matters

The `name` is used throughout WarpX — in diagnostics output, in collision-process configuration (`species=[electrons, ions]` references species by name when serialized), in field-name conventions (`rho_electrons` is the electron charge density), and in callbacks (`multi_pc.get_particle_container_from_name('electrons')`). Choose names that are short, descriptive, and won't collide.

Convention: lowercase, singular or plural form depending on the species. `electrons`, `ions`, `argon_ions`, `beam_electrons`, `photons`, `He_ions`, `Kr_ions`.

### `particle_type` vs explicit `charge`/`mass`

The clean way:

```python
electrons = picmi.Species(particle_type='electron', name='electrons', ...)
```

The strings come from the [openPMD species type extension](https://github.com/openPMD/openPMD-standard/blob/latest/EXT_SpeciesType.md). For ions, the type encodes the element and charge state:

```python
argon_ions = picmi.Species(
    particle_type='Ar',
    charge_state=1,                  # singly-ionized
    name='argon_ions',
    initial_distribution=...,
)
```

For unusual species or to override the standard charge/mass:

```python
heavy_e = picmi.Species(
    name='heavy_electron',
    charge=-picmi.constants.q_e,
    mass=10 * picmi.constants.m_e,    # 10x electron mass
    initial_distribution=...,
)
```

This pattern is useful for testing: heavier "electrons" lower the plasma frequency and let you take larger time steps for sanity checks.

## Particle types and the openPMD species extension

Common values for `particle_type`:

| Type string | What it is |
|---|---|
| `'electron'` | electron, charge = -q_e, mass = m_e |
| `'positron'` | positron |
| `'proton'` | proton |
| `'photon'` | photon (massless, charge 0) |
| `'H'`, `'He'`, `'Ar'`, `'Kr'`, `'Xe'`, `'N'`, `'O'`, ... | atoms (need `charge_state=N` for ions; charge_state=0 for neutrals) |
| `'D'` | deuterium |
| `'T'` | tritium |
| `'alpha'` | alpha particle |

For atoms, the underlying mass is the atomic mass for the most common isotope. If you need a specific isotope, set `mass` explicitly.

If `particle_type` is not in the list above, you can omit it and supply `charge` and `mass` directly. The cost is loss of metadata in the openPMD output.

### Common pitfall: ions without `charge_state`

```python
# WRONG: no charge state, defaults to neutral
bad_ions = picmi.Species(particle_type='Ar', name='ions', ...)

# RIGHT: explicitly Ar+1
good_ions = picmi.Species(particle_type='Ar', charge_state=1, name='ions', ...)
```

Without `charge_state`, the species is neutral (charge = 0), which is rarely what you want for a "plasma ion" species. The simulation will run without error but the ions won't respond to electric fields.

## Initial distributions

The distribution describes the *physical* particle state at t=0 — density profile and velocity statistics. Several options:

### `UniformDistribution`

```python
dist = picmi.UniformDistribution(
    density=1e18,                              # m^-3
    lower_bound=[-1e-3, -1e-3, 0.0],          # optional spatial limits
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th_e, v_th_e, v_th_e],    # thermal spread, m/s
    directed_velocity=[0, 0, c*0.1],          # drift velocity, m/s
    fill_in=False,                             # for moving window: refill behind window
)
```

The simplest case: uniform density inside the bounds, Maxwellian-like velocity distribution. `rms_velocity` is the thermal spread (1-sigma per component); `directed_velocity` is the bulk flow.

For a thermal species at temperature T:
```python
v_th = sqrt(k_B * T / m)  # in m/s, with T in K and m in kg
# or for T in eV:
v_th = sqrt(T * q_e / m)
```

### `AnalyticDistribution`

For spatially non-uniform density or velocity:

```python
dist = picmi.AnalyticDistribution(
    density_expression='n0 * exp(-(x**2 + y**2) / (2*r0**2))',
    n0=1e18,                                  # additional params used in expression
    r0=5e-4,
    momentum_expressions=[
        '0',
        '0',
        'uz_drift * (1 - (x**2 + y**2)/r0**2)',
    ],
    uz_drift=1e7,
    lower_bound=[-1e-3, -1e-3, 0],
    upper_bound=[1e-3, 1e-3, 1e-2],
    rms_velocity=[v_th, v_th, v_th],
    fill_in=False,
)
```

The string expressions are evaluated by WarpX's parser at runtime. Variables `x`, `y`, `z` are positions in meters. Additional parameters are passed as keyword arguments and substituted in.

**Always set `lower_bound` / `upper_bound`.** The default is the entire simulation box, and the particle initializer creates particles everywhere inside those bounds before filtering by the density expression. For a small high-density region in a large box, this wastes huge amounts of memory and time. Restrict the bounds to where density is non-negligible.

There's also `warpx_momentum_spread_expressions` for spatially varying thermal spread (a feature beyond the PICMI standard).

### `GaussianBunchDistribution`

For relativistic beams (laser-plasma, beam-plasma):

```python
beam = picmi.GaussianBunchDistribution(
    n_physical_particles=1e10,                # total real particle count
    rms_bunch_size=[1e-6, 1e-6, 1e-5],        # in meters
    rms_velocity=[c*1e-4, c*1e-4, c*1e-3],
    centroid_position=[0, 0, 0],
    centroid_velocity=[0, 0, c*0.99],          # γ*v
    velocity_divergence=[0, 0, 0],
)
```

`centroid_velocity` is the *normalized momentum* γv (in m/s with γ implied), not just v. For ultra-relativistic beams (γ=10), the third component is `c * sqrt(1 - 1/γ²) * γ ≈ c * γ` for forward motion.

### `ParticleListDistribution`

Explicit list of particles:

```python
import numpy as np

# 1000 test particles at specific positions
N = 1000
x_vals = np.random.uniform(-1e-3, 1e-3, N)
y_vals = np.random.uniform(-1e-3, 1e-3, N)
z_vals = np.zeros(N)
uz_vals = np.full(N, c * 0.1)
weights = np.full(N, 1e10 / N)

dist = picmi.ParticleListDistribution(
    x=x_vals,
    y=y_vals,
    z=z_vals,
    ux=np.zeros(N),
    uy=np.zeros(N),
    uz=uz_vals,
    weight=weights,
)
```

Useful for: prescribed particle distributions from a previous simulation, importing from external codes, debugging with specific initial conditions, single-particle dynamics tests.

When you need *thousands or more* prescribed particles, prefer reading from an openPMD file via the `*_external_particle_field` mechanism rather than embedding the data in the Python script.

## Layouts: how particles fill a cell

The layout determines the *numerical* count and placement of macroparticles, decoupled from the physical distribution.

### `GriddedLayout`

Deterministic placement on a per-cell sub-grid:

```python
layout = picmi.GriddedLayout(
    n_macroparticle_per_cell=[2, 2, 4],   # 16 particles per cell in 3D
    grid=grid,                             # optional, defaults to simulation grid
)
```

The numbers are per-dimension, so `[2, 2, 4]` means 2 in x, 2 in y, 4 in z. Total = 16 macroparticles per cell. Particles are placed at fixed sub-cell positions, the same in every cell.

Use for: deterministic initial conditions, reproducibility, low-noise simulations.

### `PseudoRandomLayout`

Random placement:

```python
layout = picmi.PseudoRandomLayout(
    n_macroparticles_per_cell=16,         # total per cell
    seed=42,                               # for reproducibility
    grid=grid,
)
```

Or specify total count instead of per-cell:

```python
layout = picmi.PseudoRandomLayout(
    n_macroparticles=10000,                # total across simulation
    seed=42,
)
```

Use for: distributions where regular placement would introduce structure (Gaussian bunches), or for simulations where reproducibility isn't required.

### Choosing the count

The "right" particles per cell depends on physics:
- **Cold/quiescent plasma**: 4-16 per cell is often enough.
- **Hot plasma with thermal noise**: 16-64+ per cell to suppress shot noise (which scales as 1/sqrt(N)).
- **Collisional plasma where collision statistics matter**: 32-128+, especially for rare events (ionization, charge exchange).
- **Beam in plasma background**: beam can be lower (4-8) since it's already localized; plasma background needs the usual count.
- **Capacitive discharge benchmarks (Turner et al. 2013)**: 200-500 per cell for the reference cases.

If you're not sure, start low (8 per cell) and increase until results stop changing.

## WarpX-specific species options

### `warpx_do_not_push` / `warpx_do_not_gather` / `warpx_do_not_deposit`

```python
neutral_background = picmi.Species(
    particle_type='Ar',
    charge_state=0,
    name='neutrals',
    initial_distribution=...,
    warpx_do_not_push=True,         # frozen background
    warpx_do_not_deposit=True,      # don't contribute to charge/current (they're neutral)
    warpx_do_not_gather=True,       # don't gather fields (no Lorentz force)
)
```

Useful for representing a frozen neutral background as simulated particles (for DSMC). The particles exist for collision pairing purposes but don't otherwise participate in the PIC loop.

For tracer particles (test particles that follow fields but don't perturb them):

```python
tracers = picmi.Species(
    particle_type='electron',
    name='tracers',
    initial_distribution=...,
    warpx_do_not_deposit=True,      # don't perturb fields
    # warpx_do_not_gather=False (default) — they DO gather fields
    # warpx_do_not_push=False (default) — they DO move
)
```

### `warpx_save_previous_position`

```python
electrons = picmi.Species(...,
    warpx_save_previous_position=True,
)
```

Adds `prev_x`, `prev_y`, `prev_z` as runtime attributes per particle. Required for some diagnostics and for trajectory reconstruction. Small memory cost.

### `warpx_self_fields_*`

For per-species electrostatic self-field calculation (used when `add_species(..., initialize_self_field=True)`):

```python
beam = picmi.Species(...,
    warpx_self_fields_required_precision=1e-11,
    warpx_self_fields_absolute_tolerance=0.0,
    warpx_self_fields_max_iters=200,
    warpx_self_fields_verbosity=2,
)
```

These control the multigrid solver convergence for the species' own field. Defaults are sensible; touch only if convergence issues arise.

## Boundary handling: reflection, absorption, thermal

The particle BCs on the grid (`lower_boundary_conditions_particles`) determine the default behavior. WarpX adds finer-grained control per species via reflection models.

### `warpx_reflection_model_*`

A string expression returning the *probability* of reflection (0 to 1), as a function of particle speed `v` in m/s:

```python
electrons = picmi.Species(...,
    warpx_reflection_model_zlo='1.0',                      # always reflect at z_lo
    warpx_reflection_model_zhi='0.5',                       # 50% reflection at z_hi
    warpx_reflection_model_xlo='if(v < 1e5, 1.0, 0.0)',    # reflect slow particles, absorb fast
)
```

Particles not reflected are subject to the particle BC. So with `bc_zlo_particles='absorbing'` and `reflection_model_zlo='1.0'`, all electrons hitting `z_lo` reflect; with `'0.5'`, half reflect and half are absorbed.

### Thermal re-injection

For particle BCs set to `'thermal'`:

```python
grid = picmi.Cartesian2DGrid(...,
    lower_boundary_conditions_particles=['thermal', 'periodic'],
    upper_boundary_conditions_particles=['thermal', 'periodic'],
    warpx_boundary_u_th={'electrons': 1e6, 'ions': 1e3},   # u_th in WarpX units (= sqrt(T*q_e/m)/c)
)
```

Wait — actually `u_th` here is the *normalized* thermal speed (β = v/c). Check the current docs; this convention has shifted. The formula in the current docs: `u_th = sqrt(T*q_e/mass)/clight` with T in eV.

## Save-at-boundary diagnostics

To record particles as they cross specific boundaries (for energy/flux/angle diagnostics):

```python
electrons = picmi.Species(...,
    warpx_save_particles_at_xlo=True,
    warpx_save_particles_at_xhi=True,
    warpx_save_particles_at_ylo=False,
    warpx_save_particles_at_yhi=False,
    warpx_save_particles_at_zlo=True,
    warpx_save_particles_at_zhi=True,
    warpx_save_particles_at_eb=True,        # for embedded boundary
)
```

This produces "boundary scraper" output in the diagnostic directory. Each crossing event is recorded with position, momentum, weight, and timestamp. Particularly useful for:
- Energy spectra of particles hitting an electrode
- Angular distributions
- Total flux as a function of time
- Particle accounting (total deposited charge per surface)

Access from openPMD-viewer:

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('./diags/boundary_scrapers_lo/')
# treat like normal particle data
```

## Particle resampling

For simulations where particle count grows (e.g., from ionization in a discharge), unbounded growth eventually exhausts memory. WarpX has a Russian-roulette-style resampling mechanism:

```python
electrons = picmi.Species(...,
    warpx_do_resampling=True,
    warpx_resampling_trigger_intervals='1000',                # every 1000 steps
    warpx_resampling_trigger_max_avg_ppc=200,                  # or when avg > 200 per cell
    warpx_resampling_algorithm='leveling_thinning',            # or 'velocity_coincidence_thinning'
    warpx_resampling_min_ppc=10,
)
```

When triggered: in cells with too many particles, some are randomly merged (their weights combined) so that the cell averages back toward the target. The algorithm preserves charge, momentum, and energy in expectation but introduces noise.

For low-temp plasma with steady ionization, resampling is essential — without it, a discharge will exhaust memory in a few thousand steps. Tune `trigger_max_avg_ppc` based on memory budget; 200-500 is typical.

The downside: resampling adds noise. Don't resample tracer species or species where you need to track individual particle histories.

## Beam injection (continuous injection)

To inject particles at a specific surface continuously (e.g., a particle source like an emission cathode):

The PICMI standard `PlasmaInjector` isn't in WarpX's PICMI implementation as a separate class — instead, continuous injection is handled by setting up a species with a small-volume initial distribution and a non-zero `lower_bound`/`upper_bound` that overlaps with `xmin_particles` boundary, combined with `fill_in` semantics.

For more explicit control, use `warpx_*` injection options on the species or add particles via callbacks at runtime (see `runtime-extension.md`).

The clean way to add particles at runtime:

```python
from pywarpx import particle_containers, callbacks

@callbacks.callfromafterstep
def inject_more_particles():
    electron_wrapper = particle_containers.ParticleContainerWrapper('electrons')
    # Inject 1000 electrons at a specific location each step
    electron_wrapper.add_particles(
        x=np.random.uniform(-1e-4, 1e-4, 1000),
        y=np.random.uniform(-1e-4, 1e-4, 1000),
        z=np.full(1000, 1e-5),
        ux=np.zeros(1000),
        uy=np.zeros(1000),
        uz=np.full(1000, 1e6),  # drift in +z
        w=np.full(1000, 1e8),    # weight
        unique_particles=False,   # equally distributed across MPI ranks
    )
```

See `runtime-extension.md` for the full callback machinery.

## Adding particles at runtime

The `ParticleContainerWrapper` provides `add_particles` for in-script injection:

```python
from pywarpx import particle_containers

pc_wrapper = particle_containers.ParticleContainerWrapper('species_name')
pc_wrapper.add_particles(
    x=array_x, y=array_y, z=array_z,
    ux=array_ux, uy=array_uy, uz=array_uz,
    w=array_w,
    unique_particles=True,           # True: each rank adds its own; False: distribute among ranks
    # **kwargs for any runtime attributes
)
```

`unique_particles` semantics:
- `True` (default): each MPI rank calls `add_particles` and adds its own particles. Total count = N * number_of_ranks.
- `False`: the call adds the *total* count distributed across ranks. Use when you want a fixed total regardless of rank count.

For low-temp plasma applications:
- **Emission from a surface**: `unique_particles=False`, sample from a surface distribution each step.
- **Bulk injection** (e.g., uniform background refilling): use `picmi.UniformDistribution` with `fill_in=True` rather than runtime injection.
- **External plasma source**: read positions/velocities from a file in a callback, inject them.

## Putting it together: example setups

### Capacitive discharge: electrons + ions in argon

```python
# Two species: electrons and Ar+ ions
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=picmi.UniformDistribution(
        density=1e15,
        lower_bound=[None, 0.0],
        upper_bound=[None, 6.7e-2],
        rms_velocity=[v_th_e, v_th_e, v_th_e],
    ),
    warpx_save_particles_at_xlo=True,
    warpx_save_particles_at_xhi=True,
)

ions = picmi.Species(
    particle_type='Ar',
    charge_state=1,
    name='ar_ions',
    initial_distribution=picmi.UniformDistribution(
        density=1e15,
        lower_bound=[None, 0.0],
        upper_bound=[None, 6.7e-2],
        rms_velocity=[v_th_ion, v_th_ion, v_th_ion],
    ),
    warpx_save_particles_at_xlo=True,
    warpx_save_particles_at_xhi=True,
    warpx_do_resampling=True,
    warpx_resampling_trigger_intervals='1000',
)

layout = picmi.GriddedLayout(n_macroparticle_per_cell=[256, 1, 1])
sim.add_species(electrons, layout=layout)
sim.add_species(ions, layout=layout)
```

### Beam-plasma interaction

```python
# Drive beam
beam_dist = picmi.GaussianBunchDistribution(
    n_physical_particles=1e10,
    rms_bunch_size=[5e-6, 5e-6, 30e-6],
    rms_velocity=[1e4, 1e4, 1e6],
    centroid_position=[0, 0, -100e-6],
    centroid_velocity=[0, 0, picmi.constants.c*0.9999],
)
beam = picmi.Species(particle_type='electron', name='beam', initial_distribution=beam_dist)
sim.add_species(beam, layout=picmi.PseudoRandomLayout(n_macroparticles=100000))

# Plasma background
plasma_dist = picmi.AnalyticDistribution(
    density_expression='n0 * (z > z_start)',
    n0=1e23,
    z_start=0,
    lower_bound=[-50e-6, -50e-6, 0],
    upper_bound=[50e-6, 50e-6, 1e-3],
)
plasma_e = picmi.Species(particle_type='electron', name='plasma_e',
                         initial_distribution=plasma_dist)
plasma_i = picmi.Species(particle_type='H', charge_state=1, name='plasma_i',
                         initial_distribution=plasma_dist)
layout = picmi.GriddedLayout(n_macroparticle_per_cell=[1, 1, 1])
sim.add_species(plasma_e, layout=layout)
sim.add_species(plasma_i, layout=layout)
```

### Frozen neutral background + simulated electrons (typical low-temp plasma)

```python
# Electrons participate fully in PIC
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=picmi.UniformDistribution(density=1e15, ...),
)
sim.add_species(electrons, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[32, 32, 4]))

# Neutrals are present for DSMC but otherwise frozen
neutrals = picmi.Species(
    particle_type='Ar',
    charge_state=0,
    name='neutrals',
    initial_distribution=picmi.UniformDistribution(density=2.5e25, ...),  # 1 atm at 300 K
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
    warpx_do_not_deposit=True,
)
sim.add_species(neutrals, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[8, 8, 2]))
```

The neutrals exist so that DSMC binary collisions can pair them with electrons or ions. They don't push, don't feel fields, and don't deposit charge — they're stationary collision partners.

For pure MCC against a *non-simulated* thermal background, you don't need a neutral species at all — see `collisions-and-reactions.md`.
