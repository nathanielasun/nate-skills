# Diagnostics and output

This reference covers WarpX's diagnostic and output system from the PICMI side: full-field and particle diagnostics, reduced diagnostics, checkpoints, output formats, particle filtering, and lab-frame (back-transformed) diagnostics. Read this when configuring what gets written to disk and when.

## Table of contents

1. [Overview: diagnostic types](#overview-diagnostic-types)
2. [FieldDiagnostic](#fielddiagnostic)
3. [ParticleDiagnostic](#particlediagnostic)
4. [Output formats](#output-formats)
5. [Reduced diagnostics](#reduced-diagnostics)
6. [Checkpoints and restarts](#checkpoints-and-restarts)
7. [Particle filtering and subsampling](#particle-filtering-and-subsampling)
8. [Lab-frame diagnostics for boosted-frame runs](#lab-frame-diagnostics-for-boosted-frame-runs)
9. [Particle-field diagnostics (derived quantities)](#particle-field-diagnostics-derived-quantities)
10. [Performance considerations](#performance-considerations)

## Overview: diagnostic types

WarpX produces several distinct diagnostic streams, each with its own configuration class in PICMI:

| Class | Purpose | Output |
|---|---|---|
| `picmi.FieldDiagnostic` | Snapshots of grid fields (E, B, J, ρ, φ) | One file per snapshot (or one openPMD series) |
| `picmi.ParticleDiagnostic` | Snapshots of particle data | One file per snapshot per species |
| `picmi.ElectrostaticFieldDiagnostic` | Alias for `FieldDiagnostic` (deprecated, kept for compatibility) | same as `FieldDiagnostic` |
| `picmi.Checkpoint` | Full simulation state for restart | Plotfile or openPMD with full state |
| `picmi.LabFrameFieldDiagnostic` | Lab-frame back-transformed fields (boosted-frame runs) | openPMD series in lab frame |
| Reduced diagnostics (no PICMI class) | Scalar/vector quantities per step (energy, particle counts, etc.) | ASCII text file |
| Boundary scrapers (via Species options) | Particles crossing boundaries | openPMD-style per-event records |

Multiple diagnostics with different periods and contents can coexist — typical setup: full snapshots every 1000 steps + reduced diagnostics every 1 step.

## FieldDiagnostic

```python
field_diag = picmi.FieldDiagnostic(
    name='diag1',                      # used in output directory name
    grid=grid,                          # the grid to diagnose
    period=1000,                        # write every N steps
    data_list=['rho', 'E', 'B', 'J', 'phi'],
    write_dir='./diags',
    step_min=0,                         # skip diagnostics before this step
    step_max=None,                      # default: unbounded
    parallelio=True,
    warpx_format='openpmd',             # 'plotfile', 'openpmd', 'checkpoint', 'ascent', 'sensei'
    warpx_openpmd_backend='h5',         # 'bp', 'h5', 'json' (for openpmd format)
    warpx_file_prefix='diag1',
    warpx_file_min_digits=6,
)
sim.add_diagnostic(field_diag)
```

### `data_list` options

The data list specifies which fields to write. Standard entries:
- `'E'` — all three components of E (Ex, Ey, Ez written separately)
- `'B'` — all three components of B
- `'J'` — all three components of J
- `'rho'` — total charge density
- `'rho_<species>'` — per-species charge density (e.g., `'rho_electrons'`)
- `'phi'` — scalar potential (electrostatic mode only)

Individual components: `'Ex'`, `'Ey'`, `'Ez'`, `'Bx'`, `'By'`, `'Bz'`, `'Jx'`, `'Jy'`, `'Jz'`.

To diagnose all standard fields:
```python
data_list = ['rho', 'E', 'B', 'J']
```

For electrostatic-only simulations:
```python
data_list = ['rho', 'E', 'phi']
```

For per-species charge densities (useful for visualizing species separately):
```python
data_list = ['rho_electrons', 'rho_ar_ions', 'E', 'phi']
```

### `warpx_format`

The output format. Major options:

- **`'openpmd'`** — the recommended modern format. Self-describing, interoperable with openPMD-viewer, openPMD-api, yt, and many other tools. Supports HDF5 (`warpx_openpmd_backend='h5'`) and ADIOS2 (`'bp'`) backends.
- **`'plotfile'`** — AMReX's native format. Directory structure with binary `Level_N/` subdirectories. Readable by yt directly.
- **`'checkpoint'`** — for restart files; full simulation state including particles.
- **`'ascent'`** — in-situ visualization via Ascent (compile-time flag required).
- **`'sensei'`** — in-situ visualization via SENSEI (compile-time flag required).

Default to `'openpmd'` for new scripts. Use `'plotfile'` only if your post-processing pipeline specifically depends on yt's plotfile loader.

### `warpx_random_fraction` and `warpx_uniform_stride` (subsampling)

For huge particle counts where writing all particles would be wasteful:

```python
particle_diag = picmi.ParticleDiagnostic(
    ...,
    warpx_random_fraction=0.01,      # write 1% of particles, randomly selected
    # or
    warpx_uniform_stride=100,        # write every 100th particle
)
```

`random_fraction` is statistically uniform — preserves distribution shape on average. `uniform_stride` is deterministic — useful when you need reproducible particle subsets across runs.

### `warpx_plot_filter_function`

A string expression filter, evaluated per particle. Only particles for which the expression is truthy are written:

```python
particle_diag = picmi.ParticleDiagnostic(
    ...,
    warpx_plot_filter_function='(uz > 1e7) and (z > 0)',   # only forward-moving in z>0
)
```

Variables: position (`x`, `y`, `z`), momentum (`ux`, `uy`, `uz`), weight (`w`). Useful for diagnostic restriction to specific phase-space regions (e.g., the high-energy tail of a distribution).

## ParticleDiagnostic

```python
particle_diag = picmi.ParticleDiagnostic(
    name='diag1',                       # often shared with field diag for openPMD consolidation
    period=1000,
    species=[electrons, ar_ions],       # which species (or None for all)
    data_list=['position', 'momentum', 'weighting', 'fields'],
    write_dir='./diags',
    parallelio=True,
    warpx_format='openpmd',
    warpx_random_fraction=None,
    warpx_uniform_stride=None,
    warpx_plot_filter_function=None,
)
sim.add_diagnostic(particle_diag)
```

### `data_list` options

- `'position'` — particle positions
- `'momentum'` — particle momenta (ux, uy, uz)
- `'weighting'` — macroparticle weight
- `'fields'` — interpolated fields at particle positions (large memory cost)
- Specific runtime attributes by name, e.g., `'prev_x'`, `'ionizationLevel'`, `'orig_x'`

If `data_list` is `None` (default), an implementation-default set is written.

### Sharing diagnostic name between Field and Particle

Setting both `FieldDiagnostic` and `ParticleDiagnostic` to `name='diag1'` (and the same `write_dir`) consolidates them into the same openPMD series — particles and fields at the same iteration appear together. This is the standard pattern for openPMD output:

```python
field_diag = picmi.FieldDiagnostic(name='diag1', ..., period=1000, write_dir='./diags', warpx_format='openpmd')
particle_diag = picmi.ParticleDiagnostic(name='diag1', ..., period=1000, write_dir='./diags', warpx_format='openpmd')
```

Both diagnostics with the same name end up in the same series, with iterations indexed by step number.

## Output formats

The choice of `warpx_format` determines what tools can read the output.

### openPMD (recommended)

```python
field_diag = picmi.FieldDiagnostic(
    ...,
    warpx_format='openpmd',
    warpx_openpmd_backend='h5',          # or 'bp' for ADIOS2
    warpx_file_prefix='diag1',
    warpx_file_min_digits=6,             # zero-pad step numbers to 6 digits
)
```

Output structure (HDF5 backend):
```
diags/
├── diag1/
│   ├── openpmd_000000.h5
│   ├── openpmd_001000.h5
│   ├── openpmd_002000.h5
│   └── ...
```

Read with openPMD-viewer:
```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('./diags/diag1/')
```

### Plotfile

```python
field_diag = picmi.FieldDiagnostic(
    ...,
    warpx_format='plotfile',
    warpx_file_prefix='plt',
)
```

Output structure:
```
diags/
├── diag1/
│   ├── plt00000/
│   │   ├── Header
│   │   ├── Level_0/
│   │   └── ...
│   ├── plt01000/
│   └── ...
```

Read with yt:
```python
import yt
ds = yt.load('./diags/diag1/plt00000/')
```

## Reduced diagnostics

Reduced diagnostics produce scalar or low-dimensional quantities at high cadence (typically every step). They're configured via `warpx_*` parameters on the Simulation rather than as separate PICMI classes:

```python
sim = picmi.Simulation(
    ...,
    warpx_reduced_diags_names=['fp', 'pe', 'be', 'np'],
    # Plus per-diagnostic configuration:
)
```

Each reduced diagnostic name corresponds to a quantity type. The configuration syntax in PICMI is awkward (requires setting attributes directly on the underlying `pywarpx.warpx` instance) — the cleanest approach is often to set them via the input file mechanism.

The reduced diagnostic types include:
- **`ParticleEnergy`** — total kinetic energy per species
- **`FieldEnergy`** — total E and B energy
- **`ParticleNumber`** — total particle count per species
- **`BeamRelevant`** — beam-physics quantities (emittance, energy spread, etc.)
- **`ParticleHistogram`** — 1D or 2D histograms of particle attributes
- **`FieldProbe`** — point-probe field values
- **`LoadBalanceCosts`** — load balancing diagnostics
- **`ChargeOnEB`** — charge on embedded boundary
- **`FieldMaximum`** — peak field values
- **`RhoMaximum`** — peak density

Each produces an ASCII text file (typically `diags/reducedfiles/<name>.txt`) with one row per time step. Cheap, dense in time, easy to plot in Python or any plotting tool.

For PICMI configuration of reduced diagnostics, the cleanest approach is to:

```python
# Use the underlying pywarpx prefix machinery
from pywarpx import warpx as wx

wx.warpx.reduced_diags_names = "fp pe np"
wx.fp.type = "FieldEnergy"
wx.fp.intervals = "10"
wx.fp.path = "./diags/reducedfiles/"
wx.pe.type = "ParticleEnergy"
wx.pe.intervals = "10"
wx.np.type = "ParticleNumber"
wx.np.intervals = "10"
```

Or write these as `warpx.*` parameters in a preprocessor-mode input file. Check the current readthedocs Reduced Diagnostics page for the exact configuration of each type.

## Checkpoints and restarts

```python
checkpoint = picmi.Checkpoint(
    period=10000,
    write_dir='./diags',
    name='chk',
    warpx_file_prefix='chk',
    warpx_file_min_digits=6,
)
sim.add_diagnostic(checkpoint)
```

This writes a full simulation state (fields + all particles + RNG state) every `period` steps. Checkpoints are stored as plotfiles (not openPMD).

### Restarting

To restart from a checkpoint:

```python
sim = picmi.Simulation(
    ...,
    warpx_amr_restart='./diags/chk001000',    # path to a checkpoint directory
)
```

The simulation resumes from the step the checkpoint was written. Diagnostics continue from there too.

Restarts preserve:
- All particle data (positions, momenta, weights, attributes)
- All field data (E, B, J, ρ, ϕ if present)
- Time, step counter, and random state

Restarts do *not* preserve:
- Python state (callbacks, ML model parameters, accumulator variables)

If your simulation has stateful Python (e.g., a ML model whose parameters update over time), save and load that state separately.

## Particle filtering and subsampling

For particle diagnostics of large populations:

```python
particle_diag = picmi.ParticleDiagnostic(
    ...,
    warpx_random_fraction=0.01,                          # ~1% of particles
    # or
    warpx_uniform_stride=100,                            # every 100th
    # or
    warpx_plot_filter_function='(z > 0.01) and (uz > 1e7)',  # forward-going only
)
```

These can combine: `random_fraction` after `plot_filter_function` applies the filter and then subsamples.

For high-energy tail diagnostics, filtering by energy:
```python
warpx_plot_filter_function='sqrt(ux**2 + uy**2 + uz**2) / c > 0.5'    # particles above v/c = 0.5
```

(where `c` is implicit — confirm in current docs if `c` is available in the parser scope or whether you need to substitute a literal.)

## Lab-frame diagnostics for boosted-frame runs

For boosted-frame simulations (`gamma_boost > 1` on the Simulation), the on-the-fly diagnostics are in the boosted frame. To get lab-frame output, use `LabFrameFieldDiagnostic`:

```python
btd = picmi.LabFrameFieldDiagnostic(
    name='btd',
    grid=grid,
    num_snapshots=10,
    dt_snapshots=1e-12,                        # lab-frame time between snapshots
    data_list=['E', 'B', 'rho', 'J'],
    z_subsampling=1,
    time_start=0.0,
    write_dir='./diags',
    parallelio=True,
    warpx_format='openpmd',
    warpx_buffer_size=256,
)
sim.add_diagnostic(btd)
```

WarpX performs the back-transformation automatically — at each step it identifies which lab-frame snapshot times intersect the current boosted-frame domain and accumulates contributions. The snapshots are written progressively.

For boosted-frame simulations of long beam-plasma propagation, BTD is essential — the boosted-frame domain may be short and travels through the physics, while lab-frame snapshots see the full propagation.

Particle BTD is also available via `LabFrameParticleDiagnostic` (similar API).

## Particle-field diagnostics (derived quantities)

To compute derived per-species fields (e.g., per-species density, mean velocity, temperature) as separate MultiFabs at diagnostic time:

```python
# Define the particle-field diagnostics:
density_diag = picmi.ParticleFieldDiagnostic(
    name='ne',
    func='1',                        # numerator: 1 = density (counts of particles per cell)
    do_average=False,                # if True, divides by particle count = average instead of sum
    filter='',                        # optional filter (string expression)
)
mean_uz_diag = picmi.ParticleFieldDiagnostic(
    name='mean_uz',
    func='uz',                       # numerator: uz
    do_average=True,                 # divides by particle count = mean uz
    filter='',
)

# Then attach to a FieldDiagnostic:
field_diag = picmi.FieldDiagnostic(
    name='diag1',
    grid=grid,
    period=1000,
    data_list=['rho', 'E', 'B'],
    warpx_particle_fields_to_plot=[density_diag, mean_uz_diag],
    warpx_particle_fields_species=['electrons'],     # which species
)
```

The result is an additional field per species (e.g., `ne_electrons`, `mean_uz_electrons`) in the diagnostic output. Useful for visualizing per-species moments without having to post-process particle dumps.

## Performance considerations

### Period choice

A few rules of thumb:
- **Full snapshots every 100-1000 steps** is typical for production runs. Smaller periods produce huge output volume.
- **Reduced diagnostics every 1-10 steps** is fine — these are scalars/vectors and trivial in size.
- **Per-step diagnostics for the entire run** is justified only when needed (e.g., capturing transient phenomena). For most science, snapshot every 100-1000 steps and reduced diagnostics every step is sufficient.

For a typical capacitive discharge benchmark (50,000 steps), 50 field+particle snapshots is plenty; that's `period=1000`.

### Output volume

- Field diagnostics: scales with grid size × number_of_fields × bytes_per_value. For a 256³ grid in double precision with 12 fields (E, B, J, rho×3 = 10 + ϕ), that's ~1.5 GB per snapshot. With 50 snapshots, 75 GB total. Plan disk space accordingly.
- Particle diagnostics: scales with particle count × attributes_per_particle × bytes. For 100M particles with 7 attributes (x, y, z, ux, uy, uz, w) in double precision, that's ~5.6 GB per species per snapshot. Use `warpx_random_fraction` or `warpx_uniform_stride` aggressively for huge species.

### openPMD vs plotfile

For modest output sizes (< 100 GB total), either format is fine. For very large outputs:
- ADIOS2 BP (set `warpx_openpmd_backend='bp'`) is more efficient than HDF5 for parallel writes.
- Plotfile native parallel I/O can be fastest if you don't need openPMD interoperability.
- openPMD with ADIOS2 + parallel HDF5 is the safest bet for HPC.

### Parallel I/O

`parallelio=True` (default) writes via MPI-IO. On HPC systems, this requires that openPMD-api was built with parallel HDF5 (or ADIOS2). If your build is serial-only, set `parallelio=False` — each rank writes its own file, slower but works.

For very large simulations on systems with weak file systems, consider striping the output directory (`lfs setstripe` on Lustre) for higher write bandwidth.

### Subsampling for huge particle counts

A simulation with 10⁹ macroparticles produces ~50 GB per particle snapshot in double precision. Two strategies:

1. **Random subsample**: `warpx_random_fraction=0.001` gives 10⁶ particles per snapshot (~50 MB) — usually statistically sufficient for distributions.
2. **Reduced diagnostics**: instead of dumping particles, dump derived moments (density, mean velocity, temperature) as `ParticleFieldDiagnostic` — orders of magnitude less data.

For most analyses, moments are what you actually need; full particle dumps are valuable for phase-space diagnostics and specific events (collision spectra, boundary fluxes).
