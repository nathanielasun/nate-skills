# Output and analysis: diagnostics + post-processing

Combined reference for configuring diagnostics in PICMI scripts and reading the resulting output files with openPMD-viewer or yt.

## Diagnostic types

| PICMI class | Purpose |
|---|---|
| `picmi.FieldDiagnostic` | E, B, J, ρ, ϕ snapshots |
| `picmi.ParticleDiagnostic` | Particle position/momentum/weight snapshots |
| `picmi.Checkpoint` | Full state for restart |
| `picmi.LabFrameFieldDiagnostic` | Back-transformed for boosted-frame runs |
| (no class) Reduced diagnostics | Scalars per step (energy, particle count, etc.) |
| Boundary scrapers (via Species options) | Particles crossing boundaries |

Multiple diagnostics can coexist with different periods.

## FieldDiagnostic

```python
field_diag = picmi.FieldDiagnostic(
    name='diag1',                              # used in output dir name
    grid=grid,
    period=1000,                                # every N steps
    data_list=['rho', 'E', 'B', 'J', 'phi'],
    write_dir='./diags',
    warpx_format='openpmd',                     # 'plotfile', 'openpmd', 'ascent', 'sensei'
    warpx_openpmd_backend='h5',                 # 'h5', 'bp', 'json'
    warpx_file_prefix='diag1',
)
sim.add_diagnostic(field_diag)
```

### data_list options

- `'E'`, `'B'`, `'J'`: all three components (Ex, Ey, Ez written separately)
- `'rho'`: total charge density
- `'rho_<species>'`: per-species (e.g., `'rho_electrons'`)
- `'phi'`: scalar potential (ES only)
- Individual: `'Ex'`, `'Ey'`, `'Bx'`, `'Jz'`, etc.

For ES discharge work:
```python
data_list = ['rho', 'rho_electrons', 'rho_ar_ions', 'E', 'phi']
```

## ParticleDiagnostic

```python
part_diag = picmi.ParticleDiagnostic(
    name='diag1',                              # same name as FieldDiagnostic = openPMD consolidation
    period=1000,
    species=[electrons, ar_ions],               # None = all species
    data_list=['position', 'momentum', 'weighting'],
    write_dir='./diags',
    warpx_format='openpmd',
    # Subsampling for huge particle counts:
    warpx_random_fraction=None,                 # e.g., 0.01 = 1%
    warpx_uniform_stride=None,                  # e.g., 100 = every 100th
    warpx_plot_filter_function=None,            # e.g., '(uz > 1e7) and (z > 0)'
)
sim.add_diagnostic(part_diag)
```

`name='diag1'` for both Field and Particle → consolidated into one openPMD series with particles+fields at matching iterations.

### data_list

- `'position'`, `'momentum'`, `'weighting'`
- `'fields'` — interpolated fields at particle positions (heavy)
- Specific runtime attributes by name

## Output formats

### openPMD (recommended)

```python
warpx_format='openpmd', warpx_openpmd_backend='h5'    # or 'bp' for ADIOS2
```

Output: `diags/diag1/openpmd_NNNNNN.h5`. Read with openPMD-viewer or openPMD-api.

### Plotfile (legacy AMReX)

```python
warpx_format='plotfile'
```

Output: `diags/diag1/pltNNNNNN/` directories. Read with yt.

**Default to openPMD for new scripts.** Use plotfile only if your post-processing pipeline depends on yt's loader.

## Checkpoints

```python
checkpoint = picmi.Checkpoint(
    period=10000,
    write_dir='./diags',
    name='chk',
    warpx_file_prefix='chk',
)
sim.add_diagnostic(checkpoint)
```

Writes full state every period. To restart:

```python
sim = picmi.Simulation(...,
    warpx_amr_restart='./diags/chk001000',
)
```

**Restarts preserve**: fields, particles, RNG, time, step counter.
**Restarts do NOT preserve**: Python state (callbacks, ML models, accumulators).

## Reduced diagnostics

High-cadence scalar/vector quantities written to ASCII per step. Configured via `pywarpx.warpx` prefix machinery (PICMI doesn't have a clean class for them):

```python
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

Types include: `ParticleEnergy`, `FieldEnergy`, `ParticleNumber`, `BeamRelevant`, `ParticleHistogram`, `FieldProbe`, `ChargeOnEB`, `LoadBalanceCosts`, `FieldMaximum`, `RhoMaximum`.

Output: ASCII text in `diags/reducedfiles/<name>.txt`, one row per step. Easy to read with `np.loadtxt` or pandas.

## Performance considerations

| Setting | Cost |
|---|---|
| `period=1` (every step) | Massive — only use for short runs or critical transients |
| `period=100-1000` | Typical production |
| Full particle dump 100M particles | ~5.6 GB per species per snapshot |
| `warpx_random_fraction=0.01` | 1% subsample → ~50 MB instead of 50 GB |
| ParticleFieldDiagnostic (moments only) | Orders of magnitude less than full dump |

For huge runs: use `warpx_random_fraction` aggressively, or compute moments via `ParticleFieldDiagnostic` instead of dumping full particles.

ADIOS2 BP backend (`warpx_openpmd_backend='bp'`) is often 2-5x faster than HDF5 on HPC for parallel writes.

---

# Post-processing

Now the analysis side. Three main tools:
- **openPMD-viewer** — best for openPMD output (the recommended format)
- **openPMD-api** — parallel reads, full openPMD control
- **yt** — for plotfile format only

## openPMD-viewer

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('./diags/diag1/')
```

### Inspection

```python
ts.iterations             # numpy array of step numbers
ts.t                       # numpy array of physical times (s)
ts.avail_fields            # ['E', 'B', 'rho', ...]
ts.avail_components        # {'E': ['x', 'y', 'z'], ...}
ts.avail_species           # species names
ts.geometry                # 'cartesian', 'thetaMode', etc.
```

### Reading fields

```python
# By iteration:
Ex, info = ts.get_field(field='E', coord='x', iteration=1000)
rho, info = ts.get_field(field='rho', iteration=1000)
phi, info = ts.get_field(field='phi', iteration=1000)

# By time (closest available):
Ex, info = ts.get_field(field='E', coord='x', t=1.5e-9)
```

Returns `(array, info)`. The `info` object has axis coordinates:

```python
info.x, info.z              # cell-center positions
info.dx                      # cell size
info.imshow_extent           # (xmin, xmax, ymin, ymax) for matplotlib
info.axes                    # axis labels
```

Plot:

```python
import matplotlib.pyplot as plt
plt.imshow(Ex.T, extent=info.imshow_extent, origin='lower',
           aspect='auto', cmap='RdBu_r')
plt.colorbar(label='E_x [V/m]')
```

### Reading particles

```python
x, y, z, ux, uy, uz, w = ts.get_particle(
    var_list=['x', 'y', 'z', 'ux', 'uy', 'uz', 'w'],
    species='electrons',
    iteration=1000,
)
```

Returns flat numpy arrays of length N_particles.

Filtering at read time:

```python
x, ux = ts.get_particle(
    var_list=['x', 'ux'],
    species='electrons', iteration=1000,
    select={'z': [0.005, 0.010]},        # only z in this range
)
```

### Particle tracking

```python
from openpmd_viewer import ParticleTracker

# Pick particles at iter 100:
pt = ParticleTracker(ts, iteration=100, species='electrons',
                      select={'z': [0.005, 0.010]})

# Same particles at iter 500:
x, ux = ts.get_particle(var_list=['x', 'ux'],
                         species='electrons', iteration=500,
                         select=pt)
```

Identifies particles by unique ID; particles that left the box are silently dropped.

### Iterating over the time series

```python
import numpy as np
mean_ne = []
for it in ts.iterations:
    rho, info = ts.get_field(field='rho', iteration=it)
    mean_ne.append(-np.mean(rho) / picmi.constants.q_e)
plt.plot(ts.t, mean_ne)
```

Long series + large files → I/O-bound. Cache intermediate results.

## openPMD-api

When you need parallel reads or full openPMD metadata:

```python
import openpmd_api as io
series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only)

for iter_num, iteration in series.iterations.items():
    Ex = iteration.meshes['E']['x']
    data = Ex.load_chunk()
    series.flush()        # actually load
    # data is numpy array
```

Parallel:

```python
from mpi4py import MPI
comm = MPI.COMM_WORLD
series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only, comm)

# Each rank reads a chunk
iteration = series.iterations[1000]
Ex = iteration.meshes['E']['x']
extent = Ex.shape
chunk = extent[0] // comm.size
offset = [comm.rank * chunk] + [0] * (len(extent) - 1)
local = Ex.load_chunk(offset, [chunk] + list(extent[1:]))
series.flush()
```

Pandas / DASK integration:

```python
df = iteration.particles['electrons'].to_df()
# columns: position_x, position_y, momentum_x, ..., weight
hot = df[df['momentum_z'] > 1e7]
```

## yt (for plotfile format)

```python
import yt

ds = yt.load('./diags/diag1/plt00100/')
ds.field_list
ds.particle_types

# Field as numpy array (uniformly gridded):
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                       dims=ds.domain_dimensions)
Ex = ad0['Ex'].to_ndarray()

# Slice plot:
slc = yt.SlicePlot(ds, 'z', 'Ex')
slc.set_cmap('Ex', 'RdBu_r')
slc.save('Ex_slice.png')

# Phase plot:
ad = ds.all_data()
phase = yt.ParticlePhasePlot(ad, 'particle_position_x',
                              'particle_velocity_x', 'particle_weight')
phase.save('phase.png')
```

Use yt only if your output is plotfile format. Otherwise openPMD-viewer is lighter.

## Reduced diagnostic files

```python
import pandas as pd
df = pd.read_csv('./diags/reducedfiles/particle_energy.txt',
                 sep=r'\s+', comment='#')
df.plot(x='time(s)', y='total(J)', logy=True)
```

Column layout depends on the diagnostic type — check the header.

## Boundary scrapers

When you set `warpx_save_particles_at_xlo=True` etc. on a Species, particles crossing that boundary are recorded in `diags/<scraper_name>/`. Format is openPMD-style:

```python
scraper = OpenPMDTimeSeries('./diags/boundary_scrapers_xlo/')
x, ux, w, t = scraper.get_particle(
    var_list=['x', 'ux', 'w', 'timestamp'],
    species='electrons',
    iteration=scraper.iterations[-1],
)

# Energy spectrum:
E_eV = 0.5 * picmi.constants.m_e * (ux**2 + uy**2 + uz**2) / picmi.constants.q_e
plt.hist(E_eV, bins=200, weights=w, log=True)
```

`timestamp` is per-particle: simulation time when the particle crossed.

## Common analysis patterns

### Density vs time

```python
ne_t = []
for it in ts.iterations:
    rho_e, info = ts.get_field(field='rho_electrons', iteration=it)
    ne_t.append(-np.mean(rho_e[10:-10]) / picmi.constants.q_e)
plt.semilogy(ts.t, ne_t)
```

### Phase space at fixed time

```python
x, ux, w = ts.get_particle(var_list=['x', 'ux', 'w'],
                            species='electrons', iteration=10000)
H, xedges, vedges = np.histogram2d(x, ux, bins=200, weights=w)
plt.imshow(H.T, origin='lower',
           extent=[xedges[0], xedges[-1], vedges[0], vedges[-1]],
           aspect='auto', cmap='viridis')
```

### Comparing two runs

```python
ts1 = OpenPMDTimeSeries('./run1/diags/diag1/')
ts2 = OpenPMDTimeSeries('./run2/diags/diag1/')

for t in np.linspace(0, 1e-8, 11):
    E1, _ = ts1.get_field(field='E', coord='x', t=t)
    E2, _ = ts2.get_field(field='E', coord='x', t=t)
    diff = np.linalg.norm(E1 - E2) / np.linalg.norm(E1)
    print(f"t={t:.2e}: rel diff = {diff:.3e}")
```

## Performance tips

- **Don't load full 3D fields into one rank**: 256³ doubles = 128 MB per field. Use openPMD-api `load_chunk`, yt `covering_grid`, or stick to 2D slices.
- **Use `select=`** for particle reads: avoids materializing all particles to subset later.
- **Cache derived quantities** to `.npy` for iterative analysis.
- **Match write format to read pattern**: HDF5 for portability; BP5 for HPC parallel reads.

## Critical convention notes

- **openPMD-viewer `ux`/`uy`/`uz` are momentum/mass = γv** (not raw v). For non-relativistic (γ≈1) they're approximately velocity. For relativistic, compute `γ = sqrt(1 + (u/c)²)` and `v = u/γ`. Always check `info.unitDimension` for the exact convention.
- **`info.imshow_extent`** is sized for matplotlib's `imshow` directly — pass it to the `extent` argument.
- **`covering_grid(level=L)`** for AMR data flattens to that level's resolution. Non-AMR runs only have level 0.
