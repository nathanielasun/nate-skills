# Post-processing WarpX output in Python

This reference covers analyzing WarpX output files: openPMD-viewer for the standard openPMD HDF5/ADIOS2 output, openPMD-api for parallel and full-control reads, yt for AMReX plotfiles, reduced diagnostic ASCII files, and boundary-scraper output. Read this when you've run a simulation and need to extract, plot, or post-process the results.

## Table of contents

1. [Which tool for which output](#which-tool-for-which-output)
2. [openPMD-viewer: the workhorse for openPMD output](#openpmd-viewer-the-workhorse-for-openpmd-output)
3. [openPMD-api: parallel and full-control](#openpmd-api-parallel-and-full-control)
4. [yt: AMReX plotfiles](#yt-amrex-plotfiles)
5. [Reduced diagnostic ASCII files](#reduced-diagnostic-ascii-files)
6. [Boundary-scraper output](#boundary-scraper-output)
7. [Common analysis patterns](#common-analysis-patterns)
8. [Parallel post-processing](#parallel-post-processing)
9. [Performance and memory tips](#performance-and-memory-tips)

## Which tool for which output

| Output type | Files | Best tool |
|---|---|---|
| `warpx_format='openpmd'` field+particle | `diag1/openpmd_NNNNNN.h5` (or `.bp`) | **openPMD-viewer** for most work |
| Same, but need parallel reads or full openPMD metadata | same | **openPMD-api** |
| `warpx_format='plotfile'` | `diag1/pltNNNNNN/` directories | **yt** |
| Reduced diagnostics | `diags/reducedfiles/*.txt` | numpy / pandas |
| Boundary scrapers (`warpx_save_particles_at_*`) | `diags/boundary_scrapers*/` (openPMD format) | **openPMD-viewer** or openPMD-api |
| Checkpoints | `chkNNNNNN/` | restart only — not for analysis |
| Lab-frame back-transformed (BTD) | `diags/btd/openpmd_NNNNNN.h5` | **openPMD-viewer** |

openPMD-viewer is the right starting tool for nearly any analysis. Drop to openPMD-api when you need parallel chunked reads or access to attributes openPMD-viewer doesn't expose. Use yt only for plotfile output.

## openPMD-viewer: the workhorse for openPMD output

```python
from openpmd_viewer import OpenPMDTimeSeries

ts = OpenPMDTimeSeries('./diags/diag1/')
```

The `OpenPMDTimeSeries` object wraps a directory of openPMD files and exposes them as a unified series.

### Inspecting available data

```python
ts.iterations          # numpy array of available step numbers
ts.t                    # numpy array of physical times (seconds)
ts.avail_fields         # list of available field names: ['E', 'B', 'rho', ...]
ts.avail_components     # dict: {'E': ['x', 'y', 'z'], 'B': ['x', 'y', 'z'], ...}
ts.avail_species        # list of species names
ts.avail_record_components  # dict per species of available attributes
ts.fields_metadata      # metadata per field (geometry, axis labels, etc.)
ts.geometry             # 'cartesian', 'thetaMode', etc.
```

The first call to any of these does a fast scan of the directory; subsequent calls are cached.

### Reading fields

```python
# By iteration number:
Ex, info = ts.get_field(field='E', coord='x', iteration=1000)
rho, info = ts.get_field(field='rho', iteration=1000)
phi, info = ts.get_field(field='phi', iteration=1000)

# By time (closest available):
Ex, info = ts.get_field(field='E', coord='x', t=1.5e-9)
```

The return is a 2-tuple: `(array, FieldMetaInformation)`. The `info` object has axes coordinates, grid spacing, geometry info:

```python
info.x                  # numpy array of x positions (cell centers)
info.z                  # numpy array of z positions
info.dx                 # cell size
info.imshow_extent      # (xmin, xmax, ymin, ymax) for matplotlib imshow()
info.axes               # ['z', 'x'] or similar
```

For 2D simulations, the field array is 2D `(nz, nx)`. For 3D, 3D arrays. Easy plotting:

```python
import matplotlib.pyplot as plt

Ex, info = ts.get_field(field='E', coord='x', iteration=1000)
plt.imshow(Ex.T, extent=info.imshow_extent, origin='lower', aspect='auto', cmap='RdBu_r')
plt.colorbar(label='E_x [V/m]')
plt.xlabel(info.axes[0])
plt.ylabel(info.axes[1])
plt.show()
```

### Reading particles

```python
# Get particle positions and momenta at a specific iteration
x, y, z, ux, uy, uz, w = ts.get_particle(
    var_list=['x', 'y', 'z', 'ux', 'uy', 'uz', 'w'],
    species='electrons',
    iteration=1000,
)
```

Each variable is returned as a flat numpy array of length `N_particles`. For 2D simulations, only `x` and `z` are meaningful as positions.

To select a subset of particles:

```python
# By boolean filter (any expression in particle variables):
x, ux = ts.get_particle(
    var_list=['x', 'ux'],
    species='electrons',
    iteration=1000,
    select={'z': [0.005, 0.010]},      # only particles with z in this range
)
```

`select` is a dict where keys are variable names and values are `[min, max]` ranges. Multiple keys produce an AND filter.

### Particle tracking via `ParticleTracker`

To reconstruct trajectories of specific particles across iterations:

```python
from openpmd_viewer import ParticleTracker

# Select particles at iteration 100 by some criterion
pt = ParticleTracker(
    ts, iteration=100, species='electrons',
    select={'z': [0.005, 0.010]},
    # Optionally constrain by attribute:
    # preserve_particle_index=True,
)

# Get the same particles' positions at iteration 500
x, ux = ts.get_particle(
    var_list=['x', 'ux'],
    species='electrons',
    iteration=500,
    select=pt,
)
```

The `ParticleTracker` identifies particles by their unique ID (assigned at creation) and follows them through subsequent dumps. Particles that left the simulation between the two iterations are silently dropped.

### Iterating over a time series

```python
import numpy as np

mean_ne = []
for iteration in ts.iterations:
    rho, info = ts.get_field(field='rho', iteration=iteration)
    mean_ne.append(np.mean(rho) / -picmi.constants.q_e)  # n_e = -rho/q_e for negative charges

plt.plot(ts.t, mean_ne)
plt.xlabel('t [s]')
plt.ylabel('mean n_e [m^-3]')
```

For long series with large files, this becomes I/O-bound. Cache intermediate results to a `.npy` file.

### Geometry-specific notes

For RZ output (`geometry='thetaMode'`), fields come as `(nr, nz, n_theta)` arrays with mode information available via `info.thetaMode`. To get a Cartesian-friendly slice:

```python
# Reconstruct the field in a (z, x) plane at theta = 0:
Ex, info = ts.get_field(field='E', coord='x', iteration=1000, theta=0.0)
```

The `theta` argument reconstructs an azimuthal slice from the modal data.

## openPMD-api: parallel and full-control

When you need parallel reads or full access to openPMD metadata:

```python
import openpmd_api as io

series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only)

for iteration_num, iteration in series.iterations.items():
    # Field access
    E = iteration.meshes['E']
    Ex = E['x']
    data = Ex.load_chunk()    # lazy
    series.flush()             # actually load
    # data is a numpy array

    # Particle access
    electrons = iteration.particles['electrons']
    x = electrons['position']['x'].load_chunk()
    series.flush()
```

For parallel reads with MPI:

```python
from mpi4py import MPI
import openpmd_api as io

comm = MPI.COMM_WORLD
series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only, comm)

# Each rank reads a chunk:
iteration = series.iterations[1000]
Ex = iteration.meshes['E']['x']

# Determine this rank's chunk
extent = Ex.shape
chunk_size = extent[0] // comm.size
offset = [comm.rank * chunk_size] + [0] * (len(extent) - 1)
local_extent = [chunk_size] + list(extent[1:])

local_data = Ex.load_chunk(offset, local_extent)
series.flush()
# Each rank now has its slice of Ex
```

For Pandas/DASK/RAPIDS integration: see [openPMD-api docs on Pandas](https://openpmd-api.readthedocs.io/en/latest/analysis/pandas.html). The general pattern is `iteration.particles[species].to_df()` returns a pandas DataFrame; `to_dask()` returns a dask DataFrame for out-of-core operations.

```python
import pandas as pd

iteration = series.iterations[1000]
df = iteration.particles['electrons'].to_df()
# df has columns: position_x, position_y, position_z, momentum_x, ..., weight, ...

# Filter and aggregate
hot = df[df['momentum_z'] > 1e7]
print(f"Hot fraction: {len(hot) / len(df):.3%}")
```

## yt: AMReX plotfiles

For simulations with `warpx_format='plotfile'`:

```python
import yt

ds = yt.load('./diags/diag1/plt00100/')
ds.field_list       # list of fields
ds.particle_types   # list of particle species

# Field access:
ad = ds.all_data()                       # full domain
Ex = ad['boxlib', 'Ex'].to_ndarray()      # to plain numpy

# Slice plot:
slc = yt.SlicePlot(ds, 'z', 'Ex')        # slice perpendicular to z
slc.set_log('Ex', False)                  # linear color
slc.set_cmap('Ex', 'RdBu_r')
slc.save('Ex_slice.png')

# Projection (integrated along axis):
prj = yt.ProjectionPlot(ds, 'z', 'rho')
prj.save('rho_proj.png')

# Particle phase plot:
phase = yt.ParticlePhasePlot(
    ad, 'particle_position_x', 'particle_velocity_x', 'particle_weight',
    deposition='cic',
)
phase.save('phase.png')
```

### Numpy from yt

To get a uniformly-gridded numpy array (with AMR refinement flattened to a single level):

```python
ad0 = ds.covering_grid(
    level=0,
    left_edge=ds.domain_left_edge,
    dims=ds.domain_dimensions,
)
Ex_3d = ad0['Ex'].to_ndarray()   # shape: (nx, ny, nz)
```

`covering_grid(level=L)` returns the simulation domain at refinement level `L`. For non-AMR runs, `level=0` is the only choice.

For particles:

```python
ad = ds.all_data()
x = ad['electrons', 'particle_position_x'].to_ndarray()
ux = ad['electrons', 'particle_momentum_x'].to_ndarray()
w = ad['electrons', 'particle_weight'].to_ndarray()
```

### yt-vs-openPMD-viewer

Use yt when:
- The output is plotfile format and converting is impractical.
- You need yt's plotting features (SlicePlot, ProjectionPlot, ParticlePhasePlot) directly.
- You're already in a yt workflow.

Use openPMD-viewer when:
- Output is openPMD format (the recommended default).
- You're doing custom analysis with numpy/scipy.
- You want a lighter dependency footprint.

## Reduced diagnostic ASCII files

Reduced diagnostics produce `diags/reducedfiles/*.txt` with columns including step number, time, and the reduced quantities. Plain ASCII, one row per write.

Read with pandas:

```python
import pandas as pd
df = pd.read_csv('./diags/reducedfiles/particle_energy.txt',
                 sep=r'\s+', comment='#')
# Inspect columns
print(df.columns)
# Plot
df.plot(x='time(s)', y='total(J)', logy=True)
```

Or numpy:

```python
import numpy as np
data = np.loadtxt('./diags/reducedfiles/particle_energy.txt')
# Columns: step, time, energy_per_species, total
```

The exact column layout depends on the reduced diagnostic type. Check the header (first commented lines) of the file.

## Boundary-scraper output

When you set `warpx_save_particles_at_xlo=True` etc. on a Species, particles crossing that boundary are recorded in `diags/<scraper_name>/`. The format is openPMD-style — read with openPMD-viewer or openPMD-api like any other particle dump.

```python
from openpmd_viewer import OpenPMDTimeSeries
scraper = OpenPMDTimeSeries('./diags/boundary_scrapers_xlo/')

# Get all electrons that crossed xlo during the simulation:
x, ux, w, t = scraper.get_particle(
    var_list=['x', 'ux', 'w', 'timestamp'],
    species='electrons',
    iteration=scraper.iterations[-1],   # the latest snapshot contains all
)
```

`timestamp` is a per-particle attribute: the simulation time when each particle crossed the boundary. Useful for:
- Energy spectra of particles arriving at an electrode.
- Time-dependent flux (histogram by `timestamp`).
- Angle distributions (compute angles from `ux/uy/uz`).

```python
import matplotlib.pyplot as plt

# Energy spectrum on xlo
E_kin = 0.5 * picmi.constants.m_e * (ux**2 + uy**2 + uz**2)   # joules
E_eV = E_kin / picmi.constants.q_e
plt.hist(E_eV, bins=200, weights=w, log=True)
plt.xlabel('E [eV]')
plt.ylabel('Counts (weighted)')

# Flux vs time on xlo
plt.hist(t, bins=200, weights=w / dt)   # dt = bin width; arbitrary normalization
plt.xlabel('t [s]')
plt.ylabel('Flux (arb)')
```

## Common analysis patterns

### Density-time plot

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('./diags/diag1/')

ne_t = []
for it in ts.iterations:
    rho_e, info = ts.get_field(field='rho_electrons', iteration=it)
    # In 1D, rho_e is 1D; take mean over interior cells
    ne_avg = np.mean(rho_e[10:-10]) / -picmi.constants.q_e
    ne_t.append(ne_avg)

plt.semilogy(ts.t, ne_t)
plt.xlabel('t [s]')
plt.ylabel('<n_e> [m^-3]')
```

### Phase space at fixed time

```python
x, ux, w = ts.get_particle(
    var_list=['x', 'ux', 'w'],
    species='electrons',
    iteration=10000,
)
# Compute v from u (assumes u = γv, but for non-relativistic, u ≈ v)
import scipy.constants as sc
v = ux  # for non-relativistic

# 2D histogram
H, xedges, vedges = np.histogram2d(x, v, bins=200, weights=w)
plt.imshow(H.T, origin='lower',
           extent=[xedges[0], xedges[-1], vedges[0], vedges[-1]],
           aspect='auto', cmap='viridis', norm=plt.matplotlib.colors.LogNorm())
plt.colorbar(label='Counts')
plt.xlabel('x [m]')
plt.ylabel('v_x [m/s]')
```

### Energy-time plot

```python
KE_t = []
for it in ts.iterations:
    ux, uy, uz, w = ts.get_particle(
        var_list=['ux', 'uy', 'uz', 'w'],
        species='electrons', iteration=it,
    )
    # Non-relativistic KE: (1/2) m u^2  [u = velocity here, not γv; check the convention]
    # In WarpX openPMD output, 'ux' is the momentum component p_x = γ m v_x.
    # For non-relativistic: u_x ≈ v_x.
    KE = 0.5 * picmi.constants.m_e * (ux**2 + uy**2 + uz**2)
    KE_total = np.sum(w * KE)
    KE_t.append(KE_total)

plt.plot(ts.t, KE_t)
plt.xlabel('t [s]')
plt.ylabel('Total electron KE [J]')
```

**Convention note**: openPMD-viewer's `ux`, `uy`, `uz` for particles is *momentum per unit mass* in SI units — i.e., γβc, equal to γv. For non-relativistic plasma `γ ≈ 1` so this is just velocity. For relativistic, compute `γ = sqrt(1 + (u/c)²)` and `v = u/γ`. Always double-check the convention in the openPMD metadata (`info.unitDimension`) for safety.

### Comparing two simulations

```python
ts1 = OpenPMDTimeSeries('./run1/diags/diag1/')
ts2 = OpenPMDTimeSeries('./run2/diags/diag1/')

# Compare a field at matching times:
for t in np.linspace(0, 1e-8, 11):
    Ex1, info = ts1.get_field(field='E', coord='x', t=t)
    Ex2, info = ts2.get_field(field='E', coord='x', t=t)
    diff = np.linalg.norm(Ex1 - Ex2) / np.linalg.norm(Ex1)
    print(f"t={t:.2e}: relative diff = {diff:.3e}")
```

For benchmarking against a reference solution: use the L2 or L_inf norm of the difference; if `< 1%` for a converged solution, results agree.

## Parallel post-processing

For datasets too large for serial reads:

### openPMD-api with MPI

As shown above — each rank reads its chunk. Good for the I/O step; subsequent computation is whatever you make it.

### DASK

```python
import openpmd_api as io
import dask.dataframe as dd

series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only)
iteration = series.iterations[10000]
ddf = iteration.particles['electrons'].to_dask()
# Operations are lazy; .compute() to execute:
ne_mean = ddf['weight'].sum().compute() / volume
```

For analyzing the entire time series in parallel, run multiple iterations through a dask cluster. Good for large parameter sweeps where each simulation produces moderate output but you have many.

### mpi4py + openPMD-api hybrid

```python
from mpi4py import MPI
import openpmd_api as io

comm = MPI.COMM_WORLD
series = io.Series('./diags/diag1/openpmd_%T.h5', io.Access.read_only, comm)

# Distribute iterations across ranks
iter_list = list(series.iterations.keys())
my_iters = iter_list[comm.rank::comm.size]
for it in my_iters:
    iteration = series.iterations[it]
    # ... process ...
```

## Performance and memory tips

### Don't load full 3D fields into one rank

A 256³ double-precision field is 128 MB; six of them (E, B) is 768 MB. Reading multiple times into memory exhausts RAM quickly. Either:
- Chunk reads with openPMD-api (`load_chunk(offset, extent)`).
- Use yt's `covering_grid` with a coarser resolution.
- Only load 2D slices: openPMD-viewer's `get_field` for a 3D dataset will return the 3D array; for slices, use openPMD-api or yt's `SlicePlot`.

### Use `select=` aggressively for particle reads

For 10⁸-particle simulations, reading all particles into one array uses 8 bytes × 10⁸ × 6 (x, y, z, ux, uy, uz) = 4.8 GB minimum. If you only need a subset, use `select=` to filter at read time.

### Cache derived quantities

If you compute a moment for plotting, save it to disk so you don't recompute next time:

```python
cache_file = './cache_ne_t.npy'
if not os.path.exists(cache_file):
    ne_t = compute_density_time_series(ts)
    np.save(cache_file, ne_t)
else:
    ne_t = np.load(cache_file)
```

For long simulations and iterative analysis, this saves hours.

### Match write format to read pattern

If your typical analysis is "read full field at a few times," HDF5 with single-file-per-step is fine.

If your typical analysis is "read one slice across many times," BP5 (ADIOS2) with stride patterns is faster.

If your typical analysis is small-N particles with high time resolution, particle dump rate matters more than chunking — see reduced diagnostics.

### ADIOS2 BP5 vs HDF5

For very large outputs (TB+), BP5 is often 2-5x faster for parallel reads on HPC systems. HDF5 is more portable (every Python install has h5py available). Default to HDF5 unless I/O becomes a bottleneck.

To use BP5 backend at write time:
```python
field_diag = picmi.FieldDiagnostic(..., warpx_openpmd_backend='bp')
```

Read with the same APIs — the backend is transparent to openPMD-viewer and openPMD-api.
