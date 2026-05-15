# Install, run, debug, read output

Combined operational reference. Install pywarpx, run with MPI, debug, and read output with openPMD-viewer or yt.

## Install

```bash
# Most common — pip:
pip install pywarpx[3D]                          # 3D only
pip install pywarpx[1D,2D,3D,RZ]                  # all dimensions

# Or conda:
conda install -c conda-forge warpx
```

Verify:
```bash
python3 -c "from pywarpx import picmi; print(picmi.__file__)"
```

## Build from source (for GPU or development)

```bash
git clone --depth 1 https://github.com/BLAST-WarpX/warpx.git
cd warpx

cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="3"
cmake --build build -j 8 --target pip_install
```

Verify:
```bash
python3 -c "from pywarpx import picmi, warpx; print('OK')"
```

### Common cmake options

| Option | Default | Meaning |
|---|---|---|
| `WarpX_PYTHON` | OFF | Build pywarpx |
| `WarpX_DIMS` | "3" | Dimensions; "1;2;3;RZ" for all |
| `WarpX_COMPUTE` | OMP | `OMP`, `CUDA`, `HIP`, `SYCL` |
| `WarpX_FFT` | OFF | For PSATD |
| `WarpX_QED` | ON | QED module |
| `WarpX_PYTHON_IPO` | ON | LTO (OFF for dev speed) |
| `AMReX_CUDA_ARCH` | autodetect | CUDA arch (80=A100, 90=H100) |

For typical low-temp plasma:
```bash
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;3;RZ" \
    -DWarpX_FFT=OFF -DWarpX_QED=OFF
cmake --build build -j 8 --target pip_install
```

For GPU:
```bash
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="3" \
    -DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=80
```

Fast rebuild during iteration:
```bash
cmake --build build -j 8 --target pip_install_nodeps
```

## pywarpx module map

| Module | Use |
|---|---|
| `pywarpx.picmi` | PICMI classes (the public interface) |
| `pywarpx.callbacks` | `installcallback`, `callfromafterstep` |
| `pywarpx.fields` | Field wrappers: `ExWrapper`, `RhoFPWrapper` |
| `pywarpx.particle_containers` | `ParticleContainerWrapper` |
| `pywarpx.libwarpx.amr` | pyAMReX (raw AMReX) |

## Run a script

```bash
# Serial:
python3 my_script.py

# MPI:
mpiexec -n 4 python3 my_script.py

# Hybrid MPI+OpenMP:
export OMP_NUM_THREADS=4
mpiexec -n 8 python3 my_script.py        # 32 cores total

# GPU (1 GPU per rank typical):
mpiexec -n 4 python3 my_script.py
```

### MPI init order matters

```python
# Correct:
from mpi4py import MPI                  # initializes MPI first
from pywarpx import picmi                # uses already-initialized MPI

comm = MPI.COMM_WORLD
```

## HPC submission (Slurm)

```bash
#!/bin/bash
#SBATCH --job-name=my_run
#SBATCH --time=2:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4

module load python/3.11 openmpi/4.1 hdf5
source /path/to/venv/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun python3 my_script.py
```

Output to scratch, not $HOME:
```python
field_diag = picmi.FieldDiagnostic(..., write_dir='/scratch/$USER/run/diags')
```

Site-specific machine profiles: `Tools/machines/<system>/` in WarpX repo (perlmutter-nersc, frontier-olcf, lumi-csc, etc.).

## Common errors

| Error | Fix |
|---|---|
| `ImportError: libwarpx_pyDIM.so cannot be opened` | Wrong dim variant. `pip install pywarpx[3D]` or rebuild with right `WarpX_DIMS` |
| `AttributeError: module 'pywarpx.picmi' has no attribute 'XxxYyy'` | Class renamed in your version. Check `https://warpx.readthedocs.io/en/<version>/usage/python.html` |
| `Could not parse expression '...'` | Bad string expression. Use `x`, `y`, `z`, `t`; `sin/cos/exp/sqrt`; ternary as `if(cond, t, f)`; pass params as kwargs |
| `MPI_Init was already called` / hangs at start | Import `mpi4py.MPI` before `pywarpx` |
| `Cross-section file not found` | Clone `BLAST-WarpX/warpx-data`. Use absolute paths |
| Simulation hangs at start | MPI ranks diverge. Test with `-n 1` first. Add rank-numbered prints |
| GPU not detected | Rebuild with `-DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=<arch>` |
| Memory exhaustion with ionization | Enable `warpx_do_resampling=True` on charged species |

Check GPU available:
```bash
python3 -c "from pywarpx import libwarpx; print('GPU:', libwarpx.Config.have_gpu)"
```

## Debugging

### Verbose

```python
sim = picmi.Simulation(..., verbose=2)         # 0=silent, 1=normal, 2=debug
```

### Inspect generated input file

```python
# Run this once before sim.step() to see the C++ params your PICMI maps to:
sim.write_input_file('inputs_debug')
sim.step()
```

Then `cat inputs_debug` and verify each `warpx_*` arg shows as expected. Many bugs become obvious.

### Trap NaN/Inf early

```python
sim = picmi.Simulation(...,
    warpx_amrex_fpe_trap_invalid=True,
    warpx_amrex_fpe_trap_zero=True,
    warpx_amrex_fpe_trap_overflow=True,
)
```

### Print debugging in MPI

```python
if MPI.COMM_WORLD.rank == 0:
    print(f"Step {step}: max E = {max_E}", flush=True)
```

### Reproducing failures

1. Run 100 steps on 128 ranks. Still fails? → issue is initial.
2. Run 100 steps on 1 rank. Still fails? → issue is non-MPI.
3. Run with `verbose=2`. What's the last printed thing?
4. Reduce dims (3D→2D), grid size (256→32), particles for fast iteration.

---

# Reading output

Choose tool by output format:

| Output format | Tool |
|---|---|
| `warpx_format='openpmd'` (default) | **openPMD-viewer** |
| Same, need parallel/full metadata | openPMD-api |
| `warpx_format='plotfile'` | **yt** |
| Boundary scrapers | openPMD-viewer (openPMD format) |
| Reduced diagnostics | numpy / pandas |

## openPMD-viewer (main tool)

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('./diags/diag1/')
```

### Inspection

```python
ts.iterations             # step numbers
ts.t                       # physical times (s)
ts.avail_fields            # ['E', 'B', 'rho', ...]
ts.avail_components        # {'E': ['x','y','z'], ...}
ts.avail_species           # species names
```

### Read fields

```python
# By iteration:
Ex, info = ts.get_field(field='E', coord='x', iteration=1000)
rho, info = ts.get_field(field='rho', iteration=1000)
phi, info = ts.get_field(field='phi', iteration=1000)

# By time (closest available):
Ex, info = ts.get_field(field='E', coord='x', t=1.5e-9)
```

The `info` object has:
```python
info.x, info.z              # cell-center positions
info.dx                      # cell size
info.imshow_extent           # for matplotlib imshow extent
```

### Plot

```python
import matplotlib.pyplot as plt
plt.imshow(Ex.T, extent=info.imshow_extent, origin='lower',
           aspect='auto', cmap='RdBu_r')
plt.colorbar(label='E_x [V/m]')
```

### Read particles

```python
x, y, z, ux, uy, uz, w = ts.get_particle(
    var_list=['x', 'y', 'z', 'ux', 'uy', 'uz', 'w'],
    species='electrons',
    iteration=1000,
)

# Filter at read time:
x, ux = ts.get_particle(
    var_list=['x', 'ux'],
    species='electrons', iteration=1000,
    select={'z': [0.005, 0.010]},        # range filter
)
```

**Note**: `ux`/`uy`/`uz` are `γv` (momentum/mass), not raw velocity. For non-relativistic γ≈1, so they're approximately velocity.

### Time series patterns

```python
import numpy as np

# Density vs time:
ne_t = []
for it in ts.iterations:
    rho_e, info = ts.get_field(field='rho_electrons', iteration=it)
    ne_t.append(-np.mean(rho_e) / picmi.constants.q_e)
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

## yt (for plotfile format only)

```python
import yt
ds = yt.load('./diags/diag1/plt00100/')
ds.field_list
ds.particle_types

# Field as numpy:
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                       dims=ds.domain_dimensions)
Ex = ad0['Ex'].to_ndarray()

# Slice plot:
slc = yt.SlicePlot(ds, 'z', 'Ex')
slc.save('Ex_slice.png')
```

Use yt only if output is plotfile format. Otherwise openPMD-viewer is lighter.

## Boundary scrapers

When `warpx_save_particles_at_xlo=True` on a Species, particles crossing that boundary are recorded:

```python
scraper = OpenPMDTimeSeries('./diags/boundary_scrapers_xlo/')
x, ux, w, t = scraper.get_particle(
    var_list=['x', 'ux', 'w', 'timestamp'],
    species='electrons',
    iteration=scraper.iterations[-1],
)

# Energy spectrum on xlo:
E_eV = 0.5 * picmi.constants.m_e * (ux**2 + uy**2 + uz**2) / picmi.constants.q_e
plt.hist(E_eV, bins=200, weights=w, log=True)
```

`timestamp` = simulation time when each particle crossed.

## Reduced diagnostics

ASCII files in `diags/reducedfiles/`:

```python
import pandas as pd
df = pd.read_csv('./diags/reducedfiles/particle_energy.txt',
                 sep=r'\s+', comment='#')
df.plot(x='time(s)', y='total(J)', logy=True)
```

## Performance tips

- **Don't load full 3D fields** into one rank — use openPMD-api `load_chunk` or yt `covering_grid` with coarser level.
- **Use `select=`** for particle reads to avoid materializing all then subsetting.
- **Cache derived quantities** to `.npy` for iterative analysis.

## Where to ask

- Bugs: `https://github.com/BLAST-WarpX/warpx/issues`
- Usage: `https://github.com/BLAST-WarpX/warpx/discussions`
- Slack: invite linked from readthedocs

Always include: WarpX version, OS, Python version, minimal reproducer, full error output.
