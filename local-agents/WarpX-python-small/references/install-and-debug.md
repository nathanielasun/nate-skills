# Install, run, and debug

## Quick install

```bash
# Most common — pip from PyPI:
pip install pywarpx[3D]                # 3D only
pip install pywarpx[1D,2D,3D,RZ]        # all dimensions

# Or conda-forge:
conda install -c conda-forge warpx
```

Verify:
```bash
python3 -c "from pywarpx import picmi; print(picmi.__file__)"
```

## Building from source

For GPU, site-specific MPI, or development on the C++ side:

```bash
git clone --depth 1 https://github.com/BLAST-WarpX/warpx.git
cd warpx

# Activate Python environment first!

# Build:
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
| `WarpX_PYTHON` | OFF | Build Python bindings |
| `WarpX_DIMS` | "3" | Dimensions; "1;2;3;RZ" for all |
| `WarpX_COMPUTE` | OMP | `OMP`, `CUDA`, `HIP`, `SYCL`, `NOACC` |
| `WarpX_FFT` | OFF | For PSATD solver |
| `WarpX_QED` | ON | QED module (turn off if not needed) |
| `WarpX_PYTHON_IPO` | ON | LTO for Python lib (slow; OFF for dev iteration) |
| `AMReX_CUDA_ARCH` | autodetect | CUDA capability, e.g., `80`, `90` |

For low-temp plasma dev (no PSATD, no QED, with Python, all dims):

```bash
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;3;RZ" \
      -DWarpX_FFT=OFF -DWarpX_QED=OFF
cmake --build build -j 8 --target pip_install
```

For GPU (A100):
```bash
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="3" \
      -DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=80
```

Faster rebuild during iteration:
```bash
cmake --build build -j 8 --target pip_install_nodeps    # skip dep resolution
```

## pywarpx module organization

| Module | Use |
|---|---|
| `pywarpx.picmi` | PICMI classes (the public interface) |
| `pywarpx.callbacks` | `installcallback`, `callfromafterstep`, etc. |
| `pywarpx.fields` | Field wrappers: `ExWrapper`, `RhoFPWrapper`, etc. |
| `pywarpx.particle_containers` | `ParticleContainerWrapper` |
| `pywarpx.warpx` | Native (non-PICMI) interface |
| `pywarpx.libwarpx` | C++ library; `libwarpx.amr` is pyAMReX |

Typical scripts use the first four. The rest are for advanced uses.

## Running a script

### Serial

```bash
python3 my_script.py
```

### MPI

```bash
mpiexec -n 4 python3 my_script.py
```

### MPI initialization order matters

If using `mpi4py`, import BEFORE pywarpx:

```python
# Correct order:
from mpi4py import MPI                   # initializes MPI
from pywarpx import picmi                 # uses already-initialized MPI

comm = MPI.COMM_WORLD
if comm.rank == 0:
    print(f"Running with {comm.size} ranks")
```

### Hybrid MPI + OpenMP

```bash
export OMP_NUM_THREADS=4
mpiexec -n 8 python3 my_script.py        # 32 cores total
```

### GPU runs

```bash
# 1 GPU per MPI rank (typical)
mpiexec -n 4 python3 my_script.py
```

Set `warpx_amrex_use_gpu_aware_mpi=True` on the Simulation only if your MPI is CUDA-aware.

## HPC submission (Slurm)

```bash
#!/bin/bash
#SBATCH --job-name=my_picmi_run
#SBATCH --time=2:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4

module load python/3.11 openmpi/4.1 hdf5
source /path/to/venv/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun python3 my_script.py
```

For GPU:
```bash
#SBATCH --gpus-per-task=1
#SBATCH --ntasks-per-node=4
export OMP_NUM_THREADS=1
srun python3 my_script.py
```

Site-specific machine profiles in `Tools/machines/<system>/` in the WarpX repo:
- `Tools/machines/perlmutter-nersc/`
- `Tools/machines/frontier-olcf/`
- `Tools/machines/lumi-csc/`

Each has `submit.sh` templates and required module loads.

### Output to scratch, not $HOME

```python
field_diag = picmi.FieldDiagnostic(...,
    write_dir='/scratch/$USER/my_run/diags',
)
```

Long runs fill $HOME quickly.

## Common errors and fixes

### `ImportError: libwarpx_pyDIM.so cannot be opened`

pywarpx not properly installed, or wrong dim variant.

**Fix**: `pip show pywarpx` to verify install. If from source: `cmake --build build --target pip_install` must have completed. Install correct extras: `pip install pywarpx[3D]`.

### `AttributeError: module 'pywarpx.picmi' has no attribute 'XxxYyy'`

PICMI class renamed or removed in your version.

**Fix**: check the PICMI reference for your version: `https://warpx.readthedocs.io/en/<version>/usage/python.html`. Class names like `ElectrostaticFieldDiagnostic` have been deprecated in favor of `FieldDiagnostic`.

### `Could not parse expression '...'`

WarpX runtime parser can't make sense of your string.

**Fix**: variables `x`, `y`, `z`, `t` are available; constant `pi`; functions `sin`, `cos`, `exp`, `sqrt`, `log`, etc.; ternary as `if(cond, t, f)` not Python syntax; custom params via kwargs:

```python
dist = picmi.AnalyticDistribution(
    density_expression='n0 * exp(-(x/r0)**2)',
    n0=1e18,                  # passed as kwarg
    r0=1e-3,
    ...
)
```

### `MPI_Init was already called` / hangs at startup

mpi4py and pywarpx conflict over MPI.

**Fix**: import `mpi4py.MPI` first (it initializes MPI), then pywarpx.

### `Cross-section file not found`

Path wrong or `warpx-data` not cloned.

**Fix**:
```bash
git clone https://github.com/BLAST-WarpX/warpx-data.git
```
Use absolute paths in PICMI scripts or run from the directory containing `warpx-data/`.

### Simulation hangs immediately

MPI ranks diverge — usually a callback making MPI calls inconsistently across ranks, or differing inputs per rank.

**Fix**: run with `-n 1` first to verify correctness. Then add MPI. Print rank-numbered checkpoints to find divergence:
```python
print(f"rank {MPI.COMM_WORLD.rank} at point X", flush=True)
```

### GPU not detected

```python
from pywarpx import libwarpx
print(libwarpx.Config.have_gpu)             # should be True
```

**Fix**: rebuild with `-DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=<sm>`. Use the right arch: `80` for A100, `90` for H100.

### Memory exhaustion in discharge runs

Ionization creates new particles each step; without resampling, count grows unbounded.

**Fix**: enable resampling on each charged species:
```python
electrons = picmi.Species(...,
    warpx_do_resampling=True,
    warpx_resampling_trigger_intervals='1000',
    warpx_resampling_trigger_max_avg_ppc=200,
)
```

### Reproducibility issues

Two identical runs give different results.

**Fix**:
```python
sim = picmi.Simulation(..., warpx_random_seed=42,
                       warpx_serialize_initial_conditions=True)
layout = picmi.PseudoRandomLayout(n_macroparticles_per_cell=16, seed=42)
```

GPU runs can have inherent non-determinism from atomic ops; expect small differences.

## Debugging

### Verbose flag

```python
sim = picmi.Simulation(..., verbose=2)         # 0=silent, 1=normal, 2=debug
```

Level 2 prints solver iterations, particle counts, MPI communication.

### Inspect generated input file

Even running interactively, generate the input file once to check `warpx_*` mappings:

```python
sim.write_input_file('inputs_debug')
sim.step()
```

Then `cat inputs_debug` — verify every `warpx_*` argument shows as the expected C++ parameter. Many bugs become obvious here.

### Catching NaN/Inf early

```python
sim = picmi.Simulation(...,
    warpx_amrex_fpe_trap_invalid=True,
    warpx_amrex_fpe_trap_zero=True,
    warpx_amrex_fpe_trap_overflow=True,
)
```

Traps FP exceptions immediately instead of propagating NaNs. Enable for debugging, disable for production (small cost).

### pdb in MPI scripts

```python
import pdb
if MPI.COMM_WORLD.rank == 0:
    pdb.set_trace()         # only rank 0 enters debugger
```

Other ranks continue — may deadlock if rank-0 code needs sync. For most cases, print debugging is easier:

```python
if MPI.COMM_WORLD.rank == 0:
    print(f"Step {step}: max E = {max_E}", flush=True)
```

### Reproducing a failure

If a sim fails at 10000 steps on 128 ranks:

1. Run 100 steps on 128 ranks. Still fails? → issue is initial.
2. Run 100 steps on 1 rank. Still fails? → issue is non-MPI.
3. Run with `verbose=2`. What's the last printed thing?
4. Reduce dims (3D→2D), grid size (256→32), particles for fast iteration.

Fix on a small case, then scale back up.

### Reading C++ stack traces

When WarpX crashes inside C++:

| Error | Likely cause |
|---|---|
| "Particle outside box" | Callback moved particles without `Redistribute()` |
| "PoissonSolver did not converge" | Tighten `required_precision`, increase `maximum_iterations`, or check BCs |
| "ParmParse: parameter not understood" | Check `inputs_debug` — a `warpx_*` arg didn't map to valid C++ param |
| Segfault in MCC | Check cross-section file paths and product species registration |

Search the WarpX issue tracker for the exact error — it's usually been seen before.

## Where to ask for help

- **Bugs**: `https://github.com/BLAST-WarpX/warpx/issues`
- **Usage questions**: `https://github.com/BLAST-WarpX/warpx/discussions`
- **Real-time chat**: WarpX Slack (invite linked from readthedocs)

Include: WarpX version, install method, OS, Python version, minimal reproducer, full error output. Smaller reproducers get faster responses.

## Cheat sheet

```bash
# Run serially
python3 my_script.py

# Run with MPI
mpiexec -n 4 python3 my_script.py

# Inspect generated input file
python3 -c "exec(open('my_script.py').read().replace('sim.step()', 'sim.write_input_file(\"inputs\")'))"
cat inputs

# Check pywarpx config
python3 -c "from pywarpx import libwarpx; print('GPU:', libwarpx.Config.have_gpu)"

# Validate PICMI API exists
python3 -c "from pywarpx import picmi; print(dir(picmi))"
```
