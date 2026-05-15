# Installation, build, and running

This reference covers getting `pywarpx` installed, building from source when needed, running PICMI scripts with MPI, HPC submission patterns, and diagnosing the most common Python-side errors. Read this for setup, deployment, or when something won't run.

## Table of contents

1. [Installation options](#installation-options)
2. [Building from source for development](#building-from-source-for-development)
3. [Module organization in pywarpx](#module-organization-in-pywarpx)
4. [Running a PICMI script](#running-a-picmi-script)
5. [MPI and parallelism](#mpi-and-parallelism)
6. [HPC submission patterns](#hpc-submission-patterns)
7. [Common errors and what they mean](#common-errors-and-what-they-mean)
8. [Debugging](#debugging)

## Installation options

### Option 1: pip (recommended for most users)

```bash
pip install pywarpx           # 3D by default
# Or, specify dimensionality:
pip install pywarpx[3D]
pip install pywarpx[2D]
pip install pywarpx[1D]
pip install pywarpx[RZ]
# Or install all dims:
pip install pywarpx[1D,2D,3D,RZ]
```

The PyPI wheels are CPU-only and serial-MPI-ready (with `mpi4py` if installed). For GPU or for HPC-specific MPI builds, install from source.

After install, verify:
```bash
python3 -c "from pywarpx import picmi; print(picmi.__file__)"
```

### Option 2: conda-forge

```bash
conda install -c conda-forge warpx
# Or for a specific dim:
conda install -c conda-forge py-warpx
```

conda-forge packages include MPI variants. To get a GPU build via conda, look for `warpx-cuda*` packages or build from source.

### Option 3: from source (developer / HPC)

For HPC systems with site-specific MPI/CUDA, GPU-enabled builds, or development on the C++ side: build from source. See the next section.

### Option 4: HPC system modules

Many HPC sites have WarpX as a module:
```bash
module load warpx                # syntax varies
which warpx.3d                    # find the binary
```

When using a system module, you still need a Python environment with `pywarpx` for PICMI scripts (or use `write_input_file()` and run the binary directly).

## Building from source for development

For Python-side development, you need WarpX built with `WarpX_PYTHON=ON`:

```bash
git clone --depth 1 https://github.com/BLAST-WarpX/warpx.git
cd warpx
# Activate your Python environment first

# Build for 3D Cartesian with Python:
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="3"

# For all dimensions (longer build):
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;3;RZ"

# For GPU (CUDA):
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="3" \
      -DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=90

# For faster iteration (skip link-time optimization):
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_PYTHON_IPO=OFF -DWarpX_DIMS="3"

cmake --build build -j 8 --target pip_install
```

The `pip_install` target builds the C++ libraries, then `pip install`s the resulting wheel into the active Python environment.

For iterative development:
```bash
cmake --build build -j 8 --target pip_install_nodeps
```

`pip_install_nodeps` skips dependency resolution — much faster on rebuilds.

After building, verify:
```bash
python3 -c "from pywarpx import picmi, warpx; print('OK')"
```

If you get `ImportError: libwarpx_pyXXX.so cannot be opened`: the build succeeded but the import path is wrong. Check that the right Python environment is active and that `pip show pywarpx` finds the installation.

### Build options summary

| Option | Default | Meaning |
|---|---|---|
| `WarpX_PYTHON` | OFF | Build Python bindings |
| `WarpX_DIMS` | "3" | Which dimensions to build; semicolon-separated: "1;2;3;RZ" |
| `WarpX_COMPUTE` | OMP | Compute backend: `OMP`, `CUDA`, `HIP`, `SYCL`, `NOACC` |
| `WarpX_MPI` | ON | MPI support |
| `WarpX_OPENPMD` | ON | openPMD output support |
| `WarpX_QED` | ON | QED module |
| `WarpX_FFT` | OFF | FFT support (for PSATD solver) |
| `WarpX_PYTHON_IPO` | ON | Link-time optimization for Python lib (slow build) |

For low-temp plasma work, the typical build:
```bash
cmake -S . -B build -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;3;RZ" \
      -DWarpX_FFT=OFF -DWarpX_QED=OFF
```

QED and FFT are usually not needed for low-temp plasma and slow the build.

## Module organization in pywarpx

The `pywarpx` Python package is organized as:

| Submodule | Purpose |
|---|---|
| `pywarpx.picmi` | PICMI standard classes + WarpX extensions (the public interface) |
| `pywarpx.warpx` | The native (non-PICMI) interface — prefix-instance machinery, lower-level |
| `pywarpx.libwarpx` | The pybind11-wrapped C++ library; `libwarpx.amr` is pyAMReX |
| `pywarpx.particle_containers` | `ParticleContainerWrapper` and friends |
| `pywarpx.callbacks` | Hook installation: `installcallback`, `callfromafterstep`, etc. |
| `pywarpx.fields` | Field wrappers: `ExWrapper`, `RhoFPWrapper`, etc. |
| `pywarpx._libwarpx` | Internal: the libwarpx import shim |
| `pywarpx.Bucket` | Internal: prefix bucket class for native interface |
| `pywarpx.WarpX` | Internal: the WarpX prefix container |

For typical scripts, only `pywarpx.picmi`, `pywarpx.callbacks`, `pywarpx.fields`, and `pywarpx.particle_containers` are used directly. The rest are internal but occasionally helpful for debugging or unusual needs.

## Running a PICMI script

### Serial

```bash
python3 my_script.py
```

That's it. Output goes to wherever your diagnostics specified (`./diags/` by default).

### Parallel (MPI)

```bash
mpiexec -n 4 python3 my_script.py
# Or:
mpirun -n 4 python3 my_script.py
```

The PICMI script must NOT initialize MPI before `pywarpx` does. If you use `mpi4py`, ensure that `mpi4py` is initialized before `pywarpx` (typically by importing `mpi4py.MPI` before `pywarpx`):

```python
# Correct order:
from mpi4py import MPI                    # initializes MPI
from pywarpx import picmi                  # uses the already-initialized MPI

# All ranks should reach this point
comm = MPI.COMM_WORLD
```

If only rank 0 writes plots, gate output:
```python
if MPI.COMM_WORLD.rank == 0:
    plt.savefig('result.png')
```

### Specifying dimensions

The choice of `Grid` class in your script selects the dimensionality. If you have `pywarpx[2D]` installed and try to use `Cartesian3DGrid`, you'll get an import error or runtime crash. Install the right dim wheel for your script.

### Output and working directory

The diagnostic directory is set per-diagnostic via `write_dir`. By default this is the current working directory. For HPC jobs, set `write_dir` to an absolute path within scratch storage to avoid filling shared $HOME directories.

## MPI and parallelism

### Domain decomposition

WarpX (via AMReX) automatically decomposes the simulation domain across MPI ranks. The decomposition respects:
- `warpx_max_grid_size` (maximum block size per dimension)
- `warpx_blocking_factor` (block sizes are multiples of this)
- `warpx_numprocs` if specified (explicit rank-per-axis assignment)

For 16 MPI ranks in a 3D simulation:
- AMReX default: decompose into approximately 16 boxes total, in some balanced shape.
- With `warpx_numprocs=[2, 2, 4]`: enforce 2x2x4 = 16 ranks.

### Hybrid MPI+OpenMP

Inside each MPI rank, OpenMP parallelism is used for some loops if WarpX was built with `WarpX_COMPUTE=OMP` (the default). Set `OMP_NUM_THREADS` to control:

```bash
export OMP_NUM_THREADS=4
mpiexec -n 8 python3 my_script.py    # 8 MPI ranks, 4 threads each = 32 cores
```

On HPC, `OMP_NUM_THREADS` is usually set by the resource manager (Slurm `--cpus-per-task`); check your site docs.

### GPU runs

With CUDA/HIP builds:

```bash
mpiexec -n 4 python3 my_script.py    # 1 GPU per MPI rank (typical)
```

WarpX expects one GPU per MPI rank by default. If you have N GPUs on a node, use N MPI ranks for that node.

For `warpx_amrex_use_gpu_aware_mpi=True`: requires CUDA-aware MPI (OpenMPI, MPICH built against UCX, MVAPICH, etc.). On HPC sites this is usually pre-configured; on workstations you may need to build MPI with the right flags.

### Load balancing

For non-uniform plasma density or moving features:

```python
sim = picmi.Simulation(...,
    warpx_load_balance_intervals='100',           # rebalance every 100 steps
    warpx_load_balance_costs_update='heuristic',  # or 'timers', 'gpuclock'
    warpx_load_balance_with_sfc=False,             # space-filling-curve (vs knapsack)
)
```

Load balancing redistributes boxes among ranks based on accumulated cost estimates. The `'timers'` strategy uses actual wall-clock measurements; `'heuristic'` uses particle and cell counts. `'gpuclock'` is for GPU runs.

For mostly-uniform low-temp plasma, default settings are fine.

## HPC submission patterns

### General Slurm submission

```bash
#!/bin/bash
#SBATCH --job-name=my_picmi_run
#SBATCH --output=logs/%j.out
#SBATCH --time=2:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0

module load python/3.11
module load openmpi/4.1
module load hdf5

source /path/to/venv/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun python3 my_script.py
```

For 4 nodes × 8 tasks × 4 threads = 128 cores total, with 32 MPI ranks.

### GPU submission

```bash
#SBATCH --gpus-per-task=1
#SBATCH --ntasks-per-node=4

export OMP_NUM_THREADS=1
srun python3 my_script.py
```

(Exact directives depend on the system; check site docs.)

### Machine profiles

WarpX provides machine-specific build/run profiles in `Tools/machines/<system>/`. Look there for examples:
- `Tools/machines/perlmutter-nersc/` — NERSC Perlmutter
- `Tools/machines/frontier-olcf/` — OLCF Frontier
- `Tools/machines/lumi-csc/` — LUMI
- `Tools/machines/leonardo-cineca/` — Leonardo

Each contains a `submit.sh` template, build instructions, and `*_warpx.profile.example` shell file with the right module loads.

### Persistent output

For long runs, write diagnostics to scratch space, not $HOME:

```python
field_diag = picmi.FieldDiagnostic(
    ...,
    write_dir='/scratch/$USER/my_run/diags',     # expand at script time, not runtime
)
```

After the run, archive needed outputs to long-term storage (HPSS, S3, etc.). Reduced diagnostics and final field snapshots are usually the priority.

## Common errors and what they mean

### `ImportError: libwarpx_pyDIM.so cannot be opened`

`pywarpx` is not properly installed, or the dim variant you're using isn't built.

**Fix**: verify with `pip show pywarpx`. If you built from source, ensure `cmake --build build --target pip_install` completed. If you installed via pip, ensure you used the right extras: `pip install pywarpx[3D]`.

### `AttributeError: module 'pywarpx.picmi' has no attribute 'XxxYyy'`

PICMI class name has changed or was removed in your version.

**Fix**: check the PICMI reference for your specific WarpX version: `https://warpx.readthedocs.io/en/<version>/usage/python.html`. Names like `ElectrostaticFieldDiagnostic` have been deprecated in favor of `FieldDiagnostic`.

### `Could not parse expression '...'`

WarpX's runtime parser can't make sense of your string expression (for `density_expression`, `Ex_expression`, applied potential, etc.).

**Fix**: check the syntax. Variables `x`, `y`, `z`, `t` are available; constants `pi` is available; functions `sin`, `cos`, `exp`, `sqrt`, etc. Custom parameters must be passed as kwargs:
```python
dist = picmi.AnalyticDistribution(
    density_expression='n0 * exp(-(x/r0)**2)',
    n0=1e18,                  # passed as kwarg
    r0=1e-3,                   # passed as kwarg
    ...
)
```

If you used `**` for exponentiation, that should work (WarpX parser supports `**`). For ternary, use `if(cond, t, f)` not Python syntax.

### `MPI_Init was already called` or hangs at startup

`mpi4py` and `pywarpx` are fighting over MPI initialization.

**Fix**: import `mpi4py.MPI` first (which initializes MPI), then import `pywarpx`. Or set:
```python
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI
from pywarpx import picmi   # pywarpx initializes MPI
```

### `Cross-section file not found`

Path to MCC cross-section file is wrong, or `warpx-data` not cloned.

**Fix**: clone `BLAST-WarpX/warpx-data` alongside your script. Use absolute paths or paths relative to where you launch the job, not relative to the script. For HPC, copy `warpx-data/` to scratch and reference there.

### Simulation hangs immediately

MPI processes are out of sync, often because of a callback that does MPI calls inconsistently across ranks, or because of differing input on different ranks.

**Fix**: run with a single rank first (`python3 script.py` without `mpiexec`) to verify correctness. Add `print(f"rank {MPI.COMM_WORLD.rank} at checkpoint X")` at key points to find where ranks diverge.

### GPU not detected

Build was CPU-only, or CUDA wasn't found at build time, or wrong CUDA arch.

**Fix**: rebuild with `-DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH=<sm>`. The arch is the compute capability of your GPU, e.g., `80` for A100, `90` for H100. Verify with:
```bash
python3 -c "from pywarpx import libwarpx; print(libwarpx.Config.have_gpu)"
```

### Reproducibility issues

Two runs with the same script give different results.

**Fix**: set `warpx_random_seed` on Simulation and `seed` on PseudoRandomLayout. For full reproducibility, also set `warpx_serialize_initial_conditions=True` (only for testing — slows initialization). Note that GPU runs can have inherent nondeterminism from atomic operations regardless of seed.

### Memory exhaustion in long discharge runs

Ionization in MCC creates new particles each step; without resampling, particle count grows unboundedly.

**Fix**: enable resampling on charged species:
```python
electrons = picmi.Species(...,
    warpx_do_resampling=True,
    warpx_resampling_trigger_intervals='1000',
    warpx_resampling_trigger_max_avg_ppc=200,
)
```

Apply to ions as well. Resampling has a small noise cost but is essential for steady-state discharge simulations.

## Debugging

### Verbose flag

```python
sim = picmi.Simulation(..., verbose=2)
```

Levels: 0 (silent), 1 (normal), 2 (debugging). Higher levels print per-step information about solver iterations, particle counts, MPI communication patterns.

### Check the generated input file

Even when running interactively, call `write_input_file()` first to see what WarpX is actually configured to do:

```python
# At the bottom of your script, before sim.step():
sim.write_input_file('inputs_debug')   # writes the C++-style input
sim.step()
```

Then `cat inputs_debug` and check that every `warpx_*` argument you set in PICMI shows up as the corresponding C++ parameter. Mismatches are the source of many "I told it X but it's doing Y" bugs.

### pdb in MPI scripts

```python
import pdb

if MPI.COMM_WORLD.rank == 0:
    pdb.set_trace()    # only rank 0 enters debugger
```

Other ranks continue, which may deadlock if the rank-0 code needs to synchronize. For genuine interactive MPI debugging, use specialized tools: `ddt`, `totalview`, or `mpi4py-fft`'s parallel pdb.

For most use cases, `print` debugging is easier — gated to rank 0 to keep output clean:

```python
if MPI.COMM_WORLD.rank == 0:
    print(f"Step {step}: max E = {max_E}")
```

### Diagnosing slow runs

```bash
# Lightweight Python-side profile:
python3 -c "import cProfile; cProfile.run('exec(open(\"my_script.py\").read())', sort='cumtime')"
```

This profiles the entire script including the C++ work; the C++ work shows as a single big call to `step`. To profile *Python overhead specifically* (callbacks, custom code), use `line_profiler`:

```python
from line_profiler import LineProfiler

@profile      # decorator from line_profiler
def my_callback():
    # ...
```

For C++-side profiling, use the WarpX-internal timers (run with `verbose=2` or set `warpx_print_timers=True`).

### Catching NaN/Inf early

```python
sim = picmi.Simulation(...,
    warpx_amrex_fpe_trap_invalid=True,
    warpx_amrex_fpe_trap_zero=True,
    warpx_amrex_fpe_trap_overflow=True,
)
```

These trap floating-point exceptions immediately (vs propagating NaNs through the simulation), making bugs much easier to localize. Enable for debugging; can be slow for production.

### Reproducing a failure on a smaller scale

If a simulation fails after 10,000 steps on 128 ranks, isolate by:
1. Running 100 steps on 128 ranks — does it still fail? If yes, the issue is initial.
2. Running 100 steps on 1 rank — does it still fail? If yes, the issue is non-MPI.
3. Running 100 steps with `verbose=2` — what's the last thing it prints?
4. Reducing dimensions (3D → 2D), grid size (256 → 32), or particle count for fast iteration.

Reproducing a bug on a single core with a small problem is often the fastest path to a fix. After fixing, scale back up.

### Reading the C++ stack trace

When WarpX crashes from inside the C++ code, you'll get a stack trace like:
```
amrex::abort: ...
SIGSEGV
==> stack: ...
```

The trace may be in C++ but tells you which subsystem failed. Common patterns:
- "Particle outside box" → position update went wrong; check that `Redistribute()` is called after any callback that moves particles.
- "PoissonSolver did not converge" → either tighten `required_precision`, increase `maximum_iterations`, or check that boundary conditions are well-posed.
- "ParmParse: parameter not understood" → check `inputs_debug` — a `warpx_*` argument didn't map to a valid C++ parameter; check spelling.

When in doubt, search the WarpX issue tracker for the exact error message: someone has likely encountered it before.

### Asking for help

For genuine head-scratchers:
- [WarpX Discussions](https://github.com/BLAST-WarpX/warpx/discussions) for usage questions
- [WarpX Issues](https://github.com/BLAST-WarpX/warpx/issues) for bugs
- Include: WarpX version, install method, OS, Python version, minimal reproducer, full error output.

When asking, write the *smallest* script that reproduces the issue. Stripping irrelevant configuration is often what reveals the bug.
