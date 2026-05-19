# Running Smilei and troubleshooting

This reference covers launching Smilei (the `smilei` executable, MPI/OpenMP, GPU mode, machine files, checkpoint/restart) and a symptom-table for the most common failure modes. The symptom table is opinionated about what to check first based on the typical priority of causes.

## Table of contents

1. The `smilei` and `smilei_test` executables
2. MPI and OpenMP basics
3. Choosing patch count
4. Running on a cluster
5. GPU mode
6. Checkpoints and restart
7. Machine files
8. Symptom table: simulation won't start
9. Symptom table: simulation crashes during run
10. Symptom table: simulation runs but results are wrong
11. Symptom table: simulation is unexpectedly slow
12. Logging and diagnostics for debugging

---

## 1. The `smilei` and `smilei_test` executables

After `make`, two binaries are produced:

- **`smilei`** — the full executable. Runs the PIC loop end-to-end.
- **`smilei_test`** — the test/debug variant. Parses the namelist, initializes data structures, but does NOT execute the time loop. Catches namelist errors, missing references, illegal patch counts, etc. Always run this before launching the real simulation:

```bash
./smilei_test my_namelist.py
```

The whole test takes seconds and reports parse errors, missing species references, invalid blocks, etc. ~80% of "my namelist won't run" issues are caught here.

To run the actual simulation:

```bash
./smilei my_namelist.py
```

This runs on a single MPI rank by default. For real production runs use MPI launching (next section).

## 2. MPI and OpenMP basics

Smilei uses hybrid MPI + OpenMP:

- **MPI ranks** own disjoint sets of patches (the simulation domain is decomposed by patches, with Hilbert-curve ordering for load balance).
- **OpenMP threads within a rank** parallelize over patches on that rank.

Typical launch:

```bash
mpirun -np 4 ./smilei my_namelist.py                  # 4 MPI ranks, default OpenMP
```

```bash
# 4 MPI ranks, each with 4 OpenMP threads, totaling 16 cores
export OMP_NUM_THREADS=4
mpirun -np 4 ./smilei my_namelist.py
```

Rule of thumb: on a node with N cores, choose `MPI_ranks × OMP_threads = N`. For a single 32-core node, common choices:
- 32 MPI × 1 OMP — high MPI cost, minimal OMP overhead
- 8 MPI × 4 OMP — balanced, often a good starting point
- 4 MPI × 8 OMP — fewer ranks, more OMP per rank, useful for memory-bound

For multi-node runs, MPI typically maps one or two ranks per socket and OpenMP fills the cores within that socket.

## 3. Choosing patch count

Patches are the fundamental decomposition unit. Two constraints:

1. **`number_of_patches` must be a power of 2 in each direction.**
2. **There must be at least one patch per MPI rank, ideally several.**

Hard guidelines:

- For N MPI ranks, total patches should be `≥ 4N`, preferably `8–16N`.
- For very large runs (1000+ ranks), `64–256N` patches helps load balance.

For a 2D simulation on 64 MPI ranks: `number_of_patches = [32, 16]` gives 512 patches, 8 per rank — good.

For 3D: `[16, 8, 4]` gives 512 patches.

If you violate the power-of-2 constraint:

```python
number_of_patches = [12, 8]              # 12 is not a power of 2 → ERROR
```

The simulation refuses to start with an opaque error message. Always check.

## 4. Running on a cluster

Smilei comes with helper scripts in `scripts/`. The typical pattern with SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=smilei_run
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=smilei_%j.out

module load gcc openmpi hdf5 python

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ./smilei my_namelist.py
```

For 8 nodes × 4 MPI/node × 4 OMP threads = 128 cores. Adjust patches accordingly (`8 × 4 = 32` MPI ranks, so `≥256` patches).

**Always start with a short test run** (`simulation_time = 5.` or so) before committing 24 hours of compute time. Verify scaling and output writing work as expected.

## 5. GPU mode

Smilei v5.1 supports GPU offload via OpenMP target or OpenACC, depending on compilation. To use:

1. **Compile with GPU support**. Use a GPU-aware machine file (NVIDIA: `nvhpc` toolchain; AMD: `rocm`).
2. **Enable in the namelist**:

```python
Main(
    # ...
    gpu_computing = True,
)
```

3. **Launch with one MPI rank per GPU**:

```bash
mpirun -np 4 ./smilei my_namelist.py     # for 4 GPUs
```

OpenMP threads inside the rank are used for CPU-side tasks (I/O, control flow). The PIC loop runs on the GPU.

GPU performance gains depend strongly on configuration. Rule of thumb: GPU wins for large simulations with high PPC and dense particle/cell distribution; CPU may win for small simulations or sparse setups (where GPU overhead dominates).

GPU support varies by feature. As of v5.1:
- Cartesian (1D, 2D, 3D) geometry: supported
- AMcylindrical: increasingly supported (check release notes)
- Collisions: partial support
- Some ionization paths: partial support

When in doubt, run CPU-only and check the docs for the specific feature/version.

## 6. Checkpoints and restart

For long simulations or wall-time-limited cluster jobs, checkpoint the state periodically and restart from the checkpoint.

```python
Checkpoints(
    dump_step = 10000,                         # dump every 10000 steps
    dump_minutes = 60.0,                       # alternative: dump every N minutes
    exit_after_dump = True,                    # exit Smilei after writing — controller resubmits
    keep_n_dumps = 2,                          # keep last 2 dumps on disk
    restart_dir = None,                        # for restart, set to previous job's dir
)
```

On the first run, set `restart_dir = None`. The job runs, dumps state, exits.

On restart:

```python
Checkpoints(
    dump_step = 10000,
    dump_minutes = 60.0,
    exit_after_dump = True,
    keep_n_dumps = 2,
    restart_dir = "/path/to/previous/dump",    # ← previous job's checkpoint dir
)
```

Smilei reads the dump and continues from where it left off. Output files are appended; diagnostics persist across checkpoint boundaries.

For chained jobs, use a small wrapper script (bash, SLURM array) that increments the restart directory on each successive submission.

## 7. Machine files

Smilei's makefile system uses "machine files" — small shell files in `scripts/compile_tools/machine/` that set compiler flags for a specific cluster or CPU family. Examples shipped: `intel-skylake`, `amd-rome`, `ibm-power9`, `nvidia-a100`.

To use:

```bash
make machine=intel-skylake
```

This sets compiler optimization flags appropriate for the architecture. If your cluster isn't covered, copy an existing file and edit:

```bash
cp scripts/compile_tools/machine/intel-skylake scripts/compile_tools/machine/my-cluster
# edit my-cluster to match your compiler / SIMD width
make machine=my-cluster
```

Performance with the right machine file can be 2–4× faster than the generic build for vectorizable physics loops.

## 8. Symptom table: simulation won't start

| Symptom | Most likely cause | Check |
|---|---|---|
| "ERROR_NAMELIST: invalid keyword" | Misspelled namelist parameter | Spell check; consult docs for parameter name |
| "Patch number not a power of 2" | `number_of_patches` violates rule 2 | Set `[2^a, 2^b]` |
| "Species 'X' not found" | Referenced species not defined | Check `species1=...` in Collisions; check `ionization_electrons=...` |
| Hangs at startup, no output | MPI deadlock during initialization | Check `mpirun` config, run with `-np 1` first |
| "HDF5: file already exists" | Previous run wasn't cleaned | `rm -rf` the output dir, or change output location |
| `smilei: command not found` | Binary not built or PATH issue | Re-run `make`; verify `./smilei` exists |
| Hangs reading namelist | Infinite loop or external file read in Python preamble | Check namelist Python code for non-terminating constructs |
| ImportError on Python module | Python module used in namelist not installed | `pip install` the missing module |
| "openMP threads = 1" but expecting more | `OMP_NUM_THREADS` not set or overridden | `echo $OMP_NUM_THREADS` before launch |

## 9. Symptom table: simulation crashes during run

| Symptom | Most likely cause | Check |
|---|---|---|
| Segfault during particle push | NaN momentum or position | Profile callable returning NaN; absurd initial momentum |
| Segfault during field solve | Memory exhausted (too many particles) | Decrease PPC, decrease box size, or increase memory |
| "MPI_Abort: rank X" | Per-rank crash; check stdout/stderr of that rank | Often hardware issue or specific particle pathology |
| Energy diverges exponentially | Numerical instability — under-resolved | Decrease timestep, finer grid, check CFL |
| Charge conservation broken | Maxwell solver / current deposition inconsistency | Verify standard solver; reduce timestep |
| Crash after restart | Checkpoint corrupt or version mismatch | Restart from earlier dump; ensure same Smilei version |
| Memory grows unboundedly | Leak from ionization or particle injection | Cap `maximum_charge_state`; check ParticleInjector rates |
| Specific patch slow → entire rank stalls | Load imbalance | Enable `LoadBalancing`; reduce patch size |

## 10. Symptom table: simulation runs but results are wrong

| Symptom | Most likely cause | Check |
|---|---|---|
| Collision rates much too low/high | `reference_angular_frequency_SI` not set | Hard rule 1 |
| Ionization rate zero everywhere | Wrong model for the regime (ADK at high `γ_K`) | Compute `γ_K` at peak; switch to PPT or `from_rate` |
| Density not where you expect | Profile callable in normalized units, user thought SI | Print profile inline; convert SI distances to `L_r` |
| Plasma "explodes" at `t = 0` | Quasi-neutrality violated; missing Poisson init | Initialize with `momentum_initialization = "cold"` first, add temperature gradually |
| Laser doesn't enter the box | Wrong `box_side`; `periodic` BC on launch face | Match `EM_boundary_conditions` and `Laser.box_side` |
| Laser intensity off by 2× | Linear vs circular polarization confusion | Check `ellipticity` and `a₀` consistency |
| Spurious EM pulse at `t = 0` | Dense relativistic beam without field init | Set `relativistic_field_initialization = True` on the beam species |
| Numerical Cherenkov radiation | Relativistic particles + standard Yee | Use `Lehe` solver or B-TIS3 interpolation |
| No ionization in envelope mode | `tunnel` model with envelope laser | Use `tunnel_envelope_averaged` |
| Charge state never advances past Ar+1 | `maximum_charge_state = 1` (default may be 1) | Set to expected peak charge state |
| Unphysical temperature spike on collisions | Too few PPC | Bump PPC to ≥ 32 |
| Energy not conserved | Under-resolution; too coarse `Δt`; wrong solver | Halve `cell_length`; check CFL margin |

## 11. Symptom table: simulation is unexpectedly slow

| Symptom | Most likely cause | Check |
|---|---|---|
| Per-step time grows with simulation time | Particle count growing (injection, ionization) | Expected; budget compute time accordingly |
| Some MPI ranks idle | Load imbalance | Enable `LoadBalancing`; check `DiagPerformances` |
| Slow file I/O | Many small `DiagFields` writes | Combine into one block; increase `every`; use `subgrid` |
| Slow on cluster but fast on workstation | Wrong machine file | Re-build with appropriate `machine=` |
| GPU mode slow | GPU overhead with small simulation | Try CPU mode; or scale up |
| Hangs at periodic intervals | Checkpoint writing | Expected; reduce `dump_step` frequency or use async I/O |
| Vectorization not engaged | Compile-time options missing | Use the right machine file with SIMD flags |
| Collision module dominates time | Many collision blocks, high PPC | Reduce non-essential collision channels; check `every=` |

## 12. Logging and diagnostics for debugging

### `print_every`

```python
Main(
    print_every = 100,                         # print step status every 100 steps
)
```

Smilei writes per-step diagnostic info to stdout. Includes: time, timestep number, kinetic energy, EM energy, total particles per species. Watching this stream is the fastest way to catch divergence.

### Verbose mode

Run with `--verbose` to get extra info from initialization and load balancing:

```bash
./smilei --verbose my_namelist.py
```

### Check the namelist's effective state

The `smilei.py` file in the output directory is the namelist *as Smilei interpreted it* — including resolved profiles, expanded defaults, and computed derived quantities. Compare to your input namelist to catch discrepancies.

### `DiagPerformances`

Always-on for any serious cluster run. Diagnoses load imbalance via per-patch timing. See `diagnostics.md` §10.

### Standalone profile testing

For complex density / temperature profiles, test in standalone Python before plugging in:

```python
# Save this to a file, run it, plot the output
import numpy as np
import matplotlib.pyplot as plt

def n_profile(x, y):
    if x < 5.: return 0.0
    if x < 15.: return 0.1
    return 0.

xs = np.linspace(0, 30, 300)
ys = np.array([n_profile(x, 10.) for x in xs])
plt.plot(xs, ys); plt.show()
```

If the profile looks wrong here, it's wrong in the simulation too. This 60-second standalone test catches half of all "where is my plasma?" mysteries.

### When to ask for help

The Smilei community is responsive on the [GitHub Discussions](https://github.com/SmileiPIC/Smilei/discussions) page. Before posting:

1. Run `smilei_test` and include the output if there's a parse error.
2. State the Smilei version (in `Smilei/README.md` or `git log`).
3. Share the full namelist (or a minimal reproducer if possible).
4. Describe the expected vs. observed behavior precisely (numbers, not "doesn't work").
5. Note your machine and launch configuration (number of ranks/threads, machine file used).

Most simulation bugs reduce to one of the issues in tables 8–11 above. The community has seen them all.
