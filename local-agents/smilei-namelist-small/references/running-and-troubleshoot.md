# Running Smilei and troubleshooting

Launching Smilei (`smilei`, `smilei_test`, MPI/OpenMP, GPU, checkpoints) and symptom tables for common failure modes.

## Table of contents

1. The `smilei` and `smilei_test` executables
2. MPI and OpenMP
3. Patch count rules
4. Cluster run template (SLURM)
5. GPU mode
6. Checkpoints
7. Symptom table — won't start
8. Symptom table — crashes during run
9. Symptom table — runs but results are wrong
10. Symptom table — unexpectedly slow
11. Debugging utilities

---

## 1. The `smilei` and `smilei_test` executables

After `make`:

- **`smilei`** — full executable, runs the PIC loop
- **`smilei_test`** — parses namelist, initializes data structures, but skips the time loop. Catches ~80% of namelist errors in seconds. **Always run before launching:**

```bash
./smilei_test my_namelist.py        # quick sanity check
./smilei my_namelist.py             # real run (single rank)
```

## 2. MPI and OpenMP

Hybrid parallelism: MPI ranks own disjoint sets of patches; OpenMP threads parallelize patches within a rank.

```bash
# 4 MPI ranks
mpirun -np 4 ./smilei my_namelist.py

# 4 MPI × 4 OpenMP = 16 cores
export OMP_NUM_THREADS=4
mpirun -np 4 ./smilei my_namelist.py
```

For an N-core node, choose `MPI_ranks × OMP_threads = N`. Common balanced choice on 32-core node: `8 MPI × 4 OMP`. For memory-bound runs: `4 MPI × 8 OMP`. For high MPI cost tolerance: `32 MPI × 1 OMP`.

## 3. Patch count rules

- `number_of_patches` must be a power of 2 per direction.
- Total patches ≥ 4 × MPI_ranks, preferably 8–16×.

| MPI ranks | Patches (2D) | Patches (3D) |
|---|---|---|
| 64 | `[32, 16]` = 512 | `[16, 8, 4]` = 512 |
| 256 | `[64, 32]` = 2048 | `[16, 16, 8]` = 2048 |
| 1024 | `[128, 64]` = 8192 | `[32, 16, 16]` = 8192 |

Violating power-of-2 → opaque startup error.

## 4. Cluster run template (SLURM)

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

8 nodes × 4 MPI/node × 4 OMP = 128 cores → 32 MPI ranks → patches ≥ 256.

**Always start with a short test** (`simulation_time = 5.`) before committing 24-hour jobs.

## 5. GPU mode

Requires GPU-aware compile (NVIDIA `nvhpc`, AMD `rocm`).

```python
Main(
    # ...
    gpu_computing = True,
)
```

```bash
mpirun -np 4 ./smilei my_namelist.py    # one MPI rank per GPU
```

GPU wins for large simulations with high PPC and dense distributions. CPU may win for small or sparse setups. As of v5.1, Cartesian geometries are supported; AMcylindrical and collisions have partial support — check release notes.

## 6. Checkpoints

For long simulations or wall-time-limited cluster jobs:

```python
Checkpoints(
    dump_step = 10000,
    dump_minutes = 60.0,
    exit_after_dump = True,                # exit after writing
    keep_n_dumps = 2,
    restart_dir = None,                    # for restart: set to previous job's dir
)
```

First run: `restart_dir = None`. On restart, set `restart_dir` to the previous dump directory. Diagnostic output appends across checkpoint boundaries.

For chained jobs, use a wrapper script that increments the restart directory.

## 7. Symptom table — won't start

| Symptom | Most likely cause | Check |
|---|---|---|
| "ERROR_NAMELIST: invalid keyword" | Misspelled parameter | Docs / examples |
| "Patch number not a power of 2" | Rule 2 violated | Use `[2^a, 2^b]` |
| "Species 'X' not found" | Referenced species not defined | `species1=`, `ionization_electrons=` |
| Hangs at startup | MPI deadlock | Run `-np 1` first |
| "HDF5: file already exists" | Previous run wasn't cleaned | `rm -rf` output dir |
| `smilei: command not found` | Not built or PATH issue | Re-run `make` |
| Hangs reading namelist | Infinite Python loop in preamble | Inspect namelist Python code |
| ImportError | Python module missing | `pip install` |
| OpenMP threads = 1 unexpectedly | `OMP_NUM_THREADS` not set | `echo $OMP_NUM_THREADS` |

## 8. Symptom table — crashes during run

| Symptom | Most likely cause | Check |
|---|---|---|
| Segfault in particle push | NaN momentum/position | Profile callable; absurd initial momentum |
| Segfault in field solve | Memory exhausted | Decrease PPC, box size, or increase memory |
| "MPI_Abort: rank X" | Per-rank crash | Check that rank's stdout |
| Energy diverges exponentially | Numerical instability — under-resolved | Decrease `Δt`, finer grid, CFL margin |
| Charge conservation broken | Maxwell/current inconsistency | Verify standard solver; reduce `Δt` |
| Crash after restart | Corrupt checkpoint or version mismatch | Restart from earlier dump |
| Memory grows unboundedly | Ionization or injector leak | Cap `maximum_charge_state`; check `ParticleInjector` rates |
| One patch stalls entire rank | Load imbalance | `LoadBalancing`; reduce patch size |

## 9. Symptom table — runs but results are wrong

| Symptom | Most likely cause | Check |
|---|---|---|
| Collision rates wrong | `reference_angular_frequency_SI` not set | Hard rule 1 |
| Ionization rate zero everywhere | Wrong model for the regime | Compute `γ_K`; PPT or `from_rate` |
| Density not where expected | Profile in normalized vs SI | Print profile; convert SI → `L_r` |
| Plasma "explodes" at t=0 | Quasi-neutrality violation | `momentum_initialization = "cold"` first, then add T |
| Laser doesn't enter box | Wrong `box_side`; periodic on launch face | Match `EM_boundary_conditions` with `Laser.box_side` |
| Laser intensity off by 2× | Linear vs circular polarization | Check `ellipticity` and `a₀` |
| Spurious EM pulse at t=0 | Dense relativistic beam, no field init | `relativistic_field_initialization = True` |
| Numerical Cherenkov radiation | Relativistic + standard Yee | Use `Lehe` solver |
| No ionization in envelope mode | `tunnel` with envelope laser | Use `tunnel_envelope_averaged` |
| Charge state stuck at +1 | `maximum_charge_state = 1` (default) | Set to expected peak |
| Unphysical temperature spike | Too few PPC | PPC ≥ 32 for collisional |
| Energy not conserved | Under-resolution; coarse `Δt`; wrong solver | Halve `cell_length`; check CFL |

## 10. Symptom table — unexpectedly slow

| Symptom | Most likely cause | Check |
|---|---|---|
| Some MPI ranks idle | Load imbalance | `LoadBalancing`; `DiagPerformances` |
| Slow file I/O | Many small `DiagFields` writes | Combine; `every >>`; `subgrid` |
| Slow on cluster, fast on workstation | Wrong machine file | Re-build with `machine=` flag |
| GPU mode slow | Overhead with small sim | Try CPU; or scale up |
| Hangs periodically | Checkpoint writing | Expected; reduce `dump_step` frequency |
| Vectorization not engaged | Compile-time SIMD flags missing | Use correct machine file |
| Collision module dominates time | Many blocks, high PPC | Reduce non-essential channels |

## 11. Debugging utilities

**`print_every`**: per-step status to stdout. Watch for divergence.

```python
Main(print_every = 100)
```

**`--verbose`**: extra init and load-balancing info.

```bash
./smilei --verbose my_namelist.py
```

**Effective namelist**: `smilei.py` in the output directory is the namelist *as Smilei interpreted it*. Compare to input to catch discrepancies.

**`DiagPerformances`**: always-on for cluster runs to diagnose load imbalance.

**Standalone profile testing**: for complex density/temperature profiles, plot in standalone Python before plugging in:

```python
import numpy as np, matplotlib.pyplot as plt
xs = np.linspace(0, 30, 300)
plt.plot(xs, [n_profile(x, 10.) for x in xs]); plt.show()
```

**Community help**: [GitHub Discussions](https://github.com/SmileiPIC/Smilei/discussions). Before posting: run `smilei_test`, include version (`git log`), share full namelist or minimal reproducer, describe expected vs observed behavior in numbers.

Most simulation bugs reduce to issues in tables 7–10. The community has seen them all.
