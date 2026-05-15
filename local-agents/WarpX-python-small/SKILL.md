---
name: warpx-python-small
description: Use whenever the user works with WarpX from Python — writing PICMI input scripts, running simulations via `pywarpx`, extending WarpX at runtime with callbacks (`installafterstep`, `callfromafterstep`), coupling to ML/AI frameworks, accessing fields/particles in situ via pyAMReX, or post-processing WarpX output with openPMD-viewer or yt. Trigger on "PICMI", `pywarpx`, `picmi.Simulation`, `picmi.Species`, `picmi.MCCCollisions`/`DSMCCollisions`/`CoulombCollisions`, `picmi.FieldIonization`, `installafterstep`/`callfromafterstep` and other callback hooks, `Simulation.extension`, `ParticleContainerWrapper`, openPMD-viewer (`OpenPMDTimeSeries`, `get_field`, `get_particle`), or "yt" with WarpX context. Optimized for low-temperature plasma. This is the compact variant tuned for mid-size models; for C++ codebase work use the warpx-cpp or warpx-cpp-small skill.
---

# WarpX Python Development (Compact)

This skill is for **PICMI scripts and runtime Python extension of WarpX**. For modifying WarpX's C++ source, use `warpx-cpp` / `warpx-cpp-small`.

PICMI looks natural to write — which means it is easy to hallucinate plausible-looking but wrong class names and arguments. **Always verify against the real docs** before relying on any specific API.

## Step 0: Check for in-repo context FIRST

1. **`AGENTS.md`** / **`CLAUDE.md`** at repo root — version-matched authoritative.
2. **`.claude/skills/`** in the repo — in-repo workflow skills override this skill.
3. **WarpX PICMI reference**: `https://warpx.readthedocs.io/en/latest/usage/python.html`
4. **WarpX Python extending guide**: `https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html`
5. **PICMI standard**: `https://picmi-standard.github.io/` (for non-WarpX-specific classes)

**Verification ladder**: actual `pywarpx.picmi` docstrings → readthedocs → this skill. Never invent an API.

## The four Python layers — know which you're in

| Layer | Module | Use for |
|---|---|---|
| **PICMI standard** | `pywarpx.picmi` | Cross-code declarative setup (same in PIConGPU, FBPIC) |
| **WarpX-specific PICMI** | same | `warpx_*` keyword args on every PICMI class — most control lives here |
| **Native pywarpx** | `pywarpx.warpx`, `pywarpx.libwarpx`, `pywarpx.callbacks`, `pywarpx.fields`, `pywarpx.particle_containers` | Callbacks, runtime data access |
| **pyAMReX** | `pywarpx.libwarpx.amr` | Raw `MultiFab` / `ParticleContainer` access |

Typical script: PICMI for setup, native pywarpx for callbacks, pyAMReX inside callbacks for data manipulation. Post-processing: openPMD-viewer or yt.

## What to read in this skill

| Task | Read |
|---|---|
| Setting up Simulation/Grid/Solver/Species/Distribution/Layout/applied fields | `references/picmi-basics.md` |
| **MCC, DSMC, Coulomb, ionization (low-temp plasma)** | `references/collisions.md` |
| Callbacks, accessing simulation state at runtime, ML coupling, custom Python solvers | `references/runtime-and-ml.md` |
| Configuring diagnostics + reading openPMD/plotfile output with viewer/yt | `references/output-and-analysis.md` |
| Installing pywarpx, running with MPI, common errors, debugging | `references/install-and-debug.md` |

Read at most one or two per task.

## Hard rules (non-negotiable)

1. **`warpx_*` prefix**: WarpX-specific features go in `warpx_*` keyword arguments on every PICMI class. Most of the actual control surface lives here, not in the PICMI base arguments.
2. **SI units everywhere**: meters, seconds, kilograms, coulombs, tesla, V/m. PICMI is SI-only. Use `picmi.constants.c`, `picmi.constants.q_e`, `picmi.constants.m_e`, etc.
3. **One dimensionality per script**: 1D, 2D, 3D, or RZ. Selected by the `Grid` class. The pywarpx binaries are per-dim — install `pywarpx[3D]` to match a `Cartesian3DGrid` script.
4. **Ions must have `charge_state`**: a species with `particle_type='Ar'` and NO `charge_state` is NEUTRAL (charge=0). Always specify explicitly.
5. **Two operating modes — don't confuse them**: `sim.step()` runs interactively (callbacks available); `sim.write_input_file('inputs')` generates a text file for separate C++ run (no callbacks).
6. **MPI init order**: import `mpi4py.MPI` BEFORE `pywarpx` if using both. pywarpx initializes MPI on import if mpi4py hasn't already.
7. **Callbacks must be MPI-consistent**: every rank must call the same collective operations in the same order, or you'll deadlock.
8. **Don't loop over particles in Python**: use vectorized numpy/cupy operations. A 10ms Python overhead × 100,000 steps = 1000s of waste.

## The two operating modes in detail

### Interactive mode: `sim.step()`
```python
sim = picmi.Simulation(...)
# ... configure ...
sim.step(nsteps=1000)         # runs inside Python; callbacks work
```
Run with `python3 script.py` or `mpiexec -n N python3 script.py`. The simulation runs in this Python process.

### Preprocessor mode: `sim.write_input_file()`
```python
sim = picmi.Simulation(...)
# ... configure ...
sim.write_input_file('inputs')   # produces C++ input file; no callbacks
```
Then run separately: `mpiexec -n N warpx.3d inputs`. No Python at runtime.

**Use interactive for**: development, ML coupling, custom physics, in-situ diagnostics. **Use preprocessor for**: large-scale HPC, sharing inputs, integration with C++-only pipelines.

Even in interactive mode, calling `write_input_file()` once to inspect the generated input is the fastest way to verify that `warpx_*` arguments mapped to the C++ parameters you intended.

## Low-temp plasma quick reference

Code paths for typical MCC/DSMC/ionization work:

| Feature | API |
|---|---|
| Background MCC vs thermal neutral gas | `picmi.MCCCollisions(name, species, background_density, background_temperature, scattering_processes, background_mass, ndt)` |
| DSMC between simulated species | `picmi.DSMCCollisions(name, species=[s1, s2], scattering_processes, ndt)` |
| Coulomb collisions | `picmi.CoulombCollisions(name, species=[s1, s2], CoulombLog=None, ndt)` |
| ADK field ionization | `picmi.FieldIonization(model='ADK', ionized_species, product_species)` |
| Electrostatic Multigrid solver (workhorse for discharges) | `picmi.ElectrostaticSolver(grid, method='Multigrid', ...)` |
| Time-varying applied potential | `warpx_potential_lo_x='V_RF * sin(2*pi*f_RF*t)'` on `Grid` |
| Particle resampling (essential with ionization runaway) | `warpx_do_resampling=True` on `Species` |
| Save particles crossing boundaries | `warpx_save_particles_at_xlo=True` etc. on `Species` |

Canonical example: `Examples/Physics_applications/capacitive_discharge/inputs_base_1d_picmi.py` — Turner et al. 2013 RF helium discharge benchmark. Start there for any capacitive/RF discharge work.

## Working pattern for a typical task

1. **Find an example to copy.** `Examples/` in the BLAST-WarpX repo has 100+ working PICMI scripts. For low-temp plasma, start from `capacitive_discharge/`. Adapting a known-good script is faster and more reliable than building from scratch.
2. **Read the relevant reference file(s).**
3. **Start small.** 16×16 cells, 100 steps, one species. Verify it runs and produces readable output before scaling up.
4. **Verify with `write_input_file()`.** Check that your PICMI parameters mapped to the expected C++ parameters.
5. **For callbacks: test serially first.** Then add MPI. Most callback bugs are MPI-related and harder to debug under MPI.

## Critical caveats

- **PICMI API evolves across releases.** Argument names (especially `warpx_*`) shift. Verify against the version you have installed.
- **The native pywarpx layer is less stable than PICMI.** `Simulation.extension`, particle iterators, callback hooks change more often. Pin pywarpx version or expect to adapt.
- **Don't try to write dimension-agnostic code.** Each script targets one dimensionality.

## Canonical links

- WarpX docs: `https://warpx.readthedocs.io/`
- PICMI reference: `https://warpx.readthedocs.io/en/latest/usage/python.html`
- Extending guide: `https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html`
- pyAMReX docs: `https://pyamrex.readthedocs.io/`
- openPMD-viewer: `https://openpmd-viewer.readthedocs.io/`
- Examples: `Examples/` in `https://github.com/BLAST-WarpX/warpx`
- Cross sections: `https://github.com/BLAST-WarpX/warpx-data`
