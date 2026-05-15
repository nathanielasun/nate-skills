---
name: warpx-python-mini
description: Use whenever the user works with WarpX from Python — writing PICMI input scripts, running simulations via `pywarpx`, extending WarpX at runtime with callbacks (`installafterstep`, `callfromafterstep`), coupling to ML/AI frameworks, accessing fields/particles via pyAMReX, or post-processing WarpX output with openPMD-viewer or yt. Trigger on "PICMI", `pywarpx`, `picmi.Simulation`, `picmi.Species`, `picmi.MCCCollisions`/`DSMCCollisions`/`CoulombCollisions`, `picmi.FieldIonization`, `installafterstep` and other callback hooks, `Simulation.extension`, `ParticleContainerWrapper`, openPMD-viewer (`OpenPMDTimeSeries`, `get_field`, `get_particle`), or "yt" with WarpX context. Optimized for low-temperature plasma. This is the ultra-compact variant for small models; for C++ work use warpx-cpp-mini.
---

# WarpX Python (Mini)

For PICMI scripts and runtime Python extension. For C++ work use `warpx-cpp-mini`.

PICMI looks natural to write — which means it is easy to invent plausible-looking but wrong class names. **Always verify against the real docs before relying on any specific API.**

## ALWAYS verify APIs before using them

1. **`AGENTS.md`** / **`CLAUDE.md`** at repo root — version-matched authoritative
2. **PICMI reference**: `https://warpx.readthedocs.io/en/latest/usage/python.html`
3. **Extending guide**: `https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html`
4. **PICMI standard** (cross-code base): `https://picmi-standard.github.io/`

If you're not sure something exists, search the docs first.

## The 8 hard rules

| # | Rule |
|---|---|
| 1 | **`warpx_*` prefix**: WarpX-specific features go in `warpx_*` kwargs on every PICMI class. Most control lives there. |
| 2 | **SI units only**: meters, seconds, kg, coulombs, tesla, V/m. Use `picmi.constants.c`, `q_e`, `m_e`, etc. |
| 3 | **One dimensionality per script**: chosen by the `Grid` class. Install matching `pywarpx[3D]`/`[2D]`/etc. |
| 4 | **Ions REQUIRE `charge_state`**: `particle_type='Ar'` alone is NEUTRAL. Always: `picmi.Species(particle_type='Ar', charge_state=1, ...)`. |
| 5 | **Two modes — don't confuse**: `sim.step()` = interactive (callbacks work). `sim.write_input_file('inputs')` = preprocessor (no callbacks). |
| 6 | **MPI order**: `from mpi4py import MPI` BEFORE `from pywarpx import picmi`. Otherwise MPI init conflicts. |
| 7 | **MPI-consistent callbacks**: every rank must do the same collective ops in the same order, or deadlock. |
| 8 | **No Python loops over particles**: use vectorized numpy/cupy. A 10ms loop × 100,000 steps = 1000s of waste. |

## The four Python layers

| Layer | Module | Use for |
|---|---|---|
| PICMI standard | `pywarpx.picmi` | Cross-code declarative setup |
| WarpX-specific PICMI | same | `warpx_*` kwargs (most control) |
| Native pywarpx | `pywarpx.callbacks`, `pywarpx.fields`, `pywarpx.particle_containers` | Runtime, callbacks |
| pyAMReX | `pywarpx.libwarpx.amr` | Raw MultiFab/ParticleContainer |

## What to read

| Task | Read |
|---|---|
| Set up Simulation/Grid/Solver/Species/Distribution + diagnostics | `references/picmi-basics.md` |
| **MCC, DSMC, Coulomb, ionization (low-temp plasma)** | `references/collisions.md` |
| Callbacks, ML coupling, accessing simulation state at runtime | `references/runtime-and-ml.md` |
| Install, MPI, read output with openPMD-viewer/yt, common errors | `references/install-and-output.md` |

## The two operating modes

### Interactive: `sim.step()`
```python
sim = picmi.Simulation(...)
# ... configure ...
sim.step(nsteps=1000)
```
Runs inside Python. Callbacks work. Run with `python3 script.py` or `mpiexec -n N python3 script.py`.

### Preprocessor: `sim.write_input_file()`
```python
sim = picmi.Simulation(...)
# ... configure ...
sim.write_input_file('inputs')        # generates C++ input file; no callbacks
```
Then run separately: `mpiexec -n N warpx.3d inputs`.

**Tip**: even when running interactively, call `write_input_file()` once first to inspect the generated C++ input — it's the fastest way to verify your `warpx_*` args mapped to expected parameters.

## Low-temp plasma quick reference

| Feature | API |
|---|---|
| Background MCC vs thermal gas | `picmi.MCCCollisions(name, species, background_density, background_temperature, scattering_processes, background_mass, ndt)` |
| DSMC between simulated species | `picmi.DSMCCollisions(name, species=[s1, s2], scattering_processes, ndt)` |
| Coulomb collisions | `picmi.CoulombCollisions(name, species=[s1, s2], CoulombLog=None, ndt)` |
| ADK field ionization | `picmi.FieldIonization(model='ADK', ionized_species, product_species)` |
| ES Multigrid solver (workhorse) | `picmi.ElectrostaticSolver(grid, method='Multigrid', required_precision=1e-7)` |
| Time-varying applied voltage | `warpx_potential_lo_x='V_RF * sin(2*pi*f_RF*t)'` on Grid |
| Particle resampling (essential w/ ionization) | `warpx_do_resampling=True` on Species |
| Save particles at boundaries | `warpx_save_particles_at_xlo=True` etc. on Species |

**Canonical example**: `Examples/Physics_applications/capacitive_discharge/inputs_base_1d_picmi.py` — Turner et al. 2013 RF helium discharge. Start by copying this for any RF discharge work.

## Working pattern

1. **Find an example to copy**. `Examples/` in `https://github.com/BLAST-WarpX/warpx` has 100+ working scripts.
2. **Read the relevant reference file**.
3. **Start small** — 16×16 cells, 100 steps, one species. Verify it runs.
4. **Verify with `write_input_file()`**.
5. **For callbacks: test serially first**, then add MPI.

## Links

- WarpX docs: `https://warpx.readthedocs.io/`
- PICMI reference: `https://warpx.readthedocs.io/en/latest/usage/python.html`
- openPMD-viewer: `https://openpmd-viewer.readthedocs.io/`
- Examples: `https://github.com/BLAST-WarpX/warpx`
- Cross sections: `https://github.com/BLAST-WarpX/warpx-data`
