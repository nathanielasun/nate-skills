---
name: warpx-python
description: Use whenever the user works with WarpX from Python — writing PICMI input scripts, running simulations via `pywarpx`, extending WarpX at runtime with callbacks (`installafterstep`, `callfromafterstep`, etc.), coupling to ML/AI frameworks, accessing fields/particles in situ via pyAMReX, or post-processing WarpX output with openPMD-viewer, openPMD-api, or yt. Trigger on "PICMI", `pywarpx`, `picmi.Simulation`, `picmi.Species`, `picmi.MCCCollisions`/`DSMCCollisions`/`CoulombCollisions`, `picmi.FieldIonization`, "WarpX input script", `installafterstep`/`callfromafterstep` and any other callback hooks, `Simulation.extension`, `ParticleContainerWrapper`, `multi_particle_container`, openPMD-viewer (`OpenPMDTimeSeries`, `get_field`, `get_particle`), or "yt" with WarpX context. Optimized for low-temperature-plasma work (MCC/DSMC/ionization, capacitive/inductive discharges) but covers the full Python surface of WarpX. For C++ internals work, defer to warpx-cpp.
---

# WarpX Python Development Skill

This skill is for **working with WarpX from Python** — writing PICMI input scripts, running simulations, extending behavior at runtime, and analyzing output. Bias toward low-temperature plasma (MCC/DSMC/discharges) but the structural advice applies across applications.

WarpX exposes Python in four overlapping layers, and the most common source of confusion is mixing them up. Knowing which layer you're in is half the battle.

| Layer | Lives in | What it does |
|---|---|---|
| **PICMI** (standard) | `pywarpx.picmi` | Cross-code, declarative simulation setup. Same class names as PIConGPU, FBPIC, Smilei, etc. |
| **WarpX-specific PICMI** | same module | `warpx_*` keyword arguments on every class; expose features outside the PICMI standard. |
| **Native pywarpx** | `pywarpx.warpx`, `pywarpx.libwarpx`, etc. | Low-level Python-facing wrapper of WarpX's C++ classes. Used inside callbacks and for direct data access. |
| **pyAMReX** | `pywarpx.libwarpx.amr` (= `amrex.space{1,2,3}d`) | The AMReX Python bindings underneath. Use for raw `MultiFab`/`ParticleContainer` access. |

A typical "low-temp-plasma simulation in Python" stays in PICMI layer for setup, dips into native pywarpx layer for runtime callbacks, and uses pyAMReX layer for direct data manipulation inside those callbacks. Post-processing is yet another stack — openPMD-viewer or yt for reading output files.

## Before doing anything else: check for in-repo context

WarpX ships with its own LLM-facing context that takes precedence over this skill:

- **`AGENTS.md`** / **`CLAUDE.md`** at the repo root. Version-matched build/test/style instructions.
- **`.claude/skills/`** in the repo. The in-repo skills (`/warpx-answer-user-question`, `/warpx-new-paper-highlight`) are workflow-specific and authoritative for the workflows they cover.
- **PICMI standard documentation** at `https://picmi-standard.github.io/` — the cross-code base. Any class in `pywarpx.picmi` that doesn't have a `warpx_` prefix is documented there.
- **WarpX readthedocs Python page** at `https://warpx.readthedocs.io/en/latest/usage/python.html` — the full PICMI reference for WarpX, including all `warpx_*` extensions.
- **pyAMReX docs** at `https://pyamrex.readthedocs.io/` — for raw AMReX bindings.
- **openPMD-viewer docs** at `https://openpmd-viewer.readthedocs.io/` — for post-processing.

For Context7 MCP server users, the WarpX, AMReX, pyAMReX, openPMD-api, openPMD-viewer, PICMI, and pybind11 docs are all queryable through Context7.

## What to read in this skill, when

| If the task involves... | Read |
|---|---|
| Starting a new PICMI input script, choosing solver/grid/geometry, applied fields, mirrors, plasma lenses | `references/picmi-fundamentals.md` |
| Species setup, initial distributions, layouts, beam injection, particle reflection/absorption | `references/species-and-injection.md` |
| **MCC, DSMC, Coulomb collisions, field ionization, low-temp plasma discharge configuration** | `references/collisions-and-reactions.md` |
| Diagnostics: ParticleDiagnostic, FieldDiagnostic, openPMD output, reduced diags, checkpoints | `references/diagnostics-and-output.md` |
| **Callbacks (`installafterstep`, decorators), accessing simulation state at runtime, ML coupling, custom Python solvers, real-time analysis** | `references/runtime-extension.md` |
| **Post-processing**: reading output with openPMD-viewer, openPMD-api, or yt; computing reduced quantities; time-series analysis | `references/post-processing.md` |
| Installing pywarpx, running with mpirun, HPC submission, debugging Python-side errors | `references/install-and-run.md` |

References are designed to be readable independently. For most tasks you'll need one or two; if you find yourself opening all seven, the request is probably actually multiple tasks.

## The two operating modes — never confuse them

A PICMI script can be used in two distinct ways. The mode determines what's available and what isn't.

### Interactive mode: `Simulation.step()`
The simulation runs inside the Python process. You can install callbacks. You can access fields and particles in situ. You can couple to ML frameworks, real-time analysis, custom Python operators.

```python
sim = picmi.Simulation(...)
# ... configure ...
sim.step(nsteps=1000)  # blocks until done
# After step(), can access final state via sim.extension
```

Requirement: `pywarpx` must be installed (built with `WarpX_PYTHON=ON`) and `python3 script.py` is the entrypoint. With MPI: `mpiexec -n N python3 script.py`.

### Preprocessor mode: `Simulation.write_input_file()`
The PICMI script generates a C++-style input file. The actual simulation is then run by the standalone WarpX binary.

```python
sim = picmi.Simulation(...)
# ... configure ...
sim.write_input_file('inputs_test')  # produces text file
# Then run separately: mpiexec -n N warpx.3d inputs_test
```

No callbacks possible. No in-situ access. But works with any WarpX binary (HPC modules, etc.) and produces a portable input file. Useful for benchmarking, large-scale HPC runs, or when you want to share an input.

**Which to use?** Interactive for: development, ML coupling, custom physics, in-situ diagnostics, single-node workflows. Preprocessor for: large-scale HPC, portability, sharing inputs, integration with C++-only pipelines.

The skill assumes interactive mode unless otherwise noted.

## The `warpx_*` prefix convention — read this carefully

PICMI is a cross-code standard. To expose WarpX-specific functionality that isn't in the standard, every PICMI class accepts `warpx_*`-prefixed keyword arguments. This is where most of the actual control surface lives.

Examples:
```python
sim = picmi.Simulation(
    solver=es_solver,
    time_step_size=1e-12,
    max_steps=1000,
    warpx_current_deposition_algo='direct',          # not in PICMI standard
    warpx_particle_pusher_algo='boris',
    warpx_random_seed=42,
    warpx_load_balance_intervals=100,
    warpx_amrex_use_gpu_aware_mpi=True,
)

grid = picmi.Cartesian2DGrid(
    nx=128, ny=128,
    xmin=0, xmax=0.02, ymin=0, ymax=0.02,
    bc_xmin='dirichlet', bc_xmax='dirichlet',
    bc_ymin='periodic', bc_ymax='periodic',
    warpx_max_grid_size=32,
    warpx_blocking_factor=8,
    warpx_potential_lo_x='100.0 * sin(2*pi*13.56e6*t)',  # parsed expression
)

electrons = picmi.Species(
    particle_type='electron', name='electrons',
    initial_distribution=...,
    warpx_save_particles_at_xlo=True,
    warpx_do_not_deposit=False,
    warpx_reflection_model_xhi='1.0',  # always reflect at upper x
    warpx_do_resampling=True,
    warpx_resampling_trigger_intervals='100',
)
```

The `warpx_*` arguments map directly to entries in the C++-style input file. To find the exact name for a feature: check the [Parameter List](https://warpx.readthedocs.io/en/latest/usage/parameters.html) on readthedocs (the C++-input documentation), find the parameter, then prepend `warpx_` and use it as a kwarg.

This indirection is the single biggest source of confusion for new users. PICMI documents the abstract feature; the implementation-specific parameters live in the C++ docs. They're cross-referenced but easy to miss.

## Decision rules when writing a PICMI script

A few high-value rules that come up repeatedly.

**On structure.** A canonical PICMI script has six sections in order: imports, constants/parameters, grid and solver, species and injection, collisions and other physics, diagnostics, finally `sim.step()` or `sim.write_input_file()`. Following this order produces scripts that other WarpX users can read at a glance.

**On `step()` granularity.** For pure simulation runs, call `sim.step(max_steps)` once. For runs with custom logic between steps, call `sim.step(1)` in a loop. The latter is slower because it interrupts the C++ loop every step; only use it when you genuinely need per-step Python control.

**On callbacks vs `step(1)` loops.** Both can intervene every step, but callbacks are faster (the C++ loop continues uninterrupted) and they have access to specific simulation lifecycle points (after deposition, after Esolve, etc.). Use callbacks unless you specifically need to run code *between* full steps.

**On expression-as-string parameters.** Many PICMI parameters accept either numeric values or string expressions (`'100.0 * sin(2*pi*f*t)'`). String expressions are parsed by WarpX's runtime parser, support `x`, `y`, `z`, `t` as variables, and can reference user-defined parameters. This is the right way to specify time-varying or spatially-varying parameters — don't try to implement them via callbacks.

**On units.** PICMI is SI: meters, seconds, kilograms, coulombs, tesla, volts/meter. Always. The constants module provides `picmi.constants.c`, `picmi.constants.q_e`, etc., also in SI. WarpX internally uses SI throughout; don't try to use other unit systems.

**On dimensionality.** Each PICMI script targets exactly one dimensionality (1D, 2D, 3D, RZ). The `pywarpx` package builds separate binaries per dimension (`pywarpx.warpx_1d`, `pywarpx.warpx_2d`, etc.). The dimensionality is selected by the `Grid` class you instantiate. If you want to run the same physics in multiple dimensions, parameterize the script and have `if args.dim == 3:` blocks; don't try to write dimension-agnostic code.

**On RZ mode.** Cylindrical symmetric (`CylindricalGrid`) with `n_azimuthal_modes >= 1`. Mode 0 is azimuthally symmetric; higher modes capture asymmetric perturbations. Some features (notably hybrid PIC) only support mode 0. The radial direction has `bc_rmin = 'none'` implicitly at r=0 (axis); avoid setting it.

**On embedded boundaries.** Use `picmi.EmbeddedBoundary(implicit_function=...)` or `stl_file=...`. The implicit function is an analytic expression (parsed at runtime); negative values are "inside the wall," positive are "in the simulation domain." EB is incompatible with RZ as of late 2025; check the current docs if you need both.

**On collisions setup.** MCC for thermal-background collisions (cheapest), DSMC for binary collisions between simulated species, Coulomb for charged-charged within or across species. They can coexist in the same simulation. See `collisions-and-reactions.md` for the detailed configuration.

**On diagnostics cadence.** Field and particle diagnostics with `period` set too small (every step) will dominate the wall-time. Typical: every 100-1000 steps for big-picture diagnostics; reduced diagnostics every step are usually cheap. Use `warpx_random_fraction` to subsample particle output for huge species.

## Low-temperature-plasma Python quick reference

If the simulation involves neutral gas, MCC, atmospheric discharges, or hybrid-PIC for ion kinetics — here are the relevant code paths in Python:

- **`picmi.MCCCollisions(name, species, background_density, background_temperature, scattering_processes, background_mass=None, ndt=1)`** — collisions against a thermal neutral background. `scattering_processes` is a dict like `{'elastic': {'cross_section': '/path/to/elastic.dat'}, 'ionization': {'cross_section': '/path/to/ion.dat', 'energy': 15.76, 'species': 'ions'}}`. `background_density`/`background_temperature` can be string expressions of `(x, y, z, t)`.
- **`picmi.DSMCCollisions(name, species, scattering_processes, ndt=1)`** — binary collisions between simulated particles. Used when neutrals are simulated as particles.
- **`picmi.CoulombCollisions(name, species, CoulombLog=None, ndt=1)`** — Perez 2012 binary Coulomb between charged species (intra- or inter-species).
- **`picmi.FieldIonization(model='ADK', ionized_species, product_species)`** — ADK ionization model; only ADK is currently implemented.
- **`warpx_resampling_trigger_intervals` / `warpx_resampling_trigger_max_avg_ppc`** on `Species` — important for steady-state low-temp plasma where ionization runaway can blow up particle count.
- **`picmi.ElectrostaticSolver(grid, method='Multigrid', ...)`** — Multigrid (MLMG) is the workhorse for low-frequency / electrostatic discharge work. FFT is faster but less flexible with boundary conditions.
- **`warpx_potential_lo_x`** etc. on `Grid` — string expressions for time-varying applied voltage on Dirichlet boundaries. The canonical RF discharge driver: `'V_RF * sin(2*pi*f_RF*t)'`.
- **Embedded boundary potential**: `picmi.EmbeddedBoundary(implicit_function=..., potential='...')` — Dirichlet BC on an internal surface.

Canonical reference example: the **capacitive discharge** examples at `Examples/Physics_applications/capacitive_discharge/`, based on the Turner et al. (Phys. Plasmas 20, 013507, 2013) benchmarks. The 1D PICMI input there is the cleanest starting point for an RF discharge in WarpX. It demonstrates: 1D Cartesian grid, electrostatic Multigrid solver, two species (electrons + ions), MCC with multiple processes, time-varying applied potential, and basic diagnostics.

When in doubt about how to configure an MCC scenario, start by copying that example and adapting.

## Working pattern for a typical task

1. **Identify the task type.** New simulation? Modification to an existing script? Adding a callback? Post-processing? Coupling? Each maps to different references and conventions.
2. **Find precedent.** WarpX's `Examples/` directory has 100+ working PICMI scripts covering everything from laser wakefield to capacitive discharge to magnetic reconnection. For low-temp plasma specifically, start from the capacitive discharge or DSMC examples.
3. **Read the relevant reference file(s).** Don't reason about PICMI from first principles — the conventions have specific quirks.
4. **Iterate small.** Start with a simple grid (16x16 cells), short run (100 steps), one species. Verify the script runs, the diagnostics produce output, the output is readable by openPMD-viewer. Then scale up.
5. **Check `Simulation.write_input_file()` output.** Even if you intend to run interactively, calling `write_input_file()` and reading the generated text file is the fastest way to verify that your PICMI parameters mapped to the C++ parameters you intended. Misunderstood `warpx_*` arguments show up here clearly.
6. **Verify dimensionality.** If you wrote a 2D grid and the binary errors out about dimension mismatch, you probably need `pip install pywarpx[2D]` or have `WARPX_DIMS=2` set at build time. The error is unhelpfully vague.
7. **For callbacks: test outside of MPI first.** Race conditions and rank-dependent state make callback debugging much harder under MPI. Develop on a single process, then add MPI.
8. **Diagnostic outputs are the API.** If a callback or analysis needs simulation data, prefer reading it from the diagnostic output (via openPMD-viewer or in-situ via `Simulation.extension`) rather than trying to compute it from inputs. The diagnostic is what other tools speak.

## Honest caveats

- **API evolves.** Class signatures and keyword argument names shift across WarpX releases. The `warpx_*` parameters in particular have been renamed several times. Verify any specific keyword against the current `pywarpx.picmi` docstrings or the readthedocs Parameter List before relying on it.
- **PICMI standard vs WarpX implementation.** Some `picmi.*` classes are pure PICMI standard (cross-code consistent); others are WarpX-only or have WarpX-specific quirks (`EmbeddedBoundary`, `CoulombCollisions`, `DSMCCollisions`, `MCCCollisions`, `FieldIonization` are WarpX-specific). For WarpX-specific classes, expect the API to differ from other PICMI codes.
- **The native pywarpx layer is less stable than PICMI.** `pywarpx.warpx`, `pywarpx.particle_containers`, `Simulation.extension` change more often than PICMI. Pin your pywarpx version or be prepared to adapt scripts.
- **Callbacks have implicit ordering.** Multiple callbacks installed at the same hook run in installation order. If two callbacks need to communicate, use module-level state or an explicit shared object, not assumption about order.
- **GPU and ML coupling has overhead.** Moving data between WarpX's GPU memory and a separate ML framework's GPU memory typically requires a copy (CPU round-trip in the worst case). Use `pti.soa().to_cupy()` to keep data on the GPU when possible, but verify that your ML framework can consume cupy arrays without copying.
- **This skill is for Python work.** For C++ internals (modifying particle pushers, adding new MCC processes at the source level, optimizing kernels), use `warpx-cpp`. PICMI exposes existing C++ machinery — it doesn't let you write new C++ machinery from Python.

## Where to find canonical references

- **WarpX docs**: `https://warpx.readthedocs.io/`
- **PICMI standard**: `https://picmi-standard.github.io/`
- **WarpX PICMI reference**: `https://warpx.readthedocs.io/en/latest/usage/python.html`
- **WarpX extending-Python guide**: `https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html`
- **pyAMReX docs**: `https://pyamrex.readthedocs.io/`
- **openPMD-viewer**: `https://openpmd-viewer.readthedocs.io/`
- **openPMD-api**: `https://openpmd-api.readthedocs.io/`
- **yt project**: `https://yt-project.org/`
- **WarpX examples**: `Examples/` directory in the [BLAST-WarpX repository](https://github.com/BLAST-WarpX/warpx). Especially `Examples/Physics_applications/capacitive_discharge/` for low-temp plasma.
- **warpx-data** (cross sections): `https://github.com/BLAST-WarpX/warpx-data`
- **Discussions**: `https://github.com/BLAST-WarpX/warpx/discussions`

For LLM-assisted work specifically, the same expectations apply as for C++ contributions: read every line of generated code before relying on it, verify against the actual source/docs (PICMI APIs in particular look very natural to hallucinate), test with a small case before committing to a large run, and don't push code through this skill that you haven't read and understood.
