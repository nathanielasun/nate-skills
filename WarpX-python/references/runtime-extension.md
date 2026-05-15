# Runtime extension: callbacks, pyAMReX, and ML/custom-code coupling

This reference covers how to extend WarpX from Python while it runs: installing callbacks at PIC-loop hook points, accessing simulation state mid-run, modifying particles and fields in situ, coupling to ML/AI frameworks, and using Python in place of part of the PIC loop. This is the distinctive capability that makes PICMI more than just an input-script standard.

## Table of contents

1. [The big picture](#the-big-picture)
2. [Callback hook points](#callback-hook-points)
3. [Installing callbacks: decorators and explicit installation](#installing-callbacks-decorators-and-explicit-installation)
4. [Accessing simulation state: Simulation.extension](#accessing-simulation-state-simulationextension)
5. [Reading and writing fields from Python](#reading-and-writing-fields-from-python)
6. [Reading and writing particles from Python](#reading-and-writing-particles-from-python)
7. [The ParticleContainerWrapper convenience](#the-particlecontainerwrapper-convenience)
8. [pyAMReX: direct access to AMReX data structures](#pyamrex-direct-access-to-amrex-data-structures)
9. [GPU vs CPU: numpy and cupy](#gpu-vs-cpu-numpy-and-cupy)
10. [Coupling to ML frameworks](#coupling-to-ml-frameworks)
11. [Replacing parts of the PIC loop with Python](#replacing-parts-of-the-pic-loop-with-python)
12. [Real-time diagnostics](#real-time-diagnostics)
13. [Patterns and pitfalls](#patterns-and-pitfalls)

## The big picture

Inside a PICMI script running in interactive mode, after `sim.step()` is called, the C++ PIC loop is in charge. To inject Python logic at specific points, install **callbacks** — Python functions that WarpX invokes at well-defined moments in each step.

The fundamental capabilities:

1. **Run Python code at hook points** — `afterstep`, `beforeEsolve`, `afterdeposition`, etc. The C++ loop continues uninterrupted between hooks.
2. **Read simulation state** — `Simulation.extension.warpx` exposes the C++ WarpX object. Particles, fields, geometry, time, step number are all accessible.
3. **Modify simulation state** — update particle positions/momenta, add new particles, modify field values, install custom Poisson solvers.
4. **Pass data in/out** — numpy arrays for CPU, cupy arrays for GPU (zero-copy when possible).

The cost model: a callback invocation is a Python function call, fast (microseconds). The actual work inside the callback dominates. As long as you use vectorized array operations (numpy/cupy) rather than Python-level loops over particles, callbacks scale to billions of particles.

## Callback hook points

WarpX exposes these PIC-loop hook points (in approximate order within a step):

| Hook | When it fires |
|---|---|
| `beforeInitEsolve` | Before the initial E-field solve at simulation start |
| `afterInitEsolve` | After the initial E-field solve at simulation start |
| `afterinit` | After all initialization (fields and particles) |
| `beforestep` | At the start of each step |
| `beforecollisions` | Before the collision operator runs |
| `aftercollisions` | After the collision operator runs |
| `beforeEsolve` | Before the E-field solver |
| `poissonsolver` | *Replaces* the standard Poisson solve (electrostatic only) — see below |
| `afterEsolve` | After the E-field solver |
| `beforedeposition` | Before particle current/charge deposition |
| `afterdeposition` | After particle current/charge deposition |
| `particleinjection` | At particle injection time, after position advance, before deposition |
| `particleloader` | When the standard particle loader runs |
| `particlescraper` | Just after particle boundary conditions, before lost particles are processed |
| `appliedfields` | When applied fields are set on the grid |
| `afterstep` | At the end of each step |
| `afterdiagnostics` | After diagnostic output has been written |
| `oncheckpointsignal` | When a checkpoint signal is received |
| `afterrestart` | After a restart from checkpoint |

`afterstep` is by far the most commonly used hook — for ML inference, custom diagnostics, parameter sweeps, anything that needs to run "once per step." For modifying particles or fields before they're used in the next step, `afterstep` is also where you want to be (after everything in step N is done, before step N+1 begins).

`poissonsolver` is special: it *replaces* the default Poisson solver. Use it to install a custom Python-level Poisson solver (machine-learning-based, multigrid variant, etc.). Only fires in electrostatic mode.

## Installing callbacks: decorators and explicit installation

Two equivalent ways to install:

### Decorator style (cleaner)

```python
from pywarpx.callbacks import callfromafterstep

@callfromafterstep
def my_diagnostic():
    # Function name is automatically installed at the 'afterstep' hook
    print("Step finished")
```

The function is installed at the hook when the decorator runs (i.e., when Python imports/executes the script).

Corresponding decorators for other hooks: `callfrombeforestep`, `callfromafterEsolve`, `callfromafterdeposition`, `callfromparticlescraper`, `callfromparticleloader`, `callfromparticleinjection`, `callfromafterstep`, `callfromafterdiagnostics`, `callfromafterinit`, `callfrombeforeEsolve`, `callfrombeforedeposition`.

### Explicit installation

```python
from pywarpx.callbacks import installcallback, installafterstep

def my_diagnostic():
    print("Step finished")

# Either generic:
installcallback('afterstep', my_diagnostic)

# Or hook-specific:
installafterstep(my_diagnostic)
```

Both styles work. Decorators are cleaner for static installation; explicit is more flexible (e.g., conditionally installing only if certain command-line options are set).

### Uninstalling

```python
from pywarpx.callbacks import uninstallafterstep, isinstalledafterstep

if isinstalledafterstep(my_diagnostic):
    uninstallafterstep(my_diagnostic)
```

Useful for callbacks that should fire only for a phase of the simulation:

```python
# Install at start
installafterstep(transient_diagnostic)

# Later in the script, after some condition:
sim.step(1000)
uninstallafterstep(transient_diagnostic)  # don't run after first 1000 steps
sim.step(10000)
```

### Multiple callbacks at one hook

Multiple functions can be installed at the same hook; they fire in installation order. Each callback runs independently; if one raises an exception, behavior is implementation-defined (often the simulation aborts).

## Accessing simulation state: Simulation.extension

Inside a callback (or after `sim.step()` returns), the C++ WarpX state is accessible through `sim.extension`:

```python
from pywarpx import picmi, callbacks

sim = picmi.Simulation(...)

@callbacks.callfromafterstep
def my_callback():
    warpx = sim.extension.warpx           # the C++ WarpX object
    Config = sim.extension.Config         # compile-time configuration

    # Current step and time
    step = warpx.getistep(0)              # step number on level 0
    t = warpx.gett_new(0)                  # current physical time, s
    dt = warpx.getdt(0)                    # time step size, s

    # GPU availability
    on_gpu = Config.have_gpu

    # Geometry
    geom = warpx.Geom(0)                   # AMReX Geometry on level 0
    prob_lo = geom.ProbLoArray()
    prob_hi = geom.ProbHiArray()
    dx = geom.CellSizeArray()
```

`sim.extension.warpx` is the Python binding of the C++ `WarpX` class. Most C++ methods are exposed; the most commonly-used:

| Method | Returns |
|---|---|
| `warpx.getistep(lev)` | Step number on level `lev` |
| `warpx.gett_new(lev)` | New (current) physical time on level |
| `warpx.gett_old(lev)` | Old (previous) physical time |
| `warpx.getdt(lev)` | Time step size on level |
| `warpx.finest_level` | Finest mesh refinement level number |
| `warpx.multi_particle_container()` | The `MultiParticleContainer` (= C++ `mypc`) |
| `warpx.Geom(lev)` | AMReX `Geometry` on level |
| `warpx.boxArray(lev)` | AMReX `BoxArray` on level |
| `warpx.DistributionMap(lev)` | AMReX `DistributionMapping` on level |

`Simulation.extension` only exists during a simulation run. Before `sim.step()` is called, it's not set. After the simulation finishes but before the Python process ends, it remains accessible. If you exit `sim.step()` and re-enter (e.g., calling `sim.step(1)` in a loop), the extension is the same object throughout.

## Reading and writing fields from Python

The cleanest way to access fields is through `pywarpx.fields`:

```python
from pywarpx import fields

# Get the Ex field on level 0 as a wrapper
Ex_wrapper = fields.ExWrapper(level=0, include_ghosts=False)

# Get a numpy view of all the data on this rank
Ex_np = Ex_wrapper[...]                    # ellipsis returns the whole local array

# Or just a slice:
Ex_slice = Ex_wrapper[64:128, :, :]

# Modify in place:
Ex_np *= 0.5                                # halve all values
# Or write through the wrapper:
Ex_wrapper[64:128, :, :] = 0.0
```

Wrappers for all standard fields: `ExWrapper`, `EyWrapper`, `EzWrapper`, `BxWrapper`, `ByWrapper`, `BzWrapper`, `JxWrapper`, `JyWrapper`, `JzWrapper`, `RhoFPWrapper`, `PhiFPWrapper` (for ES), etc. Check `pywarpx.fields` for the current list.

For per-species charge densities (set up when you have multiple species):

```python
rho_e = fields.RhoFPWrapper(level=0, species='electrons')
```

### Access patterns

```python
@callbacks.callfromafterstep
def field_diagnostic():
    Ex = fields.ExWrapper()[...]            # full local array as numpy
    max_Ex = np.max(np.abs(Ex))
    print(f"Step {sim.extension.warpx.getistep(0)}: max|Ex| = {max_Ex}")
```

For data that spans multiple MPI ranks, each rank's wrapper returns only its local portion. Use MPI reductions (`mpi4py`) to combine:

```python
from mpi4py import MPI
comm = MPI.COMM_WORLD

@callbacks.callfromafterstep
def global_field_max():
    Ex = fields.ExWrapper()[...]
    local_max = np.max(np.abs(Ex))
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    if comm.rank == 0:
        print(f"Global max |Ex| = {global_max}")
```

### Modifying fields between steps

```python
@callbacks.callfromafterstep
def apply_external_pulse():
    """Add a transient external Ex pulse at z = z0."""
    t = sim.extension.warpx.gett_new(0)
    if 1e-9 < t < 1.5e-9:                   # only during this time window
        Ex = fields.ExWrapper()
        Ex_np = Ex[...]
        # ... compute the pulse profile ...
        Ex_np[:, :, z0_index] += pulse_value
```

**Caveat**: modifying fields directly is fine for diagnostics and applied perturbations, but it bypasses WarpX's internal consistency (e.g., divergence cleaning, NCI filters). For physically meaningful field modifications, prefer `picmi.AnalyticAppliedField` or `picmi.AnalyticInitialField` over manual modification.

## Reading and writing particles from Python

Particles are more complex than fields because they're partitioned per MPI rank and per box. There are two routes:

### High-level: `ParticleContainerWrapper`

```python
from pywarpx import particle_containers

elec_wrapper = particle_containers.ParticleContainerWrapper('electrons')

# Get all positions, momenta, weights as numpy arrays (concatenated across local boxes):
x_arr = elec_wrapper.get_particle_x()
y_arr = elec_wrapper.get_particle_y()
z_arr = elec_wrapper.get_particle_z()
ux_arr = elec_wrapper.get_particle_ux()
weights = elec_wrapper.get_particle_weight()

# Adding particles:
elec_wrapper.add_particles(
    x=np.array([1e-3, 2e-3]),
    y=np.array([0.0, 0.0]),
    z=np.array([0.0, 0.0]),
    ux=np.array([1e6, -1e6]),
    uy=np.array([0.0, 0.0]),
    uz=np.array([0.0, 0.0]),
    w=np.array([1e8, 1e8]),
    unique_particles=False,
)

# Removing/marking particles:
# (typically by setting an attribute or via the per-box iteration below)
```

The wrapper handles MPI: getters return only the local rank's particles; `add_particles` with `unique_particles=False` distributes the total across ranks. For analysis spanning all particles globally, use MPI reductions or `comm.gather`.

### Low-level: per-box iteration

For better performance (no concatenation overhead) and for direct access to per-particle SoA arrays:

```python
@callbacks.callfromafterstep
def per_box_processing():
    warpx = sim.extension.warpx
    multi_pc = warpx.multi_particle_container()
    pc = multi_pc.get_particle_container_from_name('electrons')

    for lvl in range(pc.finest_level + 1):
        for pti in pc.iterator(pc, level=lvl):
            # pti is a particle tile iterator
            soa = pti.soa()

            if sim.extension.Config.have_gpu:
                arrays = soa.to_cupy()
            else:
                arrays = soa.to_numpy()

            # arrays is a structure with .real and .int subfields
            # Real components, by SoA index:
            x_local = arrays.real[0]        # x
            y_local = arrays.real[1]        # y (in 3D) or z (in 2D)
            # ... etc.

            # Modify in place (this is a view, not a copy):
            ux = arrays.real[3]              # ux for typical SoA layout
            ux[:] *= 0.99                    # apply a damping

            # Iterating over particles within the tile is slow in Python
            # — use vectorized array operations
```

The SoA index mapping (which `.real[i]` is which quantity) is fixed by WarpX and depends on the dimensionality and the features enabled. Common layout for 3D Cartesian + no QED:
- `arrays.real[0]` = x
- `arrays.real[1]` = y
- `arrays.real[2]` = z
- `arrays.real[3]` = w (weight)
- `arrays.real[4]` = ux
- `arrays.real[5]` = uy
- `arrays.real[6]` = uz

Plus runtime attributes (ionization level, optical depth, etc.) appended after. Check the current source for the precise order; this varies across versions.

### Redistribution after position changes

If you modify particle positions (not just momenta), particles can end up in the wrong box. After your callback, call:

```python
pc.Redistribute()
```

This moves particles to the correct box owners and removes any marked-invalid particles. Skipping this can produce wrong results in the next step.

## The ParticleContainerWrapper convenience

`particle_containers.ParticleContainerWrapper` exposes a higher-level API:

```python
pc = particle_containers.ParticleContainerWrapper(species_name='electrons')

# Get arrays (concatenated across local boxes):
pc.get_particle_x()
pc.get_particle_y()
pc.get_particle_z()
pc.get_particle_ux()
pc.get_particle_uy()
pc.get_particle_uz()
pc.get_particle_weight()
pc.get_particle_arrays(comp_name, level=0)    # generic accessor by name

# Add particles:
pc.add_particles(x=..., y=..., z=..., ux=..., uy=..., uz=..., w=...,
                 unique_particles=True, **kwargs)
# kwargs sets runtime attributes; missing ones default to 0

# Add/remove runtime attributes:
pc.add_real_comp(name)
pc.add_int_comp(name)

# Total number of particles on this rank:
pc.get_particle_count(local=True)        # local
pc.get_particle_count(local=False)       # globally reduced (requires MPI)

# Charge density deposition:
pc.deposit_charge_density(level=0)        # forces a rho deposition update
```

This is convenient for occasional access. For per-step hot-loop work, the lower-level per-box iteration is faster.

## pyAMReX: direct access to AMReX data structures

The underlying AMReX bindings live in `pywarpx.libwarpx.amr` (equivalently `amrex.space{1,2,3}d` for the appropriate dimension):

```python
from pywarpx import libwarpx
amr = libwarpx.amr            # equivalent to: import amrex.space3d as amr (for 3D)
```

Useful APIs:
- `amr.ParallelDescriptor` — MPI rank info: `MyProc()`, `NProcs()`, `IOProcessor()`.
- `amr.MultiFab` — the underlying grid container.
- `amr.ParticleContainer_*` — underlying particle containers.
- `amr.Geometry` — domain geometry.
- `amr.Box`, `amr.BoxArray`, `amr.DistributionMapping` — index space and partitioning.

For typical PICMI work, you rarely need pyAMReX directly — `pywarpx.fields` and `pywarpx.particle_containers` wrap most use cases. Drop to pyAMReX when:
- You need an operation not exposed by the wrappers (e.g., `Redistribute` on a specific container).
- You're writing a custom solver and need raw `MultiFab` access.
- You're working with a level or component that the wrappers don't expose.

The pyAMReX docs at `https://pyamrex.readthedocs.io/` are the canonical reference.

## GPU vs CPU: numpy and cupy

For GPU-built WarpX:
- Field wrapper `[...]` returns a **cupy** array (zero-copy view into GPU memory).
- `pti.soa().to_cupy()` returns cupy arrays.
- Numpy operations on these arrays trigger CPU↔GPU transfers (slow).
- Cupy operations (`cupy.sum`, `cupy.sqrt`, etc.) execute on the GPU.

For CPU-built WarpX:
- All wrappers return numpy arrays directly.
- `to_numpy()` is a no-op view.

For portable code:

```python
import numpy as np
from pywarpx import callbacks

@callbacks.callfromafterstep
def portable_callback():
    have_gpu = sim.extension.Config.have_gpu
    xp = cupy if have_gpu else np    # array module dispatch

    # ...
    pti_arrays = pti.soa().to_cupy() if have_gpu else pti.soa().to_numpy()
    ux = pti_arrays.real[4]
    max_ux = xp.max(xp.abs(ux))     # works on either
```

This dispatch pattern is the idiomatic way to write a single callback that runs on both backends.

**Don't** transfer GPU arrays to CPU for printing/logging in every step — that's a CPU↔GPU roundtrip and will dominate runtime. Reduce on the GPU first, then transfer the scalar.

## Coupling to ML frameworks

A common use case: use a trained ML model to drive a custom term in the simulation, or use the simulation to generate training data.

### Pattern A: ML model as a custom Poisson solver

For an electrostatic simulation, install a callback at `poissonsolver` to replace the standard solve:

```python
import torch
model = load_my_potential_predictor()       # trained PyTorch model

@callbacks.installcallback('poissonsolver')
def ml_poisson_solver():
    """Replace the standard Poisson solve with a neural-net prediction."""
    # Read the charge density
    rho = fields.RhoFPWrapper()[...]

    # Convert to model input (on GPU if available)
    if sim.extension.Config.have_gpu:
        rho_torch = torch.as_tensor(rho)    # zero-copy from cupy via __cuda_array_interface__
    else:
        rho_torch = torch.from_numpy(rho)

    # Predict potential
    with torch.no_grad():
        phi_pred = model(rho_torch)

    # Write back to WarpX
    phi = fields.PhiFPWrapper()
    if sim.extension.Config.have_gpu:
        phi[...] = phi_pred.cpu().numpy()    # or zero-copy if framework supports
    else:
        phi[...] = phi_pred.numpy()
```

This pattern lets you trade physical accuracy for speed — useful when the Poisson solve dominates cost and an approximate solution is acceptable for the use case (e.g., evolutionary algorithm, design optimization, real-time control).

### Pattern B: ML model as an extra force or source

Install at `beforeEsolve` or `afterstep`:

```python
@callbacks.callfromafterstep
def apply_ml_correction():
    """Use an ML model to add a correction to particle momenta."""
    elec_wrapper = particle_containers.ParticleContainerWrapper('electrons')
    x = elec_wrapper.get_particle_x()
    y = elec_wrapper.get_particle_y()
    z = elec_wrapper.get_particle_z()

    features = build_features(x, y, z)
    corrections = model(features)            # shape: (N, 3) for du = (dux, duy, duz)

    # Apply correction to all particles (this is per-rank; corrections must be local)
    # Lower-level access needed because ParticleContainerWrapper doesn't expose write methods
    for pti in elec_wrapper.iterator():
        arrays = pti.soa().to_cupy() if have_gpu else pti.soa().to_numpy()
        # Find slices corresponding to this tile's particles in the corrections array
        # ... apply ...
```

This is more involved because writing requires the per-box iteration; the `ParticleContainerWrapper` getters return concatenated arrays but writing back to specific particles requires care.

### Pattern C: Simulation as training data generator

```python
@callbacks.callfromafterstep
def collect_training_data():
    """Every 100 steps, sample a snapshot of (state, action, reward) for RL training."""
    step = sim.extension.warpx.getistep(0)
    if step % 100 == 0:
        state = build_state_vector(sim)        # e.g., density field, particle moments
        action = current_action(sim)           # the control we're testing
        reward = compute_reward(sim)
        training_buffer.append((state, action, reward))
```

Pair with periodic checkpoints so you can resume training across multiple simulation runs.

### Caveats for ML coupling

- **Data transfer overhead**: GPU↔GPU within the same device is fast; CPU↔GPU is slow. Ideally, ML model and WarpX are both on the same GPU. Use cupy ↔ PyTorch/JAX zero-copy interop (via `__cuda_array_interface__` or `__dlpack__`).
- **Synchronization**: when WarpX is running asynchronous GPU ops and your callback queries the data, you may need to synchronize. `cupy.cuda.runtime.deviceSynchronize()` is the heavy hammer.
- **Reproducibility**: ML models often introduce nondeterminism via dropout, batch norm in training mode, atomic operations, etc. Switch to eval mode (`model.eval()`) and use deterministic algorithms if reproducibility matters.
- **Performance**: a Python callback at every step that takes 10 ms adds 10 ms × N_steps to runtime. For 100,000-step simulations, that's 1000 seconds of overhead. Profile early.

## Replacing parts of the PIC loop with Python

Beyond ML, callbacks can implement custom physics directly in Python. Some examples:

### Custom collision model

```python
@callbacks.callfromaftercollisions
def my_custom_attachment():
    """Apply a custom electron-attachment process between simulated species."""
    # Per-cell attachment rate computed from density:
    rho_e = fields.RhoFPWrapper(species='electrons')[...]
    rho_n = fields.RhoFPWrapper(species='neutrals')[...]    # if present
    rate = k_attach * rho_e * rho_n
    # ... apply by marking some electrons for deletion ...
```

### Custom particle injection

```python
@callbacks.callfromparticleinjection
def inject_at_emitter():
    """Inject electrons from an emitter surface with field-dependent rate."""
    E_at_surface = compute_E_at_surface()
    rate = fowler_nordheim_current(E_at_surface, work_function=4.5)  # A/m^2
    n_to_inject = rate * dt * area / q_e

    elec_wrapper = particle_containers.ParticleContainerWrapper('electrons')
    elec_wrapper.add_particles(
        x=sample_emitter_x(n_to_inject),
        y=sample_emitter_y(n_to_inject),
        z=np.full(n_to_inject, surface_z),
        ux=np.zeros(n_to_inject),
        uy=np.zeros(n_to_inject),
        uz=sample_emitter_uz(n_to_inject),
        w=np.full(n_to_inject, weight),
        unique_particles=False,
    )
```

### Custom external field

For fields that depend on accumulated simulation state (not just an analytic expression), use a callback:

```python
phase = 0.0

@callbacks.callfromafterstep
def update_feedback_field():
    """Apply an external field whose amplitude depends on accumulated charge at boundary."""
    global phase
    boundary_charge = compute_boundary_charge()
    phase += k_feedback * boundary_charge * dt
    # Apply via direct field write or via setting a parser parameter:
    sim.extension.warpx.set_external_E_phase(phase)    # hypothetical, depends on API
```

The cleanest mechanism is via `picmi.AnalyticAppliedField` with parser expressions, but those are static (can't depend on simulation state). For state-dependent fields, callbacks are the right tool.

## Real-time diagnostics

Beyond what `picmi.FieldDiagnostic` and `picmi.ParticleDiagnostic` produce on disk, callbacks let you compute and visualize diagnostics during the run:

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
line, = ax.plot([], [])

step_times = []
mean_ne = []

@callbacks.callfromafterstep
def live_diagnostic():
    if sim.extension.warpx.getistep(0) % 100 != 0:
        return

    # Sample only on rank 0
    if MPI.COMM_WORLD.rank != 0:
        return

    rho_e = fields.RhoFPWrapper(species='electrons')[...]
    n_e_avg = np.mean(rho_e) / picmi.constants.q_e

    step_times.append(sim.extension.warpx.gett_new(0))
    mean_ne.append(n_e_avg)

    line.set_data(step_times, mean_ne)
    ax.relim()
    ax.autoscale_view()
    plt.pause(0.001)
```

For HPC runs without an interactive display, the equivalent pattern is to write derived quantities to a small ASCII file each step (or every N steps) and tail it during the run.

## Patterns and pitfalls

### Globals and closures

Callback functions can use globals or be defined with closures over outer-scope variables:

```python
my_state = {'count': 0}

@callbacks.callfromafterstep
def use_outer_state():
    my_state['count'] += 1
    # ... do things with my_state ...
```

This is idiomatic and fine. Just be aware that the callback persists for the simulation lifetime — closures hold their captured variables in memory.

### MPI rank-dependent callbacks

Each callback fires on every MPI rank. If you want a callback to do work only on rank 0:

```python
from mpi4py import MPI
comm = MPI.COMM_WORLD

@callbacks.callfromafterstep
def rank0_only():
    if comm.rank != 0:
        return
    # Rank-0 logic here
```

For callbacks that need to collect data from all ranks, gather/reduce before doing the work.

### MPI deadlock from inconsistent callbacks

If the callback does MPI operations on some ranks but not others, you'll deadlock. The two patterns to avoid:

```python
# BAD: rank 0 does an allreduce, others don't
@callbacks.callfromafterstep
def bad_rank0():
    if comm.rank != 0:
        return
    global_max = comm.allreduce(local_max, op=MPI.MAX)   # other ranks never call this
    # → deadlock
```

```python
# GOOD: all ranks call the allreduce, only rank 0 prints
@callbacks.callfromafterstep
def good_pattern():
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    if comm.rank == 0:
        print(global_max)
```

### Callback order matters

Multiple callbacks at one hook fire in installation order. If one callback modifies particle data and another reads it, the order determines what the second sees. Either:
- Install callbacks in a deliberate order.
- Make callbacks independent of each other.
- Use explicit synchronization (one callback that orchestrates several).

### Don't mutate fields during the PIC loop

Modifying fields between `beforeEsolve` and `afterEsolve` causes WarpX to overwrite your changes during the solve. The safe windows for field modification:
- Between steps (`afterstep` of step N, before `beforestep` of step N+1).
- Before the initial solve (`beforeInitEsolve`) for setting initial state.
- For applied fields, prefer `picmi.AnalyticAppliedField` over manual modification.

### Particles modified out-of-box

If you change particle positions in a callback such that particles end up outside their current box, call `pc.Redistribute()` before returning. Without this, the next step's field-gather and current-deposition will be wrong for those particles.

### Closures over `sim` may not pickle for checkpointing

If you want to restart-resume a simulation that uses callbacks, the callbacks themselves are not part of the WarpX checkpoint — the Python state isn't saved. After restart, your script re-runs from the top, re-installs the callbacks, and continues. Callbacks that depend on accumulated Python state (e.g., a counter, a learned ML model) need to save and load that state separately.

### Performance: profile your callbacks

A 10 ms callback at every step is 1000 s of overhead for 100,000 steps. Common culprits:
- CPU↔GPU data transfers for every print/log.
- Python-level loops over particles (always use vectorized numpy/cupy).
- Allocating big arrays in the callback (allocate once outside, reuse).
- `comm.allreduce` for scalar reductions when a less synchronous primitive would do.

Use `cProfile` or `line_profiler` on a representative short run before scaling up.
