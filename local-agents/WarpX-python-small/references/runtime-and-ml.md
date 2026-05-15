# Runtime extension: callbacks, pyAMReX, ML coupling

How to extend WarpX from Python while it runs: install callbacks, access simulation state, modify fields/particles in situ, couple to ML frameworks. Only works in **interactive mode** (`sim.step()`), not preprocessor mode (`write_input_file`).

## The mental model

A PICMI script in interactive mode runs the C++ PIC loop. To inject Python logic, install **callbacks** — Python functions WarpX invokes at well-defined points in each step:

```python
from pywarpx import picmi, callbacks

sim = picmi.Simulation(...)
# ... configure ...

@callbacks.callfromafterstep
def my_callback():
    # runs after every step
    ...

sim.step()
```

The C++ loop continues uninterrupted between hooks. Callback overhead is just a Python function call (microseconds) — the work inside dominates. Use vectorized numpy/cupy to scale to billions of particles.

## Callback hook points

| Hook | Fires |
|---|---|
| `beforeInitEsolve` | Before initial E-field solve at start |
| `afterInitEsolve` | After initial E-field solve at start |
| `afterinit` | After all initialization |
| `beforestep` | Start of each step |
| `beforecollisions` / `aftercollisions` | Around collision operator |
| `beforeEsolve` / `afterEsolve` | Around E-field solve |
| `poissonsolver` | **Replaces** default Poisson solve (ES only) |
| `beforedeposition` / `afterdeposition` | Around current/charge deposition |
| `particleinjection` | At particle injection time |
| `particleloader` | When standard loader runs |
| `particlescraper` | After particle BCs, before lost particles processed |
| `afterstep` | End of each step — **most common hook** |
| `afterdiagnostics` | After diagnostic output |
| `afterrestart` | After restart from checkpoint |
| `oncheckpointsignal` | On checkpoint signal |

## Installing callbacks

Two equivalent styles:

```python
from pywarpx.callbacks import callfromafterstep, installafterstep

# Decorator style (cleaner):
@callfromafterstep
def my_diagnostic():
    ...

# Explicit:
def my_other():
    ...
installafterstep(my_other)
```

Decorators for other hooks: `callfrombeforestep`, `callfromafterEsolve`, `callfromafterdeposition`, `callfromparticlescraper`, `callfromparticleloader`, `callfromparticleinjection`, `callfromafterdiagnostics`, `callfromafterinit`, etc.

Generic form: `installcallback('afterstep', my_func)` or `uninstallcallback('afterstep', my_func)`.

### Multiple callbacks at one hook

Fire in installation order. Don't rely on order across callbacks — make them independent.

### Uninstall mid-simulation

```python
installafterstep(transient_diag)
sim.step(1000)
uninstallafterstep(transient_diag)    # off for the rest of the run
sim.step(10000)
```

## Accessing simulation state

Inside a callback (or after `sim.step()` returns):

```python
@callbacks.callfromafterstep
def my_callback():
    warpx = sim.extension.warpx          # C++ WarpX object
    Config = sim.extension.Config         # compile-time config

    step = warpx.getistep(0)              # step number on level 0
    t = warpx.gett_new(0)                  # physical time, s
    dt = warpx.getdt(0)                    # time step size
    on_gpu = Config.have_gpu               # bool

    geom = warpx.Geom(0)                   # AMReX Geometry
    prob_lo = geom.ProbLoArray()
    prob_hi = geom.ProbHiArray()
    dx = geom.CellSizeArray()
```

`sim.extension` only exists during/after a simulation run.

## Reading and writing fields

Use `pywarpx.fields`:

```python
from pywarpx import fields

# Get Ex field on level 0:
Ex_wrapper = fields.ExWrapper(level=0, include_ghosts=False)

# Full local array as numpy/cupy view:
Ex = Ex_wrapper[...]

# Slice:
Ex_slice = Ex_wrapper[64:128, :, :]

# Modify in place:
Ex[...] *= 0.5

# Or write through wrapper:
Ex_wrapper[64:128, :, :] = 0.0
```

Wrappers: `ExWrapper`, `EyWrapper`, `EzWrapper`, `BxWrapper`, `ByWrapper`, `BzWrapper`, `JxWrapper`, `JyWrapper`, `JzWrapper`, `RhoFPWrapper`, `PhiFPWrapper` (ES only). Per-species: `fields.RhoFPWrapper(level=0, species='electrons')`.

Each rank gets its local portion. For global stats:

```python
from mpi4py import MPI
comm = MPI.COMM_WORLD

@callbacks.callfromafterstep
def global_max():
    Ex = fields.ExWrapper()[...]
    local_max = np.max(np.abs(Ex))
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    if comm.rank == 0:
        print(f"Global max |Ex| = {global_max}")
```

**Caveat**: directly modifying fields bypasses internal consistency (divergence cleaning, NCI filters). For physically-motivated external fields, prefer `picmi.AnalyticAppliedField` over manual modification.

## Reading and writing particles

### High-level: ParticleContainerWrapper

```python
from pywarpx import particle_containers

pc = particle_containers.ParticleContainerWrapper('electrons')

# Get arrays (concatenated across local boxes):
x = pc.get_particle_x()
y = pc.get_particle_y()
z = pc.get_particle_z()
ux = pc.get_particle_ux()
uy = pc.get_particle_uy()
uz = pc.get_particle_uz()
w = pc.get_particle_weight()

# Add particles:
pc.add_particles(
    x=np.array([1e-3, 2e-3]),
    y=np.array([0.0, 0.0]),
    z=np.array([0.0, 0.0]),
    ux=np.array([1e6, -1e6]),
    uy=np.zeros(2),
    uz=np.zeros(2),
    w=np.array([1e8, 1e8]),
    unique_particles=False,             # distribute total across MPI ranks
)

# Count:
n_local = pc.get_particle_count(local=True)
n_global = pc.get_particle_count(local=False)   # requires MPI reduction
```

`unique_particles`:
- `True` (default): each rank adds its own. Total = N × n_ranks.
- `False`: total N distributed across ranks.

### Low-level: per-box iteration (faster for hot loops)

```python
@callbacks.callfromafterstep
def per_box_processing():
    warpx = sim.extension.warpx
    multi_pc = warpx.multi_particle_container()
    pc = multi_pc.get_particle_container_from_name('electrons')

    for lvl in range(pc.finest_level + 1):
        for pti in pc.iterator(pc, level=lvl):
            if sim.extension.Config.have_gpu:
                arrays = pti.soa().to_cupy()
            else:
                arrays = pti.soa().to_numpy()

            # SoA index layout for 3D Cartesian (varies — check current source):
            # arrays.real[0] = x, [1] = y, [2] = z, [3] = w
            # arrays.real[4] = ux, [5] = uy, [6] = uz

            ux = arrays.real[4]
            ux[:] *= 0.99    # in-place damping (this is a view, not a copy)
```

After modifying positions, call `pc.Redistribute()` to move particles to correct box owners.

## GPU vs CPU dispatch

```python
import numpy as np

@callbacks.callfromafterstep
def portable_callback():
    have_gpu = sim.extension.Config.have_gpu

    if have_gpu:
        import cupy as cp
        xp = cp
    else:
        xp = np

    Ex = fields.ExWrapper()[...]    # numpy on CPU, cupy on GPU
    max_E = xp.max(xp.abs(Ex))      # works on either backend
```

For GPU-built WarpX:
- Field wrappers return **cupy** arrays (zero-copy GPU view).
- `pti.soa().to_cupy()` returns cupy.
- numpy operations on these → CPU↔GPU transfer (slow).
- cupy operations stay on GPU.

For CPU-built WarpX:
- All wrappers return numpy.
- `to_numpy()` is a no-op view.

**Don't transfer GPU arrays to CPU every step** for printing — reduce on GPU first, transfer the scalar.

## ML coupling patterns

### Pattern A: ML model as Poisson solver

```python
import torch
model = load_my_potential_predictor()

@callbacks.installcallback('poissonsolver')
def ml_poisson_solver():
    """Replace standard Poisson solve with neural-net prediction."""
    rho = fields.RhoFPWrapper()[...]

    if sim.extension.Config.have_gpu:
        rho_torch = torch.as_tensor(rho)         # zero-copy from cupy
    else:
        rho_torch = torch.from_numpy(rho)

    with torch.no_grad():
        phi_pred = model(rho_torch)

    phi = fields.PhiFPWrapper()
    if sim.extension.Config.have_gpu:
        phi[...] = phi_pred.cpu().numpy()
    else:
        phi[...] = phi_pred.numpy()
```

Trades physical accuracy for speed. Useful for evolutionary algorithms, real-time control, design optimization.

### Pattern B: ML correction to particle momenta

```python
@callbacks.callfromafterstep
def apply_ml_correction():
    pc = particle_containers.ParticleContainerWrapper('electrons')
    x = pc.get_particle_x()
    y = pc.get_particle_y()
    z = pc.get_particle_z()

    features = build_features(x, y, z)
    corrections = model(features)        # shape (N, 3)

    # Apply via per-box iteration (writing back requires it)
    # ... see per-box pattern above
```

### Pattern C: Simulation as RL training data generator

```python
training_buffer = []

@callbacks.callfromafterstep
def collect_data():
    step = sim.extension.warpx.getistep(0)
    if step % 100 == 0:
        state = build_state_vector(sim)
        action = current_action(sim)
        reward = compute_reward(sim)
        training_buffer.append((state, action, reward))
```

### Caveats for ML coupling

- **Data transfer**: GPU↔GPU same device is fast; CPU↔GPU slow. Put WarpX and ML model on the same GPU.
- **Zero-copy interop**: cupy ↔ PyTorch via `__cuda_array_interface__`; cupy ↔ JAX via `__dlpack__`.
- **Synchronization**: GPU ops may be async; `cupy.cuda.runtime.deviceSynchronize()` if needed.
- **Reproducibility**: switch model to `eval()` mode; ML training-mode ops (dropout, batchnorm) introduce nondeterminism.

## Real-time diagnostics

```python
import matplotlib.pyplot as plt

step_times = []
mean_ne = []

@callbacks.callfromafterstep
def live_plot():
    step = sim.extension.warpx.getistep(0)
    if step % 100 != 0 or MPI.COMM_WORLD.rank != 0:
        return

    rho_e = fields.RhoFPWrapper(species='electrons')[...]
    n_e = -np.mean(rho_e) / picmi.constants.q_e
    step_times.append(sim.extension.warpx.gett_new(0))
    mean_ne.append(n_e)
    # update plot...
```

For HPC without display: write derived quantities to a small ASCII file and tail during the run.

## MPI pitfalls

### MPI deadlock from inconsistent callbacks

```python
# BAD — rank 0 calls allreduce, others don't → deadlock
@callbacks.callfromafterstep
def bad():
    if comm.rank != 0:
        return
    global_max = comm.allreduce(local_max, op=MPI.MAX)
```

```python
# GOOD — all ranks call allreduce, only rank 0 prints
@callbacks.callfromafterstep
def good():
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    if comm.rank == 0:
        print(global_max)
```

### MPI initialization order

Import `mpi4py.MPI` **before** `pywarpx` if using both. pywarpx initializes MPI on import otherwise.

```python
# Correct:
from mpi4py import MPI               # initializes MPI
from pywarpx import picmi             # uses already-initialized MPI
```

## Restart and callback state

WarpX checkpoints save fields + particles + RNG. They do NOT save Python state — callbacks, accumulator variables, ML model parameters.

After restart, your script re-runs from the top and re-installs callbacks. Save and load any stateful Python objects (ML models, buffers) separately.

## Common pitfalls

- **Capturing `m_*` in lambdas (C++ context)** — doesn't apply here; Python closures are fine.
- **Forgetting to call `pc.Redistribute()`** after modifying particle positions → particles in wrong box → wrong gather/deposit next step.
- **Writing fields between `beforeEsolve` and `afterEsolve`** → solver overwrites your changes.
- **Python-level loop over particles** → 100x slower than vectorized numpy. Always vectorize.
- **Inconsistent collective MPI calls across ranks in callbacks** → deadlock.
- **Modifying particle count in a hot callback** without batching → expensive memory reallocations.

## Performance

A 10ms callback × 100,000 steps = 1000s of overhead. Profile callbacks before scaling:

```python
import cProfile
cProfile.run('sim.step(100)', sort='cumtime')
```

Common culprits:
- CPU↔GPU transfers every step
- Python loops instead of vectorized operations
- Allocating big arrays inside the callback (allocate once, reuse)
- Synchronous `comm.allreduce` for scalars when a less-blocking primitive would do
