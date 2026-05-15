# Runtime extension and ML coupling

How to inject Python logic into WarpX while it runs. Only works in **interactive mode** (`sim.step()`), not preprocessor mode.

## The mental model

A PICMI script in interactive mode runs the C++ PIC loop. To run Python at specific points, install **callbacks**:

```python
from pywarpx import picmi, callbacks

sim = picmi.Simulation(...)

@callbacks.callfromafterstep
def my_callback():
    # runs after every step
    ...

sim.step()
```

The C++ loop continues uninterrupted between hooks. Callback overhead is just a Python call. Use vectorized numpy/cupy.

## Callback hook points

| Hook | Fires |
|---|---|
| `afterinit` | After all initialization |
| `beforestep` / `afterstep` | Around each step |
| `beforeEsolve` / `afterEsolve` | Around E-field solve |
| `beforedeposition` / `afterdeposition` | Around current/charge deposition |
| `beforecollisions` / `aftercollisions` | Around collisions |
| `poissonsolver` | **Replaces** default Poisson solve (ES only) |
| `particleinjection` | At particle injection time |
| `particlescraper` | After particle BCs |
| `afterdiagnostics` | After diagnostic output |

`afterstep` is by far the most common.

## Install callbacks: decorator or explicit

```python
from pywarpx.callbacks import callfromafterstep, installafterstep, uninstallafterstep

# Decorator (cleaner):
@callfromafterstep
def my_diag():
    ...

# Explicit:
def my_other():
    ...
installafterstep(my_other)

# Generic:
installcallback('afterstep', my_func)

# Uninstall:
uninstallafterstep(my_diag)
```

Other decorators: `callfrombeforestep`, `callfromafterEsolve`, `callfromafterdeposition`, `callfromparticleinjection`, `callfromafterdiagnostics`, etc.

## Access simulation state

Inside any callback (or after `sim.step()` returns):

```python
@callbacks.callfromafterstep
def my_callback():
    warpx = sim.extension.warpx          # C++ WarpX object
    Config = sim.extension.Config         # compile-time config

    step = warpx.getistep(0)              # step number
    t = warpx.gett_new(0)                  # physical time (s)
    dt = warpx.getdt(0)                    # time step size
    on_gpu = Config.have_gpu               # bool
```

## Read/write fields

```python
from pywarpx import fields

# Get Ex on level 0:
Ex_wrapper = fields.ExWrapper(level=0, include_ghosts=False)

# Full local array (numpy on CPU, cupy on GPU):
Ex = Ex_wrapper[...]

# Modify in place:
Ex[...] *= 0.5

# Write through wrapper:
Ex_wrapper[64:128, :, :] = 0.0
```

Wrappers: `ExWrapper`, `EyWrapper`, `EzWrapper`, `BxWrapper`, `ByWrapper`, `BzWrapper`, `JxWrapper`, `JyWrapper`, `JzWrapper`, `RhoFPWrapper`, `PhiFPWrapper` (ES only).

Per-species: `fields.RhoFPWrapper(level=0, species='electrons')`.

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
        print(f"Max |Ex| = {global_max}")
```

## Read/write particles

### High-level: ParticleContainerWrapper

```python
from pywarpx import particle_containers

pc = particle_containers.ParticleContainerWrapper('electrons')

# Get arrays (local to this rank):
x = pc.get_particle_x()
y = pc.get_particle_y()
z = pc.get_particle_z()
ux = pc.get_particle_ux()
w = pc.get_particle_weight()

# Add particles:
pc.add_particles(
    x=np.array([1e-3, 2e-3]),
    y=np.array([0.0, 0.0]),
    z=np.array([0.0, 0.0]),
    ux=np.array([1e6, -1e6]),
    uy=np.zeros(2), uz=np.zeros(2),
    w=np.array([1e8, 1e8]),
    unique_particles=False,     # distribute total across MPI ranks
)
```

`unique_particles`:
- `True` (default): each rank adds its own → total = N × n_ranks
- `False`: total N distributed across ranks

### After modifying positions: Redistribute

```python
@callbacks.callfromafterstep
def move_particles():
    # ... modify particle positions ...
    pc.Redistribute()          # fix box ownership
```

Skip this → particles in wrong box → next step's gather/deposit wrong.

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

    Ex = fields.ExWrapper()[...]      # numpy on CPU, cupy on GPU
    max_E = xp.max(xp.abs(Ex))         # works on either
```

For GPU-built WarpX:
- Field wrappers return **cupy** arrays (zero-copy GPU view)
- numpy ops on these → CPU↔GPU transfer (slow)
- Use cupy ops to stay on GPU

For CPU-built WarpX:
- All wrappers return numpy

**Don't transfer GPU→CPU every step** for printing. Reduce on GPU, then transfer the scalar.

## ML coupling patterns

### Pattern A: ML model replaces Poisson solver

```python
import torch
model = load_my_potential_predictor()

@callbacks.installcallback('poissonsolver')
def ml_poisson():
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

### Pattern B: ML correction to particle momenta

```python
@callbacks.callfromafterstep
def apply_ml_correction():
    pc = particle_containers.ParticleContainerWrapper('electrons')
    x = pc.get_particle_x()
    y = pc.get_particle_y()
    z = pc.get_particle_z()

    features = build_features(x, y, z)
    corrections = model(features)            # (N, 3)
    # ... apply via per-box iteration to modify in place ...
```

### Pattern C: Simulation generates RL training data

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

### ML coupling caveats

| Concern | Fix |
|---|---|
| CPU↔GPU transfer overhead | Put WarpX + ML model on same GPU |
| cupy ↔ PyTorch interop | Use `torch.as_tensor()` (zero-copy via `__cuda_array_interface__`) |
| cupy ↔ JAX interop | Use `__dlpack__` |
| GPU async ops | `cupy.cuda.runtime.deviceSynchronize()` if needed |
| ML non-determinism | Use `model.eval()` mode |

## MPI pitfalls

### Inconsistent callbacks deadlock

```python
# BAD — only rank 0 calls allreduce → other ranks hang
@callbacks.callfromafterstep
def bad():
    if comm.rank != 0:
        return
    global_max = comm.allreduce(local_max, op=MPI.MAX)   # DEADLOCK
```

```python
# GOOD — all ranks call allreduce, only rank 0 uses result
@callbacks.callfromafterstep
def good():
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    if comm.rank == 0:
        print(global_max)
```

### MPI init order

```python
# Correct order:
from mpi4py import MPI                   # initializes MPI first
from pywarpx import picmi                 # uses already-initialized MPI
```

If pywarpx imports first, MPI may be initialized in a way that conflicts with mpi4py.

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

For HPC without display: write derived quantities to an ASCII file and tail during the run.

## Restart and callback state

WarpX checkpoints save fields + particles + RNG. They do NOT save Python state (callbacks, accumulators, ML model parameters).

After restart, your script re-runs from the top and re-installs callbacks. Save Python state separately:

```python
import pickle

# Save:
@callbacks.callfromafterdiagnostics
def save_python_state():
    with open(f'./checkpoints/python_state.pkl', 'wb') as f:
        pickle.dump({'buffer': training_buffer, ...}, f)

# Load (at script top):
if os.path.exists('./checkpoints/python_state.pkl'):
    with open('./checkpoints/python_state.pkl', 'rb') as f:
        state = pickle.load(f)
```

## Common pitfalls

| Pitfall | Fix |
|---|---|
| Modifying positions without `Redistribute()` | Always call `pc.Redistribute()` after position changes |
| Writing fields between `beforeEsolve` and `afterEsolve` | Solver overwrites — use `afterstep` instead |
| Python loop over particles | Always vectorize with numpy/cupy |
| Inconsistent MPI calls across ranks | Every rank must do same collective ops |
| Allocating big arrays each call | Allocate once outside callback, reuse |
| GPU sync overhead every step | Reduce on GPU, transfer only scalars |

## Performance

A 10ms callback × 100,000 steps = 1000s overhead. Profile early:

```python
import cProfile
cProfile.run('sim.step(100)', sort='cumtime')
```
