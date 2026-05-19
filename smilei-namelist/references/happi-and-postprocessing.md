# happi and post-processing

This reference covers `happi`, the Python module bundled with Smilei for opening simulations and post-processing diagnostic output. Every diagnostic from `diagnostics.md` is read via happi here. For unit conversion semantics see `basics-and-units.md` §6.

## Table of contents

1. Installing and importing happi
2. Opening a simulation
3. Accessing namelist parameters
4. Reading diagnostics — overview
5. `Scalar` access
6. `Field` access
7. `Probe` access
8. `ParticleBinning` access
9. `Screen`, `RadiationSpectrum` access
10. `TrackParticles` access
11. `NewParticles` access
12. `Performances` access
13. Unit conversion with Pint
14. Plotting, animating, exporting
15. Multi-simulation comparison
16. Common mistakes

---

## 1. Installing and importing happi

happi is built alongside Smilei:

```bash
cd <smilei-directory>
make happi
```

This installs the `happi` module into the user-site Python directory. Once installed:

```python
import happi
```

If `make happi` was never run or Smilei was moved, the import fails. The fix is to re-run `make happi` in the current Smilei directory or set `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/smilei/happi:$PYTHONPATH
```

happi has no dependencies beyond NumPy, h5py, and matplotlib. Optional: `pint` for SI unit handling (recommended).

## 2. Opening a simulation

```python
S = happi.Open("/path/to/simulation_directory")
```

The path should be the directory containing the `smilei.py` (auto-generated input dump), `Scalars.txt`, and the `*.h5` files. Common patterns:

```python
# Open relative to current directory
S = happi.Open("../results")

# Open from a specific directory
S = happi.Open("/scratch/runs/exp_5/")

# Open with explicit reference frequency (if namelist didn't set it)
S = happi.Open("/path/to/sim", reference_angular_frequency_SI = 2*3.14159*3e8/248e-9)

# Pass pint=True if you want Pint-based unit handling
S = happi.Open("/path/to/sim", pint = True)

# Multiple simulations at once
S = happi.Open(["sim_a", "sim_b", "sim_c"])
```

If the simulation didn't finish, happi still opens it and gives access to whatever output was written. Useful for monitoring runs in progress.

## 3. Accessing namelist parameters

The full namelist is preserved as an attribute:

```python
S = happi.Open("path/to/sim")

# Direct access to anything in the namelist
print(S.namelist.Main.timestep)               # the timestep
print(S.namelist.Main.geometry)               # "2Dcartesian", etc.
print(S.namelist.Main.cell_length)            # [0.39, 0.39]
print(S.namelist.Main.reference_angular_frequency_SI)
print(S.namelist.Species[0].name)             # first defined species
print([s.name for s in S.namelist.Species])   # all species

# Custom variables you defined in the namelist also survive
lambda_SI = S.namelist.lambda_SI              # if you defined it
```

This is invaluable for parameter scans and reproducible analysis — your post-processing script doesn't need to be hand-edited per simulation.

## 4. Reading diagnostics — overview

Each diagnostic type has a corresponding accessor:

| Diagnostic block | Accessor |
|---|---|
| `DiagScalar` | `S.Scalar(name)` |
| `DiagFields` | `S.Field(index, field_name)` |
| `DiagProbe` | `S.Probe(index, field_name)` |
| `DiagParticleBinning` | `S.ParticleBinning(index)` |
| `DiagScreen` | `S.Screen(index)` |
| `DiagRadiationSpectrum` | `S.RadiationSpectrum(index)` |
| `DiagTrackParticles` | `S.TrackParticles(species_name, attributes)` |
| `DiagNewParticles` | `S.NewParticles(species_name)` |
| `DiagPerformances` | `S.Performances()` |

Indices match the order of blocks in the namelist (0-indexed).

All accessors return a "result" object with common methods:

| Method | Purpose |
|---|---|
| `.getData()` | List of arrays (one per timestep) |
| `.getTimes()` / `.getTimesteps()` | Time/step values |
| `.getAxis(name)` | Coordinates along a named axis |
| `.plot()` | Static plot of first timestep (or specified) |
| `.slide()` | Interactive slider through timesteps |
| `.animate()` | Animated GIF/MP4 |
| `.streak()` | 2D heatmap with time on one axis |
| `.toVTK()` | Export to VTK for ParaView/VisIt |
| `.set(...)` | Modify plot parameters in place |

## 5. `Scalar` access

```python
# List available scalars
print(S.getScalars())

# Get the kinetic energy time series
Ukin = S.Scalar("Ukin")
times = Ukin.getTimes()                       # array of times
data  = Ukin.getData()                        # array of values

# Plot
Ukin.plot()

# Multiple on one axis (matplotlib)
import matplotlib.pyplot as plt
S.Scalar("Ukin").plot(figure=1)
S.Scalar("Uelm").plot(figure=1)
S.Scalar("Utot").plot(figure=1)
plt.show()
```

`DiagScalar` data lives in `Scalars.txt` and is parsed by happi on demand. Cheap.

## 6. `Field` access

```python
# What's available?
print(S.fieldInfo())

# Get Ex from the first DiagFields block, at all timesteps
Ex = S.Field(0, "Ex")
data = Ex.getData()                           # list of 2D arrays
times = Ex.getTimes()
x = Ex.getAxis("x")                           # 1D array of x positions
y = Ex.getAxis("y")

# Visualize the first frame
Ex.plot(timestep=0)

# Slide through frames interactively
Ex.slide()

# Apply a slice
Ex_slice = S.Field(0, "Ex", subset={"y": [20., 20.]})       # at fixed y=20
Ex_slice.slide()

# Convert to SI on the fly
Ex_si = S.Field(0, "Ex", units=["um", "fs", "GV/m"])
Ex_si.slide()
```

The `subset` argument projects to a lower-dimensional slice — useful for picking a 1D line out of a 2D field.

## 7. `Probe` access

```python
# Inspect probe geometry
print(S.probeInfo())

# Get the first probe — a 1D line along x
Ex_probe = S.Probe(0, "Ex")
Ex_probe.slide()

# Time-space "streak" (x on horizontal, t on vertical)
Ex_probe.streak()

# 2D probe — gives 2D heatmap
Bz = S.Probe(1, "Bz")
Bz.slide()
```

`DiagProbe` outputs are typically cheaper than `DiagFields` (sparse sampling), so streak plots over many timesteps are practical.

## 8. `ParticleBinning` access

```python
# What ParticleBinning blocks exist?
print(S.getDiags())                          # general diagnostic listing

# Get the first ParticleBinning
pb = S.ParticleBinning(0)
data = pb.getData()                           # list of arrays
axes_info = pb.getAxis(0)                    # bin centers of axis 0

# Plot
pb.plot()                                    # 1D → line, 2D → heatmap
pb.slide()                                   # over time
pb.streak()                                  # 1D-time
```

For an energy spectrum:

```python
# axes were [["ekin", 1e-4, 1., 100, "logscale"]]
spec = S.ParticleBinning(0)
spec.set(xlog=True, ylog=True)               # log axes for spectrum
spec.plot()
```

For combining multiple binning diagnostics (e.g., dividing one by another to get a ratio):

```python
# Energy per particle = weight_ekin / weight
ekin = S.ParticleBinning(0)                  # weight_ekin
nparticles = S.ParticleBinning(1)            # weight
ratio = ekin / nparticles                    # happi supports arithmetic between accessors
ratio.plot()
```

## 9. `Screen`, `RadiationSpectrum` access

```python
scr = S.Screen(0)
scr.plot()                                   # screen output is typically 2D (angle, energy)
scr.streak()
```

Same interface as `ParticleBinning`.

## 10. `TrackParticles` access

```python
# Get trajectory data
tr = S.TrackParticles("electron", axes=["x", "y", "px", "py"])
data = tr.getData()                          # dict of arrays per attribute, per timestep

# Plot trajectories
tr.plot()                                    # one line per tracked particle
tr.slide()                                   # advance over time
```

Memory note: `getData()` loads all particle trajectories into RAM. For tracking many particles over many timesteps this can be large. Use `select=` to restrict:

```python
tr = S.TrackParticles("electron", axes=["x", "y"], select={"px": [0.5, 10.]})
```

This filters at load time — only particles whose px is between 0.5 and 10 are loaded.

## 11. `NewParticles` access

```python
# All ionization / pair-creation events
new = S.NewParticles("electron", axes=["x", "y", "px"])
times = new.getTimes()                        # creation time per event
data = new.getData()
```

For ionization-front diagnostics, plot creation time vs position:

```python
import numpy as np
import matplotlib.pyplot as plt
events = new.getData()
xs = events["x"]
ts = new.getTimes()
plt.scatter(xs, ts, s=1)                      # ionization front in (x, t) space
```

## 12. `Performances` access

```python
perf = S.Performances()
perf.plot()                                   # per-patch compute time, particle count, etc.
```

Useful for diagnosing load imbalance. Look for outlier patches with much higher computational cost than median.

## 13. Unit conversion with Pint

When happi is opened with `pint = True` or when units are specified per call, conversion uses Pint:

```python
# Open with pint
S = happi.Open("/path/to/sim", pint=True)

# Specify desired units per accessor
Ex = S.Probe(0, "Ex", units=["um", "fs", "GV/m"])
Ex.slide()

# Densities in cm⁻³
Rho = S.Field(0, "Rho_electron", units=["um", "fs", "cm**-3"])

# Energy
ekin = S.ParticleBinning(0, units=["um", "fs", "keV"])
```

Pint syntax (`"cm**-3"`, `"GV/m"`) is standard. If the units don't make dimensional sense for the quantity, Pint raises.

### Manual conversion

If you don't want Pint, work directly:

```python
omega_r = S.namelist.Main.reference_angular_frequency_SI
c = 2.99792458e8

# Length: c/omega_r in meters
L_r_SI = c / omega_r

# Time: 1/omega_r in seconds
T_r_SI = 1.0 / omega_r

# Field: m_e c omega_r / e in V/m
import scipy.constants as sc
E_r_SI = sc.m_e * sc.c * omega_r / sc.e

# Get raw normalized data
Ex_normalized = S.Field(0, "Ex").getData()
# Convert to SI
Ex_SI = [arr * E_r_SI for arr in Ex_normalized]
```

## 14. Plotting, animating, exporting

```python
diag = S.Field(0, "Rho_electron")

# Static plot of one frame
diag.plot(timestep=500)

# Interactive slider
diag.slide()

# Animated GIF
diag.animate(filename="density.gif", fps=10)

# 2D space-time "streak" (useful for 1D probes)
S.Probe(0, "Ex").streak()

# Export to VTK for ParaView/VisIt
diag.toVTK(folder="vtk_output")

# Modify plot parameters
diag.set(vmin=0, vmax=2.0, cmap="viridis", xlabel="x [μm]")
diag.plot()
```

`.set()` modifies the accessor's plot config in place; subsequent calls use the new settings.

## 15. Multi-simulation comparison

```python
# Open three simulations together
S_list = happi.Open(["sim_a", "sim_b", "sim_c"])

# Compare scalar evolution
for s in S_list:
    s.Scalar("Ukin").plot(figure=1)

# Or with custom labels
import matplotlib.pyplot as plt
for s, label in zip(S_list, ["a","b","c"]):
    times = s.Scalar("Ukin").getTimes()
    data = s.Scalar("Ukin").getData()
    plt.plot(times, data, label=label)
plt.legend()
plt.show()
```

For parameter scans, a common pattern is iterating over a directory of runs:

```python
import glob
runs = sorted(glob.glob("/scratch/scan_*/"))
for run in runs:
    S = happi.Open(run)
    # ... analyze ...
```

## 16. Common mistakes

### Mistake 1: opening a directory that's not a Smilei simulation

```python
S = happi.Open("/wrong/path")    # WARNING: no valid Smilei simulation
```
happi issues a warning and returns an invalid object. Method calls then fail confusingly. Verify the directory contains `smilei.py` before opening.

### Mistake 2: requesting units when `reference_angular_frequency_SI` wasn't set

```python
# Namelist had no reference_angular_frequency_SI
S = happi.Open("path")
S.Field(0, "Ex", units=["um","fs","V/m"])    # raises
```
Either re-run with `reference_angular_frequency_SI` set, or pass it to `happi.Open()`:

```python
S = happi.Open("path", reference_angular_frequency_SI = 2*3.14159*3e8/800e-9)
```

### Mistake 3: trying to plot a diagnostic that wasn't enabled

```python
S.ParticleBinning(0)    # but the namelist had no DiagParticleBinning blocks
```
Returns `None` or an empty object. Use `S.getDiags()` to see what's available.

### Mistake 4: forgetting `pint=True` and then using Pint syntax

```python
S = happi.Open("path")
S.Field(0, "Ex", units=["um", "fs", "GV/m"])    # works without pint, simple unit conversion
```

Most simple conversions work without `pint=True`. Pint is needed for complex compound units (e.g., `"V*m/A"`) or non-standard expressions.

### Mistake 5: `getData()` returns nested arrays

```python
data = S.Field(0, "Ex").getData()
print(data[0])    # this is a 2D NumPy array
print(data.shape) # AttributeError — data is a Python list, not array
```
`getData()` returns a Python list of arrays, one per output timestep. To stack them:

```python
import numpy as np
full = np.array(S.Field(0, "Ex").getData())    # shape: (n_times, nx, ny)
```

### Mistake 6: too many `.plot()` calls without `plt.show()`

In interactive Python without `%matplotlib`, each `.plot()` creates a new figure but doesn't display. Add `import matplotlib.pyplot as plt; plt.show()` at the end of your script.

### Mistake 7: opening simulations in a loop without cleanup

```python
results = []
for run in many_runs:
    S = happi.Open(run)
    results.append(S.Scalar("Ukin").getData())
```
For ~100 runs this is fine. For 1000+ runs, accumulated open HDF5 file handles can exhaust the OS limit. Close explicitly or extract data and discard the `S` object.

### Mistake 8: assuming HDF5 layout is stable across versions

`getData()` and friends are the supported interface. Reading the HDF5 files directly with h5py works but the layout has changed between Smilei versions. Use happi for portability.

### Mistake 9: passing wrong type to `subset`

```python
Ex = S.Field(0, "Ex", subset={"y": 20.})        # WRONG — wants a list/range, not a scalar
Ex = S.Field(0, "Ex", subset={"y": [20.,20.]})  # CORRECT — narrow slice at y=20
```

### Mistake 10: confusing `every` (namelist) and getTimesteps (post-processing)

`getTimesteps()` returns the timesteps at which output was actually written (i.e., multiples of the namelist's `every`). It is not the simulation's total timestep count.
