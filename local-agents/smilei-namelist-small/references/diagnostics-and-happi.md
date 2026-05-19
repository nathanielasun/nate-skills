# Diagnostics and happi

Every diagnostic block in Smilei and how to read the output via `happi`. Merged from frontier's two output references.

## Table of contents

1. Diagnostic philosophy and output structure
2. `DiagScalar` — global integrated quantities
3. `DiagFields` — grid-based snapshots
4. `DiagProbe` — sampled lines/surfaces
5. `DiagParticleBinning` — distribution-function projections
6. `DiagTrackParticles` — selected-particle trajectories
7. `DiagScreen`, `DiagRadiationSpectrum`, `DiagNewParticles`, `DiagPerformances`
8. Filter functions
9. Cadence and file size
10. happi — opening a simulation
11. Reading diagnostics with happi
12. Unit conversion
13. Plotting, animating, exporting
14. Multi-simulation comparison
15. Common mistakes

---

## 1. Diagnostic philosophy and output structure

Smilei writes all diagnostics to HDF5 files in the simulation output directory. Each diagnostic block produces files indexed by appearance order: `Fields0.h5`, `Fields1.h5`, `Probes0.h5`, `ParticleBinning0.h5`, etc. Data is always in normalized units; conversion happens in happi.

Common timing controls:
- `every = N` — output cadence in timesteps. `every = 0` disables.
- `time_average = N` — average over N timesteps centered on output time (default 1 = none).

## 2. `DiagScalar` — global integrated quantities

```python
DiagScalar(
    every = 50,
    precision = 10,
    vars = ["Utot", "Ukin", "Uelm"],          # optional restriction
)
```

Writes flat text (`Scalars.txt`). Cheap. **Always enable** — first diagnostic in any namelist.

| Name | Meaning |
|---|---|
| `Utot` | Total energy (kinetic + EM + bound + radiated) |
| `Ukin` | Total kinetic across all species |
| `Uelm` | EM-field energy |
| `Uexp` | Energy expected from absorbed laser (sanity check) |
| `Ubal` | `Utot - Uexp` (should be ~0) |
| `Ukin_<species>` | Per-species kinetic |
| `Ntot_<species>` | Total macro-particle count |
| `Zavg_<species>` | Mean charge (tracks ionization progress) |

Quick health check after a run: `Utot` should be conserved to within 1%.

## 3. `DiagFields` — grid-based snapshots

```python
DiagFields(
    every = 200,
    fields = ["Ex", "Ey", "Ez", "Bz", "Rho_electron", "Rho_ion"],
    time_average = 1,
    subgrid = [slice(None), slice(0, 100)],   # output only a slice
)
```

Most expensive diagnostic. Use `every >> 1` and `subgrid` to control file size.

Available fields: `Ex/Ey/Ez`, `Bx/By/Bz`, `Jx/Jy/Jz`, `Jx_<species>`, `Rho`, `Rho_<species>`, `Env_A`, `Env_E` (envelope mode).

**`time_average > 1`** smooths physics. Use only for envelope-like quantities, never for instantaneous fields (laser snapshots, polarization).

**`subgrid`** with `slice(None, None, 4)` downsamples by 4× per axis — 16× file-size reduction in 2D for slowly-varying fields.

**File size estimate**: 2D grid `1000 × 500`, four fields, `every = 200`, 1000 output frames → 16 GB. Plan ahead.

## 4. `DiagProbe` — sampled lines/surfaces

```python
DiagProbe(
    every = 100,
    origin = [0., 20.],
    corners = [[60., 20.]],                   # 1D line: one corner; 2D plane: two
    number = [300],
    fields = ["Ex", "Ey", "Bz", "Rho"],
)
```

Geometry:
- 0 corners → 0D point at `origin`
- 1 corner → 1D line from `origin` to `corners[0]`
- 2 corners → 2D parallelogram
- 3 corners → 3D box

A 1D probe along the laser axis is the cheapest way to see laser propagation — small output, can use `every = 10` aggressively.

## 5. `DiagParticleBinning` — distribution-function projections

```python
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 200,
    species = ["electron"],
    axes = [
        ["x", 0., 60., 200],
        ["px", -0.5, 0.5, 100],
    ],
)
```

**`deposited_quantity`** options:

| Value | Deposits |
|---|---|
| `"weight"` | macro-particle weight (number density) |
| `"weight_charge"` | weight × charge (charge density) |
| `"weight_charge_vx"` | weight × charge × vx (current component) |
| `"weight_p"` | weight × |p| |
| `"weight_ekin"` | weight × kinetic energy (energy spectrum) |
| `"weight_chi"` | quantum parameter χ |

**Axes**: `[name, min, max, nbins]`; append `"logscale"` for log-spaced bins:

```python
axes = [["ekin", 1e-4, 1e2, 100, "logscale"]]
```

**Common patterns:**

```python
# Density map
DiagParticleBinning(
    deposited_quantity = "weight_charge", every = 200, species = ["electron"],
    axes = [["x", 0., Lx, 400], ["y", 0., Ly, 200]],
)

# Energy spectrum
DiagParticleBinning(
    deposited_quantity = "weight", every = 100, species = ["electron"],
    axes = [["ekin", 1e-4, 10., 200, "logscale"]],
)

# Phase-space (x, px)
DiagParticleBinning(
    deposited_quantity = "weight", every = 200, species = ["electron"],
    axes = [["x", 0., Lx, 300], ["px", -1., 1., 200]],
)

# Charge-state distribution
DiagParticleBinning(
    deposited_quantity = "weight", every = 200, species = ["Ar0"],
    axes = [["charge", -0.5, 18.5, 19]],
)
```

## 6. `DiagTrackParticles` — selected-particle trajectories

```python
DiagTrackParticles(
    species = "electron",
    every = 100,
    flush_every = 1000,
    filter = lambda p: p.px > 1.0,
    attributes = ["x", "y", "px", "py", "pz", "w"],
)
```

**`filter`** selects which particles to track — array-vectorized callable returning a boolean mask. Always use a narrow filter; tracking all particles can produce terabytes.

`flush_every >> every` reduces I/O frequency.

## 7. `DiagScreen`, `DiagRadiationSpectrum`, `DiagNewParticles`, `DiagPerformances`

```python
# Flux through a surface
DiagScreen(
    shape = "plane", point = [Lx, 20.], vector = [1., 0.],
    direction = "forward",                    # or "backward", "both", "canceling"
    every = 100, species = ["electron"],
    deposited_quantity = "weight_ekin",
    axes = [["a", -10., 10., 100], ["ekin", 1e-4, 10., 100, "logscale"]],
)

# Synchrotron-like emission spectrum (requires reference_angular_frequency_SI)
DiagRadiationSpectrum(
    every = 500, species = ["electron"],
    axes = [["x", 0., Lx, 100], ["gamma", 1., 1000., 100, "logscale"]],
    photon_energy_axis = [1., 1e6, 500, "logscale"],
)

# Ionization / pair creation events
DiagNewParticles(
    species = "electron",
    every = 100,
    attributes = ["x", "y", "px", "py", "pz", "w"],
)

# Per-patch timing (essential when load balancing)
DiagPerformances(every = 200, patch_information = True)
```

## 8. Filter functions

Several diagnostics accept `filter` — a vectorized boolean callable:

```python
# Forward-moving only
filter = lambda p: p.px > 0

# In a region
filter = lambda p: (p.x > 20.) & (p.x < 30.) & (p.y > 15.) & (p.y < 25.)

# High-energy
filter = lambda p: p.px**2 + p.py**2 > 1.0

# Specific charge state
filter = lambda p: p.charge > 5
```

Must use array operators (`&`, `|`, comparisons) — never `if` (treats arrays as scalars).

## 9. Cadence and file size

| Diagnostic | Sane cadence |
|---|---|
| `DiagScalar` | `every = 10–100` (cheap, always on) |
| `DiagFields` full grid | `every = 200–500` |
| `DiagFields` with `subgrid` | `every = 100` |
| `DiagProbe` 1D | `every = 10–50` |
| `DiagProbe` 2D | `every = 100–500` |
| `DiagParticleBinning` | `every = 100–500` |
| `DiagTrackParticles` | `every = 50–200`, narrow filter |
| `DiagScreen` | `every = 50–200` |

Reasonable starting config for a 2D excimer run:

```python
DiagScalar(every = 50)
DiagFields(every = 500, fields = ["Ex","Ey","Bz","Rho_electron","Rho_ion"])
DiagProbe(every = 50, origin=[0., Ly/2], corners=[[Lx, Ly/2]],
          number=[400], fields=["Ex","Ey","Bz"])
DiagParticleBinning(deposited_quantity="weight_ekin", every=200,
                    species=["electron"],
                    axes=[["ekin", 1e-4, 1., 100, "logscale"]])
```

## 10. happi — opening a simulation

```python
import happi

# Single simulation
S = happi.Open("/path/to/simulation_directory")

# Override reference frequency if namelist didn't set it
S = happi.Open("/path", reference_angular_frequency_SI = 2*3.14159*3e8/248e-9)

# Multiple simulations
S = happi.Open(["sim_a", "sim_b", "sim_c"])

# Access namelist
print(S.namelist.Main.timestep)
print(S.namelist.Main.reference_angular_frequency_SI)
print([s.name for s in S.namelist.Species])
```

Install: `cd <smilei-dir> && make happi`. If import fails: `export PYTHONPATH=/path/to/smilei/happi:$PYTHONPATH`.

## 11. Reading diagnostics with happi

| Diagnostic block | Accessor |
|---|---|
| `DiagScalar` | `S.Scalar(name)` |
| `DiagFields` | `S.Field(index, field_name)` |
| `DiagProbe` | `S.Probe(index, field_name)` |
| `DiagParticleBinning` | `S.ParticleBinning(index)` |
| `DiagScreen` | `S.Screen(index)` |
| `DiagRadiationSpectrum` | `S.RadiationSpectrum(index)` |
| `DiagTrackParticles` | `S.TrackParticles(species, attributes)` |
| `DiagNewParticles` | `S.NewParticles(species)` |
| `DiagPerformances` | `S.Performances()` |

All accessors return objects with common methods: `.getData()`, `.getTimes()`, `.getAxis(name)`, `.plot()`, `.slide()`, `.animate()`, `.streak()`, `.toVTK()`, `.set(...)`.

```python
# Scalar time series
Ukin = S.Scalar("Ukin")
Ukin.plot()

# Field snapshots
Ex = S.Field(0, "Ex")
Ex.slide()                                    # interactive slider
Ex.streak()                                   # 2D time-space

# Field with slice projection
Ex_at_y20 = S.Field(0, "Ex", subset={"y": [20., 20.]})

# Probe streak (1D line over time)
S.Probe(0, "Ex").streak()

# ParticleBinning
pb = S.ParticleBinning(0)
pb.set(xlog=True, ylog=True)                  # log axes
pb.plot()

# TrackParticles with selection
tr = S.TrackParticles("electron",
                      axes=["x","y","px"],
                      select={"px": [0.5, 10.]})
```

`getData()` returns a Python list of arrays (one per output timestep). To stack: `np.array(diag.getData())`.

## 12. Unit conversion

The `units` argument converts on the fly via Pint-style syntax:

```python
S.Field(0, "Ex", units=["um","fs","GV/m"]).slide()
S.Field(0, "Rho_electron", units=["um","fs","cm**-3"]).slide()
S.ParticleBinning(0, units=["um","fs","keV"]).plot()
```

Requires `reference_angular_frequency_SI` in namelist (or passed to `happi.Open()`). Without it, conversion raises.

**Manual conversion factors** if you want raw numbers:

```python
omega_r = S.namelist.Main.reference_angular_frequency_SI
c = 2.99792458e8
L_r = c / omega_r                             # meters
T_r = 1.0 / omega_r                           # seconds
import scipy.constants as sc
E_r = sc.m_e * sc.c * omega_r / sc.e          # V/m

Ex_normalized = S.Field(0, "Ex").getData()
Ex_SI = [arr * E_r for arr in Ex_normalized]
```

## 13. Plotting, animating, exporting

```python
diag = S.Field(0, "Rho_electron")

diag.plot(timestep=500)                       # static
diag.slide()                                  # interactive slider
diag.animate(filename="density.gif", fps=10)  # animated GIF/MP4
diag.streak()                                 # 2D time-space (for 1D probes)
diag.toVTK(folder="vtk_output")               # ParaView/VisIt export

diag.set(vmin=0, vmax=2.0, cmap="viridis", xlabel="x [μm]")
diag.plot()
```

`.set()` modifies plot config in place.

## 14. Multi-simulation comparison

```python
# Open three simulations
S_list = happi.Open(["sim_a", "sim_b", "sim_c"])

# Compare scalars
import matplotlib.pyplot as plt
for s, label in zip(S_list, ["a","b","c"]):
    times = s.Scalar("Ukin").getTimes()
    data = s.Scalar("Ukin").getData()
    plt.plot(times, data, label=label)
plt.legend(); plt.show()

# Parameter scan over a directory tree
import glob
for run in sorted(glob.glob("/scratch/scan_*/")):
    S = happi.Open(run)
    # ... analyze ...
```

## 15. Common mistakes

**`every = 1` on full `DiagFields`**: fills the disk in minutes. Use `every >> 1` or `DiagProbe` for high-cadence time series.

**`time_average > 1` on instantaneous quantity**: smears laser snapshots. Use only for envelope or averaged quantities.

**Phase-space axis range too narrow**: bins clip silently; check `Scalar.txt` max momenta and set ranges with margin.

**Filter using `if`**: must be vectorized (`&`, `|`, comparisons). `if p.x > 20: return True` fails on arrays.

**Tracking too many particles**: tracking all of an electron species can produce terabytes. Always narrow with `filter`.

**Many `DiagFields` with overlapping fields**: separate blocks → separate I/O passes. Consolidate.

**Forgetting `DiagPerformances` with `LoadBalancing`**: can't diagnose imbalance without it.

**Requesting units without `reference_angular_frequency_SI`**: raises. Either re-run with rule 1 satisfied, or pass to `happi.Open(..., reference_angular_frequency_SI=...)`.

**`subset` as scalar**: `subset={"y": 20.}` fails. Use `subset={"y": [20., 20.]}` for narrow slice.

**Reading HDF5 directly with h5py**: works but layout has changed across Smilei versions. Use happi for portability.
