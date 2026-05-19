# Diagnostics

This reference covers every diagnostic block in Smilei: `DiagScalar`, `DiagFields`, `DiagProbe`, `DiagParticleBinning`, `DiagTrackParticles`, `DiagScreen`, `DiagRadiationSpectrum`, `DiagNewParticles`, `DiagPerformances`. Each section gives syntax, common patterns, output structure, and file-size considerations. For accessing the output via Python, see `happi-and-postprocessing.md`.

## Table of contents

1. Diagnostic philosophy and output structure
2. `DiagScalar` — global integrated quantities
3. `DiagFields` — grid-based field snapshots
4. `DiagProbe` — sampled lines/surfaces in the box
5. `DiagParticleBinning` — distribution-function projections
6. `DiagTrackParticles` — selected-particle trajectories
7. `DiagScreen` — flux through a surface
8. `DiagRadiationSpectrum` — photon emission spectra
9. `DiagNewParticles` — ionization / pair creation events
10. `DiagPerformances` — per-patch timing
11. Filter functions (the `filter` argument)
12. Output cadence and file size
13. Common mistakes

---

## 1. Diagnostic philosophy and output structure

Smilei writes all diagnostics to HDF5 files in the simulation output directory. Each diagnostic block produces one or more files (`Scalars.txt`, `Fields0.h5`, `Probes0.h5`, `ParticleBinning0.h5`, etc.). The filename suffix is the index — multiple `DiagFields` blocks become `Fields0.h5`, `Fields1.h5`, etc. in the order they appear in the namelist.

Diagnostic data is always written in normalized (Smilei) units. Conversion to SI happens at post-processing time via `happi` (see `happi-and-postprocessing.md`).

All diagnostics share two timing-control patterns:

- `every = N` — output cadence in timesteps. `every = 100` means every 100 timesteps. Setting `every = 0` disables the diagnostic.
- `time_range` or `timesteps` — start/stop or explicit list of times for output.

Output cadence is the single biggest control on disk usage. A 3D simulation with full-grid `DiagFields` at `every = 10` can easily produce terabytes.

## 2. `DiagScalar` — global integrated quantities

```python
DiagScalar(
    every = 50,                                 # cadence
    precision = 10,                             # decimal digits of output
    vars = ["Utot", "Ukin", "Uelm", "Uexp"],    # restrict to subset (optional)
)
```

Writes a flat text file (`Scalars.txt`) with one row per output time. Each column is a globally-integrated quantity. Cheap to compute, cheap to write, always enable. The first diagnostic you should add to any namelist.

Available scalars (subset):

| Name | Meaning |
|---|---|
| `Utot` | Total energy (kinetic + EM + bound + radiated) |
| `Ukin` | Total kinetic energy across all species |
| `Uelm` | Total EM-field energy |
| `Uexp` | Energy expected from absorbed laser (sanity check) |
| `Ubal` | Energy balance `Utot - Uexp` (should be ~0) |
| `Ukin_<species>` | Per-species kinetic energy |
| `Uelm_Ex`, `Uelm_Ey`, etc. | Per-component EM energy |
| `Ntot_<species>` | Total macro-particle count |
| `Zavg_<species>` | Mean charge (useful for tracking ionization progress) |

A common quick check after a run: load `Scalars.txt` and verify `Utot` is conserved to within 1% (or whatever your run can tolerate). If it isn't, the simulation is under-resolved or there's a numerical instability.

## 3. `DiagFields` — grid-based field snapshots

```python
DiagFields(
    every = 200,                                # cadence
    fields = ["Ex", "Ey", "Ez", "Bz", "Rho_electron", "Rho_ion"],
    time_average = 1,                           # time-average over N steps centered on output time
    subgrid = [slice(None), slice(0, 100)],     # output only a slice
)
```

Writes full-grid HDF5 snapshots of the named fields. The most general (and most expensive) diagnostic.

Available fields:

| Field | Description |
|---|---|
| `Ex`, `Ey`, `Ez` | Electric field components |
| `Bx`, `By`, `Bz` | Magnetic field components |
| `Jx`, `Jy`, `Jz` | Total current density |
| `Jx_<species>`, ... | Per-species current density |
| `Rho` | Total charge density |
| `Rho_<species>` | Per-species charge density |
| `Env_A`, `Env_E` | Envelope laser amplitude (envelope mode only) |

### `time_average`

If `time_average > 1`, the output is averaged over `time_average` timesteps centered on the output time. Use for slowly-varying fields (envelope-like density profiles, average current) or to smooth out laser-cycle oscillation. **Do not use `time_average > 1` for instantaneous quantities** (laser snapshots, instantaneous polarization) — it smears physics.

### `subgrid`

Save only a sub-region using Python slicing syntax:

```python
subgrid = [slice(0, 500), slice(None)]          # x: cells 0-500, y: full range
subgrid = [slice(None, None, 4), slice(None)]   # x: every 4th cell (downsample)
```

The slice format is `[slice for axis 0, slice for axis 1, ...]`. Downsampling with `slice(None, None, 4)` is a common space-saver when you want spatial maps but don't need full resolution.

### File size estimate

For a 2D simulation with grid `1000 × 500`, four fields, `every = 200`, simulation_time covering 1000 output frames: `1000 × 500 × 4 × 8 bytes × 1000 frames = 16 GB`. Multiply by `time_average` overhead. Plan ahead.

## 4. `DiagProbe` — sampled lines/surfaces in the box

```python
DiagProbe(
    every = 100,
    origin = [0., 20.],                         # starting point in normalized units
    corners = [[60., 20.]],                     # 1D line: one corner; 2D plane: two corners
    number = [300],                             # number of probe points along each dimension
    fields = ["Ex", "Ey", "Bz", "Rho"],
    time_average = 1,
)
```

`DiagProbe` samples fields at a regularly-spaced set of points within the box. The `corners` list defines the geometry:

- 0 corners → 0D probe (single point at `origin`)
- 1 corner → 1D probe (line from `origin` to `corners[0]`)
- 2 corners → 2D probe (parallelogram with `origin` + the two corner vectors)
- 3 corners → 3D probe (box)

`number` is the number of sample points along each probe axis.

A 1D probe along the laser axis is the cheapest way to see laser propagation through plasma. Output is small (1D × N_times), so cadence can be aggressive (`every = 10` is fine).

```python
# 1D probe along x at y = 20 (the box midline)
DiagProbe(
    every = 50,
    origin = [0., 20.],
    corners = [[Lx, 20.]],
    number = [400],
    fields = ["Ex", "Ey", "Bz"],
)

# 2D probe in the xy plane at z = 0 (3D simulation)
DiagProbe(
    every = 200,
    origin = [0., 0., 0.],
    corners = [[Lx, 0., 0.], [0., Ly, 0.]],
    number = [200, 100],
    fields = ["Rho", "Bz"],
)
```

## 5. `DiagParticleBinning` — distribution-function projections

```python
DiagParticleBinning(
    deposited_quantity = "weight",              # what to deposit
    every = 200,                                # cadence
    species = ["electron"],                     # one or more species
    axes = [
        ["x", 0., 60., 200],                    # axis 0: x from 0 to 60, 200 bins
        ["px", -0.5, 0.5, 100],                 # axis 1: px from -0.5 to 0.5, 100 bins
    ],
    time_average = 1,
)
```

Projects the macro-particle distribution onto a chosen subspace. The most flexible diagnostic — handles density maps, phase-space slices, energy spectra, charge-state distributions, etc.

### `deposited_quantity`

Controls what is summed into each bin:

| Value | What's deposited |
|---|---|
| `"weight"` | Macro-particle weight (number) |
| `"weight_charge"` | Weight × charge (gives charge density) |
| `"weight_charge_vx"` | Weight × charge × vx (current density component) |
| `"weight_p"` | Weight × |p| (momentum magnitude sum) |
| `"weight_ekin"` | Weight × kinetic energy (energy spectrum) |
| `"weight_vx_px"` etc. | Various products useful for fluxes |
| `"weight_chi"` | Quantum parameter χ (relativistic) |

A user-defined Python callable can also be deposited, e.g., for custom weighting.

### `axes`

List of `[name, min, max, nbins]` for each binning dimension. Up to 3 axes typical. Each entry can also append `"logscale"` for a log-spaced axis or `"edge_inclusive"` for boundary handling:

```python
axes = [
    ["ekin", 1e-4, 1e2, 100, "logscale"],       # log-spaced energy bins
    ["x", 0., 60., 200],
]
```

### Common patterns

**Density map (2D)**:
```python
DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = 200, species = ["electron"],
    axes = [["x", 0., Lx, 400], ["y", 0., Ly, 200]],
)
```

**Energy spectrum**:
```python
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 100, species = ["electron"],
    axes = [["ekin", 1e-4, 10., 200, "logscale"]],
)
```

**Phase-space slice (x, px)**:
```python
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 200, species = ["electron"],
    axes = [["x", 0., Lx, 300], ["px", -1., 1., 200]],
)
```

**Charge-state distribution**:
```python
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 200, species = ["Ar0"],
    axes = [["charge", -0.5, 18.5, 19]],
)
```

**Angular distribution (with `chi` and direction)**:
```python
import math
DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 500, species = ["electron"],
    axes = [["x", 0., Lx, 200], ["theta_xy", -math.pi, math.pi, 90]],
)
```

## 6. `DiagTrackParticles` — selected-particle trajectories

```python
DiagTrackParticles(
    species = "electron",
    every = 100,
    flush_every = 1000,
    filter = lambda particles: particles.px > 1.0,    # select particles to track
    attributes = ["x", "y", "px", "py", "pz", "w"],
)
```

Tracks individual macro-particles' state over time. Output is per-particle time series.

The `filter` argument selects which particles to track. It's a callable that receives a struct-of-arrays of all particle attributes and returns a boolean mask. The example above tracks any electron with `px > 1` (relativistic forward-moving).

`every` controls the cadence of state recording. `flush_every` controls how often the buffer is written to disk. For long simulations, `flush_every` should be much larger than `every` to reduce I/O overhead.

This diagnostic is expensive in memory (one record per tracked particle per output) and disk. For tracking 10,000 particles over 1000 output frames with 7 attributes, that's `10⁴ × 10³ × 7 × 8 bytes ≈ 0.6 GB`. Manageable. Tracking 10⁷ particles is not.

## 7. `DiagScreen` — flux through a surface

```python
DiagScreen(
    shape = "plane",                            # or "sphere"
    point = [Lx, 20.],                          # point on the screen
    vector = [1., 0.],                          # normal to the plane
    direction = "forward",                      # or "backward", "both", "canceling"
    every = 100,
    species = ["electron"],
    deposited_quantity = "weight_ekin",
    axes = [
        ["a", -10., 10., 100],                  # transverse coordinate
        ["ekin", 1e-4, 10., 100, "logscale"],
    ],
)
```

Records the flux of particles through a virtual surface. The "plane" shape is the most common — at a fixed location in the box, what's the time-integrated flux of particles crossing it?

Useful for:
- Beam spectra at a sampling point (e.g., what comes out of the back of the target)
- Forward-emission angular distributions
- Charged-particle flux through an aperture

`direction = "canceling"` records both directions but subtracts back-crossing particles — useful when you want net forward flux.

The output `axes` typically include a transverse coordinate (`"a"`, the transverse position on the screen) and an energy or momentum axis.

## 8. `DiagRadiationSpectrum` — photon emission spectra

```python
DiagRadiationSpectrum(
    every = 500,
    species = ["electron"],
    axes = [["x", 0., Lx, 100], ["gamma", 1., 1000., 100, "logscale"]],
    photon_energy_axis = [1., 1e6, 500, "logscale"],    # in keV
)
```

For high-energy electrons (γ > 1), records the synchrotron-like radiation spectrum each electron would emit, binned over the chosen axes. The output is a spectrum, not actual emitted photons (for actual photons see Monte-Carlo radiation reaction).

Requires `reference_angular_frequency_SI` to express photon energies in keV.

For excimer regimes this diagnostic is rarely active — electrons are sub-relativistic and synchrotron emission is negligible. Used in UHI laser-plasma where bremsstrahlung-style emission is detectable.

## 9. `DiagNewParticles` — ionization / pair creation events

```python
DiagNewParticles(
    species = "electron",
    every = 100,
    attributes = ["x", "y", "px", "py", "pz", "w"],
)
```

Records every event where a new particle is created in the named species — useful for tracking ionization in time and space, or photon-driven pair creation in QED-mode runs.

For excimer collisional runs with ionization on, this writes one record per ionization event. With many events per timestep, file size grows quickly; use `every` to subsample, or apply a `filter` to record only events in a region of interest.

## 10. `DiagPerformances` — per-patch timing

```python
DiagPerformances(
    every = 200,
    flush_every = 1000,
    patch_information = True,
)
```

Records per-patch compute time, particle count, and memory. The diagnostic for load balancing — if `LoadBalancing` is enabled but performance is still poor, this is how to find which patches are bottlenecks. Useful when running on large clusters.

## 11. Filter functions (the `filter` argument)

Several diagnostics accept a `filter` callable that selects which particles to include:

```python
# DiagParticleBinning — bin only electrons with px > 0
DiagParticleBinning(
    # ...
    species = ["electron"],
    axes = [["x", 0., Lx, 200]],
    filter = lambda particles: particles.px > 0,
)

# DiagTrackParticles — track only particles in a region
def select(p):
    return (p.x > 20.) & (p.x < 30.) & (p.y > 15.) & (p.y < 25.)
DiagTrackParticles(species="electron", every=100, filter=select, ...)

# DiagScreen — only forward-moving particles
DiagScreen(..., filter = lambda p: p.px > 0)
```

The callable receives a Smilei `Particles` proxy that exposes all particle attributes as NumPy-like arrays (vectorized comparison expected). The returned boolean array selects which particles to include.

Common filter patterns:
- `p.px > p_threshold` — forward-moving / energetic
- `(p.x > xmin) & (p.x < xmax)` — within a region
- `p.charge > 0` — only positive charge states (in ladders)
- `p.weight > 1e-10` — filter out zero-weight macro-particles

## 12. Output cadence and file size

Some rules of thumb for keeping output manageable:

| Diagnostic | Sane cadence |
|---|---|
| `DiagScalar` | `every = 10` to `100` — always on, cheap |
| `DiagFields` (full grid) | `every = 200` to `500` — expensive; use subgrid or sparse cadence |
| `DiagFields` (subset of fields) | `every = 100` |
| `DiagProbe` (1D line) | `every = 10` to `50` — cheap |
| `DiagProbe` (2D plane) | `every = 100` to `500` |
| `DiagParticleBinning` | `every = 100` to `500` |
| `DiagTrackParticles` | `every = 50` to `200`, narrow filter |
| `DiagScreen` | `every = 50` to `200` |

For a typical excimer 2D simulation (~10,000 timesteps, ~few-GB output budget), a reasonable starting configuration:

```python
DiagScalar(every = 50)
DiagFields(every = 500, fields = ["Ex","Ey","Bz","Rho_electron","Rho_ion"])
DiagProbe(every = 50, origin=[0., Ly/2], corners=[[Lx, Ly/2]], number=[400], fields=["Ex","Ey","Bz"])
DiagParticleBinning(deposited_quantity="weight_ekin", every=200, species=["electron"], axes=[["ekin", 1e-4, 1., 100, "logscale"]])
```

Then add more diagnostics as needed — but only after a test run shows the basics are working.

## 13. Common mistakes

### Mistake 1: `every = 1` on `DiagFields` for a full grid

This is the fastest way to fill a disk. A 2D simulation at full resolution writing every timestep can produce 100 GB in minutes. Use `every >> 1` and rely on `DiagProbe` for high-cadence time series.

### Mistake 2: forgetting `subgrid` for downsampling

If you want spatial maps but don't need full resolution, `subgrid = [slice(None, None, 4), slice(None, None, 4)]` cuts file size by 16× with minimal physics loss for slowly-varying fields.

### Mistake 3: `time_average` on an instantaneous quantity

Laser snapshots, polarization, anything oscillating at the laser frequency — these need `time_average = 1`. Setting `time_average = 5` smears the laser. Use `time_average > 1` only for envelope quantities, averaged densities, or thermal noise reduction.

### Mistake 4: phase-space axis range too narrow

```python
axes = [["px", -1., 1., 100]]            # but max px = 5 in your problem
```
Particles outside the range are dropped from the bin total. The diagnostic doesn't error — you just see a clipped distribution. Check `Scalar.txt` for max momenta and set bin ranges with margin.

### Mistake 5: filter callable that's not vectorized

```python
def filter(p):
    if p.x > 20:    # treats p.x as scalar; fails on arrays
        return True
    return False
```
Smilei passes arrays; filter must use array operations. Use `&`, `|`, comparison operators directly, not `if` statements.

### Mistake 6: tracking too many particles

```python
DiagTrackParticles(species="electron", filter = lambda p: True, ...)
```
Tracks every electron. For a typical run that's `~10⁹` particle states recorded — terabytes of data. Always use a narrow filter that selects a small fraction.

### Mistake 7: defining many `DiagFields` with overlapping field lists

```python
DiagFields(every = 100, fields = ["Ex", "Ey"])
DiagFields(every = 100, fields = ["Bz"])
DiagFields(every = 100, fields = ["Rho_electron"])
```
Three separate files, three I/O passes. More efficient as one block:

```python
DiagFields(every = 100, fields = ["Ex","Ey","Bz","Rho_electron"])
```

Unless you want different cadences for different fields, consolidate.

### Mistake 8: forgetting to enable `DiagPerformances` when load-balancing

If a run is slow and you've enabled `LoadBalancing`, add `DiagPerformances` and inspect — without it, you can't diagnose whether load balancing is working.
