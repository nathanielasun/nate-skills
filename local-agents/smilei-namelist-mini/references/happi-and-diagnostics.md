# Diagnostics and happi

Every diagnostic block in Smilei and how to read output via `happi`. Merged for the mini tier.

## Diagnostic blocks — at a glance

All diagnostics share `every = N` (cadence in timesteps; 0 disables). Output is in normalized units; conversion happens in happi.

### `DiagScalar` — global quantities (always enable)

```python
DiagScalar(every = 50, vars = ["Utot", "Ukin", "Uelm"])    # vars optional
```

Writes flat text (`Scalars.txt`). Cheap. Available: `Utot`, `Ukin`, `Uelm`, `Uexp` (expected absorbed laser energy), `Ubal` (`Utot - Uexp`), `Ukin_<species>`, `Ntot_<species>`, `Zavg_<species>`.

**Health check**: `Utot` should be conserved to within ~1%.

### `DiagFields` — grid snapshots (expensive)

```python
DiagFields(
    every = 200,
    fields = ["Ex", "Ey", "Ez", "Bz", "Rho_electron", "Rho_ion"],
    time_average = 1,
    subgrid = [slice(None), slice(0, 100)],     # optional sub-region
)
```

Available: `Ex/Ey/Ez`, `Bx/By/Bz`, `Jx/Jy/Jz`, `Jx_<species>`, `Rho`, `Rho_<species>`, `Env_A`, `Env_E` (envelope mode).

**Don't use `time_average > 1` for instantaneous quantities** (laser, polarization). Smears physics.

**File size**: 2D grid 1000×500, 4 fields, `every=200`, 1000 frames → 16 GB. Plan ahead.

### `DiagProbe` — sampled lines/surfaces (cheap)

```python
DiagProbe(
    every = 100,
    origin = [0., 20.],
    corners = [[60., 20.]],                     # 1 corner → 1D line; 2 → 2D plane; 3 → 3D box
    number = [300],
    fields = ["Ex", "Ey", "Bz", "Rho"],
)
```

Best diagnostic for high-cadence time series. A 1D probe along the laser axis is the cheapest way to see laser propagation.

### `DiagParticleBinning` — distribution projections

```python
DiagParticleBinning(
    deposited_quantity = "weight",              # or weight_charge, weight_ekin, weight_p, etc.
    every = 200,
    species = ["electron"],
    axes = [
        ["x", 0., 60., 200],
        ["px", -0.5, 0.5, 100],
    ],
)
```

Common `deposited_quantity` values:

| Value | Deposits |
|---|---|
| `"weight"` | number density |
| `"weight_charge"` | charge density |
| `"weight_charge_vx"` | current density component |
| `"weight_ekin"` | energy density (spectrum) |

Append `"logscale"` to an axis for log-spaced bins:

```python
axes = [["ekin", 1e-4, 1e2, 100, "logscale"]]
```

### `DiagTrackParticles` — particle trajectories

```python
DiagTrackParticles(
    species = "electron",
    every = 100,
    flush_every = 1000,
    filter = lambda p: p.px > 1.0,              # narrow filter REQUIRED
    attributes = ["x", "y", "px", "py", "pz", "w"],
)
```

Tracking all particles can produce terabytes. **Always use a narrow filter.** Filter must be vectorized (`&`, `|`, comparisons — never `if`).

### Other diagnostics (briefly)

```python
# Flux through a surface
DiagScreen(shape="plane", point=[Lx, 20.], vector=[1.,0.], direction="forward",
           every=100, species=["electron"], deposited_quantity="weight_ekin",
           axes=[["a", -10., 10., 100], ["ekin", 1e-4, 10., 100, "logscale"]])

# Per-event particle creation (ionization, pair generation)
DiagNewParticles(species="electron", every=100,
                 attributes=["x","y","px","py","pz","w"])

# Per-patch timing (essential for cluster runs with LoadBalancing)
DiagPerformances(every=200, patch_information=True)
```

## Diagnostic cadence guidelines

| Diagnostic | Sane cadence |
|---|---|
| `DiagScalar` | 10–100 |
| `DiagFields` full grid | 200–500 |
| `DiagFields` with subgrid | 100 |
| `DiagProbe` 1D | 10–50 |
| `DiagProbe` 2D | 100–500 |
| `DiagParticleBinning` | 100–500 |
| `DiagTrackParticles` | 50–200 |

Reasonable starting set for a 2D excimer run:

```python
DiagScalar(every = 50)
DiagFields(every = 500, fields = ["Ex","Ey","Bz","Rho_electron","Rho_ion"])
DiagProbe(every = 50, origin=[0., Ly/2], corners=[[Lx, Ly/2]],
          number=[400], fields=["Ex","Ey","Bz"])
DiagParticleBinning(deposited_quantity="weight_ekin", every=200,
                    species=["electron"],
                    axes=[["ekin", 1e-4, 1., 100, "logscale"]])
```

## happi — opening a simulation

Install via `make happi` in the Smilei directory.

```python
import happi

S = happi.Open("/path/to/simulation_directory")

# Multiple simulations
S = happi.Open(["sim_a", "sim_b", "sim_c"])

# Override reference frequency if missing in namelist
S = happi.Open("/path", reference_angular_frequency_SI = 2*3.14159*3e8/248e-9)

# Access namelist
print(S.namelist.Main.timestep)
print(S.namelist.Main.reference_angular_frequency_SI)
print([s.name for s in S.namelist.Species])
```

## happi — reading diagnostics

| Diagnostic block | Accessor |
|---|---|
| `DiagScalar` | `S.Scalar(name)` |
| `DiagFields` | `S.Field(index, field_name)` |
| `DiagProbe` | `S.Probe(index, field_name)` |
| `DiagParticleBinning` | `S.ParticleBinning(index)` |
| `DiagScreen` | `S.Screen(index)` |
| `DiagTrackParticles` | `S.TrackParticles(species, attributes)` |
| `DiagNewParticles` | `S.NewParticles(species)` |
| `DiagPerformances` | `S.Performances()` |

All accessors share: `.getData()`, `.getTimes()`, `.getAxis()`, `.plot()`, `.slide()`, `.animate()`, `.streak()`, `.toVTK()`, `.set()`.

```python
# Scalar time series
S.Scalar("Ukin").plot()

# Field snapshots
Ex = S.Field(0, "Ex")
Ex.slide()                                       # interactive slider
Ex.streak()                                      # 2D space-time

# Field slice
S.Field(0, "Ex", subset={"y": [20., 20.]}).slide()

# Probe streak
S.Probe(0, "Ex").streak()

# ParticleBinning with log axes
pb = S.ParticleBinning(0)
pb.set(xlog=True, ylog=True)
pb.plot()

# Selective particle tracking
tr = S.TrackParticles("electron", axes=["x","y","px"],
                      select={"px": [0.5, 10.]})
```

`getData()` returns a Python list of arrays (one per output timestep). To stack: `np.array(diag.getData())`.

## happi — unit conversion

Pint-style units, requires `reference_angular_frequency_SI`:

```python
S.Field(0, "Ex", units=["um","fs","GV/m"]).slide()
S.Field(0, "Rho_electron", units=["um","fs","cm**-3"]).slide()
S.ParticleBinning(0, units=["um","fs","keV"]).plot()
```

Without `reference_angular_frequency_SI`, conversion raises. Either re-run with rule 1 satisfied, or pass it explicitly to `happi.Open()`.

Manual conversion factors:

```python
omega_r = S.namelist.Main.reference_angular_frequency_SI
L_r = 2.99792458e8 / omega_r                     # meters
T_r = 1.0 / omega_r                              # seconds

Ex_norm = S.Field(0, "Ex").getData()
# E_r in V/m: m_e c omega_r / e
import scipy.constants as sc
E_r = sc.m_e * sc.c * omega_r / sc.e
Ex_SI = [arr * E_r for arr in Ex_norm]
```

## happi — plotting and export

```python
diag = S.Field(0, "Rho_electron")

diag.plot(timestep=500)                          # static
diag.slide()                                     # interactive slider
diag.animate(filename="density.gif", fps=10)     # GIF/MP4
diag.streak()                                    # time-space (for 1D probes)
diag.toVTK(folder="vtk_output")                  # ParaView/VisIt

diag.set(vmin=0, vmax=2.0, cmap="viridis", xlabel="x [μm]")
diag.plot()
```

## happi — parameter scan

```python
import glob
import matplotlib.pyplot as plt

for run in sorted(glob.glob("/scratch/scan_*/")):
    S = happi.Open(run)
    times = S.Scalar("Ukin").getTimes()
    data = S.Scalar("Ukin").getData()
    plt.plot(times, data, label=run.split('/')[-2])
plt.legend(); plt.show()
```

## Common mistakes

- `every = 1` on `DiagFields` — fills disk in minutes. Use `>>1` or `DiagProbe`.
- `time_average > 1` on laser fields — smears oscillations.
- Phase-space axis too narrow — bins clip silently. Check `Scalar.txt` max momenta.
- Filter using `if` — must use vectorized `&`/`|`/comparisons.
- Tracking too many particles — terabytes. Always narrow with `filter`.
- Requesting `units=...` without `reference_angular_frequency_SI` — raises.
- `subset={"y": 20.}` (scalar) — fails. Use `subset={"y": [20., 20.]}` for narrow slice.
- Reading HDF5 directly with h5py — layout changes across Smilei versions. Use happi.
