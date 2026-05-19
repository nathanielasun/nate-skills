# Solvers and fields — source layer

Maxwell solvers (Yee, Lehe, Cowan, Bouchard), the `ElectroMagn` class hierarchy, PML implementation, field interpolation, and current deposition. The non-particle half of the PIC loop.

## Table of contents

1. `ElectroMagn` — the field container
2. Geometry specializations
3. The Yee staggered grid
4. Maxwell solvers — overview
5. `MaxwellSolver_Yee` — the default
6. `MaxwellSolver_Lehe` — relativistic-LWFA optimized
7. `MaxwellSolver_Cowan`, `MaxwellSolver_Bouchard`
8. Ghost cells
9. EM boundary conditions
10. PML implementation
11. Current deposition (Esirkepov)
12. Adding a new Maxwell solver
13. Common pitfalls

---

## 1. `ElectroMagn` — the field container

`src/ElectroMagn/ElectroMagn.h` is the base class. It owns the EM field arrays for one patch:

```cpp
class ElectroMagn {
public:
    Field* Ex;
    Field* Ey;
    Field* Ez;
    Field* Bx;
    Field* By;
    Field* Bz;
    Field* Jx;
    Field* Jy;
    Field* Jz;
    Field* rho;

    // Per-species currents and densities (for diagnostics)
    std::map<std::string, Field*> rho_s;
    std::map<std::string, Field*> Jx_s, Jy_s, Jz_s;

    // PML auxiliary fields (allocated if PML boundaries)
    std::vector<Field*> PML_fields;

    // Methods
    virtual void initElectroMagn(...) = 0;
    virtual void solveMaxwell(...) = 0;
    virtual void centerMagneticFields() = 0;
    virtual void saveMagneticFields() = 0;             // for time-staggered output
    virtual void applyEMBoundaryConditions(...) = 0;
};
```

**Key points:**

- One `ElectroMagn` instance per patch.
- `Field*` is `src/Field/Field.h` — a wrapper around a `std::vector<double>` with shape metadata.
- The `Jx_s` etc. maps allocate per-species current arrays only if `DiagFields` requests them.
- PML fields are allocated lazily — only on patches whose boundary touches a PML face.

## 2. Geometry specializations

`ElectroMagn` is abstract; concrete derived classes per geometry:

| File | Geometry |
|---|---|
| `ElectroMagn1D.{h,cpp}` | 1D cartesian |
| `ElectroMagn2D.{h,cpp}` | 2D cartesian |
| `ElectroMagn3D.{h,cpp}` | 3D cartesian |
| `ElectroMagnAM.{h,cpp}` | AM cylindrical (modes 0..M) |

The geometry classes implement `initElectroMagn` (allocates `Field*` instances with the right shape), `applyEMBoundaryConditions` (calls into geometry-specific BCs), and per-geometry helpers.

**AM cylindrical** is the most complex: instead of one `Ex` field, it stores `Ex_m` (complex field) per azimuthal mode `m = 0..M-1`. Look at `ElectroMagnAM.cpp` carefully before modifying — the indexing is non-trivial.

The Maxwell *solver* (next section) is geometry-dispatched at factory time. The same physics (Yee, Lehe, etc.) has separate implementations per geometry.

## 3. The Yee staggered grid

Smilei uses the **Yee grid** — fields are staggered by half a cell:

```
2D layout (cell at origin):

  Ey(0, 0.5) ──── Ey(1, 0.5) ────
  │              │
  Bz(0.5, 0.5)   Bz(1.5, 0.5)
  │              │
  Ey(0, 0) ───── Ex(0.5, 0) ── Ey(1, 0) ──── Ex(1.5, 0) ──
```

- `Ex` lives at integer-y, half-x positions
- `Ey` lives at half-y, integer-x positions
- `Bz` lives at half-x, half-y (cell centers in 2D)
- Density `rho` lives at integer-integer (cell corners)
- Currents `Jx`, `Jy` mirror their corresponding E-field locations

**Why staggered**: the second-order centered finite difference of `B` on the Yee grid gives `curl B` at the E-field location with second-order accuracy and no DC mode error. This is the foundation of Yee's algorithm.

**Indexing convention in code**: `Ex(i, j)` is stored at array index `[i + nx*j]` (row-major) — but the *physical* location is `(i*dx + dx/2, j*dy)` for `Ex` (offset by half a cell in `x`). The offset is implicit; the array index doesn't carry it. When implementing physics that needs the physical location of a field point, compute it from the cell index plus the convention table.

`src/Field/Field2D.h` and `Field3D.h` have helpful comments documenting the offsets.

## 4. Maxwell solvers — overview

Smilei provides four Maxwell solvers selected by `Main.maxwell_solver`:

| String | Solver | Files | Use |
|---|---|---|---|
| `"Yee"` | Standard Yee | `MaxwellSolver_Yee_*.cpp` | Default — robust, second-order |
| `"Lehe"` | Lehe (cross-direction term) | `MaxwellSolver_Lehe_*.cpp` | Relativistic LWFA — kills numerical Cherenkov |
| `"Cowan"` | Cowan (improved dispersion) | `MaxwellSolver_Cowan_*.cpp` | Reduced grid anisotropy |
| `"Bouchard"` | Bouchard (pseudo-spectral-like) | `MaxwellSolver_Bouchard_*.cpp` | Higher-order, more expensive |

Each is a class deriving from `MaxwellSolver`:

```cpp
class MaxwellSolver {
public:
    virtual void operator()(ElectroMagn* EMfields) = 0;
};
```

The solver is invoked once per timestep per patch in `Patch::solveMaxwell()`.

## 5. `MaxwellSolver_Yee` — the default

`src/ElectroMagn/MaxwellSolver/MaxwellSolver_Yee_2D.cpp` (sketch):

```cpp
void MaxwellSolver_Yee_2D::operator()(ElectroMagn* EMfields) {
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz);
    Field2D* Jx2D = static_cast<Field2D*>(EMfields->Jx);
    Field2D* Jy2D = static_cast<Field2D*>(EMfields->Jy);

    // Bz update from curl E
    for (i = 0; i < nx-1; ++i)
        for (j = 0; j < ny-1; ++j) {
            (*Bz2D)(i, j) -= dt_over_dx * ((*Ey2D)(i+1, j) - (*Ey2D)(i, j))
                            - dt_over_dy * ((*Ex2D)(i, j+1) - (*Ex2D)(i, j));
        }

    // Ex, Ey update from curl B and J
    for (i = 0; i < nx; ++i)
        for (j = 1; j < ny; ++j) {
            (*Ex2D)(i, j) += dt_over_dy * ((*Bz2D)(i, j) - (*Bz2D)(i, j-1))
                           - dt * (*Jx2D)(i, j);
        }
    // ... similar for Ey ...
}
```

Two stencils, each second-order. The leapfrog scheme advances B by `dt/2`, then E by `dt`, then B by another `dt/2` — but in Smilei this is rolled into one call per step using a half-step convention. `centerMagneticFields()` interpolates `B` to integer times for diagnostics.

**CFL**: `dt² × (1/dx² + 1/dy² + 1/dz²) ≤ 1` for Yee. With `timestep_over_CFL = 0.95`, the actual `dt = 0.95 × dt_CFL_max`.

**Ghost cells**: standard Yee needs 2 ghost cells per face (one for E, one for B, but in practice 2 suffices for all stencils).

## 6. `MaxwellSolver_Lehe` — relativistic-LWFA optimized

The Lehe solver adds a cross-direction term to the standard Yee that cancels numerical Cherenkov radiation from relativistic particles. The math:

```
∂_t E = c² (∂_x B_z - β_lehe ∂_y² ∂_x B_z) - J
```

The extra term `β_lehe ∂_y² ∂_x B_z` exactly cancels the numerical-Cherenkov instability when `β_lehe = 1/8` (in 2D).

`src/ElectroMagn/MaxwellSolver/MaxwellSolver_Lehe_2D.cpp` implements this. Code structure is similar to Yee but with one extra finite-difference term per E-component update.

**Cost**: ~30% more arithmetic per cell than Yee.
**Ghost cells**: 3 per face (the extra finite-difference reaches further).
**CFL**: same as Yee (the extra term is a perturbation).

Use for LWFA simulations with relativistic injected beams. Avoid for excimer or low-`a₀` work — the extra cost isn't justified and the dispersion is slightly different.

## 7. `MaxwellSolver_Cowan`, `MaxwellSolver_Bouchard`

**Cowan**: a higher-order finite-difference that reduces grid anisotropy. Useful when waves propagate at 45° to the grid; the standard Yee has visible angular dispersion. Cost: ~2× Yee. Ghost cells: 3.

**Bouchard**: a pseudo-spectral-like solver in the spatial directions. Higher-order spatial accuracy at the cost of larger stencil. Useful for accurate dispersion relations. Cost: ~4× Yee. Ghost cells: 4–5.

Both are niche — most production work uses Yee or Lehe. Look at the source files for the exact stencils if you need to modify these.

## 8. Ghost cells

Each patch has a halo of "ghost cells" outside its interior. These hold the values from neighboring patches needed for the Maxwell stencil. After each Maxwell step, ghost cells are stale and must be refreshed via inter-patch communication.

**Per-solver ghost-cell count:**

| Solver | Ghost cells per face |
|---|---|
| Yee | 2 |
| Lehe | 3 |
| Cowan | 3 |
| Bouchard | 4–5 |

Set in `Params::PicParams` from the chosen solver name. The number is checked at startup — if your patch is smaller than the ghost extent, Smilei refuses to start.

**Ghost-cell exchange**: handled by `SmileiMPI::exchangeFields()`. Adjacent patches send their interior edges to the neighbor's ghost region. The communication is non-blocking + waitall.

**Patches at the simulation boundary** have ghost cells filled by boundary conditions, not communication.

## 9. EM boundary conditions

`src/ElectroMagnBC/` per-geometry:

| BC | File | Effect |
|---|---|---|
| `"silver-muller"` | `ElectroMagnBC*_Refl.cpp`, `*SM.cpp` | First-order absorbing |
| `"PML"` | `ElectroMagnBC*_PML.cpp` | Perfectly Matched Layer (see §10) |
| `"periodic"` | (handled at patch-exchange level) | Wraps |
| `"reflective"` | `ElectroMagnBC*_Refl.cpp` | Perfect conductor |

The base class:

```cpp
class ElectroMagnBC {
public:
    virtual void apply(ElectroMagn* EMfields, double time, ...) = 0;
};
```

Called from `Patch::applyEMBoundaryConditions()` after the Maxwell solver. The boundary condition writes to ghost cells (for absorbing) or modifies the boundary face (for reflective/conductor).

**Laser injection** is implemented as a special boundary condition on top of Silver-Müller or PML. The injected `B` field at the boundary is the laser profile; the absorbing BC handles the reflected/outgoing wave separately.

## 10. PML implementation

The Perfectly Matched Layer is a region of extra cells beyond the simulation interior where the wave equation has artificial damping that attenuates outgoing waves. Reflections from a well-tuned PML are typically below `-60 dB`.

`src/ElectroMagnBC/ElectroMagnBC2D_PML.cpp` (2D Cartesian) shows the pattern:

- PML allocates auxiliary fields `E_pml`, `B_pml`, and stretched-coordinate parameters.
- The Maxwell update in the PML region uses stretched coordinates `(σ_x, σ_y)` that grow inside the layer.
- The damping profile is typically polynomial: `σ(d) = σ_max (d/L_PML)^n` for some order `n`.

**Configuration in the namelist:**

```python
Main(
    EM_boundary_conditions = [["PML"], ["PML"]],
    number_of_pml_cells = [[12, 12], [12, 12]],   # per face per direction
)
```

12 cells is typical. More cells reduce reflection at the cost of compute and memory.

**When to use PML over Silver-Müller:**

- Long-time simulations where reflections accumulate
- Oblique-incidence lasers (Silver-Müller reflects strongly off-normal)
- Dense plasmas where overdense reflections must be cleanly absorbed

For straightforward UHI or sub-relativistic excimer work with normal-incidence lasers, Silver-Müller suffices.

## 11. Current deposition (Esirkepov)

`src/Projector/` implements the Esirkepov current-deposition scheme. The key property: `∂_t ρ + ∇·J = 0` is satisfied exactly on the Yee grid, so no charge-correction step (no Poisson solve) is needed during the time loop.

The Esirkepov shape function for a particle is the difference of the new and old particle position's shape functions, integrated over the timestep. For a particle moving from `(x₀, y₀)` to `(x₁, y₁)` in `dt`:

```
W_x[i,j] = (S_new[i,j] - S_old[i,j])           (in x-direction support)
J_x[i+1/2, j] += charge × W_x[i,j] / dt
```

Concrete files: `Projector2D2Order.cpp`, `Projector2D4Order.cpp`, `Projector3D2Order.cpp`, etc. The order corresponds to the spline order of the particle shape function.

**Per-particle cost**: ~50 floating-point operations in 2D, ~100 in 3D. The hottest loop in many simulations.

**SIMD**: the loop over particles is `#pragma omp simd`. Vectorization is essential. For complex stencils (high-order Esirkepov), the compiler may not auto-vectorize cleanly — check assembly if performance is suspicious.

**GPU**: `Projector*` has `_GPU_KERNEL_` paths. Per-particle deposition on GPU requires atomic adds to the grid (multiple particles in the same cell collide), which is the main performance bottleneck.

## 12. Adding a new Maxwell solver

1. **Pick** a name (e.g., `"MySolver"`).
2. **Create** `src/ElectroMagn/MaxwellSolver/MaxwellSolver_MySolver_2D.{h,cpp}` and 3D, AM as needed.
3. **Derive** from `MaxwellSolver`.
4. **Implement** `operator()(ElectroMagn*)` with the geometry-specific update stencil.
5. **Register** in `src/ElectroMagn/MaxwellSolver/MaxwellSolverFactory.h`:
   ```cpp
   else if (params.maxwell_solver == "MySolver") {
       return new MaxwellSolver_MySolver_2D(params);
   }
   ```
6. **Declare** the ghost-cell count in `Params::PicParams::setMaxwellSolver(...)`.
7. **Update** CFL in `PicParams::computeTimeStep()` if your solver's CFL differs.
8. **Document** in `doc/Sphinx/Source/namelist.rst`.
9. **Test** — add a `validation/tst_*.py` with the new solver and a reference output.

If your solver uses fewer ghost cells, you may improve compute efficiency. If more, patches must be at least as large as the ghost extent.

## 13. Common pitfalls

**Wrong field location**: physical position of `Ex(i,j)` is `(i*dx + dx/2, j*dy)`, not `(i*dx, j*dy)`. The half-cell offset bites when implementing new physics that needs absolute positions.

**Forgetting to update ghost cells**: any code that reads from neighbor patches must do so after `exchangeFields()`. Reading stale ghost values produces silent errors that grow with patch boundary length.

**Patch smaller than ghost extent**: a patch with 4 interior cells per direction and a Lehe solver needing 3 ghost cells per face has only 4 - 6 = -2 effective cells. Smilei refuses to start. Fix: larger patches or fewer patches.

**Atomic operations in current deposition on GPU**: omitting atomic adds gives wrong physics when particles in the same cell race. Always atomic on GPU; CPU can use sequential (per-particle in serial inside a patch).

**Mixing solver order with interpolation order**: `Main.interpolation_order = 2` and `maxwell_solver = "Bouchard"` (which is high-order spatially) is inconsistent — the spatial accuracy of the field-particle coupling is limited by the lower order. Use matching orders for clean convergence studies.

**Adding J directly without sign convention**: `Jx` in Smilei is in units where the Maxwell update is `∂_t E = c² curl B - J / ε₀_norm`. The sign and prefactor are baked into the update; new deposition code must match.

**PML on a periodic boundary**: meaningless and causes parse error. PML implies open boundary, periodic implies closed.

**Modifying `Field*` size at runtime**: field arrays are sized at patch construction. Don't `resize()` during the simulation — many invariants depend on stable sizes.

**AM mode coupling**: in cylindrical AM, modes evolve independently in vacuum but couple through nonlinear sources (current density depends on particle motion, which depends on all modes). When adding a new physics module for AM, think carefully about how to project to/from modes.
