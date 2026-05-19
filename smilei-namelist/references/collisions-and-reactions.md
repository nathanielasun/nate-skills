# Collisions and reactions

This reference covers Smilei's binary collision module (Pérez 2012 scheme), collisional ionization, nuclear reactions, and bound-electron screening. For field ionization see `ionization.md`. For when to enable collisions in excimer-laser plasmas see `excimer-and-LTP-recipes.md`.

## Table of contents

1. The Pérez 2012 algorithm in two paragraphs
2. The `Collisions` block — basic syntax
3. Intra-species vs inter-species
4. The Coulomb logarithm
5. `CollisionalIonization`
6. Bound-electron screening
7. Nuclear reactions
8. Multiple `Collisions` blocks
9. Particles-per-cell requirements
10. Validation against test cases
11. Performance and load balancing
12. Diagnostics for collisions
13. Common mistakes

---

## 1. The Pérez 2012 algorithm in two paragraphs

Smilei implements binary collisions via the algorithm of Pérez, Gremillet, Decoster, Drouin & Lefebvre (Phys. Plasmas 2012). Within each cell at each timestep, macro-particles of the two specified species are paired and a single binary collision is computed per pair using a Monte-Carlo scattering angle drawn from the relativistic Coulomb cross-section. The pairing is randomized each timestep to ensure all macro-particles eventually interact with all others on the cell scale. The deflection angle is sampled to give the correct cumulative mean-square deflection over the timestep.

Two subtleties matter. **First**, the scheme handles different macro-particle weights by adjusting the effective collision probability — heavier (in weight) macro-particles collide less often per pairing but produce larger deflections, conserving the underlying physical collision frequency. **Second**, the scheme is statistical, not deterministic: collision rates have variance that scales as `1/sqrt(PPC)`. Low PPC produces noisy and biased collision rates. The practical floor for collisional runs is `PPC ≥ 32`; `PPC ≥ 64` is the safe default for excimer/LTP regimes.

## 2. The `Collisions` block — basic syntax

```python
Collisions(
    species1 = ["electron"],          # list of species names colliding from group 1
    species2 = ["ion"],               # list of species names colliding from group 2
    coulomb_log = 0.,                 # 0 → compute Coulomb log automatically; >0 → fix to this value
    debug_every = 0,                  # diagnostic output cadence (0 = off)
    ionizing = False,                 # enable collisional ionization (see §5)
)
```

This single block enables electron-ion collisions for the named pairs. To cover all relevant pairs in a plasma you typically need multiple blocks (see §8).

`species1` and `species2` are lists, not scalars, so you can group multiple species:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1", "Ar2", "Ar3"],     # electron collides against all three ion charge states
)
```

This pairs each electron with each ion species. Each pair is independently processed. For mixed-species plasmas (e.g., D-T or hydrogen-helium), grouping is the right pattern.

## 3. Intra-species vs inter-species

For collisions of a species with itself, list the same species in both groups:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["electron"],
    coulomb_log = 0.,
)
```

Smilei automatically detects this and uses a single-list pairing scheme (each particle is paired with another in the same species, no double-counting). The same block syntax handles both cases.

For excimer-class plasmas, the relevant collision channels are typically:
- e-e (electron heating, thermalization)
- e-i (energy loss, momentum transfer)
- i-i (ion thermalization, ion-acoustic-wave damping)
- e-n (electron-neutral, if neutrals are present)

Each gets its own `Collisions` block.

## 4. The Coulomb logarithm

`coulomb_log` controls the logarithm `ln Λ` that appears in the Coulomb cross-section.

```python
coulomb_log = 0.        # automatic (computed from local plasma parameters)
coulomb_log = 5.        # fixed value
```

### Automatic mode (`coulomb_log = 0.`)

Smilei computes `ln Λ` at each timestep using the local electron temperature, electron density, and the species masses/charges. The formulas follow standard plasma-physics references (NRL formulary, Lee-More for dense regimes). This is the right choice for most physics — it handles temperature and density evolution correctly.

The automatic mode requires `reference_angular_frequency_SI` to be set, because the Coulomb logarithm depends on absolute densities and temperatures. This is hard-rule 1 from `SKILL.md` — without it, the formula uses normalized units that don't correspond to physical Λ.

### Fixed mode (`coulomb_log > 0.`)

When set to a positive value, Smilei skips the automatic computation and uses your fixed `ln Λ` everywhere. Use cases:
- Benchmarking against analytical solutions that assume constant `ln Λ`
- Reproducing legacy results that used a fixed value
- Studying sensitivity to `ln Λ` choice
- Cases where the automatic formula is suspected to misfire (very low temperature, very high density)

Typical values are `ln Λ ≈ 3–10` for laser-plasma conditions, `ln Λ ≈ 10–20` for fusion-relevant conditions.

### `coulomb_log_factor`

A multiplicative factor on the computed `ln Λ`, useful for parametric studies:

```python
coulomb_log = 0.            # automatic
coulomb_log_factor = 2.0    # scale the result by 2 — non-physical but useful for sensitivity studies
```

Defaults to 1.0. Setting it to a value other than 1 is not physical and should only be used for testing.

## 5. `CollisionalIonization`

When a `Collisions` block has `ionizing = True` and `species2` is a neutral or partially-ionized species with `atomic_number > 0`, electron impact ionization is enabled. The mechanism uses Lotz-Drawin cross-sections internally.

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1"],
    coulomb_log = 0.,
    ionizing = "electron",                # electrons go into this species
)
```

Wait — `ionizing` takes a string, not a Boolean. The value is the *name of the electron species* that receives the impact-ionized electrons. (Different from field ionization, where `ionization_electrons` is on the `Species` block.)

For a full ionization ladder, define each ion charge state as a separate `Species` and chain the collisions:

```python
# Define species
Species(name = "electron", ...)
Species(name = "Ar0", atomic_number = 18, charge = 0., ...)
Species(name = "Ar1", atomic_number = 18, charge = 1., ...)
Species(name = "Ar2", atomic_number = 18, charge = 2., ...)
# ...

# Cross-link the chain: each ionization step needs its own Collisions block
Collisions(species1=["electron"], species2=["Ar0"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar1"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar2"], ionizing="electron")
```

This is verbose but explicit. The macro-particle that gets ionized doesn't move between species — the `charge` field on the macro-particle is incremented. The neutral species (`Ar0`) gradually depopulates as macro-particles climb the ladder, but they remain in the same data structure.

**Important behavioral note**: Collisional ionization typically requires `Ionization` (field) to also be configured if the field is strong enough to ionize independently. They are independent channels; neither subsumes the other. For excimer-class plasmas at `I ~ 10¹²–10¹³ W/cm²`, both should be on.

## 6. Bound-electron screening

Bound electrons partially screen the nuclear charge seen by a colliding incident electron. For partially-ionized plasmas this materially affects cross-sections. Smilei provides an option:

```python
Collisions(
    species1 = ["electron"],
    species2 = ["Ar1"],
    coulomb_log = 0.,
    every = 1,                          # apply collisions every N timesteps (default 1)
)
```

Bound-electron screening is automatic for `CollisionalIonization`-relevant pairs. For partially-ionized plasmas the effect is a few percent on the Coulomb log — important when comparing to experiments but rarely a leading-order effect.

For more nuanced control, consult the current Smilei documentation; the implementation has been refined in recent versions and the exact knobs vary.

## 7. Nuclear reactions

For DD, DT, or other fusion-relevant reactions, Smilei provides a nuclear-reaction option built on top of the binary-collision pairing:

```python
Collisions(
    species1 = ["deuteron"],
    species2 = ["deuteron"],
    coulomb_log = 0.,
    nuclear_reaction = ["alpha", "neutron"],    # reaction products
    every = 100,                                # rate of pairing for reactions
)
```

The `nuclear_reaction` argument lists the product species names; corresponding `Species` blocks must be defined elsewhere (typically as `mass = 0`-like daughter products, or with appropriate masses for alpha particles).

Smilei supports a small set of reactions: D(d,n)³He, D(d,p)T, D(t,n)⁴He, p-B fusion (with appropriate isotope species), and similar standard fusion channels. Cross-sections are tabulated from ENDF or similar sources.

This is specialized — most excimer-laser plasma work won't use it. For ICF-relevant simulations involving fusion burn, consult the latest documentation for the precise reaction set.

## 8. Multiple `Collisions` blocks

To enable all relevant collision channels in a plasma, you typically need several blocks:

```python
# Pure electron-ion plasma (e.g., H+ at high temperature)
Collisions(species1=["electron"], species2=["electron"])
Collisions(species1=["electron"], species2=["ion"])
Collisions(species1=["ion"],      species2=["ion"])
```

```python
# Argon plasma with multiple charge states
Collisions(species1=["electron"], species2=["electron"])
Collisions(species1=["electron"], species2=["Ar1","Ar2","Ar3"])
Collisions(species1=["Ar1","Ar2","Ar3"], species2=["Ar1","Ar2","Ar3"])

# Collisional ionization ladder
Collisions(species1=["electron"], species2=["Ar1"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar2"], ionizing="electron")
Collisions(species1=["electron"], species2=["Ar3"], ionizing="electron")
```

Each block adds compute cost roughly proportional to `PPC × number_of_cells`. For a typical excimer setup this can be a 10–30% overhead on the total simulation cost.

### `every`

```python
Collisions(species1=["electron"], species2=["ion"], every = 10)
```

Apply this collision block only every N timesteps. The cross-section is internally adjusted so the time-averaged rate matches a "every-step" calculation. Useful for:
- Reducing cost on cheap-to-skip channels (e.g., heavy-heavy collisions whose effect builds slowly)
- Matching collision-time scales to timestep when collisions are fast compared to `Δt`

Default is 1 (apply every timestep). Setting `every > 1` is an optimization; never use it unless you've verified the underlying timescale separation.

## 9. Particles-per-cell requirements

Collisional accuracy depends on having enough macro-particles per cell to sample the collision distribution. Minimum guidelines:

| Channel | Minimum PPC | Recommended PPC |
|---|---|---|
| Pure collisionless | — | 8–16 |
| e-i collisions only | 16 | 32 |
| Full collision set (e-e, e-i, i-i) | 32 | 64 |
| Collisional ionization | 32 | 64 |
| Multi-charge-state ladder | 32 per state | 64+ per state |
| Nuclear reactions | 128+ | 256+ (low cross-section) |

If you see noisy or oscillatory collision-driven energy transfer in diagnostics, the most common cause is too few macro-particles. Doubling PPC and rerunning is the diagnostic.

## 10. Validation against test cases

Smilei's collision module has been validated against:
- Spitzer-Härm electron-ion temperature equilibration
- Thomas-Fermi ionization equilibrium for various Z (see Fig. 25 of the docs)
- Electron-neutral scattering against ELSEPA (see Fig. 26)
- Coulomb-log-dependent stopping power
- Bremsstrahlung electron-electron energy loss

If your physics depends critically on accurate collisions (e.g., excimer-class collisional absorption coefficients), run a benchmark of your collision configuration against analytical results before trusting full simulations. The simplest benchmark is two-temperature relaxation in a uniform plasma: initialize electrons and ions at different temperatures, evolve with only collisions enabled, and check the relaxation rate against Spitzer's formula.

## 11. Performance and load balancing

Binary collisions are an embarrassingly local operation — they happen within each cell with no neighbor communication. Performance impact:

- **Compute cost**: roughly proportional to `PPC²` per cell per block per timestep. For PPC = 64 and 3 collision blocks, this can be 20–40% of total step cost.
- **Memory**: minimal extra memory beyond the macro-particles themselves.
- **Load balance**: cells with more particles have proportionally more collision work. **Enable `LoadBalancing`** for any collisional run — without it, cells with high density (typical in laser-plasma fronts) will be 10× slower than vacuum cells, idling many MPI ranks.

```python
LoadBalancing(every = 50, cell_load = 1.0, frequency_load = 0.1)
```

This block redistributes patches across MPI ranks every 50 steps based on per-cell work weighted by `cell_load` (per-particle cost) and `frequency_load` (per-cell overhead). For collisional runs the default values are usually fine.

## 12. Diagnostics for collisions

```python
Collisions(
    species1 = ["electron"], species2 = ["ion"],
    coulomb_log = 0.,
    debug_every = 100,                  # write diagnostic data every 100 steps
)
```

When `debug_every > 0`, Smilei writes per-cell collision diagnostics (collision frequency, total energy transferred, computed `ln Λ`) to a separate HDF5 file. Available via happi as a `DiagPerformances`-like accessor in some versions; in others as a custom diagnostic. Useful for verifying that your collision configuration produces expected rates.

For overall energy conservation, watch `DiagScalar` outputs `Ukin` (kinetic energy by species) and `Uelm` (EM energy). In a closed system without external sources, `Ukin + Uelm` should be conserved; if it isn't, the collision algorithm is either too coarsely resolved (low PPC) or there's a bug.

## 13. Common mistakes

### Mistake 1: forgetting `reference_angular_frequency_SI`

This is hard-rule 1. Without it, `coulomb_log = 0.` (automatic) computes nonsense because the formula needs absolute physical units. Fix: set `Main.reference_angular_frequency_SI = 2*pi*c/lambda`.

### Mistake 2: only enabling e-i collisions, not e-e or i-i

```python
Collisions(species1=["electron"], species2=["ion"])
# Missing: e-e and i-i
```
In a dense plasma, e-e collisions drive electron thermalization faster than e-i. Omitting them means the electron distribution stays artificially non-Maxwellian. For physically meaningful collisional plasma, enable all three.

### Mistake 3: low PPC with collisions enabled

```python
Species(..., particles_per_cell = 8)
Collisions(species1=["electron"], species2=["ion"])
```
PPC = 8 is too few for reliable collision statistics. The collision rate will have ~35% variance per cell per step. Symptoms: spurious temperature spikes, oscillatory energy transfer. Fix: bump PPC to ≥ 32.

### Mistake 4: `ionizing = True` (Boolean) instead of `ionizing = "electron"` (species name)

```python
Collisions(..., ionizing = True)        # historical syntax; depending on version, may or may not work
```
Current versions expect `ionizing` to be the *name* of the electron species. Use the explicit string.

### Mistake 5: missing ion species in the ionization ladder

```python
Species(name = "Ar0", atomic_number = 18, charge = 0., ...)
# Defined Ar0 but not Ar1, Ar2, ...
Collisions(species1=["electron"], species2=["Ar0"], ionizing="electron")
```
Smilei doesn't auto-create the next charge state. You must define `Ar1`, `Ar2`, ..., explicitly. Without them, ionization from `Ar0` raises at the first ionizing collision.

### Mistake 6: `coulomb_log_factor` left at non-1 from a sensitivity study

A common debugging pattern is to set `coulomb_log_factor = 0.1` to weaken collisions and check whether collisions are responsible for an observed effect. Forgetting to reset to 1.0 before a production run silently produces wrong physics. Always check this value is 1.0 (or absent) in production namelists.

### Mistake 7: applying intra-species collisions without recognizing the cost

```python
Collisions(species1=["Ar1","Ar2","Ar3"], species2=["Ar1","Ar2","Ar3"])
```
This expands to 9 pair combinations (Ar1-Ar1, Ar1-Ar2, Ar1-Ar3, Ar2-Ar1, etc.). For dense ladders this becomes prohibitive. For ion-ion collisions among heavy species, consider whether all pairs are physically relevant — often only same-species or adjacent-state pairs matter.

### Mistake 8: forgetting `LoadBalancing`

```python
# Collisional run without LoadBalancing
Main(..., load_balancing = ... )                # not set
```
Without dynamic load balancing, MPI ranks idle waiting for slow ranks holding the dense plasma. Easy to detect with `DiagPerformances`. Always enable `LoadBalancing` for collisional or ionizing runs.

### Mistake 9: `every` set too high

```python
Collisions(..., every = 1000)
```
The internal rate adjustment is correct only when `every × Δt` is short compared to the actual collision time. For `every = 1000` and `Δt ~ 1/ω_r`, the gap is ~10³ optical periods — often longer than the collision time, breaking the assumption. Generally stick with `every = 1` unless you've explicitly verified the timescale separation.

### Mistake 10: mixing `Collisions` and an envelope-laser run

The collision algorithm operates on the actual macro-particle momenta, which for envelope-laser runs are the slow (cycle-averaged) momenta — not the fast quiver. Collision rates computed from these can be slightly off when ponderomotive heating matters. For envelope + collisions, results are usually still qualitatively right but quantitative comparisons should be cross-checked against full-EM benchmarks.
