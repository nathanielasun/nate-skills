# Collisions — source layer

The Pérez 2012 binary-collision implementation in Smilei: `Collisions` class, `BinaryProcesses` wrapper, particle pairing, intra/inter-species detection, collisional ionization, nuclear reactions, and how to add a new pair process.

## Table of contents

1. `Collisions` and `BinaryProcesses` — class layout
2. The Pérez 2012 algorithm in code
3. Particle pairing — intra vs inter-species
4. Coulomb logarithm — auto vs fixed
5. `CollisionalIonization` — implementation
6. `CollisionalNuclearReaction` — implementation
7. The `every` argument and rate rescaling
8. Per-patch state and thread safety
9. Diagnostics and debug output
10. Adding a new pair process
11. Performance considerations
12. Common pitfalls

---

## 1. `Collisions` and `BinaryProcesses` — class layout

`src/Collisions/` contains the collision module. The key classes:

```cpp
class Collisions {
public:
    std::vector<unsigned int> species_group1;     // indices into Patch::vecSpecies
    std::vector<unsigned int> species_group2;
    double clog;                                   // 0 = auto, >0 = fixed
    int every;                                     // apply every N steps
    bool intra_collisions;                         // species1 ∩ species2 nonempty?

    // For collisional ionization
    bool ionizing;
    int ionization_electrons_index;                // species index for receiving electrons

    // For nuclear reactions
    CollisionalNuclearReaction* nuclear_reaction;

    void apply(Patch* patch, int itime);
};

class BinaryProcesses {
public:
    std::vector<Collisions*> coll;

    void apply(Patch* patch, int itime);           // dispatches to each Collisions
};
```

`BinaryProcesses` is a thin wrapper that owns a list of `Collisions` instances (one per namelist `Collisions(...)` block) and dispatches them in order. Each `Collisions` instance handles one species-pair specification.

**Per-patch invocation**: `BinaryProcesses::apply` is called from `Patch::dynamics`, inside the per-patch OpenMP-parallel loop. Each patch processes its own particles independently — no cross-patch synchronization within the collision module.

## 2. The Pérez 2012 algorithm in code

The implementation lives primarily in `Collisions::apply(Patch*, int itime)`. Pseudocode:

```cpp
void Collisions::apply(Patch* patch, int itime) {
    if (itime % every != 0) return;

    // For each cell in the patch:
    for (each_cell c in patch) {
        // 1. Gather particle indices in this cell
        std::vector<size_t> idx1 = particles_in_cell(species_group1, c);
        std::vector<size_t> idx2 = particles_in_cell(species_group2, c);

        // 2. Random shuffle for pairing
        patch->rand->shuffle(idx1);
        patch->rand->shuffle(idx2);

        // 3. Pair particles (intra vs inter detection)
        if (intra_collisions) {
            // Pair idx1[2k] with idx1[2k+1]
            n_pairs = idx1.size() / 2;
        } else {
            // Pair idx1[k] with idx2[k]
            n_pairs = min(idx1.size(), idx2.size());
        }

        // 4. For each pair, compute scattering angle and update momenta
        for (each pair p) {
            // Pérez formula: relativistic Coulomb scattering angle
            // depending on local plasma parameters and Coulomb log
            double scattering_angle = sample_perez_angle(...);
            update_momenta_post_collision(p, scattering_angle);
        }
    }
}
```

The exact source is in `Collisions.cpp::apply()`. The function is long (~500 lines) because of:
- Geometry-specific cell indexing (1D, 2D, 3D, AM)
- Different code paths for inter vs intra collisions
- Auto vs fixed Coulomb log
- Optional collisional ionization (`ionizing`)
- Optional nuclear reactions (`nuclear_reaction != nullptr`)

When modifying, find the section relevant to your physics and treat the rest as orthogonal.

**Energy and momentum conservation**: each binary collision conserves total energy and total momentum of the pair (by construction — the scattering angle is sampled in the COM frame and momenta are rotated, preserving magnitudes). The algorithm is exactly conservative pair-by-pair, modulo floating-point round-off.

## 3. Particle pairing — intra vs inter-species

`intra_collisions` is set at construction by checking whether `species_group1` and `species_group2` have any common indices:

```cpp
// In Collisions::Collisions(...) constructor:
intra_collisions = false;
for (int s1 : species_group1)
    for (int s2 : species_group2)
        if (s1 == s2) { intra_collisions = true; break; }
```

This is purely a code-path optimization — the physics is the same. Intra-collisions avoid double-counting and use a single-list random pairing.

**Multi-species groups**: when `species1 = ["Ar1", "Ar2", "Ar3"]`, all three are merged into one virtual group. The cell-level pairing draws from the combined particle list. Each pair may be Ar1-Ar1, Ar1-Ar2, etc. — the physics correctly handles different masses and charges per particle.

**Weight handling**: macro-particle weights differ when species have different `particles_per_cell` or non-uniform density profiles. Pérez's algorithm handles this by adjusting collision probability inversely with the heavier macro-particle's weight — a heavy macro-particle collides less often per pairing but produces larger deflections.

```cpp
// Sketch from Collisions.cpp
double w1 = particles1->weight[i1];
double w2 = particles2->weight[i2];
double w_max = std::max(w1, w2);
double w_min = std::min(w1, w2);

double prob = base_collision_prob * w_min / w_max;
if (patch->rand->uniform() < prob) {
    // Apply collision to both particles
}
```

This preserves the underlying physical collision frequency.

## 4. Coulomb logarithm — auto vs fixed

In auto mode (`clog == 0`), the formula uses local plasma parameters:

```cpp
double computeCoulombLog(double T_norm, double n_norm,
                         double mass1, double charge1,
                         double mass2, double charge2,
                         double omega_r_SI) {
    // Convert to SI
    double n_SI = n_norm * N_r_SI(omega_r_SI);
    double T_SI = T_norm * 511e3 * 1.602e-19;     // T in Joules

    // Standard NRL formulary formula
    double lambda_D = std::sqrt(eps0 * T_SI / (n_SI * e2));
    double b_min   = ...;                          // min impact parameter (quantum or classical)

    return std::log(lambda_D / b_min);
}
```

The actual implementation in `Collisions.cpp::computeCoulombLog` handles relativistic regimes, quantum-vs-classical regime selection, and the modifications from Lee & More for dense plasmas. **Always set `reference_angular_frequency_SI`** — the formula requires SI units.

In fixed mode (`clog > 0`), the function returns `clog` directly without computation.

## 5. `CollisionalIonization` — implementation

When `ionizing = "electron"` (a species name string), an extra step runs after the standard collision:

```cpp
// In Collisions::apply, after the standard binary collision:
if (ionizing && (particles2 is an ion species)) {
    double E_inc = compute_incident_kinetic_energy(...);
    double U_i   = lookup_ionization_potential(atomic_number, current_charge_state);

    if (E_inc > U_i) {
        // Lotz-Drawin cross-section
        double sigma = lotz_drawin_cross_section(E_inc, U_i, ...);
        double prob_ionize = sigma * v_rel * dt;

        if (patch->rand->uniform() < prob_ionize) {
            // Increment ion's charge state
            particles2->charge[i2] += 1;

            // Create new electron at ion's position
            int new_idx = electron_species->particles->size();
            electron_species->particles->resize(new_idx + 1);
            electron_species->particles->position[0][new_idx] = particles2->position[0][i2];
            // ... (momentum drawn from energy budget) ...

            // Energy bookkeeping: subtract U_i from incident electron, add to freed electron's KE
            update_incident_kinetic_energy(particles1, i1, -U_i);
        }
    }
}
```

The ionization potential `U_i(atomic_number, charge_state)` is looked up from an internal table (in `src/Ionization/IonizationTables.cpp`). Same table used by field ionization.

**The freed electron's species** is determined at construction time by parsing the `ionizing` string and storing `ionization_electrons_index`. The lookup `electron_species = patch->vecSpecies[ionization_electrons_index]` happens once per collision call.

## 6. `CollisionalNuclearReaction` — implementation

When `nuclear_reaction = ["alpha", "neutron"]` etc., a `CollisionalNuclearReaction*` is allocated and called from `Collisions::apply` for each pair:

```cpp
if (nuclear_reaction) {
    double E_com = compute_center_of_mass_energy(...);
    double sigma = nuclear_reaction->crossSection(E_com);
    double prob_react = sigma * v_rel * dt;

    if (patch->rand->uniform() < prob_react) {
        // Consume the two parent particles
        particles1->eraseParticle(i1);
        particles2->eraseParticle(i2);

        // Create products with appropriate momenta
        nuclear_reaction->createProducts(patch, COM_position, COM_momentum, ...);
    }
}
```

Cross-sections come from ENDF or similar tabulated data, loaded at startup. Smilei supports D-D, D-T, D-He³, p-B fusion (with appropriate isotope species), and a small set of other reactions. Adding a new nuclear reaction requires:

1. Cross-section data (in `src/Collisions/NuclearReactionData/`)
2. A new `CollisionalNuclearReaction` subclass
3. Registration in `NuclearReactionFactory`
4. Validation against published data

This is specialized work; most contributors won't touch it.

## 7. The `every` argument and rate rescaling

When `every > 1`, the collision is applied every `every` timesteps. To preserve the physical rate, the cross-section is internally multiplied by `every`:

```cpp
double effective_prob = sigma * v_rel * (every * dt);
```

This is valid only when `every × dt` is much shorter than the actual collision time. For typical fast collisions in dense plasma, `every = 1` is required for correct physics. The `every > 1` mode is an optimization for slowly-colliding channels (heavy-heavy in dilute plasma).

## 8. Per-patch state and thread safety

Each patch has its own:

- Random number generator (`patch->rand`)
- Particle arrays for each species (`patch->vecSpecies[s]->particles`)
- Diagnostic buffers

**Within a patch**, the collision loop is single-threaded — it modifies particle momenta directly, which would race if parallelized over cells. The OpenMP parallelism is at the patch level: each thread processes one patch independently.

**Across patches**, the modifier of `vecSpecies[s]->particles` for patch P doesn't touch patch Q's particles, so no race. Inter-patch particle exchange (after the collision step but before the next step) is handled by `SmileiMPI::exchangeParticles`.

When adding new state (e.g., per-cell collision counters for diagnostics), use per-patch storage on `Patch` rather than global static, to keep the parallelism clean.

## 9. Diagnostics and debug output

`Collisions(debug_every = N)` enables per-cell collision-rate output:

```cpp
// In Collisions::apply, when debug_every > 0:
if (itime % debug_every == 0) {
    write_collision_diagnostic_data(...);          // per-patch HDF5
}
```

Output goes to a dedicated HDF5 file per `Collisions` block. The schema includes: per-cell collision counts, per-cell mean Coulomb log, per-cell total momentum and energy transferred. happi has a corresponding accessor (`S.Collisions(...)` in some versions; this varies).

For ad-hoc debugging, add `printf` or `LOG()` calls temporarily — but remove before committing, as they break OpenMP performance.

## 10. Adding a new pair process

To add a new physics interaction between particle pairs (e.g., recombination, charge exchange, photoionization-by-collision):

1. **Decide** where it fits. Is it a modification of the existing `Collisions` block (new optional argument) or a new namelist construct (new `Collisions`-like block)?

2. **Implement** the cross-section function. Most processes have a velocity-dependent cross-section `σ(v)`. Implement it in a static or inline function.

3. **Add the namelist hook** in `src/Params/PicParams.cpp` — parse the new argument, validate.

4. **Modify `Collisions::apply`** to call your new process in the per-pair loop:
   ```cpp
   if (my_new_process_enabled) {
       my_new_process->applyToPair(patch, i1, i2, ...);
   }
   ```

5. **Test** — add a validation test exercising the new process. Where possible, benchmark against an analytical limit.

For a wholly separate algorithm (not a per-pair Pérez-style scheme), consider whether it should live outside `Collisions/`. For instance, fluid-style recombination would be in a `Recombination/` directory with its own per-patch loop.

## 11. Performance considerations

Binary collisions are an O(N²) algorithm within a cell (technically O(N) when pairing is random, but practically O(N × cost-per-pair)). Cost per cell scales linearly with particle count.

**For dense plasmas with high PPC** (e.g., 128 PPC in excimer-class runs), collisions can dominate the per-step cost — 30–50% of total time is not unusual.

**Optimization strategies:**

- **`every > 1` for slow channels** (where physically justified).
- **`LoadBalancing`** — without it, cells with the most particles stall MPI ranks waiting.
- **Avoid intra-species collisions** within multi-species groups when not needed. Specifying `species1=["Ar1","Ar2","Ar3"], species2=["Ar1","Ar2","Ar3"]` is 9 pair combinations; sometimes you only need same-state collisions.
- **GPU paths** for collisions are partial; the CPU path is currently the main target for performance work.

**Vectorization**: the inner per-pair loop is not vectorized — each pair has data dependencies on previous pairs (particle momenta are modified in place). Don't waste effort trying to SIMD-ize it.

## 12. Common pitfalls

**`reference_angular_frequency_SI` not set with `clog=0`**: silently produces nonsense Coulomb log. The user-facing rule 1. From the C++ side, the formula in `computeCoulombLog` returns infinity or NaN if `omega_r_SI` is zero — check for this guard.

**Modifying particle indices in the middle of the collision loop**: `eraseParticle` may invalidate indices in the per-cell list. If you add new physics that consumes particles, do the deletion at the end of the per-pair loop.

**Race on `electron_species->particles` during ionization**: each patch has its own electron species — no race across patches. But within a patch, if you add code that pushes new electrons in parallel, you need synchronization.

**Forgetting weight handling**: when different species have different weights, the collision probability must scale with `min/max(w1, w2)`. Code that ignores this gives wrong rates for mixed-PPC simulations.

**Cross-species collisions when one is a test particle**: should typically skip (test particles don't affect dynamics). The `is_test` check is in `Collisions::apply`.

**Modifying `Collisions::apply` without updating both CPU and GPU paths**: GPU paths are partial as noted. When the CPU path is updated, the GPU path may diverge. Document any divergence in the source.

**Adding state to `Collisions` without per-patch independence**: a `Collisions` object is shared across patches on a rank (one per namelist block). Don't store per-patch data on the `Collisions` itself; use the `Patch` it's called from.

**Assuming particle order is stable across calls**: `Collisions::apply` doesn't reorder, but other modules (sorting, migration) do. Don't cache particle indices across timesteps.

**Forgetting energy conservation in custom processes**: pair scattering must conserve total energy and momentum. New processes that don't (e.g., processes with energy thresholds or inelastic losses) must explicitly debit the lost energy to a tracked quantity, otherwise `DiagScalar.Ubal` drifts.
