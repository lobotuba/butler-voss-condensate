# Butler-Voss Condensate — Cheat Sheet

Reference for `simulation.py` (the Butler-Voss Condensate research engine,
formerly the "Grain Fabric Model").
This is a physics-INSPIRED toy model, not a claim about the real universe.
A **node** here is a physical locus (it carries `u`, `v`, and a local tension),
not a bare graph vertex; edges are real tension-bearing relations, not mere adjacency.

---

## 1. Governing equations

**Per-step update** (semi-implicit / symplectic Euler) for every node *i*:

```
a_i   = c²·tension_i·L(u)_i  +  F_pot(u_i)          # acceleration
v_i  += a_i·dt ;  v_i *= damping ;  u_i += v_i·dt   # integrate
# gravity drift (advection up the tension gradient):
w     = clip( gravity_strength · ∇tension , max_drift )
u_i  -= dt·(w · ∇u)_i ;   v_i -= dt·(w · ∇v)_i
```

- **Lattice operator** (symmetric graph Laplacian):
  `L(u)_i = [ Σ_j (u_j − u_i) ] / deg_norm`
- **Potential force** (oscillon):
  `F_pot(u) = −m²u + focus·u³ − sat·u⁵`
- **Energy density**:
  `e_i = ½v_i² + V(u_i) + ½·tension_i·[Σ_j (u_j−u_i)²]/deg_norm`
  (kinetic + potential + spring/gradient)
- **Tension update**: if `e_i > energy_threshold`: `tension_i += tightening_rate·e_i`;
  else relax toward `base_tension` at `relaxation_rate`; then clip to `[min_tension, max_tension]`.

Force and energy are derived from the same potential and the same symmetric
Laplacian, so at `damping = 1` (and frozen tension) total energy is conserved
up to integrator error — a built-in correctness check.

---

## 2. State arrays (computed, not set by you)

| Name | Shape | Meaning |
|---|---|---|
| `u` | (N,) | Fabric **displacement** at each node (the field) |
| `v` | (N,) | **Velocity** (rate of change of `u`) |
| `tension` | (N,) | Local **string tension** (the gravity field) |
| `pos` | (N, D) | Node coordinates; `D` = 2 or 3 |
| `neighbor_idx`, `valid` | (N, K) | Who each node's neighbors are, and a 0/1 mask |
| `deg_norm` | scalar | Max neighbor count; normalizes the Laplacian symmetrically |
| `dvec`, `Minv` | — | Geometry for the least-squares spatial gradient `∇` |
| `N`, `D` | scalars | Number of nodes; number of dimensions |
| `step_count`, `time` | scalars | Steps taken; simulated time = `step_count·dt` |

---

## 3. Config parameters (the "givens" — all tunable)

### Integration / medium
| Param | Default | Description | Effect of raising it |
|---|---|---|---|
| `dt` | 0.06 | Time step | Faster but less stable; too big → blow-up |
| `base_wave_speed` (c) | 1.0 | Ripple propagation speed | Ripples travel faster |
| `base_tension` | 1.0 | Resting string tension (empty fabric) | Stiffer baseline medium |
| `damping` | 0.999 | Energy kept per step | `1.0` = energy-conserving (clean tests); `<1` bleeds energy |

### Self-interaction potential V(u) — what makes particles
| Param | Default | Description |
|---|---|---|
| `potential` | "oscillon" | Which `V(u)` to use (see §6) |
| `mass2` (m²) | 1.0 | Linear restoring term → oscillation frequency / "mass" scale |
| `focus` (β) | 1.0 | `+u³` **self-trapping** term. **Set to 0 for the linear control run** |
| `saturation` (γ) | 0.5 | `−u⁵` term that prevents collapse / blow-up |
| `nonlinear_strength` | 0.06 | Only used by the old `potential="soft"` model |

### Gravity (tension mechanism)
| Param | Default | Description |
|---|---|---|
| `gravity_strength` (G) | 0.3 | How strongly excitations drift up the tension gradient. `0` = pre-gravity model |
| `max_drift` | 0.5 | Safety cap on drift speed (stability) |
| `tightening_rate` | 0.004 | How fast high-energy regions tighten strings |
| `relaxation_rate` | 0.001 | How fast tension relaxes back to baseline |
| `min_tension`, `max_tension` | 0.6, 4.0 | Clamp on tension |

### Particle detection & tracking
| Param | Default | Description |
|---|---|---|
| `energy_threshold` | 0.20 | Energy above which a node counts as "in a particle" |
| `min_blob_cells` | 3 | Ignore blobs smaller than this (kills flicker noise) |
| `match_radius` | 4.0 | Max centroid movement to keep the same particle ID across frames |
| `steps_per_measure` | 5 | Simulate this many steps between each measurement / log |
| `seed` | 0 | RNG seed → reproducible runs |
| `lattice`, `size` | "hex2d", 0 | Geometry and grid size (`0` = lattice default) |

---

## 4. Lattices (`--lattice`)

| Name | Dim | Neighbors | Notes |
|---|---|---|---|
| `hex2d` | 2D | 6 | Isotropic (default) |
| `square2d` | 2D | 4 | **Anisotropic** — waves prefer the axes (artifact test) |
| `cubic3d` | 3D | 6 | Simple cubic |
| `fcc3d` | 3D | 12 | Close-packed; true 3D analog of hex |

---

## 5. Experiments (`--experiment`)

| Name | Setup | Tests |
|---|---|---|
| `lump` | One centered lump | **H1** — does a single oscillon persist? |
| `collide` | Three overlapping ripples | H1 — does overlap birth particles? |
| `gravity` | Two lumps placed apart | **H2** — does separation shrink? |
| `travel` | One lump + sideways push | **H3** — pattern moves while nodes don't |

---

## 6. Potentials (`--potential`)

| Name | Force `F(u)` | `u=0` vacuum | Produces |
|---|---|---|---|
| `oscillon` | `−m²u + β·u³ − γ·u⁵` | stable | Long-lived breathing particles |
| `soft` | `−nonlinear_strength·u³` | stable | Dispersing ripples (old model) |
| `wave` | `0` | — | Massless linear waves (pure control) |

---

## 7. CLI flags

```
--lattice {hex2d,square2d,cubic3d,fcc3d}   --size N
--experiment {lump,collide,gravity,travel} --steps N   --seed N
--potential {oscillon,soft,wave}
--mass m²   --focus β   --saturation γ   --nonlinear k
--gravity G   --damping D   --dt T   --wave-speed c
--live | --headless   --out DIR
```

### Common recipes
```bash
# Watch one oscillon breathe (live animation)
python simulation.py --live --experiment lump

# H1 control: turn focusing off -> a massive linear wave that disperses
python simulation.py --headless --experiment lump --focus 0 --damping 1

# 3D oscillon (shown as a mid-plane slice)
python simulation.py --live --lattice cubic3d --experiment lump

# Lattice-independence test: same experiment on each lattice, then diff CSVs
python simulation.py --headless --lattice hex2d    --experiment collide --steps 1500
python simulation.py --headless --lattice square2d --experiment collide --steps 1500

# Energy-conservation correctness check (expect small drift)
python simulation.py --headless --experiment lump --damping 1 --gravity 0
```

---

## 8. Output columns

Written to `--out` dir (default `condensate_runs/<experiment>_<lattice>_seed<seed>/`).

**summary.csv** — one row per measurement
`step, time, total_energy, max_energy, mean_tension, max_tension,
 n_particles, total_particle_mass, max_lifetime, separation`

**particles.csv** — one row per tracked particle per measurement
`step, time, id, mass, peak, x, y, z, speed, lifetime`
- `mass` = summed energy in the blob, `peak` = max energy in it
- `speed` = centroid speed, `lifetime` = how long that ID has existed
- `z` is blank for 2D lattices

---

## 9. Hypotheses

| | Claim | Status | Key column |
|---|---|---|---|
| **H1** | Overlap/focusing makes persistent particles | passes | `max_lifetime` |
| **H2** | Tightened strings attract masses (gravity) | advection stable but disruptive | `separation` |
| **H3** | A particle is a pattern that moves through nodes | testable | particle `speed` |
| **H4** | Results are the same across lattices | testable | compare runs |
| **H5** | 2D vs 3D differ | oscillons survive both | compare runs |

---

## 10. Quick interpretation guide

- **Particle = blob** of nodes with energy above `energy_threshold`, tracked across frames.
- **Mass** = total energy bound in the blob (no separate mass variable exists).
- **A real effect should survive changing `--lattice`.** If a result flips between
  hex / square / cubic / fcc, suspect a lattice artifact, not physics.
- **2D cannot test realistic gravity** (force falls off as ~1/r in 2D, 1/r² only in 3D).
- **Oscillons are time-dependent**, which is why they persist in 3D even though
  *static* solitons cannot (Derrick's theorem).
