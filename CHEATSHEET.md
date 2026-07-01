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

---

## 11. Layer-below program — H6–H10 (complex-field extension)

These target the **sub-quantum layer**: the structures here are *not* quantum
particles but the material a layer (or more) below them — and the medium itself
(the "space between nodes") is taken as **active**, its texture carrying the
elementary quantum numbers. Tested in `prototype_complex.py` (complex order
parameter `psi`), not yet in `simulation.py`. Requires a symmetry-broken
(`mexicanhat`) vacuum `|psi| = v0 > 0` so phase windings are protected.

| | Claim | Key observable | Status |
|---|---|---|---|
| **H6** | **Charge = topological winding** of the field through the between-space (charge from the *shape* of space) | integer winding `n = (1/2π)·Σ plaquette phase` — quantized & conserved | prototype: `n=0,±1,±2` conserved exactly |
| **H7** | **Spin / charge from symmetry** — internal field-space rotation gives a continuous conserved charge (a Q-ball; precursor of spin/angular momentum) | Noether charge `Q = Im Σ conj(psi)·psi_dot` | prototype: Q conserved (Q-ball) |
| **H8** | **The particle zoo is finite (≈3–5)** — distinct fundamental species = distinct stable sectors `(n, Q)` | census count of long-lived `(charge, energy)` types | prototype shows `n=0, ±1, ±2` + Q-ball family |
| **H9** | **Sub-quantum hierarchy** — these solitons sit below quantum particles; higher layers form by **confluence/binding** | binding energy `E_bound < ΣE_free`; **scale gap** (size/energy ratio) between base solitons and bound states | prototype: `+/-` bind (ΔE≈4) & annihilate, `n=2` splits→composite, gap ≈6.6× size |
| **H10** | **Self-cohering medium** — nodes attract (with short-range repulsion); this *sets* the inter-node spacing, **self-selects close-packed/isotropic** geometry, and stabilizes structures | equilibrium spacing; emergent lattice = fcc/hcp; longer lifetimes | prototype: cloud self-assembles to hex (psi6 0.10→1.00, coord→6) at spacing≈r0; hex energy < square |

**Why H6/H7 need a complex field:** a real scalar `u` has no orientation to
rotate (no spin) and no symmetry to conserve (no charge). Promoting `u → psi`
(complex / director) gives both — exactly the **condensate order parameter** the
project is named for. This is the Skyrme-model lineage: a "particle" one layer
down is a *texture in a medium*, not a point.

**Why H10 matters:** it answers "how much space exists between nodes?" — the
spacing is no longer a free knob but the equilibrium of attraction vs. short-range
repulsion. A self-assembling medium relaxes to close packing (the isotropic
`fcc3d` we already found stable), so the `cubic3d` anisotropy artifact would
**solve itself**.

### Prototype CLI
```bash
python prototype_complex.py     # H8 census: conserved charges + persistence
python h9_binding.py            # H9: binding curve, annihilation, n=2 splitting, scale gap
python prototype_mobile_nodes.py # H10: self-assembly, spacing, isotropy selection
```
`prototype_complex.py` columns: `topo0/topoF` (winding, start/end) ·
`Noether0/F` (U(1) charge) · `E0/EF` (energy) · `loc%F` (excess energy still localized).

### H9 first results (`h9_binding.py`)
- **Opposite charges bind, like charges repel.** Static interaction energy *rises*
  with separation for `+/-` (attractive, binding energy ≈ 4) and *falls* for `+/+`
  (repulsive). Released at rest, a `+/-` pair spirals in and **annihilates**
  (core depth `min|psi|/v0`: 0.21 → 0.95).
- **`n=2` is a composite, not a fundamental.** Two like vortices fly apart
  (charge spread 4 → 8), so higher windings are *not* new species — refining the
  H8 zoo to `n=0, ±1` plus the Noether (Q-ball) family.
- **Scale gap.** A bound `+/-` "molecule" is ≈6.6× the size and ≈1.6× the energy
  of a single base vortex — a first quantitative handle on the layer separation.
- **Caveat:** single-charge (`+/+`, `n=2`) configs have net winding, so their
  absolute energy is box-dependent; the quantitative binding uses the neutral
  `+/-` pair, whose far field cancels.

### H10 first results (`prototype_mobile_nodes.py`)
Mobile nodes interacting by a Lennard-Jones pair potential (attract at mid-range,
repel at short range), velocity-Verlet with cooling, on a free 2D cluster.
- **Spacing is set, not chosen.** A disordered cloud settles to mean spacing
  ≈ `r0 = 2^(1/6)σ` (the potential minimum) — answering "how much space exists
  between nodes" dynamically.
- **Isotropy self-selects.** Bond order `psi6` rises 0.10 → 1.00 and coordination
  → 6: the medium freezes into a **hexagonal** lattice. Static check: energy/node
  is lower for triangular/hex (−6.64) than square (−5.27), so the anisotropic
  arrangement is energetically rejected. **In 3D the same logic selects fcc/hcp —
  i.e. the isotropic `fcc3d` we found stable — so a self-organizing medium would
  never adopt the `cubic3d` arrangement that blew up.** The artifact dissolves at
  the level of the medium itself.
- **Cohesion = stability.** Mutual attraction binds the cloud into a stable
  droplet (kinetic energy cools to ~0) instead of dispersing — confirming the
  long-run-stability motivation.
- **Caveat:** finite droplet, so edge nodes lower the *average* coordination
  below 6 (interior is 6); a periodic box would give a clean 6 throughout.

### Integration Phases 0–1 (`integration_field_medium.py`)
First de-risking of "complex field **on** the self-assembled medium" (field
attached to nodes, Lagrangian).
- **Phase 0 — the operator is the crux.** A plain weighted-graph Laplacian is
  fine on a perfect lattice (≈4% error) but **degrades ~16× on the irregular
  self-assembled cloud (≈60%)**; a **least-squares (LSQ) meshfree Laplacian stays
  robust (~13%)**. So accurate field dynamics on a self-organised medium *require*
  the LSQ operator — the cheap graph operator re-introduces the artifact class.
- **Phase 1 — charge survives, shape degrades.** A seeded `n=+1` vortex keeps its
  winding **conserved start→end on the self-assembled medium** (it is
  topologically protected — the integer can't change unless the zero exits or
  pair-annihilates). What the irregular mesh degrades is the vortex *core* (it
  wanders and fills in / is under-resolved); LSQ preserves it better than graph.
- **Takeaway:** feasibility is **GO for the conserved quantum number**; the open
  quality issue is core resolution (denser medium near cores, or relax the vortex
  to equilibrium on the actual mesh before testing). Next: Phase 2 (let the medium
  move) then Phase 3 (two-way coupling = gravity-by-density).

### Integration Phase 2 (`integration_phase2.py`)
One-way coupling: the medium MOVES (Lennard-Jones MD + thermostat), the field is
Lagrangian (`psi` attached to nodes), and the LSQ operator is rebuilt from the
current positions every few steps. Question: does winding survive when nodes
**reconnect** (swap neighbors)?
- **Solid but rearranging medium → charge survives.** At moderate temperature the
  medium underwent **~308 neighbor reconnections** (15× the frozen control) while
  staying cohesive, and the `n=+1` winding stayed **conserved throughout**. The
  topological charge is robust to the lattice moving and reconnecting under it.
- **Melting is the limit.** Near melting (~1587 reconnections) the winding is
  lost — destroying the medium's order destroys the hosted particle. So there is
  effectively a **melting threshold above which particles cease to exist.**
- **Verdict: Phase 2 = GO** (charge robust to medium motion while the medium
  stays solid). Remaining: Phase 3 — two-way coupling (field ↔ node density) for
  **gravity-by-density**.

### Integration Phase 3a/3b (`integration_phase3ab.py`)
The two halves of the gravity-by-density loop, validated in isolation with an
**imposed** field (no feedback) — cheap, stable, and they fix the sign.
- **3a — Coupling B (medium → field) works, sign confirmed.** A massless test
  wave-packet flown past an imposed dense region (wave speed `c²·g(ρ)`,
  `g=1/(1+βρ)`) **refracts toward the mass when denser = slower** (impact
  parameter 5.0 → −1.75, i.e. bends across the mass line; control stays +5.0), and
  **away** when the sign is flipped (→ +8.0). So the lensing mechanism and the
  attractive sign (**denser ⇒ slower**) are pinned.
- **3b — Coupling A (field → medium) works.** Pushing mobile nodes up the gradient
  of an imposed field-energy blob (`F += α·∇e`) **compresses the medium**: central
  node count rises +22–24% (spacing 1.105 → 1.05) as `α` increases; `α=0` is flat.
- **Takeaway:** both halves are sound with the correct sign → the closed loop
  (3c: self-focusing of one lump; 3d: two lumps drifting together = gravity) has a
  validated foundation. The open risk is the feedback **stability** when the loop
  is closed.

### Integration Phase 3c (`integration_phase3c.py`) — partial
Both couplings ON at once, with a real evolving field, on one lump. Honest result:
- **Stable run achieved (milestone 1).** The naive loop blew up — `∇e` of the
  rough mesh field ejected nodes. **Smoothing the field energy before
  differentiating it** fixed that: the coupled run stays bounded, medium intact.
- **Directional effect, correct sign.** The lump's energy compresses the medium
  *persistently* (central ρ 1.06 → 1.10 vs flat ~1.05 control) and the slowed
  waves make it **spread ~60% less** (Δwidth 0.93 vs 2.26 over the run).
- **But NOT self-focusing.** Width still grows (4.9 vs 6.23) — the density well is
  too shallow (ρ~1.1 ⇒ g≈0.94) to halt dispersion. Pushing harder is capped by a
  **field-operator instability**: the LSQ Laplacian has spurious positive
  eigenvalues that blow up weakly-restored fields.
- **Verdict: 3c partial** — stable, correctly signed, sub-critical. A true
  self-trapped state (and 3d) needs a *stabilized/symmetric* field operator and a
  *deeper/more-sensitive* density→speed coupling — ideally both couplings derived
  from one Lagrangian so energy is conserved by construction.

### Symmetric field-operator scorecard (`operator_scorecard.py`)
Quantifies field operators on: symmetry `‖A−Aᵀ‖`, stability `max Re(eig)` (≤0
required — it *is* the blow-up growth rate), accuracy RMS vs an analytic
Laplacian, and an end-to-end energy-drift run.
- **The LSQ operator is the 3c blocker, measured:** asymmetric (0.81), indefinite
  (`max Re(eig)=+0.31`) → in a `damping=1` run it blows up (|u|→1011, energy drift
  −13000%).
- **Stability is solvable.** A symmetric SPH **Brookshaw** Laplacian and a
  spectrally-clipped LSQ both have `max Re(eig)=0`, a real spectrum, and run
  **bounded with ~0% energy drift** — the field no longer explodes.
- **Accuracy is the cost.** Symmetric operators are ~2× less accurate than LSQ
  (Brookshaw 43%, clipped-LSQ 66% vs LSQ 24% shape-RMS): removing the spurious
  modes sacrifices the first-derivative cancellation LSQ gets from asymmetry.
- **Path:** Brookshaw is the better symmetric candidate; being stable, it lets the
  gravity coupling be pushed far harder than LSQ allowed (good enough to retry 3c
  self-trapping, at a smaller CFL `dt`). Matching LSQ accuracy *and* staying
  symmetric needs a consistency-corrected operator (renormalized SPH / constrained
  symmetric MLS).

### 3c revisited with the Brookshaw operator — self-trapping signature
Swapping the field Laplacian to **Brookshaw** (stable) and keeping the LSQ
*gradient* for the node force unlocked the regime LSQ couldn't reach: deep wells
(β up to 60) run **without blow-up**, and the field can run **near-conservatively**
(`field_damping≈1`).
- **The lump now HOLDS its width** instead of dispersing: width 3.97 → **4.13**
  (β=60, fd≈1) vs the uncoupled control's 3.97 → **7.96** — ~96% less spreading,
  with a **persistent density well** under it (ρ 1.06 → 1.10, held). A clear
  gravity-by-density **self-focusing signature**.
- **Not yet a clean stationary soliton:** the peak amplitude still decays faster
  than damping → residual energy leakage (radiation; the ad-hoc couplings aren't
  exactly energy-conserving on a moving mesh).
- **Net:** the operator fix moved 3c from "disperses ~60% less, sub-critical" to
  "stops dispersing (self-traps) with a stable persistent well." Remaining for a
  true bound state: a variational (single-Lagrangian) coupling + radiation
  handling — then **3d** (two lumps → drift together) becomes meaningful.

### Variational coupling (`integration_phase3_variational.py`)
Both gravity-by-density couplings derived from ONE energy functional
`U = ¼ Σ_ij (γ_i+γ_j) W(r_ij)(u_i−u_j)²`, `γ=g(ρ)` — so they are exact gradients
of the same E and energy is conserved by construction (the 3c leakage fix).
- **Construction verified.** The analytic node force matches `−∇E` to **1e-9**
  (finite-difference), and at β=0 energy drifts **−0.02%** (machine precision).
  The field force is automatically a *symmetric* Laplacian (stable), with
  density-dependent weights = Coupling B; the node force = Coupling A, and the
  radial sign is **inward (compression) = gravitational**, not anti-gravity.
- **Singularity found & fixed.** `g=1/(1+β(ρ/ρ0−1))` diverges for rarefied nodes
  (denominator→0) → +10⁹% blow-up; replaced with bounded `g=exp(−β(ρ/ρ0−1))`.
- **Residual obstacle = stiffness.** With coupling on (β=10) the *explicit*
  integrator shows ~7% energy drift at `dt=0.003` (grows; shrinks with `dt`), and
  the bounded coupling is too gentle to self-focus.
- **State of the program — two partial results that don't yet combine:**
  Brookshaw-3c *self-traps but leaks energy*; the variational form *conserves
  energy but is too gentle/stiff to self-trap*. Getting **both** (energy-conserving
  self-trapping) needs a **structure-preserving / implicit integrator** so β can be
  pushed into the trapping regime while keeping conservation. Then **3d** (two
  lumps → mutual drift = gravity) becomes the decisive test.

**`dt`-check + β-sweep diagnostics (deciding the integrator):**
- **The residual drift is not `dt` error.** Shrinking `dt` 6× barely moved it
  (3.13%→2.58%, plateauing ~2.5%); both forces are FD-verified gradients of
  `energy()` to ~1e-9. So the system *is* Hamiltonian — the residual is **symplectic
  shadow-energy error of a very stiff system**, which an *energy-conserving*
  (discrete-gradient / AVF) integrator removes but an implicit *symplectic* one
  would not. (The LJ force-cap was ruled out: clipped == unclipped.)
- **The coupling is more than strong enough.** β-sweep (damping=1): β=10 disperses
  (sub-critical), but **β≥20 strongly compresses the lump** (width ~4 → minW
  1.1–1.3) then **blows up**. So a *critical* β giving a stable bound state lies
  between dispersal and collapse — the mechanism works; the regime just needs
  controlling.
- **Bounded `g(ρ)` (done) — logistic `g∈[g_min,g_max]`** caps the field-operator
  stiffness so β only sets the *sharpness* of the density→speed response, not its
  magnitude. This **cured the blow-up AND restored energy conservation** (dE ≈
  0.00% even at β=40, dt=0.002) — confirming the dt-independent ~2.5% drift was the
  rarefied-region stiffness, so the **AVF integrator is no longer needed for
  conservation**. (FD-verified the bounded `g` keeps forces = exact gradients.)
- **But the bounded coupling is too gentle to trap from a dispersing seed:**
  `minW` ≈ initial for all β ≤ 80 (even deep wells `g_min=0.02` and slow fields
  `m2=6`). Trapping happened *only* in the violent unbounded-collapse regime
  (which blew up). So a stable bound state — if it exists — must be reached by
  (a) **seeding near it** (pre-compressed medium + concentrated field, then test
  *stability*) or (b) **controlled collapse** (let it collapse; the LJ floor halts
  it at finite width). The seed-near-a-bound-state *existence test* is the cheaper
  next probe.
- **Existence test result: NO self-bound soliton (and why).** Pre-compressing the
  medium under a field lump (relax nodes with the field frozen) then releasing
  conservatively: the lump still **disperses**. The cause is decisive — the LJ
  medium is **nearly incompressible**: the central density well **saturates at
  ~1.09** (≈9%) *regardless of source strength* (amp 1→3.5, β 30→60 all give the
  same shallow well), because the steep LJ repulsive core balances the field's pull
  at modest compression. So **gravity-by-density is real but weak** — too weak for a
  single lump to self-bind against the rigid medium. (Physically apt: gravity is the
  weakest force; self-binding needs enormous mass.)
- **Redirect to 3d.** Self-binding needs *strong* gravity; mutual attraction needs
  only the *weak* force already confirmed (3a refraction). So **3d — do two lumps
  drift together? — is the right decisive gravity test**, and it does NOT require a
  bound state. With energy now conserved (bounded `g`), a clean 3d run is feasible:
  seed two lumps, measure whether their separation shrinks.

### 3d result — GRAVITY (`integration_phase3d.py`)
Two field lumps seeded at ±4.5 on the medium, energy-conserving variational
coupling, separation tracked to t=8:
| β | separation | net | dE |
|---|---|---|---|
| **0** (control) | 9.33 → 10.85 | **+1.51** (drift apart) | −0.12% |
| **40** | 9.33 → 7.78 | **−1.56** (drift together) | −0.17% |
| **60** | 9.33 → 7.48 | **−1.85** (stronger) | −0.08% |

**Two masses attract via the emergent medium force** — the gravitational
signature is unambiguous: *no* attraction without the coupling (β=0 drifts apart),
attraction *with* it, **stronger at higher β**, and **energy conserved** (~0.1%,
so it's real dynamics, not drift). This is the decisive positive result of the
whole program: **gravity emerges from the medium**, as a genuine (weak) attractive
force between mass-energy concentrations — not the leaky advection hack of H2.

### Force law F(d) — not cleanly resolvable; tentatively SHORT-RANGE
Three probes tried to extract how the attraction scales with separation:
(1) two dispersing lumps (separation) — noisy, non-monotonic; (2) core-tracking
over long time — *increased* with d (merging/dispersion artifact, wrong sign);
(3) cleanest — a test lump released at distance d on a **frozen dug well**.
- All three are confounded by weak gravity + field dispersion, and (3) also by the
  well itself: a steep heavy source *rarefies* its center (ρ/ρ0≈0.84) rather than
  compressing — the medium's response is **sign-sensitive to the source profile**.
- The test-lump drift is strong at small d (~0.95 at d=4) and falls to the **noise
  floor by d≈6**; the density perturbation is **local** (within r≈5). A `d^-2` fit
  appears but is an artifact of fitting a fast drop + noise — **not trustworthy**.
- **Defensible conclusion:** the emergent gravity is **short-range / screened**
  (dies within ~5–6 lattice units), *not* a long-range power law — consistent with a
  nearly-incompressible medium localizing the density perturbation (no long-range
  strain tail). A definitive F(d) needs a dedicated build (clean controllable well +
  proper force probe); no exponent is claimed.

### Self-binding on a compressible medium — NOT achieved (artifact caught)
Since the LJ medium is nearly incompressible, we tried a softer (Morse) medium to
see if a single lump could self-bind. **A cautionary result:** with a naive
Morse cutoff (force truncated, energy not shifted) the lump appeared to **self-trap**
(width 3.9→2.8) — but that was a **numerical artifact**: the inconsistent cutoff
leaked **~80% of the energy**, letting the system cool into a compressed state.
Fixing the cutoff (shifted-force Morse) **removed the trapping** — the lump
disperses at every β (β>0 disperses *less* than β=0, consistent with weak gravity,
but no bound state). Takeaway: **no self-bound soliton even on a compressible
medium** in the energy-conserving regime; the gravity force is too weak to overcome
dispersion. (Caveat: the soft-Morse+coupling energy accounting is itself imperfect
here (±15–34%), so a *definitive* compressible-medium test needs a properly
conservative soft potential — a real build. The clean, conserved result remains 3d:
mutual attraction on the LJ medium.)
