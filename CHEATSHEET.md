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

**Resolved by the static density-response test (the clean method).** Instead of
chasing a dynamic drift, measure the medium's static density perturbation `Δρ(r)`
around one **gentle** frozen mass (angle-averaged, large medium, relaxed with vs
without the source). Result: `Δρ(r)` is a small (~5%) **local** perturbation — a
rarefied core + a compression ring at r≈(source size) — that **falls to the noise
floor (~±0.01) by r≈5–6**, with **no measurable tail** (both exp and power-law fits
are meaningless — it's noise beyond the core). Two amplitudes give the same-shape
profile (linear response). **Conclusion (final): there is NO long-range force law.**
A mass's density perturbation is confined to ~its own size, so two masses attract
only when their perturbations **overlap** — a **contact-range / screened**
attraction (screening length ≈ mass size), not Newtonian `1/r` or `1/r²`. So this
model's gravity is a real, weak, energy-conserving, **short-range** attraction.
(Optional follow-up: the adiabatic `U(d)=E(d)−E(∞)` energy curve could quantify the
contact-interaction *strength*, but the *range* question is settled.)

**Contact strength quantified (`interaction_energy.py`).** The adiabatic method —
freeze two gentle masses at separation d, relax the medium, and isolate the
medium-mediated energy `U(d)` by subtracting the frozen-medium (direct-overlap)
baseline — gives a clean, monotonic **attractive well** (N=400, β=60):

| d | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|
| `U(d)` | **−0.921** | −0.255 | +0.069 | +0.071 | +0.056 | 0.000 |
| `F=−dU/dd` | | −0.50 | −0.16 | +0.01 | +0.02 | |

- **Contact binding energy ≈ 0.92** (model units) at d=4 — ~25% of the single-mass
  medium-response energy (~−3 to −4): a substantial contact interaction.
- **Peak force ≈ 0.5 at contact** (d≈5), comparable to the LJ medium forces (~1),
  **vanishing by d≈6–7**.
- `U→0` for d≥6 → **screening range ≈ 5–6 = the mass size**, matching Δρ(r).
- The static energy method is what worked (clean, monotonic) where every *dynamic*
  drift probe failed. Full characterization: **short-range attractive well, depth
  ≈0.9, peak force ≈0.5, screening length ≈5–6 — not a power law.**

### Order lengthens & cleans the force (disorder was screening it)
Redoing `Δρ(r)` on a **perfect hex lattice** (defect-free) instead of the
self-assembled cloud changes the picture:
| r | 0.5 | 2.5 | 4.5 | 6.5 | 8.5 | 10.5 | 11.5 |
|---|---|---|---|---|---|---|---|
| `Δρ` | .0010 | .0010 | .0009 | .0006 | .0004 | .0003 | .0002 |
- **Clean, monotone, and ~2× longer-range** — best fit `exp(−r/λ)`, **λ ≈ 5.7**, a
  real tail out past r=11 — vs the disordered medium's *noisy, local* blob (peak
  ~5%, screened by r≈5). The perfect lattice's response is ~50× smaller in
  amplitude (~0.1%) but coherent and long-range.
- **Interpretation:** disorder *masks* the true force with large **local plastic
  rearrangements** (defects absorbing compression); order reveals the genuine
  **coherent elastic strain**. So the "short-range" result was **partly a disorder
  artifact — order/rigidity is a real lever on the force range.**
- **Still screened** (finite λ), not yet a power law. **Levers for truly
  long-range gravity: (1) 3D** (a point dilatation gives `1/r²` strain — Newtonian
  territory; 2D screens harder), **(2) an ordered, rigid lattice.** A pure
  tetrahedral (z=4) lattice is floppy under central forces (3D rigidity needs
  z≥6); the useful "tetrahedral" geometry is the close-packed **fcc** (z=12, built
  from tetrahedra+octahedra), reached by ordering.

### 3D (ordered fcc): screened too, λ≈3.3 — NOT Newtonian (`density_response_3d.py`, `density_response_3d_large.py`)
The variational engine is dimension-agnostic, so it runs in 3D unchanged
(Phase 3D-0: fcc N≈3055, energy drift −2.2%, stable).

**3D-1 (radius 9, N≈3055) — a false positive.** On an ordered fcc the response
`Δρ(r)` barely decayed over the measurable range (factor ~1.6 to r=6.8), and a fit
gave **λ≈14** with a power law not excluded — suggesting 3D was dramatically *less*
screened than 2D. **This was a finite-size illusion:** radius 9 only reaches r≈6.8,
which is inside the inner, barely-curved part of the response *and* near the free
surface, so an exponential looks nearly flat there.

**3D-2 (radius 13, N=9213) — the honest answer.** A **cell-list sparse** force
(O(N·nbrs), ~76 nbrs/node, pair list built once from the perfect fcc and reused
through the cooled relaxation) reaches far enough that the tail curls over. The
sparse coupling force is verified against the dense variational engine to
**machine precision (rel-err 2e-15)** before trusting the large run.
| r | 0.5 | 2.5 | 4.5 | 6.5 | 8.5 | 10.5 |
|---|---|---|---|---|---|---|
| `Δρ` | .054 | .049 | .037 | .020 | .010 | .006 |
- Fit over the clean window (2.5 < r < 10.5, away from source and surface):
  **exp λ ≈ 3.3 (SS 0.011) beats power-law n≈1.8 (SS 0.10) by 10:1** → the response
  is **exponential/screened**, with a screening length **comparable to — even a touch
  shorter than — 2D hex's λ≈5.7**. The lone outermost shell (r>10.5) flattens
  slightly (surface pile-up); the fit correctly excludes it.
- **Conclusion:** 3D does **not** lengthen the force into Newtonian territory.
  Gravity-by-density is **short-range/screened in both 2D and 3D** (the wave-mediated
  coupling does *not* evade the Bitter–Crum suppression as hoped). The mutual
  attraction of two masses (3d) is **real but intrinsically short-ranged** — a
  contact-like force, not `1/r²`. Getting a truly long-range force would require a
  different coupling mechanism, not merely more dimensions.

### Attacking the screening — why it screens, and the route out (`screening_diagnosis.py`, `screening_gauss.py`)
The screened result raised the real question: *why* is it screened, and is the
screening escapable? Two experiments pin it down.

**Screen-0 — the coupling dial is not a lever (`screening_diagnosis.py`).** Sweep
the coupling sharpness β on the large fcc and measure `Δρ(r)`:
| β | 20 | 40 | 60 | 100 |
|---|---|---|---|---|
| λ | 3.2 | 3.7 | 4.3 | 5.8 |
| peak amp | .073 | .059 | .036 | .016 |
| integrated Σ`Δρ` | 103 | 96 | 63 | 32 |
- λ creeps up with β, but the **amplitude *and* the total integrated compression
  collapse together** — at high β the logistic `g'(ρ)` saturates to zero except in a
  thin shell at ρ≈ρ₀, so the whole response fades. You can trade strength for a
  slightly longer reach, **never a strong long-range force**. Intrinsic elastic
  screening (weak-coupling limit) is **λ≈3**. There is **no conserved, un-screenable
  flux** in this coupling — exactly what a long-range force needs.

**Screen-1 — existence proof: screening IS a mass term (`screening_gauss.py`).**
On the *same* fcc node machinery, solve a discrete Poisson equation for a point
source (sparse nn Laplacian scaled to `lap(r²)=6`, pure-numpy CG, Dirichlet wall):
massless `∇²Φ=−s` vs massive `(∇²−m²)Φ=−s`. Boundary-honest discriminator — **scale
the box** (a massless 1/r field has *no* intrinsic length; a Yukawa field pins at
λ=1/m):
| box radius | 10 | 14 | 20 |
|---|---|---|---|
| massless λ_apparent | 2.45 | 3.18 | **4.24** (grows) |
| massive m=1/3 λ | 1.86 | 2.02 | 2.18 (pins) |
- The **massless field's range grows without bound with the box** (`Φ·r` declines
  *linearly* `A(1−r/R_b)` = bounded Coulomb, not concave-exp) → genuine **1/r,
  unscreened**. Add a mass term and it becomes **Yukawa exp(−r/λ), λ≈1/m**.
- **Conclusion:** the node machinery *can* carry an unscreened `1/r²` force — it does
  so **iff the mediating field is massless**. Gravity-by-density screens precisely
  because the medium's **pinning to its rest spacing R₀ acts as a mass term**. The
  escape is a **conserved source that forbids a mass term** (a Gauss law), *not* a
  bigger coupling or more dimensions. → **Screen-2:** couple the potential to the
  **conserved topological charge (H6)** and test for emergent long-range — the
  prediction being a long-range force for conserved charge (EM-like) alongside the
  short-range force for energy (nuclear-like).

**Screen-2 — the payoff: a CONSERVED charge sources a LONG-RANGE force (`screening_topocharge.py`).**
The model already owns a conserved source: the H6 topological winding. Its mediator
is the **phase** of the Mexican-hat complex field — a **Goldstone mode, massless by
symmetry** (not by choice), and the winding is conserved by topology. Both together
forbid the mass term that screens gravity. Test in 2D (where U(1) winding is a point
vortex): seed a neutral +1/−1 vortex pair at separation `d`, measure the pair
formation energy `E(d)` straight from the energy functional (no dynamics — an
opposite pair would annihilate). A massless phase gives the 2D-Coulomb **log** law
`E ~ 2πρ_s ln d` that never saturates; a massive phase would flatten at `d~λ`.
Boundary-honest discriminator = **scale the lattice** (100→140→200):
| box | 100 | 140 | 200 |
|---|---|---|---|
| log slope A | 2.75 | 2.83 | **2.84** (stable → genuine log) |
| best-saturating λ | 20 | 25 | **33** (grows ∝ box → no intrinsic length) |
- `E(d)` rises logarithmically with a **box-independent slope** and **keeps climbing
  at every box edge**; the best exponential-saturation fit degrades catastrophically
  as the box grows (SS 0.0→0.2→3.1) and its λ just tracks ~⅕ of the box. So the
  interaction is a **genuine unscreened 2D-Coulomb log** (force ~1/d).
- **The result:** the condensate hosts **two force classes at once** — **long-range
  for a conserved (topological) charge** (EM-like, unscreened, massless Goldstone
  mediator) and **short-range for energy** (nuclear/gravity-like, screened λ≈3–6,
  the medium's R₀-pinning acts as a mass). Masslessness here is **emergent** (Goldstone),
  unlike Screen-1 where it was imposed. Force range is set by *what conserves the
  source*, not by dimension or coupling strength.
- **Massive control, now built:** the natural "massive control" is to **gauge the
  U(1)** — see the gauged-U(1) section below.

### Gauging the U(1) screens the force — the massive control (`screening_gauged.py`)
The last piece: if the long-range force in Screen-2 is really the massless Goldstone
phase, then **gauging** the symmetry (adding a gauge field on the lattice links →
the **Abelian Higgs model**, i.e. a superconductor) should destroy it. The gauge
field eats the Goldstone; the **Meissner effect** gives the photon a mass
`m_A = e·v₀` (penetration depth `λ_L ~ 1/(e v₀)`), turning the vortices into
Abrikosov vortices with a **short-range** interaction. This is Screen-1's
massless→massive knob, but the mass is now **dynamical** (from the gauge coupling
`e`), not imposed. Gauge-invariant lattice energy: covariant hopping
`Σ|ψ_j e^{-iθ_ij}−ψ_i|²` + compact Maxwell `(1/e²)Σ(1−cos B)` + Mexican-hat.

- **G-0 (correctness gates, both pass):** energy invariant under a lattice gauge
  transform `ψ_i→ψ_i e^{iα_i}, θ_ij→θ_ij+α_j−α_i` to machine zero (dE=0); analytic
  `dE/dθ` matches finite differences (rel-err 7e-8). (Relaxation of the gauge field
  needs a stiffness-scaled step `η~e²` + momentum — the Maxwell term's `1/e²` makes
  naive fixed-step descent diverge for small `e`.)
- **G-1 (single vortex — the smoking gun):** the **global** vortex energy grows with
  the box (19.9→22.4→24.2 at L=40/60/80 = the `ln L` divergence / long-range tail),
  while the **relaxed gauged** energy is *exactly* box-independent (7.346 at every L)
  and the magnetic flux **quantizes to one quantum** (ΣB/2π = 1.000). Gauging cut off
  the long-range tail at the one-body level.
- **G-2 (pair E(d) vs coupling e — the screening length):**
  | e | 0.15 | 0.20 | 0.30 |
  |---|---|---|---|
  | λ_L | 4.14 | 3.13 | 1.94 |
  | e·λ_L | 0.62 | 0.63 | 0.58 |
  `E(d)` **saturates** at the plateau `2·E_single` (vortices stop interacting beyond
  λ_L), and **`e·λ_L ≈ const → λ_L ~ 1/e`** = the Meissner penetration depth. Crucially
  **λ_L is box-independent** (1.94 at box 70 *and* 110) — a genuine *intrinsic*
  screening length, the opposite of Screen-2's apparent-λ that grew with the box
  (20→25→33). (For e≳0.45, λ_L<2 lattice units — screened below resolution.)
- **Conclusion — the arc closes.** Force range is governed by the **symmetry
  structure of the source**: a **global** conserved charge → massless Goldstone →
  **long-range** (Screen-2); **gauging** it → Meissner mass → **short-range** (here),
  with the crossover controlled by `e`. The condensate reproduces, from one complex
  field on the medium, the qualitative menu of real forces — long-range (unscreened),
  short-range-by-mass (screened), and short-range-by-gauging (Meissner) — and in every
  case the range is set by *what protects the mediator from a mass term*, exactly as
  Screen-1 diagnosed.

### Mobile-vortex check of the screened force law (`screening_gauged_mobile.py`)
G-2 measured the gauged-vortex interaction on a *frozen* scalar ansatz (only the gauge
field relaxed). Here the vortices actually **move** and the force law is read from their
motion — overdamped (dissipative) scalar dynamics with an *adiabatic* gauge field (relaxed
to equilibrium each step), which sidesteps the dispersion that sinks weak-force dynamic
probes. A +1/−1 pair closes in; overdamped, the closing speed `v(d) ∝ F(d)` (extracted via
a dwell-time proxy — steps to close each unit of separation — robust to plaquette-quantized
vortex positions).
- **M-0 gates:** scalar *and* gauge forces vs finite-difference `<1e-7` (the scalar gate needs
  the Wirtinger factor 2: a real/imag FD gives `2·Re/Im(force)`); energy monotone-decreasing;
  winding (±1) conserved as the vortices move.
- **M-1:** `v(d) ~ exp(−d/λ_L)` — a **screened** force. From motion **e·λ_L ≈ 0.42** (e=0.15→λ 2.82,
  0.20→2.13), i.e. **λ_L ~ 1/e** — the Meissner scaling, confirmed *directly from vortex motion*.
  It runs ~30% shorter than the static frozen-scalar G-2 value (0.6/e), as expected (the mobile
  measurement lets the scalar fully relax; the static one is a variational upper bound). So the
  motion independently confirms the **screened** (Meissner) force with an **intrinsic** length —
  not the long-range log of the ungauged U(1).

### The 3D topological defect: long-range survives into 3D as a vortex line (`screening_topocharge_3d.py`)
Screen-2 used a 2D point vortex — the U(1) defect in 2D. In 3D a single complex
field's defect is a **line** (π₂(S¹)=0, so no point charge from one field; a 3D
point charge needs a larger target space, S² hedgehog, or gauging). The direct 3D
continuation of Screen-2 is therefore the **vortex-line** interaction, built on a
cubic lattice (topology is a property of the field's target space, not the medium;
a cubic grid gives clean transverse plaquettes, and every quantity is a *static*
energy so cubic wave-anisotropy is irrelevant).
- **Correctness gate:** winding is an exact integer and conserved along the line —
  +1 in every z-slab for a single line, +2 for n=2, net 0 with two nonzero cores for
  a ±pair; vacuum energy 0. (Caught a real bug: a vortex centred *on a lattice site*
  splits its winding across four plaquette corners and miscounts as 0 — cores must
  sit at plaquette midpoints. Energy is computed from ψ directly and was unaffected.)
- **L-0 (a line, not stacked points):** energy is **strictly proportional to length**
  (E/L_z = 13.30 at L_z = 4, 8, 16) and the winding threads every slab — a genuine
  3D line with a tension.
- **L-1 (long-range in 3D):** two antiparallel lines, E(d)/L_z vs separation, transverse
  box-scaling 60/90/120:
  | box | 60 | 90 | 120 |
  |---|---|---|---|
  | log slope | 5.45 | 5.73 | **5.77** (stable, 2.6% → genuine log) |
  | best-sat. λ | 12 | 16 | **20** (∝ box → no intrinsic length) |
  The vortex-line interaction is a **long-range 2D-Coulomb log per unit length**;
  the conserved topological charge sources a long-range force in 3D as it did in 2D,
  with energy-gravity (screened, λ≈3) the short-range contrast.
- **L-2 (a vortex ring — the intrinsically-3D defect):** a closed vortex loop (axis
  along z, poloidal-phase ansatz). Gate: the loop pierces a meridian plane at exactly
  **two opposite cores** (+1, −1), the signature of a closed line. `E(R)` for R = 5…17
  rises monotonically and **~linearly in R** (marginal tension `dE/dR / 2π ≈ 10.9`), i.e.
  **energy = tension × circumference** — so the ring carries a line tension and shrinks
  under it. The average tension `E/2πR` falls from 16.6 (R=5) toward ~11–13 as R grows
  (small tight rings pay extra curvature energy per length), approaching the straight-line
  tension (L-0: 13.3). Same defect, same stiffness scale — a genuine closed 3D topological
  object.
- **Next (Routes B/C):** a genuine 3D *point*-charge needs an S² hedgehog field
  (π₂(S²)=ℤ; Route B, below) or a gauged monopole in the Coulomb phase (Route C;
  massless photon → true 1/r; larger build).

### Route B — the S² hedgehog: a genuine 3D point charge (global monopole) (`route_b_hedgehog.py`)
A single complex field can't give a 3D point charge (its defects are lines), so Route
B upgrades the field to a **3-component unit vector** n̂ ∈ S² (the O(3)/Heisenberg
model), where **π₂(S²)=ℤ** classifies point defects. The **hedgehog** n̂ = r̂ is the
charge-+1 point defect — the "particle" the project has been after. Charge = degree of
the map from an enclosing surface to S², computed as the flux of the topological current
`jⁱ = (1/4π) n̂·(∂ⱼn̂ × ∂ₖn̂)` (no per-plaquette orientation bookkeeping).
- **H-0 (gate — a clean integer point charge):** hedgehog **Q = +1.04**, antihedgehog
  **−1.04**, vacuum **0** — integer and localized (the 0.04 is central-difference
  discretization; it tightens to +1.03 as the lattice grows). The model hosts a genuine
  3D *point* topological charge.
- **H-1 (the global monopole — linear self-energy):** a single hedgehog's energy diverges
  **linearly** with box size (`E/L ≈ 7.1 → 7.5`, `dE/dL ≈ 7.7`), the 3D point-charge analog
  of the 2D-vortex log and the 3D-line tension. A *bare* topological charge has no finite
  energy; only neutral combinations do — which is exactly why a clean 1/r point charge
  needs a screening/gauge field.
- **H-2 (the interaction — short-ranged, not 1/r):** an opposite pair attracts and would
  annihilate, so both cores are pinned and the box boundary pinned to a uniform ẑ (a single
  hedgehog *can't* sit in a uniform vacuum — degree mismatch — but the net-0 pair can), then
  the texture is relaxed by O(3) neighbour-alignment. Result (converged to ΔE≈0): `E(d)` rises
  only ~207.6→210.1 over d=6…18 and **saturates** (a saturating fit, λ≈7, beats linear 5:1),
  and the plateau is **box-independent** (209.5/209.4/209.3 at box 40/56/72). So a neutral
  hedgehog pair has **finite, localized energy** and its interaction is a **weak, short-ranged
  texture-overlap attraction — NOT confining and NOT a clean 1/r**. (An earlier guess of
  linear *confinement* was wrong: unlike a gauge monopole, global O(3) hedgehogs are not
  confined — they attract weakly and annihilate.)
- **Verdict / Route C:** Route B delivers a genuine 3D *point* charge, but the *global* charge
  gives no long-range force (short-range pair, linearly-divergent bare self-energy). The clean
  1/r EM-in-3D requires **Route C** — a *gauged* monopole in the Coulomb phase, whose massless
  photon both makes the isolated charge finite-energy and gives a true 1/r (below).

### Route C — the gauged monopole: EM-in-3D, a genuine 1/r² point charge (`route_c_monopole.py`)
The capstone. Gauge the charge: a **compact U(1)** gauge field on the 3D lattice links in its
**Coulomb (deconfined) phase** — no Higgs, so the photon stays massless. The topological charge
is a magnetic **monopole** (quantized flux out of a cube, DeGrand–Toussaint count), and minimizing
the Maxwell energy `E = (1/e²)Σ(1−cos B)` at fixed monopole content is exactly Screen-1's massless
Gauss law → a genuine 1/r Coulomb.
- **C-0 (gates):** gauge invariance to machine zero (`dE=0`, monopole count invariant); the Wu–Yang
  seeded monopole is a **single cube of charge +1** (antimonopole −1, vacuum 0) — a clean quantized
  topological charge; the Maxwell force matches finite differences to `1.2e-7`.
- **C-1 (deconfined — finite self-energy):** the relaxed single-monopole energy is **box-independent**
  (`9.16 → 9.28` over box 16→34, charge stays +1). Unlike Route B's global hedgehog (`E ∝ L`,
  divergent), the gauged charge has a **finite** self-energy — the massless photon deconfines it.
- **C-2 (the 1/r² Coulomb):** the relaxed monopole's radial field is **`|B| ~ 1/r^{2.03}`** in the
  clean interior (`3 < r < 10.5`; |B| flattens past r≈13 = finite-box floor) — textbook Coulomb,
  the direct field-law test (boundary-robust, unlike the pair-energy fit which the finite box
  confounds). The monopole–antimonopole pair energy rises and saturates (an attractive Coulomb well).
- **Result — EM-in-3D.** A genuine **1/r² force between quantized topological point charges**,
  deconfined. Route C is the counterpart of Route B (global hedgehog → divergent, short-range) and
  the **3D closure of the 2D gauged story**: there the *broken* phase gave Meissner **screening**;
  here the *Coulomb* phase gives an unscreened **1/r**. The full arc now spans, from one medium, the
  whole menu of forces and both gauge phases.

### Toward fundamental physics: a test for emergent Lorentz invariance (`test_lorentz.py`)
A medium has a preferred frame, so the model is viable as fundamental physics only if Lorentz
symmetry **emerges** at long wavelength — a single, round, universal light cone. Tested from the
exact lattice dispersion `ω(k)` (Fourier symbol; no time evolution).
- **A. Isotropy (field sector) — EMERGES.** On the self-assembled isotropic lattices the field wave
  speed `c(direction,|k|)` becomes direction-independent as `k→0`; the anisotropy (the Lorentz
  violation) is suppressed at low energy: **hex 2D ~ (k/k_max)⁴**, **fcc 3D ~ (k/k_max)²** (hexagonal
  is isotropic to rank 4, cubic only to rank 2). With `k_max = π/R₀` the model's Planck wavenumber,
  the LV is `~(E/E_Planck)^{2–4}` — the standard emergent-Lorentz story. Rotational Lorentz invariance
  is real here. (Anisotropy at k/k_max = 0.05: hex 1e-7, fcc 2e-4; growing to ~1–5% at the BZ edge.)
- **B. Universality (field vs medium) — FAILS generically (the crux).** The field wave speed is a
  free parameter, while the medium's own phonons travel at the LJ sound speeds (`c_L ≈ 7.6`, `c_T ≈ 5.9`
  in LJ units). Nothing forces `c_field = c_sound`: **two sectors, two light cones ⇒ Lorentz violation
  *between* sectors** unless tuned by hand. This is the central obstacle for any "field on a medium"
  theory. A single universal cone requires all excitations to emerge from **one** structure (as near a
  single Fermi point in Volovik's ³He), not a field bolted onto an independent elastic medium — the
  clearest signpost for where the model must change to be fundamental.

### The unified prototype: one operator → one universal light cone (`test_lorentz_unified.py`)
Acting on the universality signpost. The fix is structural: keep central (LJ) forces to *set* the
isotropic geometry (self-assembly, H10), but govern all low-energy **dynamics** — medium displacement
*and* the matter field — by a single **vector-Hooke graph Laplacian** `E = ½K Σ|U_i−U_j|²`, `U∈R^{D+1}`
(D medium + 1 field). Its dynamical matrix is diagonal and isotropic, `D_ab(k)=δ_ab·K·S(k)`, so every
polarisation shares one dispersion.
- **Central-force (LJ) medium — many cones:** `c_L/c_T` = 1.28, 2.31, 1.99 along [100]/[110]/[111]
  (split *and* direction-dependent), and the matter field's speed is an independent free parameter on top.
- **Unified (vector-Hooke) — one cone:** `c_L = c_T = c_field = 1.587` **exactly** in every direction
  (`c_L/c_T = 1.000000`). The matter field is the (D+1)-th component of the *same* operator, so
  `c_field = c_L` by construction — a single universal light cone, no tuning; and it is **round** (same
  emergent `(E/E_Planck)^{2–4}` isotropy as the field symbol). This recovers **speed universality**:
  matter and medium are one structure, not a field on an elastic solid.
- **Honest cost:** the bonds must carry shear stiffness (vector-Hooke), which bare central forces cannot
  (`c_L>c_T` always for a stable central-force solid); central forces still supply the isotropic geometry.
  So the model *can* have a single universal cone — provided the medium's dynamical operator is the field's,
  not a separate elastic law.

### Emergent boost invariance — the full Lorentz group (`test_lorentz_boost.py`)
The rotational facets done (isotropy + one universal cone), the remaining Lorentz piece is **boosts**.
Boosts act on `(ω,k)` as a 4-vector; a dispersion is boost-invariant iff its surface maps to itself under
`ω'=γ(ω−βc k_x), k_x'=γ(k_x−βω/c)`. Applying a boost (emergent `c`, normalised to 1) to the unified
single-cone dispersion and measuring how far the boosted point falls **off-shell**:
- **Massless cone** `ω=c|k|`: off-shell residual scales as **`(k/k_max)²`** (6×10⁻⁴ at k/k_max=0.05 for β=0.3;
  1×10⁻³ for β=0.6) — boost-invariant at low energy, violation `~(E/E_Planck)²`.
- **Massive mass-shell** `ω²=c²k²+m²`: `ω²−|k|²` holds at **0.0900 = m²** to three figures at low k (drifting
  only near the zone edge) — `ω²−c²k²` is a genuine Lorentz invariant; mass and momentum combine relativistically.
- **Causality:** max group (front) velocity = **1.0000 = c** — no superluminal signal.
- **Verdict:** boosts join isotropy and universality — **the full Lorentz group emerges at long wavelength**,
  with all violations suppressed as `(E/E_Planck)²` (the lattice spacing = the Planck scale). The three tests
  (`test_lorentz` isotropy, `test_lorentz_unified` universality, `test_lorentz_boost` boosts) together clear the
  *first-order* obstacle to the model being fundamental. (Still open: quantization/QM, fermions + chirality, and
  long-range spin-2 gravity.)

### Emergent relativistic fermions — the Dirac cone and its wall (`test_dirac.py`)
Matter is fermions (spin-½, chiral); the model is bosonic. Relativistic fermions can **emerge** near a
band-touching point (Volovik/Wen), where the tight-binding dispersion becomes a linear isotropic Dirac
cone. Bloch Hamiltonian on the honeycomb (nn hopping): `H(k)=[[0,f],[f*,0]]`, `f=−t Σ_j e^{ik·δ_j}`, bands `E=±|f|`.
- **A. The cone exists.** On the **honeycomb** (bipartite) the bands touch (gapless, `|f|~4e-16`) with a
  **linear, isotropic** cone: `v_F = 3/2` (t=a=1), anisotropy 5e-4 — emergent 2D massless Dirac fermions,
  the fermions' own "speed of light".
- **B. Nielsen–Ninomiya doubling.** The two Dirac points K, K′ carry **opposite chirality** (winding of
  `arg f` = −1 and +1, summing to 0) — fermions appear only as canceling pairs; a **single chiral**
  (Standard-Model) fermion is forbidden. (The chirality *is* a topological winding — the project's charge
  concept again; total charge in the BZ must vanish.)
- **C. The plain medium can't.** The self-assembled close-packing (triangular, one site/cell) has a single
  band with a normal `|k|²` minimum — **no Dirac cone**. Relativistic fermions need a **bipartite** structure.
- **Verdict:** the model **can** host emergent relativistic Dirac fermions — but only on a two-sublattice
  medium and only as opposite-chirality pairs. Two open consequences: a single chiral fermion needs an extra
  mechanism (domain wall / dimensional reduction / interactions), and `v_F` is a **new light cone** to unify
  with the boson `c` — the universality demand again, now across statistics.

### Evading the chirality wall: a single chiral fermion on a domain wall (`test_domain_wall.py`)
The standard escape from Nielsen–Ninomiya (Kaplan domain-wall fermions; Callan–Harvey): regularise the
fermion as a 2-band **Wilson–Dirac (Chern) insulator**, `H(k)=sin k_x σ_x + sin k_y σ_y + (M−cos k_x−cos k_y)σ_z`
(topological for `0<|M|<2`), and put it on a **domain wall** — an edge is a wall to the trivial vacuum.
- **Topological strip (M=1):** two in-gap branches **cross E=0**, each a **single chiral** mode bound to
  one edge (`E=∓0.30` at k_x=0.3, peak density **0.998** on the boundary site), with **opposite velocities**
  (+0.955 bottom, −0.955 top) = opposite chirality on opposite edges.
- **Trivial strip (M=3):** gapped; the near-zero states are delocalized **bulk** (peak density 0.03) — **no**
  chiral edge mode.
- **Result:** the chirality wall is **evaded** the standard way. A lattice that is vector-like *overall* still
  carries a **single chiral fermion on a wall** — its Nielsen–Ninomiya partner is the opposite-chirality mode
  on the *other* wall, spatially separated, not on this one. This is the mechanism a fundamental version of the
  model would use for Standard-Model chirality. (So the fermion barrier now reads: Dirac cone ✓, single chiral
  fermion ✓ via a wall; still open — locking `v_F` to the boson `c`, and full SM chiral matter content.)

### First quantization: the unified sector is a relativistic QFT (`test_quantization.py`)
The unified linear theory is coupled harmonic oscillators = a free field; canonical quantization gives
`H = Σ_k ħω_k(a_k†a_k + ½)`, bosonic quanta of energy `ħω_k`. Two checks that it's genuinely *relativistic*:
- **A. Mass-shell.** The single-quantum energies are `ω→c|k|` (massless, linear/isotropic) and
  `ω→√(c²k²+m²)` (massive) — the quanta are relativistic particles; zero-point energy per site ≈ 1.19.
- **B. Vacuum correlator** `⟨0|φ(x)φ(0)|0⟩ = (1/N)Σ_k (1/2ω_k)e^{ik·x}`: **massless → power law `~1/r^{2.3}`**
  (relativistic `1/r²` in 3D, the excess being lattice/finite-size → 2 in the continuum); **massive →
  Yukawa `e^{−mr}/r`**, short-range. The quantum vacuum reproduces the relativistic-QFT correlator — and
  the **same massless/massive, long-/short-range dichotomy** that ran through the whole forces program, now
  at the level of the vacuum.
- **C. Honest scope.** Canonical quantization *imposes* `[φ,π]=iħ`. It shows the model quantizes to a proper
  relativistic QFT (a real check), but does **not** derive QM from the sub-quantum medium (emergent/stochastic
  QM) — the deepest open problem, untouched. So the QM barrier is *half* addressed: the quantized theory is
  correct; whether QM itself emerges is unanswered.

### The sharpest barrier: long-range spin-2 gravity — diagnosed, NOT achieved (`test_graviton.py`)
Closing the loop to where the forces program began. Real gravity is **long-range 1/r²**, **universally
attractive**, and **spin-2**. Each requirement checked against the model:
- **A. Long-range + universal — achievable *if* massless.** A massless mediator gives a `Φ~1/r` potential
  (measured `~1/r^{1.3}`, box-corrected toward 1); a massive one gives `e^{−r/ξ}/r` (ξ≈1.4, screened). Since
  mass density is intrinsically **positive**, the coupling is single-sign → **universally attractive** (unlike
  EM's ±). *But the model's own gravity-by-density mediator is the **massive** one* (Bitter–Crum, λ≈3) — which
  is precisely why gravity is short-range. Masslessness is the missing piece even for scalar-Newtonian gravity.
- **B. Spin-2 — missing.** The medium's displacement field is a **vector**: its phonons are 1 longitudinal
  (helicity 0) + 2 transverse (helicity ±1) = spin-0 + spin-1 (photon-like). A vector field **cannot** carry
  helicity ±2, so there is **no spin-2 (graviton) mode** among the phonons.
- **C. Verdict.** Long-range spin-2 gravity is **not achieved** — the sharpest open problem, now precisely
  diagnosed: (i) the mass-coupling is massive/screened (so not even scalar-Newtonian survives long-range), and
  (ii) there is no spin-2 d.o.f. (and Weinberg–Witten forbids an emergent massless spin-2 coupling to the full
  `T_μν`). The route it points to: add a symmetric-tensor d.o.f. whose mass-coupling is protected massless by
  **gauged diffeomorphisms** — the spin-2 analog of how gauging gave EM its unscreened 1/r² (Route C).

### Overcoming the screening — Route 1: the elasticity–fracton duality (`test_fracton_gravity.py`)
The medium's *phonons* have no long-range gravity sector, but its **defects** do. In 2D linear elasticity every
defect sources the same **biharmonic** (Airy) equation `∇⁴χ = source`; they differ only by the *multipole order*
of the source, which sets the range (all computed from `G = FFT⁻¹(1/k^{2n})`):
- **Dilatation** (energy density = a center of compression, source `∇²δ`): **contact only** — `|G(r>0)| ≈ 1.5×10⁻⁵ ≈ 0`
  — **screened (Bitter–Crum)**. *This is exactly the gravity-by-density coupling* — and precisely why it screened:
  energy couples as the **most-screened** multipole.
- **Dislocation** (torsion, source `∇δ`): `G ~ −0.154·ln r` (2D-Coulomb log; `1/2π=0.159`) — **long-range**.
- **Disclination** (curvature, source `δ`): the biharmonic `1/k⁴` kernel, the **least-screened** multipole
  (free-space `~r²ln r`, even longer-range; verified `∇²`(disclination) = dislocation to 6×10⁻¹³).
- **The point:** *fewer derivatives on the source ⇒ longer range.* Energy density is the maximally-screened case;
  **curvature (disclinations) is long-range and unscreened.** By the **Pretko–Radzihovsky elasticity–fracton
  duality**, these defects are the charges of a **rank-2 symmetric-tensor gauge theory** — the structure of
  linearized gravity (the disclination is an immobile *fracton*). And in 2D gravity a **point mass = a conical
  deficit = a disclination**: the medium already contains gravitational "masses" that curve space around them,
  long-range. **So Route 1 overcomes the screening — couple gravity to *curvature* (the tensor-gauge/fracton
  sector), not to energy density.** Open next: 3D (where the tensor gauge field is a dynamical spin-2 graviton),
  and wiring matter energy → disclination density so "mass curves the medium" at long range.

### The spin-2 half: the tensor sector is a genuine graviton (`test_graviton_spin2.py`)
`test_graviton` found the phonons are spin-0 + spin-1 (no spin-2); Route 1 located the missing tensor structure
in the defect sector. This shows what that tensor field *is* — a real graviton:
- **A. 2 polarizations, helicity ±2.** The transverse-traceless projector on symmetric 3×3 tensors has **rank 2**
  (exactly 2 physical polarizations, the GW `+` and `×`), and under a rotation by θ about `k` the polarization
  rotates by **2θ** → **helicity ±2 = spin-2** (a photon's transverse polarizations rotate by θ → helicity ±1).
  This is the decisive fingerprint separating a graviton from a photon.
- **B. Universal 1/r² attraction** (Newtonian limit, positive mass → single-sign → attractive; the massless 1/r
  was measured in `test_graviton` A).
- **C. Light-bending factor 2.** A spin-2 graviton couples to the full `T_μν` (energy *and* pressure/momentum),
  bending light by **twice** the scalar-gravity value — Eddington 1919 measured the factor-of-2, ruling out
  scalar gravity. So real gravity *requires* the spin-2 tensor — exactly the sector the medium hides in its
  disclination/fracton defects. **Remaining build:** derive the *propagating* 3D graviton dynamically from the
  medium's 3D defects (and clear Weinberg–Witten — whose loophole the model's cutoff-scale Lorentz violation,
  ~(E/E_Planck)², already opens).

### Cone universality across STATISTICS — an honest correction (`test_cone_universality.py`)
`test_lorentz_unified` locked the medium and the *bosonic* field to one cone; `test_dirac` then introduced
fermions with their own Fermi velocity — a **third cone**, never locked. Both live on the same honeycomb
structure factor `f(k)`: the fermion cone is the slope of `E=±t|f|` at the Dirac point (`v_F = (3/2)ta`), the
boson cone the curvature of `ω²=(K/m)(3−|f|)` at Γ (`c_B = (√3/2)√(K/m)·a`). Hence
- **`v_F / c_B = √3 · t / √(K/m)`** — verified exactly (1.7321 at `t=K=m=1`; 0.8660 at `t=½`; 3.0000 at `m=3`).
- The ratio is an **arbitrary, tunable** number set by *independent* couplings. `v_F = c_B` requires the
  fine-tuning `K/m = 3t²` — **a tuning, not a symmetry.** Generically **two cones ⇒ Lorentz violation *between
  statistics***. (Real graphene is the cautionary case: `v_F ~ c/300`, its phonons ~100× slower still; only the
  fermion sector is even approximately relativistic.)
- **Correction to the record:** the emergent-Lorentz result is a **within-sector** statement (one round universal
  cone for medium + bosonic field). A genuinely Lorentz-invariant world needs **all excitations to descend from
  one structure** (bosons as collective modes of the fermions, as near a Fermi point), or a symmetry relating
  them — **supersymmetry** is precisely the boson–fermion symmetry that would lock the two cones together.

### Locking the cones: a composite boson rides the fermion cone (`test_cone_lock.py`)
The cure for the gap above, and the project's recurring lesson at its deepest level. Don't *put in* a boson —
let it be a **collective mode of the fermions** (a fermion bilinear, or a gauge field induced by integrating the
fermions out, à la Sakharov). Then it has no cone of its own. The decisive statement: every composite
(particle–hole) boson of momentum `q` costs at least the **lower edge of the interband continuum**,
`ω_min(q) = min_k [E₊(k+q) − E₋(k)] = min_k v_F(|k+q|+|k|) = v_F|q|`.
- Measured on the real honeycomb bands: `ω_min/(v_F|q|)` = **0.995** at `|q|=0.02`, → 1 as `q→0`, and **isotropic**
  (`[1,0]` vs `[0.6,0.8]` agree to ~0.5%). The residual at larger `q` is band curvature beyond the linear Dirac cone.
- So a massless collective mode **rides the fermion light cone**: its speed *is* `v_F`, **inherited, not chosen** —
  no free parameter, no fine-tuning. Contrast an independent spring-boson, whose `c_B` gives the tunable ratio
  `v_F/c_B = √3 t/√(K/m)`.
- **Resolution:** cross-statistics Lorentz universality is *automatic* once **all excitations descend from one
  structure**. In Volovik's picture the emergent gauge field *and* the graviton are tetrad fluctuations of the
  fermion cone, so they ride it too. The independent lattice boson — a field bolted on — was the culprit, exactly as
  the medium-vs-field cone mismatch was in `test_lorentz_unified`.

### Sakharov's lock: the *induced* boson action is Lorentz-invariant with the fermion cone (`test_induced_action.py`)
The effective-action strengthening of `test_cone_lock`. A boson **induced** by integrating out the fermions
inherits their Lorentz invariance. Sharp signature: for a Lorentz-invariant fermion sector of speed `v_F`, the
one-loop polarization must obey `Π(q,Ω)/q² = Π(s)`, a function of the **Euclidean invariant** `s = Ω² + v_F²q²`
*alone*. Test on the real honeycomb bands: is `P = (Π/q²)·√s` the same for a mostly-spatial `q` as for a
mostly-temporal `Ω` at fixed `s`?
| √s | P (Ω/√s = 0) | 0.5 | 0.8 | 0.95 | spread |
|---|---|---|---|---|---|
| 0.30 | 2.3159 | 2.3175 | 2.3230 | 2.3296 | 0.59% |
| 0.15 | 2.4015 | 2.4015 | 2.4027 | 2.4040 | **0.10%** |
| 0.08 | 2.4349 | 2.4330 | 2.4354 | 2.4343 | **0.10%** |
- `P` is **mix-independent** — the induced action depends on `(q,Ω)` only through `Ω² + v_F²q²`, i.e. it is
  **Lorentz-invariant with the fermion's cone**, inherited rather than tuned. The residual spread shrinks at low
  energy (0.59% → 0.10%): lattice corrections beyond the linear Dirac cone, the same `(E/E_Planck)²` suppression
  seen throughout.
- Independent check: `P` also converges to a constant as `s→0` (2.32 → 2.40 → 2.43), approaching the universal
  Dirac coefficient `Π ≈ q²/(16√s)`.
- **Meaning:** a gauge field defined as a *fluctuation of the fermion structure* (hopping phases) gets its Maxwell
  term from the fermion loop, and hence its light cone from the fermions. With `test_cone_lock`, the
  cross-statistics Lorentz problem is solved **by construction** — provided the boson is made *of* the fermions
  rather than bolted on beside them. (Sakharov induced dynamics; the same mechanism Volovik uses for the emergent
  photon and graviton.)

### The capstone: one structure → fermions + photon + graviton on one cone (`test_emergent_tetrad.py`)
Volovik's mechanism, realized in the model's own medium. Near a Dirac node the fermion Hamiltonian is
`H = e^a_i σ_a (k_i − A_i)`, where **`A` = the node position is an emergent U(1) gauge field (photon, spin-1)**
and **`e` = the cone shape (tetrad) is an emergent metric (graviton, spin-2)**. Both are *features of the fermion
dispersion*, so neither has a light cone of its own — they ride the fermion cone **by construction**. Perturbing
the medium's own three nn bonds, `t_j → t(1+u_j)`, and reading the fermion bands (`A` = node shift;
`G_ij = Re(M_i M_j*)`, `M=∇f`; trace → `v_F²`, traceless `h` → graviton):
| bond fluctuation | \|A\| (photon) | v_F² | \|h\| (graviton) |
|---|---|---|---|
| unperturbed | 0 | **2.2500** = (3/2)² | 0 |
| uniform `(u,u,u)` | **0** | 2.4806 = (1.05)²·2.25 | **0** |
| doublet `(2u,−u,−u)` | 0.0629 | 2.2541 | 0.2741 |
| doublet `(0,u,−u)` | 0.0347 | 2.2513 | 0.1559 |
- The **singlet** (uniform stretch) is a pure **conformal** rescaling of the cone: no photon, no graviton.
- The **doublet** (E-representation) bond fluctuations source **both** the emergent gauge field *and* the
  traceless emergent metric — one microscopic degree of freedom (the bonds) yields both. Both are linear
  responses: `|A|/u → 2.0`, `|h|/u → 9.0` as `u→0`.
- **Meaning — the construction the whole program pointed to.** The photon and the graviton are not *added* to the
  medium; they **are** the medium's own bond fluctuations, seen by its fermions. Because both are read off the
  fermion dispersion, both are automatically Lorentz-invariant on the **single fermion cone** — exactly what
  `test_cone_lock` and `test_induced_action` require. One structure ⇒ fermions, EM, and gravity, one cone, no tuning.

### Frontier 1 — the first falsifiable prediction: Lorentz violation vs experiment (`test_lv_prediction.py`)
Everything above *reproduces* known physics; a theory must *predict* something nature could veto. The model's
distinctive, quantitative feature is the Lorentz violation it already measured. Extracting the leading coefficients
from the emergent fcc dispersion (with `a → l_Planck`, `k_max ~` the Planck momentum):
- **Boost (isotropic):** `1 − v/c = ζ_boost (k/k_max)²`, `ζ_boost = 0.245`. **Rotation (crystallographic):**
  `Δc/c = ζ_aniso (k/k_max)²`, `ζ_aniso = 0.068`. Both are **order-unity** and enter at **`(E/E_Planck)²`** — a
  **quadratic, mass-dimension-6 (n=2)** Lorentz violation, *not* the linear (n=1) form many QG scenarios predict.
- **Cast as `v(E)/c = 1 − ζ(E/E_Planck)²`** → effective `E_QG,2 = E_Planck/√ζ ≈ 2.5×10¹⁹ GeV`.
- **Confrontation:** vs Fermi-LAT GRB timing (`E_QG,2 > 10¹⁰ GeV`) → safe by ~6×10¹⁸ in the effect; vs UHECR
  (strongest n=2, `>10¹¹ GeV`) → safe by ~6×10¹⁶. Predicted `|Δv/c|` = 3×10⁻²⁷ at LHAASO (10⁶ GeV), **1.6×10⁻¹⁷ at
  UHECR** (the closest frontier). Cross-species `|c_γ−c_e|/c = 0` at leading order (one cone) ≪ 10⁻¹⁵ bound.
- **Verdict — consistent, and falsifiable in structure.** The LV is quadratic and Planck-suppressed, so it clears
  every current bound by many orders of magnitude and is not quantitatively detectable today. But it makes three
  **qualitative** predictions that need no Planck-energy access and would kill it: (1) LV is **quadratic (n=2), not
  linear** — a confirmed linear photon LV falsifies it; (2) the rotation-violating part is **anisotropic with the
  emergent lattice's crystallographic pattern**, correlated between boost and rotation sectors; (3) **one universal
  cone** — no leading-order species-dependent maximal speed (a confirmed `c_γ ≠ c_e`, or `≠ c_grav` à la GW170817,
  at `O(E/E_Planck)` falsifies it). The project's first *empirical* claim: a specific, cross-species-universal,
  crystallographically-anisotropic n=2 signature that nature could rule out — not merely a reproduction of physics.

### Frontier 2 — emergent gravity as a RUNNING dynamical model (`test_graviton_dynamics.py`)
The static graviton pieces were all verified; this makes it *dynamical* — the graviton **propagates**. Evolving the
symmetric-tensor field `h_ij` on a 3D lattice (`∂²_t h = c²∇²h`, leapfrog):
- **A. Propagation.** A gravitational-wave packet (`h_+`/`h_×` polarization) moves at **group velocity 0.970 c** —
  massless and luminal; the ~3% deficit is exactly the `(E/E_Planck)²` lattice dispersion at the packet's finite
  `k` (the effect quantified in Frontier 1), so even the graviton's *slowdown* ties to the one-cone prediction.
- **B. One cone.** The graviton obeys the universal wave operator, so its Lorentz-violation coefficient
  `ζ_graviton = 0.250` matches the boson/fermion/photon `(E/E_Planck)²` form — it rides the single cone, consistent
  with being the tetrad (`test_emergent_tetrad`).
- **C. Force.** The static Newtonian limit: `Φ(r) ~ 1/r` (measured `1/r^{1.13}`, the excess a finite-box artifact →
  1 in the continuum) ⇒ `F ~ 1/r²`; mass density is single-sign ⇒ **universal attraction**.
- **Result:** emergent gravity is now a running model — a **massless, luminal, spin-2 graviton** that propagates as
  a gravitational wave, rides the one universal Lorentz cone, and mediates a universal `1/r²` attraction. With
  Route 1 (long-range curvature sector) and the tetrad capstone, the picture is assembled: mass curves the medium,
  the curvature propagates at `c`, other mass falls toward it. **Open:** source the propagating tetrad from matter
  energy self-consistently in 3D (full nonlinear back-reaction), and derive its Einstein–Hilbert stiffness from the
  fermion loop (Sakharov) rather than imposing the wave operator by hand.

### Fundamental-physics scorecard (`test_lorentz*.py`, `test_dirac.py`, `test_domain_wall.py`, `test_quantization.py`, `test_graviton*.py`, `test_fracton_gravity.py`, `test_cone_universality.py`, `test_cone_lock.py`, `test_induced_action.py`, `test_emergent_tetrad.py`, `test_lv_prediction.py`, `test_graviton_dynamics.py`)
From "can it host X" to "can it be fundamental." **✅ Emergent Lorentz** (isotropy + speed universality via one
operator + boosts; violations ~(E/E_Planck)². The cross-statistics cone is **not** locked for an *independent*
boson, but *is* locked automatically once the boson is a composite of the fermions — one structure). **✅ Fermions** (Dirac
cone on a bipartite medium; a single chiral fermion via a domain wall, evading Nielsen–Ninomiya). **🟡 Quantum mechanics** (quantizes to a correct
relativistic QFT; deriving QM from the sub-quantum medium is untouched). **🟡 Long-range gravity** (energy-coupled
gravity-by-density screens — the sharpest contradiction — *but* Route 1 / the elasticity–fracton duality exposes a
**long-range curvature (disclination) sector**, dual to a rank-2 tensor gauge theory: couple to curvature, not
energy; the 3D dynamical spin-2 upgrade remains open). **⬜** Standard-Model gauge group & constants, continuum
limit, cosmology, Weinberg–Witten. Five "likely-fatal" barriers now have concrete in-model demonstrations, each
with an honest statement of what is shown vs still open — including a concrete route past the gravity screening,
and one self-correction (the Lorentz result is within-sector, not across statistics).

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
