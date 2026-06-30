# Butler-Voss Condensate

A physics-**inspired** research model (Robert Voss; formerly the "Grain Fabric
Model"). Space is treated as an **active medium** sampled by a network of
**nodes** — a node is a physical locus that carries field state, not a bare graph
vertex, and its edges are real tension-bearing relations. Disturbances propagate
as ripples; nonlinearity self-traps them into localized "particles"; and the
medium's own structure is the candidate source of charge, spin, and gravity.

> A toy model, not a claim about the real universe. Every idea here is turned
> into a **measurement** so it can be confirmed or refuted from the numbers.

The conceptual stance (built up across the work below): the medium — the
"condensate" — is the substance; nodes sample/excite it; **charge is a
topological winding**, **spin is a Noether charge**; the localized structures sit
**a layer below quantum particles** (which would be their bound states); and
**gravity is meant to emerge from the medium**, not be inserted.

## Status at a glance

| | Hypothesis | Result |
|---|---|---|
| **H1** | focusing makes persistent particles | ✅ supported (oscillon lives full run; control disperses) |
| **H2** | gravity via tension-advection | ❌ not supported (leaky / disruptive) |
| **H3** | a particle is a moving pattern | 🟡 weak (moves in the push direction) |
| **H4** | results independent of lattice | ⚠️ artifact-sensitive (`square2d` CFL; `cubic3d` anisotropy) |
| **H5** | 2D vs 3D differ | ✅ both host oscillons; differ as expected |
| **H6** | charge = topological winding | ✅ confirmed (integer, conserved) |
| **H7** | spin = Noether charge | ✅ confirmed (Q-ball) |
| **H8** | only ~3–5 fundamentals | ✅ *derived*: `n=0, ±1` + Noether |
| **H9** | higher layers by confluence | ✅ opposite bind, `n=2` splits, scale gap ≈6.6× |
| **H10** | self-cohering medium | ✅ self-assembles to isotropic close packing, sets spacing |
| **P0–P1** | field on a self-assembled medium | ✅ charge survives (needs LSQ operator) |
| **P2** | field on a *moving* medium | ✅ charge survives reconnections; lost at melting |
| **P3** | gravity-by-density | 🟡 halves validated (3a/3b); closed loop (3c) stable + correctly-signed but sub-critical (no self-trapping yet) |

## The three bodies of work

### 1. Base engine — H1–H5 (`simulation.py`)
The original scalar-field engine: a displacement field on fixed lattices
(`hex2d`, `square2d`, `cubic3d`, `fcc3d`), with oscillon particles, tension, a
tension-advection gravity attempt, and particle tracking → CSV.

```bash
python simulation.py --live --experiment lump                 # watch an oscillon breathe
python simulation.py --headless --experiment collide --steps 3000
python simulation.py --headless --experiment lump --damping 1 --gravity 0   # energy-conservation check
```
Key findings: focusing produces persistent particles (H1); the tension-advection
gravity term mostly drains energy (H2); `cubic3d` is artifact-prone while
isotropic lattices (`hex2d`, `fcc3d`) are stable (H4).

### 2. Layer-below program — H6–H10 (three prototypes)
Promotes the field to a complex **order parameter**, so the medium can carry the
quantum numbers, and makes the lattice itself dynamical.

```bash
python prototype_complex.py       # H6–H8: charge = winding, spin = Noether; species census
python h9_binding.py              # H9: opposite charges bind, n=2 splits, scale gap
python prototype_mobile_nodes.py  # H10: mutually-attracting nodes self-assemble to a lattice
```
Findings: charge is an integer topological winding (H6) and spin a Noether charge
(H7); the stable species reduce to `n=0, ±1` plus the Noether family (H8);
opposite charges bind into composites while like charges repel (H9); and a cloud
of mutually-attracting nodes self-selects an isotropic close-packed lattice,
*setting* its own spacing (H10) — which is why a self-organizing medium avoids the
`cubic3d` artifact entirely.

### 3. Integration — Phases 0–2 (field on the self-assembled medium)
Puts the complex field **on** the mobile medium (field attached to nodes).

```bash
python integration_field_medium.py   # Phase 0–1: meshfree operator accuracy; vortex on a frozen medium
python integration_phase2.py         # Phase 2: vortex on a MOVING, rearranging medium
```
Findings: an irregular mesh **requires a least-squares meshfree Laplacian** (the
plain graph operator degrades ~16×) (P0); a vortex's **topological charge is
conserved** on a self-assembled medium (P1) and survives the medium *moving and
reconnecting* under it (~308 reconnections), failing only when the medium
**melts** (P2). The core shape is under-resolved on the irregular mesh — an open
quality issue, not a topological one.

### Phase 3 — gravity-by-density (in progress)
Close the loop two-way: field energy compresses the medium → denser nodes slow
the local waves → other excitations refract toward the mass = attraction. This is
the non-leaky replacement for H2's tension-advection gravity, and the project's
headline open problem.

```bash
python integration_phase3ab.py   # 3a refraction (sign check) + 3b compression, in isolation
```
```bash
python integration_phase3c.py    # close the loop on one lump (both couplings on)
```
Sub-phases **3a/3b** are done: a test wave-packet refracts *toward* an imposed
dense region (confirming the sign **denser ⇒ slower**), and an imposed field-energy
blob compresses the node medium. **3c** closes the loop on one lump: it now runs
**stably** (after smoothing the node force) and shows the **correctly-signed**
effect — the lump compresses the medium and spreads ~60% less than uncoupled — but
is **sub-critical**: it does not yet self-trap. Reaching a bound state (and **3d**,
two lumps drifting together) needs a stabilized field operator and a deeper
density→speed coupling, ideally from a single Lagrangian.

## Install

```bash
pip install -r requirements.txt    # numpy, matplotlib
```
Requires Python 3.11+.

## File map

| File | Purpose |
|---|---|
| `bvc_core.py` | Shared primitives: LJ medium + self-assembly, meshfree operators, lattice builders |
| `simulation.py` | Base scalar engine — lattices, dynamics, particle tracking, CLI (H1–H5) |
| `prototype_complex.py` | Complex order parameter; topological + Noether charge; species census (H6–H8) |
| `h9_binding.py` | Binding / confluence study (H9) |
| `prototype_mobile_nodes.py` | Self-cohering medium of mobile nodes (H10) |
| `integration_field_medium.py` | Field on a self-assembled medium; meshfree operators (Phases 0–1) |
| `integration_phase2.py` | Field on a moving, rearranging medium (Phase 2) |
| `integration_phase3ab.py` | Gravity-by-density, the two halves in isolation (Phases 3a/3b) |
| `integration_phase3c.py` | Gravity-by-density, closed loop on one lump (Phase 3c) |
| `CHEATSHEET.md` | Full reference: equations, parameters, all hypotheses (H1–H10) + integration |

`Simulation - Cheat Sheet.docx` is a personal copy — left untouched.
Simulation output goes to `condensate_runs/` (git-ignored).
