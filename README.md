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
| **P3** | gravity-by-density | ✅ **two masses attract** (3d): energy-conserving variational coupling makes two lumps drift together (sep 9.3→7.5), no attraction at β=0, stronger at higher β. Halves validated (3a/3b), coupling energy-conserving (bounded `g`). No single-lump self-bound soliton (medium nearly incompressible → gravity is real but weak) |

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

### Phase 3 — gravity-by-density (✅ two masses attract)
Close the loop two-way: field energy compresses the medium → denser nodes slow
the local waves → other excitations refract toward the mass = attraction. The
non-leaky replacement for H2's tension-advection gravity.

```bash
python integration_phase3ab.py          # 3a refraction (sign) + 3b compression, in isolation
python integration_phase3c.py           # close the loop on one lump
python integration_phase3_variational.py # both couplings from one energy functional (conserving)
python integration_phase3d.py           # THE gravity test: two masses drift together
```
The arc: **3a/3b** validated the two halves and the sign (waves refract *toward* a
dense region). **3c** closed the loop on one lump (self-focuses but leaks with the
hand-tuned coupling). The **variational** form derives both couplings from one
energy functional (energy-conserving by construction); a bounded `g(ρ)` then cured
a blow-up and fixed conservation. A single lump does **not** self-bind — the medium
is nearly incompressible, so gravity-by-density is real but *weak*. But self-binding
isn't needed: **3d shows two masses drift together** (separation 9.3 → 7.5),
with **no attraction at β=0** and stronger attraction at higher β, energy conserved
— gravitational attraction emerging from the medium.

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
| `integration_phase3_variational.py` | Gravity-by-density from one energy functional (energy-conserving) |
| `integration_phase3d.py` | The gravity test: two masses drift together (Phase 3d) |
| `operator_scorecard.py` | Symmetric field-operator candidates scored (symmetry/stability/accuracy) |
| `figures.py` | Renders the headline results to `figures/` (self-assembly, vortex/charge, gravity) |
| `interaction_energy.py` | Quantifies the gravitational contact strength: the interaction well `U(d)` |
| `density_response_3d.py` | 3D (ordered fcc) density response, small medium (3D-0/3D-1) |
| `density_response_3d_large.py` | 3D sparse scale-up (3D-2): gravity is screened (λ≈3.3), not Newtonian |
| `screening_diagnosis.py` | Screen-0: coupling strength β is not a lever — screening is intrinsic |
| `screening_gauss.py` | Screen-1: existence proof — the medium carries unscreened 1/r iff the mediating field is massless (screening = a mass term) |
| `screening_topocharge.py` | Screen-2: a conserved topological charge sources a long-range (2D-Coulomb log) force — the model hosts EM-like *and* nuclear-like forces |
| `screening_topocharge_3d.py` | The 3D topological defect (Route A): long-range survives into 3D as a vortex-line interaction (log per unit length); vortex ring is a closed 3D line (E = tension × circumference) |
| `route_b_hedgehog.py` | Route B: the S² hedgehog — a genuine 3D **point** topological charge (integer); bare self-energy diverges linearly (global monopole); the pair interaction is short-ranged, not 1/r |
| `route_c_monopole.py` | Route C: the gauged monopole (compact U(1), Coulomb phase) — **EM-in-3D**: finite deconfined self-energy and a genuine 1/r² Coulomb field between quantized topological point charges |
| `screening_gauged.py` | Gauged U(1) / Abelian Higgs (G-0/G-1/G-2): gauging the symmetry screens the vortex force (Meissner) — box-independent λ_L ~ 1/e, vs Screen-2's box-growing log |
| `screening_gauged_mobile.py` | Mobile-vortex check: vortices move under the screened force (overdamped, adiabatic gauge); the force law from motion confirms λ_L ~ 1/e |
| `CHEATSHEET.md` | Full reference: equations, parameters, all hypotheses (H1–H10) + integration |

Simulation output goes to `condensate_runs/` (git-ignored); rendered figures to
`figures/` (git-ignored). Personal documents (`*.docx`) are git-ignored.
