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
| `test_lorentz.py` | Test for emergent Lorentz invariance: field-sector isotropy emerges (LV ~ (E/E_Planck)²⁻⁴), but field and medium have different light cones (speed non-universality — the central obstacle) |
| `test_lorentz_unified.py` | Unified prototype: governing medium + field by one vector-Hooke operator gives a single universal light cone (c_L = c_T = c_field exactly), recovering speed universality |
| `test_lorentz_boost.py` | Emergent boost invariance: the massless cone and massive mass-shell are boost-invariant at low energy (violation ~ (E/E_Planck)²), signals sub-luminal — the full Lorentz group emerges |
| `test_dirac.py` | Emergent relativistic fermions: a honeycomb medium gives a linear isotropic Dirac cone (v_F=3/2), but fermions come as opposite-chirality pairs (Nielsen–Ninomiya); the triangular close-packing has no cone |
| `test_domain_wall.py` | Evading the chirality wall: a Wilson–Dirac (Chern) strip binds a single chiral fermion to each edge (opposite chirality, spatially separated) — the domain-wall route to Standard-Model chirality |
| `test_quantization.py` | First quantization: the unified sector quantizes to a relativistic QFT (quanta on the mass-shell; massless→power-law, massive→Yukawa vacuum correlator). Honest: it imposes QM, doesn't derive it |
| `test_graviton.py` | The sharpest barrier: long-range spin-2 gravity — diagnosed, not achieved. Massless mediator would give universal 1/r attraction, but the model's mass-coupling is screened, and the medium has no spin-2 (only spin-0/1 phonons) |
| `test_fracton_gravity.py` | Route 1 (elasticity–fracton duality): the medium's defects source a biharmonic tensor-gauge theory — curvature (disclinations) is long-range while energy/dilatation is screened. ⚠️ **Its conclusion that this "overcomes the screening" is RETRACTED** — see `test_disclination_force.py` |
| `test_graviton_spin2.py` | The spin-2 half: the tensor-gauge field is a genuine graviton — exactly 2 transverse-traceless polarizations carrying helicity ±2 (vs a photon's ±1), universal 1/r² attraction, light-bending factor 2 |
| `test_cone_universality.py` | Honest correction: the fermion cone `v_F` is **not** locked to the boson cone `c_B` (`v_F/c_B = √3·t/√(K/m)`, a tuning not a symmetry) — emergent Lorentz is a within-sector result |
| `test_cone_lock.py` | The cure: a *composite* boson (particle–hole / induced gauge field) rides the fermion light cone — `ω_min(q) → v_F\|q\|` — so cross-statistics Lorentz universality is automatic when all excitations descend from one structure |
| `test_induced_action.py` | Sakharov's lock: integrating out the fermions induces a Lorentz-invariant boson action — the polarization depends on `(q,Ω)` only through `Ω²+v_F²q²` (0.10% spread) — the emergent gauge field inherits the fermion cone |
| `test_emergent_tetrad.py` | **Capstone:** the medium's own bond fluctuations *are* the emergent photon (Dirac-node shift, spin-1) and graviton (cone deformation / tetrad, spin-2) — both read off the fermion dispersion, so both ride one cone |
| `test_lv_prediction.py` | **Frontier 1 — first falsifiable prediction:** quadratic (n=2), crystallographically-anisotropic, cross-species-universal Lorentz violation (`ζ~0.2`, `E_QG,2~2.5×10¹⁹ GeV`) — consistent with all bounds, falsifiable in structure |
| `test_graviton_dynamics.py` | **Frontier 2 — emergent gravity, running:** the graviton propagates as a gravitational wave (group velocity 0.97c, massless), rides the one universal cone (`ζ=0.25`), and mediates a universal `1/r²` attraction |
| `test_induced_gravity.py` | **Frontier 2 completion — Sakharov:** the graviton's Einstein–Hilbert kinetic term is *generated* by the fermion stress-tensor loop `⟨TT⟩` (Lorentz-invariant, spread → 0.37% as s→0) — its dynamics inherited from the fermion cone, not imposed |
| `test_emergent_qm.py` | **Frontier 3, step 1:** the wave half of QM emerges — the field's NR envelope spreads at the exact Schrödinger rate (`D` measured 0.830 vs predicted 0.833), with `ħ/2m = c²/2Ω` a material property of the medium. The Born rule / measurement remain the open moonshot |
| `test_born_rule.py` | **Frontier 3, step 2:** the Born rule as a stochastic equilibrium (Nelson/Valentini) — an ensemble started non-Born relaxes to `\|ψ\|²` (KL `4.4→1e-4`) under the medium's diffusion `ν=ħ/2m`; `\|ψ\|²` is an attractor, not a postulate |
| `test_double_solution.py` | **Frontier 3, step 3:** de Broglie's double solution — deriving guidance `v=∇S/m`. A soliton carrying its *own* phase drifts at exactly `k` (slope **1.000**): `λ=h/p` is a theorem of the medium. But a resting soliton is *not* steered by a *separate* pilot wave (slope **≈0**): guidance-by-a-distinct-wave and definite outcomes stay a postulate — the honest F3 boundary |
| `test_shielding.py` | **Why gravity is unshieldable.** Screening is *not* loss (a superconductor screens `B` with zero dissipation) — it is **neutralization**, and mass is unipolar. Gauss-law probe: a **dilatation** charge is cancelled exactly by the medium (gate `div u ≡ C·s` to `7.5e-16`) ⇒ shieldable ⇒ short-range; a **topological** charge is **exactly invariant** (`+1.000000000`) under any smooth response, stiffness contrast, or relaxation. Topological quantization does what "no negative mass" does in nature |
| `test_tetrad_shielding.py` | The tetrad graviton is **long-range but shieldable**. Bitter–Crum kills the strain's **trace**, not its **shear** — the deviator falls as `r^-2.01` (box-gated), and the tetrad *is* the traceless part, so it reads the one unscreened sector. But an intervening shell **attenuates it ~4×** — real gravity shows none. Diagnosis: the medium's moduli are a *free* background; in GR stiffness *is* mass and can only add |
| `test_disclination_force.py` | **Route 1's force law — it fails.** Real-space clamped disc (no periodic IR saturation), box-gated. Control reproduces the known dislocation log-repulsion. Measurement: two **like disclinations REPEL**, `\|dE\| ~ R^1.97` — wrong **sign** (gravity attracts) and wrong **direction** (force *grows* with R). A disclination is a **charge**, not a **mass** |
| `test_light_bending.py` | **First assault on the spin-2 wall: `γ = 0`, derived.** Does the model bend light? The PPN `γ` decides it (0 = scalar, 1 = GR's factor of two). Measured `γ = 1.4×10⁻⁵ ≈ 0`, and the *reason* is the result: a mass can only make the medium move its nodes (a displacement `u`), and displacement-derived strain is **compatible** → a **flat** metric → the coordinate change `x→x+u` → **pure gauge**, which bends no light (`η/|strain| = 1.8×10⁻¹⁵`). Genuine curvature exists **only** as *incompatible* strain — a disclination (`η = 0.13`) — but that sector **repels** (`test_disclination_force`). So the wall is structural: a fixed-background medium responds to mass by a diffeomorphism, and cannot make GR's spatial curvature. The route through is a metric dynamical in its own right, not a finer lattice |
| `test_incompatible_gravity.py` | **The door through the wall: the graviton lives in incompatible bond DOF.** The wall tested only node displacements. But each triangular-lattice site has 3 bond lengths vs a displacement's 2 DOF → bond fluctuations = 2 (displacement, gauge) + **1 incompatible (curvature)**. Shown: the compatible sector is flat (`η=1.8×10⁻¹⁵`) and bends nothing, but the incompatible sector (mass → curvature `η~ρ`, via an Airy potential, gate `η==ρ` to `6×10⁻¹⁰`) **does deflect light**, long-range, vs the compatible sector's exact zero. Reframes the frontier from "can it bend light" (yes, in this sector) to "is the coupling **Einstein** (`γ=1`)" |
| `test_einstein_source.py` | **Is the coupling Einstein? γ isolated to one number.** `γ=Ψ/Φ` decides light bending (γ=1 = GR's factor of two). [A] Mass sources *no* curvature elastically (`η/ρ = 1.2×10⁻¹⁵`) — equilibrium elasticity returns a displacement, so the model's native gravity (the scalar amplitude mode) is γ=0 (Nordström). [B] Only a topological disclination sources curvature — and it repels. [C] An induced `η=κρ` gives `γ=κ/(4πG)`: γ=1 → bend ×2. So the crossing is **one missing ingredient** — the Einstein source coupling (mass→curvature at γ=1), which the fermion loop must induce by coupling the metric to the *conserved* stress tensor (Weinberg) |
| `test_scale_fixing.py` | **The medium's absolute scale, from data.** Matching the measured `G` with an order-unity coupling forces `a₀ = l_Planck` (an *assumption* in `test_lv_prediction`, now a **consistency result**). And since the model predicts gravity's range is `1/m_A`, every graviton-mass bound measures the medium's distance from criticality: the strongest gives **1 part in 10¹²²** — essentially the cosmological-constant number. Falsifiable: the model predicts a **Yukawa**, not a pure `1/r²`; a measured graviton mass would *measure* the tuning, not kill the model |
| `test_collapse.py` | **Minimum mass to destroy the condensate — and why it's not a black hole.** P1 (mass = confluence of loci) is real: where `2gρ > a` the condensate is locally destroyed. But there is **no critical nucleus and no runaway** — `N(σ)` is monotonic and `N_min` is gap-independent (ratio 0.96 when `a` halves, not the predicted √2), because the condensate is the *true* vacuum so a normal bubble is never favorable. Matter *pins* a suppressed region; remove it and it heals shut. **No gravitational collapse exists in this model** — the third independent reason P2's black hole isn't here (with `ρ_horizon ∝ 1/M²` and scalar gravity not bending light). *An honest self-correction: my critical-nucleus prediction was refuted by the measurement* |
| `test_critical_gravity.py` | **✅ WORKING long-range gravity.** The mediator is the condensate's **amplitude mode**: the phase is a Goldstone (shift symmetry ⇒ derivative couplings only ⇒ can never mediate a monopole force), but the amplitude is *unprotected*, so it is gapped — **which is why gravity always looked screened** — and it couples monopolarly to positive-definite energy, so **scalar exchange attracts**. Measured in 3D: **`λ·m_A = 1.00`** (gap read off the potential, not fitted), and dividing out the exponential leaves a pure **`R^-1.00`** core ⇒ `E_int = −C·e^{−R/λ}/R` **exactly**. At criticality: **Newton's law, `1/r²`, universally attractive.** 🟡 Scalar (Nordström) gravity — no light bending, no spin-2 waves |
| `test_tetrad_force.py` | **A field falloff is not a force.** Eshelby/**Crum**: two dilatation centres in an infinite *isotropic* medium have **zero** interaction. Measured: the isotropic force **collapses 69×** from `R≈10→22` (only a short-range contact term — the screened Phase-3 attraction). Control (isotropy broken) shows a real long-range force appears, **46×** larger — and **repulsive**. So the tetrad's `1/r²` field exerts **no long-range force** |
| `screening_gauged.py` | Gauged U(1) / Abelian Higgs (G-0/G-1/G-2): gauging the symmetry screens the vortex force (Meissner) — box-independent λ_L ~ 1/e, vs Screen-2's box-growing log |
| `screening_gauged_mobile.py` | Mobile-vortex check: vortices move under the screened force (overdamped, adiabatic gauge); the force law from motion confirms λ_L ~ 1/e |
| `CHEATSHEET.md` | Full reference: equations, parameters, all hypotheses (H1–H10) + integration |

Simulation output goes to `condensate_runs/` (git-ignored); rendered figures to
`figures/` (git-ignored). Personal documents (`*.docx`) are git-ignored.
