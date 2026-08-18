# Literature audit — prior art for every result, and what the simulation actually did

A companion to `EPISTEMIC_AUDIT.md` and `EXTERNAL_VALIDATION.md`. The epistemic
audit calibrated each result on the internal axis of *how much could it have come
out otherwise* (Output / Existence / Known-physics / Construction / Negative).
This document does the complementary, external check: **for every result cluster,
what is the prior art in the published literature, and is anything in the project
novel to physics?** It then records, honestly, the thing that survives that check
— **the project as a validated simulation that reproduces real laboratory and
astrophysical results** — and states precisely what that is and is not worth.

The audit was prompted by a direct question ("is there anything worth publishing?")
and carried out by a systematic literature search over every major result. The
short answer it reached is stated up front, because the honesty is the point.

## Bottom line

**No result in the project clears the novelty bar for a physics-discovery
publication.** Every quantitative result reproduces, instantiates, or independently
rediscovers physics that already exists in the literature — in two cases (the
world-crystal γ = 1 arc and sin²θ_W = 3/8) rediscovering specific papers almost
move-for-move. This is consistent with the epistemic audit's own **K** (known-
physics) tags; the literature search closes the last two gaps where the project
had hoped for an **O** (novel output).

**What genuinely survives is not a discovery but an instrument.** The project is a
single, self-contained, honestly-audited simulation platform that — fed known
microphysics — independently returns results matching real measurements, *including
ones it was never tuned to* (the GR quadrupole luminosity to sub-percent, the
Peters–Mathews inspiral rate, the GW170817 speed bound). Reproducing a real
measurement with a model built to contain the generating physics is a **consistency
check**, the strongest kind of validation — but it is not new knowledge about
nature, and it is not a prediction. The distinction is the subject of Part II.

---

# Part I — Prior-art concordance

For each cluster: the project's section(s), the canonical prior work, and the
honest novelty verdict. "Rediscovery" means the project reached it independently,
by measurement, without building from the prior paper — a real sign the framework
is internally sound, but *first* is what counts as novel, and none of these are
first.

## Foundations — solitons, charge, spin, the self-assembled medium

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| H1 (`simulation.py`) | Focusing makes persistent localized particles (oscillons) | **Bogolyubsky & Makhankov (1976)**, "pulsons"; **Gleiser (1994)** oscillons in φ⁴ — long-lived localized breathers are a 50-year-old, well-characterized phenomenon | Instantiation of known physics |
| H6 (`prototype_complex.py`) | Charge = integer topological winding, conserved | Standard **topological solitons / vortices** (winding number as a topological charge) — textbook | Instantiation |
| H7 | Spin = conserved Noether charge (Q-ball) | **Coleman (1985)**, "Q-balls," *Nucl. Phys.* B262 — non-topological soliton stabilized by a conserved Noether charge | Instantiation |
| H8–H10 | Species census (n = 0, ±1); confluence binding; self-assembled isotropic close packing setting its own spacing | Standard soliton spectra + **classical close-packing / self-assembly** condensed matter | Instantiation |
| P3 | Gravity-by-density: two lumps attract, energy-conserving | Emergent/analogue attraction in a medium — see gravity cluster below | Existence demo |

Nothing here is novel; the project's own README already labels these as
supported/confirmed *hypotheses about the toy model*, not claims about nature.

## Emergent Lorentz, fermions, quantum field theory

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| §8.1 | One emergent Lorentz cone at the elastic point; LV ∼ (E/E_P)² | **Chadha & Nielsen (1983)**, "Lorentz invariance as a low-energy phenomenon," *Nucl. Phys.* B217; **Volovik** Fermi-point scenario | Rediscovery of a known mechanism |
| §8.2 | Emergent relativistic Dirac fermions (cone + chiral domain-wall mode) | **Semenoff (1984)** graphene Dirac cone; **Jackiw–Rebbi (1976)** domain-wall zero mode; **Nielsen–Ninomiya** doubling | Textbook instantiation |
| §8.3 | Canonical quantization → relativistic QFT | Standard QFT | Known |
| §8.5, §8.14 | Induced/composite gauge bosons riding v_F; emergent Yang–Mills | **Bjorken (1963); Terazawa–Akama; Volovik–Sakharov** induced gauge fields | Known mechanism; group is an input |
| §8.15 | Chirality consistent by anomaly inflow (Chern = modes = pumped charge) | **Callan–Harvey (1985)** anomaly inflow | Instantiation; integers measured |

## Quantum mechanics from the medium

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| §8.6 | Schrödinger equation, ℏ scale, Born rule as a stochastic attractor | **Nelson (1966)**, *Phys. Rev.* 150 — Schrödinger from a configuration-space diffusion; **Valentini (1991)** sub-quantum H-theorem — Born rule as quantum-equilibrium relaxation | Rediscovery; the project's own WP-38 catch already re-scoped it to "the medium fixes ℏ, given the wave" — exactly the Nelson/Valentini boundary |
| §8.50 | Einselection + Born branches from the medium's phonon bath; single outcome stays a postulate | **Zurek (2003)**, *Rev. Mod. Phys.* 75, 715 — decoherence, einselection, pointer basis | Instantiation; the honest "single outcome is still a postulate" is the standard limitation |

## Gravity — the whole arc

This is the cluster where the project invested most, and where the literature match
is most complete. See `EPISTEMIC_AUDIT.md` §8.64–§8.71 for the internal story; here
is the external one.

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| §8.7, §8.12, §8.57 | Induced dynamical spin-2 graviton; γ decided by the substrate, no finite continuum G | **Sakharov (1967)** induced gravity — the Einstein–Hilbert term as vacuum elasticity | Instantiation; §8.57's "it's all Sakharov, no scheme-independent finite part" is the modern understanding (Visser, *Sakharov's induced gravity: a modern perspective*) |
| §8.5, §8.10, §8.13 | Emergent gravity from the medium; cosmological constant dissolved (vacuum gravitates −P, zero at equilibrium) | **Volovik**, *The Universe in a Helium Droplet* (2003); Fermi-point emergent gravity — **CC = 0 from thermodynamic equilibrium / zero pressure is Volovik's argument**, and appears verbatim in Kleinert–Zaanen (see below) | Rediscovery |
| §8.20–§8.24 | Radiation is spin-2, monopole forbidden; GR quadrupole luminosity sub-percent; Peters–Mathews inspiral | **Peters & Mathews (1963); Peters (1964)** — the GR radiation and orbital-decay laws | Reproduction of GR (see Part II — these are un-tuned matches to real inspiral data) |
| §8.31 | Spin-2 self-coupling fixed by the bootstrap | **Deser (1970)** — GR from consistent self-coupling of a massless spin-2 field | Theorem instantiated |
| §8.64 | Why γ = 0 is *forced*: fracton-dual diagnosis (elastic solid ⇒ scalar-charge ⇒ Nordström) | **Pretko & Radzihovsky (2018)**, "Fracton-Elasticity Duality," *PRL* 120, 195301 (arXiv:1711.11044) — elasticity is dual to a symmetric-tensor gauge theory; disclinations = immobile fractons | The duality is Pretko–Radzihovsky; applying it to diagnose emergent-gravity's γ is the project's framing, but the physics is theirs |
| §8.65–§8.71 | The world-crystal escape: massless graviton = second-gradient "floppy" crystal; graviton is the dual/disclination field; matter sources disclinations; gauge invariance ⇒ Einstein–Hilbert ⇒ γ = 1 | **Kleinert (1987)**, *Ann. Physik* 44, 117; **Kleinert & Zaanen (2004)**, "Nematic world crystal model of gravity," *Phys. Lett. A* 324, 361; **Kleinert (2005)**, "Emerging Gravity from Defects in World Crystal," *Braz. J. Phys.* 35, 359 | **Rediscovery, near-complete.** See the concordance below — the entire arc reproduces Kleinert's program |
| §8.71 | Massless + two-derivative + diffeo-invariant ⇒ uniquely linearized EH | **Fierz–Pauli (1939); Weinberg (1964–65); Deser (1970)** — the uniqueness theorem the project invokes by name | Known theorem, correctly applied |

### The Kleinert concordance (§8.64–§8.71 ↔ Kleinert–Zaanen 2004 + Kleinert 2005)

The gravity resolution the project built as an eight-section arc is, equation for
equation, Kleinert's world-crystal gravity:

| Project section | Kleinert equation / statement |
|---|---|
| §8.64 elastic solid ⇒ γ = 0 (first-gradient elasticity is not Einstein) | Kleinert 2005, Eq. (17): first-gradient crystal gives interaction ∝ 1/(∂²)², **"This is not the Einstein action"** |
| §8.65/§8.67 world crystal = second-gradient "floppy" massless graviton, ω² ∝ q⁴ | Kleinert 2005, Eq. (30): higher-gradient elastic energy `A = μ∫[∂(u−uᵖ)]²`; "leading elastic terms vanish" |
| §8.68 graviton = dual/disclination field; disclination = Einstein tensor | Kleinert 2005: "ηµν becomes the **Einstein tensor**"; Kleinert–Zaanen 2004: "the disclination density θ_ij **is** the Einstein tensor R_ji − ½g_ji R" |
| §8.69 matter sources curvature iff disclinations deconfine (Y = 0) | Kleinert–Zaanen 2004: dislocation **condensation ⇒ Meissner screening ⇒ confining forces become Newtonian 1/R** (Eqs. 1.19–1.20) |
| §8.70 matter is a genuine source via its stress tensor | Kleinert–Zaanen 2004, Eq. (1.21): matter couples by "the usual Einstein interaction `E_int = ∫h_ij T^ij`," giving Newton's law at `c = 8πG` |
| §8.71 gauge invariance ⇒ Einstein–Hilbert ⇒ γ = 1 | Kleinert 2005, Eq. (29): recovers the linearized Einstein–Hilbert action `(1/4κ)∫hG`; γ = 1 is immediate |
| §8.13 CC = 0 | Kleinert–Zaanen 2004: "our model has automatically a **vanishing cosmological constant** … the pressure is zero" |
| §8.1 emergent Lorentz | Kleinert–Zaanen 2004: **fluctuation-induced isotropy** (Heisenberg fixed point) restores rotational/Lorentz invariance |

The one thing Kleinert does not do in these papers is phrase the output in PPN
language or ray-trace light deflection to read γ = 1.0000 numerically. But
"linearized Einstein–Hilbert ⇒ γ = 1" is a textbook one-line corollary, so this is
a restatement, not an increment. The fracton-duality *language* wrapping the same
physics is itself published (Pretko–Radzihovsky 2018; and the "dual gravity" line,
Beekman–Zaanen and the 2023 Lifshitz-to-dual-gravity work).

## Standard-Model sector

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| §8.63 | **sin²θ_W = 3/8** from induced couplings (1/g² ∝ Tr T², equal traces) + anomaly hypercharges, no GUT group | **Terazawa, Akama & Chikashige (1976)**, "What Are the Gauge Bosons Made of?", *Prog. Theor. Phys.* 56 — composite gauge bosons from an NJL Lagrangian give **sin²θ_W = 3/8 for fractionally-charged quarks, no GUT assumed**; also **Georgi–Glashow (1974)** for the SU(5) value by the standard route | **Rediscovery, 49 years old.** The mechanism and the number are Terazawa–Akama–Chikashige's. Two-loop running on top is standard SU(5)-normalization RG |
| §8.42 | Hypercharges fixed uniquely by anomaly cancellation | Known SM property — anomaly cancellation quantizes/fixes hypercharge up to a discrete choice (Geng–Marshak; Minahan–Ramond–Warner) | Verification of a known SM fact |
| §8.43, §8.48 | Gauge group + three generations as band topology / winding | **Topological-defect zero-mode** generation models (e.g. three families from a defect of topological number three, hep-ph/0011095); the project honestly flags the *category* is fixed, the group and the number three are **not derived** | Reframing; the open problem (why three) is acknowledged unsolved in the literature too |
| §8.44 | Fermion lattice needs O(1) angular rigidity (λ* ≈ 0.57) | Project-specific lattice measurement; not a physics claim about nature | Internal measurement |

## The one prediction — Lorentz violation and the UHE frontier

| Project | Result | Prior art | Verdict |
|---|---|---|---|
| §8.8, §8.39, §8.51 | Falsifiable prediction: n = 2 (∝ (E/E_P)²), subluminal, species-universal LV, coefficient ζ ≈ 0.245 | **Quantum-gravity phenomenology**: Amelino-Camelia; **Jacobson, Liberati & Mattingly (2003–06)** modified dispersion / threshold bounds (astro-ph/0505267) | The *framework* (modified dispersion from discrete/emergent spacetime) is standard; the specific coefficient is the project's — see below |
| §8.52 | GW170817 one-cone speed bound reproduced; model's Δv/c = 0 sits inside | **Abbott et al. (2017)**, *ApJL* 848, L13 — bound −3×10⁻¹⁵ … +7×10⁻¹⁶ | Reproduction + a genuine pass of the universality half (Part II) |
| §8.54–§8.61 | n = 2 GZK threshold; proton-compositeness suppression; DGLAP running to UV | **Coleman–Glashow (1999)** LV thresholds; **Greisen–Zatsepin–Kuzmin** cutoff; standard **DGLAP** evolution | Standard framework; the model's coefficient straddles the GZK bound, then recedes below reach when run to the UV (documented in `EXTERNAL_VALIDATION.md`) |

The LV prediction is the closest the project comes to something *its own* — not the
mechanism (standard QG phenomenology) but the specific structural claim (n = 2,
subluminal, one universal cone, near-Planckian ζ). It is not novel physics in kind,
and §8.61 pushed its coefficient ~1–2 orders below current UHECR sensitivity, so it
is hard to test now. It remains the only place the project offers a number that was
not already in a prior paper — which is why Part III points to it as the seam to
mine for genuine novelty.

---

# Part II — Simulation → laboratory results

This is the part worth stating carefully, because it is the project's real and
defensible strength, and because it is easy to overclaim.

## The catalogue

Where does the simulation touch reality? Three tiers, by how much the match could
have failed:

### Tier 1 — Un-tuned matches to real measurements (the instrument working)

Results the simulation was **not** built or tuned to produce, that nonetheless land
on real measured values:

| Result | Real-world target | Match |
|---|---|---|
| §8.23 GR quadrupole luminosity | The GR radiation law (measured in binary-pulsar spin-down) | Sub-percent, nothing fitted |
| §8.24 Orbital decay | **Hulse–Taylor / Peters–Mathews** inspiral rate | Reproduces the GR rate |
| §8.52 One-cone speed | **GW170817/GRB 170817A**, Δv/c ∈ [−3×10⁻¹⁵, +7×10⁻¹⁶] | Model predicts exactly 0 — a genuine pass against a 2017 measurement |

These validate the platform: put in the microphysics, and the machine returns
numbers that agree with the sky. That a toy medium reproduces the quadrupole
formula and passes GW170817 without tuning is the strongest evidence the framework
is internally sound.

### Tier 2 — Reproductions of textbook constants (consistency, not measurement)

Values the simulation reproduces that are known constants/relations rather than
direct measurements: sin²θ_W = 3/8 (§8.63), the SM hypercharges (§8.42), the
anomaly-inflow integers (§8.15), m_W/m_Z = cos θ_W (§8.41). Real consistency
checks; each reproduces a known number by a known-or-rediscovered mechanism.

### Tier 3 — The one prediction touching live data

The n = 2 / subluminal / one-cone LV signature, confronted with GW170817,
IceCube TXS 0506+056, and the GZK frontier. Documented in full — including its two
self-corrected errors and its straddle of the GZK bound — in
`EXTERNAL_VALIDATION.md`. This is the only tier that is a *prediction* rather than
a reproduction, and it is half-confirmed (universality) / half-contested
(coefficient, now receded below reach).

## What this is, and what it is not

**It is** a validated, honest, integrative simulation: a single codebase that
ingests known microphysics and returns results consistent with both textbook theory
and real astrophysical measurements, with a fully documented honest-negative record
(seven independent γ = 0 measurements refuting the author's own hope before the
world-crystal escape; multiple retractions; two GZK errors caught and corrected).
Building this cleanly is genuinely uncommon, and the reproductions above make it
trustworthy.

**It is not** new knowledge about nature. Reproducing a real measurement with a
model constructed to contain the physics that generates it is a consistency check —
the measurement already existed, and the sim matched it; it did not predict it. The
distinction between *matched known data* and *predicted unknown data* is the line
between validation and discovery, and everything in Tiers 1–2 is validation. The
project has exactly one item on the discovery side of that line (Tier 3), and it is
not yet decidable.

So the accurate one-sentence statement of the project's standing is:

> **A validated emergent-medium simulation that independently reproduces the
> world-crystal emergent-gravity program (Kleinert–Zaanen), a composite-model
> Weinberg angle (Terazawa–Akama–Chikashige), and un-tuned matches to the GR
> radiation laws and the GW170817 speed bound — with one live, near-Planckian
> Lorentz-violation prediction poised at the edge of testability.**

That is honest, and it is citable — as a *methods / reproduction* contribution, not
a discovery.

---

# Part III — What a genuinely new result would require

The reproductions are the instrument's calibration, not its output. The value of a
validated instrument is the *next* thing it says. To move from "excellent validated
sandbox" to "novel result" needs one of:

1. **A new, near-term-testable prediction** the platform implies that Kleinert,
   Terazawa, Volovik and the QG-phenomenology literature did **not** already
   compute — and that current or near-future experiment can reach. The natural seam
   is the LV sector (the one place the project produces an own number), but §8.61
   pushed it below reach; a *different* observable of the two-sector world crystal
   (e.g. a specific, non-Kleinert consequence of the matter/gravity sector split)
   would be more promising.
2. **A genuine derivation of an input** the project currently assumes — the gauge
   group, the number three, the specific lattice — which the literature also lists
   as open. Deriving one of these (not reframing it) would be new.
3. **A publishable methods contribution**: the simulation platform itself, as a
   reproducible, auditable environment that recovers these known results from one
   substrate, submitted honestly as such (arXiv gr-qc or physics.ed-ph), citing
   Kleinert–Zaanen and Terazawa as targets it reproduces.

Route 3 is available now and is honest. Routes 1–2 are where any actual novelty
lives, and both are hard for the same reason the field finds them hard.

---

# Part IV — The novelty hunt (blind computations, then literature)

After the audit above, we tested whether the validated platform could produce a *new* prediction — one the
priors it reproduces did not already compute. Method (deliberate): compute the model's own answer **blind**
in a focused test file, lock it, *then* search the literature. Three hunts were run.

| Hunt | Blind result | Literature verdict |
|---|---|---|
| **Graviton–photon cone split** (`test_two_sector_dispersion.py`) | The two-sector world crystal splits the cones: photon is strain-elastic (ω²∝S), graviton is curvature-elastic (physical field h=∂u, ω²∝S₂/S), so ζ_g ≠ ζ_m **on the same lattice** — a signed split, graviton ~1.2× faster, Δv/c=+0.049(E/E_P)². Also resolved the §8.67 puzzle: the graviton is luminal in h, the q⁴ is only the strain u. | **Prior art.** Bimetric emergent gravity generically gives different cones for different fields (Volovik ³He-A); emergent photon+graviton LV is published (arXiv:1811.09578, 1709.02736; "Lorentz violation in Goldstone gravity", PRD 80); quartic graviton-dispersion LV is Hořava–Lifshitz + the GW-LV constraint literature; species-dependent LV coefficients are SME (gravitational Čerenkov). Not novel; also Planck-suppressed. |
| **Universal decoherence scale** (`test_medium_decoherence_scale.py`) | A universal, unshieldable Γ ∝ m²·Δx² localization from the medium's own phonon bath (slopes 2.000, 2.000), with a nonzero T→0 vacuum floor; at gravitational coupling strength it lands in the **seconds-to-hours window** current matter-wave/levitated experiments probe — *not* Planck-buried. | **Prior art on mechanism.** Diósi–Penrose (gravitational m² collapse); Gambini–Porto–Pullin "fundamental decoherence from quantum spacetime" (universal, unshieldable — our exact framing); CSL (m²Δx² form). But — unlike every earlier result — it is at a **testable scale**, and the parameter-free DP version is already *falsified* (Gran Sasso 2020–21). |
| **Evade-or-die: does it radiate?** (`test_collapse_radiation.py`) | **EVADE.** The medium's decoherence is dissipationless pure dephasing: the bath energy *saturates* (a one-time reorganization, dE/dt→0), because the medium is a T=0 equilibrium vacuum (§8.13). No continuous heating → no charge acceleration → **no spontaneous X-rays**. The model is *not* excluded by the Gran Sasso test that killed parameter-free DP. | Survival is **firm**; the "dissipationless/non-Markovian collapse" idea has cousins (dissipative CSL, Kafri–Taylor–Milburn), so it is not claimed novel — but the model **survives**, and offers a distinguishing signature. |
| **Residual dissipation of a dynamical particle** (`test_residual_dissipation.py`) | The medium is a superfluid condensate, so the **Landau criterion** applies: single-phonon emission is kinematically forbidden below v_c ≈ 0.72c, so lab matter (v/c ~ 10⁻⁶) is ~6 orders below threshold — the residual dissipation is strongly suppressed and the decoherence proceeds without measurable radiation. **Honest correction:** the Landau criterion is neither necessary nor sufficient — sub-critical channels (vortex nucleation, rotons, thermal normal fraction) leave a small, *uncomputed* residual, so this is **not** a proven exact zero. | Textbook / prior art (Landau 1941; Volovik's superfluid-vacuum emergent gravity; the decoherence-without-dissipation tension is Kafri–Taylor–Milburn + Diósi–Tilloy). Strengthens survival; the exact sub-critical rate is the honest open number. |

**What the hunt actually produced.** Not a novel mechanism — the pattern held: every checkable claim (Kleinert,
Terazawa, bimetric-LV, DP/GPP decoherence) is already in the literature, which is strong evidence the
framework is *sound*. But the decoherence direction produced something the project had lacked since the LV
coefficient receded (§8.61): a **genuine near-term falsification handle at a testable scale** — the model
predicts universal *m²Δx²* decoherence **without** the spontaneous radiation that CSL/DP predict, a signature
experiments can separate (see collapse-like decoherence *without* the X-rays → favours the model; find that
even dissipationless collapse is excluded → the model dies with it). That is the healthy, honest place a
validated toy model can reach: not a discovery, but a way to be tested and possibly killed by real apparatus.

**Honest boundary.** The decoherence coefficient's absolute scale (hence "testable now") assumes gravitational
coupling strength — motivated (the same compression coupling as §8.10) but not rigorously fixed; it is soft,
like the LV ζ. The firm content is the *m²Δx²* law, the unshieldable vacuum floor, and the no-radiation
(dissipationless) character. The remaining open number is the residual dissipation of a fully dynamical
particle (beyond the exact pure-dephasing limit).

## Provenance of this audit

Prior art established by systematic literature search (August 2026) over every
result cluster; the two most load-bearing priors (Kleinert–Zaanen 2004 / Kleinert
2005 for gravity; Terazawa–Akama–Chikashige 1976 for sin²θ_W) were read in full and
matched at the equation level (see the Kleinert concordance above). Key sources:

- H. Kleinert & J. Zaanen, *Phys. Lett. A* **324** (2004) 361.
- H. Kleinert, *Braz. J. Phys.* **35** (2005) 359; *Ann. Physik* **44** (1987) 117.
- H. Terazawa, K. Akama & Y. Chikashige, *Prog. Theor. Phys.* **56** (1976) 1935.
- M. Pretko & L. Radzihovsky, *Phys. Rev. Lett.* **120** (2018) 195301.
- G. E. Volovik, *The Universe in a Helium Droplet* (OUP, 2003).
- A. D. Sakharov, *Dokl. Akad. Nauk SSSR* **177** (1967) 70.
- E. Nelson, *Phys. Rev.* **150** (1966) 1079; A. Valentini, *Phys. Lett. A* **156** (1991) 5.
- W. H. Zurek, *Rev. Mod. Phys.* **75** (2003) 715.
- S. Coleman, *Nucl. Phys. B* **262** (1985) 263.
- S. Chadha & H. B. Nielsen, *Nucl. Phys. B* **217** (1983) 125.
- P. C. Peters & J. Mathews, *Phys. Rev.* **131** (1963) 435.
- B. P. Abbott et al. (LIGO/Virgo), *ApJL* **848** (2017) L13.
- T. Jacobson, S. Liberati & D. Mattingly, review astro-ph/0505267 (2006).
