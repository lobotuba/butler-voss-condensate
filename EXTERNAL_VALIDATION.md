# External validation: the one prediction against real data (§8.52–§8.56)

A companion to `EPISTEMIC_AUDIT.md`. The audit's stated ceiling was **external
validation** — every result in the report was internal consistency or the
reproduction of known physics, with no contact with an actual measurement. This
document records the arc that changed that: the confrontation of the model's one
genuinely falsifiable prediction with real astrophysical data, including the two
errors made and corrected along the way. It is deliberately honest about where
the model now stands, which is *on the edge*.

Sections `test_*.py` and whitepaper §8.52–§8.56; committed WP-42 → WP-46.

## The one prediction

After the gravity arc closed (γ = 0 in every channel) and the "scale-dependent γ"
gravitational prediction was retracted (§8.11/§8.16, WP-37), the model was left
with a **single live falsifiable prediction**: a Lorentz violation that is

- **quadratic** (n = 2, dimension-six, ∝ (E/E_Planck)²) — not the linear (n = 1)
  form of many quantum-gravity scenarios;
- **subluminal** (velocities fall below c at high E);
- **species-universal** — one emergent cone shared by photon, fermion, graviton
  (§8.4–§8.5); no leading-order species-dependent maximal speed;
- of coefficient **ζ ≈ 0.245**, i.e. an effective scale E_QG,2 ≈ 2 × the Planck
  energy, shown lattice-independent in §8.51 (so the falsifiable content is
  structural, not a tuned number).

Two halves are testable separately: the **universality** (do the species share a
speed?) and the **coefficient** (is ζ the size claimed?).

## The confrontations, in order

| § | Test | What it confronts | Outcome |
|---|------|-------------------|---------|
| 8.52 | `test_gw170817_onecone` | universality, at MeV | **passes** — the GW170817/GRB 170817A 1.74 s timing bounds \|Δv\|/c to [−3, +0.7]×10⁻¹⁵; the model predicts exactly 0 |
| 8.53 | `test_lv_uhe_reach` | the coefficient, UHE frontier | LHAASO (photon-decay ⇒ superluminal-only ⇒ N/A) and IceCube TXS 0506+056 (a **second** universality pass, at 10⁵ GeV) stand; **the GZK gate was wrong** |
| 8.54 | `test_lv_gzk_threshold` | the coefficient, via GZK | **corrects 8.53**: the n=2 GZK threshold does *not* cancel under universality; the bare coefficient is in ~12× tension with the observed cutoff |
| 8.55 | `test_lv_proton_compositeness` | the composite-proton suppression | the deciding factor: ξ_eff = ξ·⟨x²⟩_P ≈ 0.8× the GZK bound (toy PDFs) — marginal |
| 8.56 | `test_lv_proton_moment_data` | the same, with real error bars | anchored to measured moments, ξ_eff ≈ **0.5–1.5×** the bound — **straddles** the boundary |

## Two errors, caught and corrected

The arc is worth reading as much for its self-correction as its result.

1. **The GZK cancellation (§8.53 → §8.54).** §8.53 claimed one-cone universality
   lets the model *evade* the GZK bound, by a Coleman–Glashow cancellation. That
   is an **n = 0** (velocity) theorem — a universal maximal speed is
   unobservable. The model's violation is **n = 2**, and a universal n = 2
   coefficient does *not* cancel in thresholds: the shift is (3/2)ξ·x(1−x)·p³/M²,
   dominated by the *absolute* proton coefficient (at the GZK optimum x ≈ 0.87 the
   proton piece ≈ 0.34 dwarfs the pion ≈ 0.002). So the coefficient is not evaded
   but constrained. What *is* correctly evaded — because they are genuine
   relative-speed effects — are vacuum Čerenkov and photon decay (LHAASO's
   bounds), which one cone disarms.

2. **The moment order (§8.56).** The measured isovector moment ⟨x²⟩_{u-d} =
   ∫x²(u−d)dx equals w_iso · ⟨x⟩ of the momentum density (a k=1 moment), not k=2;
   the suppression factor Σz³ = ⟨x²⟩_P is the k=2 moment. Getting this right is
   what let the toy shapes be calibrated to the real second moment.

## Where the model stands

The chain: the fundamental coefficient ξ ≈ ζ/2 ≈ 0.12 is 1 order of magnitude
above the dimension-six GZK proton bound |η₄| ≲ 10⁻² (Jacobson–Liberati–Mattingly,
astro-ph/0505267). Proton compositeness suppresses it by Σz³ = ⟨x²⟩_P — the mean
of x² under the proton's momentum density — because the momentum is spread over
many low-x partons and x² weights them down. Anchored to the world momentum
fractions and the physical-point lattice second moment ⟨x²⟩_{u-d} = 0.083(14)
(arXiv:2605.02808, vs global-fit ~0.055), this gives

> **ξ_eff ≈ 0.5–1.5 × the GZK bound** — just inside on global-fit moments, just
> outside on the lattice moment. The model **straddles the GZK exclusion boundary.**

This is the program's **first genuine contact between its coefficient and real
data, and its first genuine tension** — external validation cutting against the
model rather than for it. It survives every test it has faced (universality
twice, the GZK threshold marginally) and is pinned to the edge of current
ultra-high-energy sensitivity: the most falsifiable, and most exposed, position a
prediction can occupy.

## What decides it (the open items)

Three concrete, near-term numbers:

1. **The lattice-vs-global-fit moment tension.** ⟨x²⟩_{u-d} = 0.083(14) (lattice)
   vs ~0.055 (fits) is a ~2× spread the answer inherits directly. Resolving it
   moves ξ_eff within [0.5×, 1.5×].
2. **A firm GZK Lorentz-violation bound.** The ~10⁻² proton bound is
   order-of-magnitude; a spectrum-level fit to current Auger/TA data would sharpen
   the boundary the model is being compared to.
3. **The LV operator's ultraviolet scale.** Σz³ = ∫x³D(x,Q²)dx *decreases* with
   Q², and the fundamental operator's natural scale is the deep UV, so the ~2 GeV
   reference value is an **upper** estimate — this one trends the model *safer*,
   and quantifying it is the operator's anomalous dimension (a real QCD
   calculation, not done here).

Beyond these, the next round of ultra-high-energy cosmic-ray data is the natural
experimental tightening. None of this requires new instruments — the model's one
prediction is decidable with existing physics and a few concrete calculations,
which is exactly what a healthy falsifiable claim should be.

## Honest status of external validation

From "no contact with data" at the start of the arc to a single prediction poised
on a real experimental frontier with a concrete way to decide it. That is a real
gain, and it is bounded: the model's reach beyond internal consistency and
known-physics reproduction remains this **one** prediction, now half-confirmed
(universality) and half-contested (coefficient). The report claims nothing more.
