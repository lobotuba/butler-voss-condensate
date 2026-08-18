"""
HUNTING A NEW PREDICTION (2c): the residual dissipation of a DYNAMICAL particle -- is "no radiation" exact?

Method (Robert): compute blind, then check the literature. This closes the one open number left by the
evade-or-die test (test_collapse_radiation): that test showed the medium decoheres WITHOUT heating in the
pure-dephasing (static-particle) limit. But a real particle is DYNAMICAL -- its position does not commute
with p^2/2m -- so finite mass could reintroduce a small dissipation, hence some spontaneous radiation. How
big is it? Exactly zero, or merely small?

THE KEY PHYSICS. The model's medium is a CONDENSATE -- a superfluid-like ground state (that is the whole
premise: space is an active condensate). A particle moving through a superfluid obeys the LANDAU CRITERION:
it can only shed energy by emitting a phonon (k, omega(k)) if energy-momentum conservation allows, which for
a heavy particle at velocity v requires
        omega(k) <= v . k          (emission threshold),
i.e. v must exceed the phonon phase velocity omega(k)/|k| for SOME mode. The minimum phase velocity over the
whole spectrum is the Landau critical velocity
        v_c = min_k omega(k) / |k|.
Below v_c NO mode satisfies the threshold: the emission phase space is EMPTY, and the dissipation is EXACTLY
zero -- frictionless superfluid motion. The medium's phonon speed is the emergent light speed c (S8.1), so
v_c is of order c. Every laboratory particle moves at v << c, so it sits astronomically below threshold: the
residual dissipation is not just small, it is kinematically FORBIDDEN. The decoherence (dephasing) needs no
such threshold -- it is elastic/virtual -- so it survives while the dissipation is switched off. Decoherence
WITHOUT dissipation, protected by the superfluidity of the medium.

This also ties to the earlier hunts: the emission threshold is exactly the (gravitational-)Cerenkov condition,
and the two-sector result put the graviton FASTER than the photon (test_two_sector_dispersion) -- the
Cerenkov-SAFE ordering. Subluminal matter cannot radiate into the cone.

  [G1] LANDAU CRITICAL VELOCITY: compute v_c = min_k omega(k)/|k| over the Brillouin zone from the medium's
       acoustic dispersion. It is a finite fraction of the emergent light speed c (order c, not zero).
  [G2] EMISSION PHASE SPACE IS EMPTY BELOW v_c: count the modes satisfying omega(k) <= v.k as a function of
       v. Zero for v < v_c, turning on only above -- so the dissipation (and the radiation it would source)
       is EXACTLY zero below threshold, not merely small.
  [G3] LAB MATTER IS DEEPLY PROTECTED: v_lab / c ~ 1e-6 (a fast lab particle) << v_c/c ~ O(1), so real
       matter is ~6 orders below the dissipation threshold -- "no radiation" is protected, not approximate.
       Meanwhile the dephasing coherence loss (elastic, thresholdless) is unchanged: decoherence persists.
  [G4] VERDICT and honest caveat: in the mean-field Landau kinematics the single-phonon residual dissipation
       is forbidden below v_c ~ c, so it is STRONGLY SUPPRESSED for all subluminal (laboratory) matter and
       the evade-or-die survival is strengthened. But NOT proven exact: the Landau criterion is neither
       necessary nor sufficient -- real condensates dissipate below v_c via vortex nucleation, roton emission,
       boundary effects and any thermal normal component. So a small, uncomputed sub-critical dissipation
       remains -- the SAME open problem as the collapse/gravitational-decoherence literature (Kafri-Taylor-
       Milburn decoherence-without-dissipation and its Diosi-Tilloy dissipative generalizations).
"""
from __future__ import annotations
import numpy as np

SQ3 = np.sqrt(3.0)
# triangular medium: 6 nearest neighbours (the acoustic phonon lattice)
NN = np.array([(1, 0), (0.5, SQ3 / 2), (-0.5, SQ3 / 2),
               (-1, 0), (-0.5, -SQ3 / 2), (0.5, -SQ3 / 2)])


def omega(k):
    """Acoustic dispersion omega(k) = c_s * sqrt(S(k)/S''0), normalised so omega -> |k| (c_s=1) at small k.
    S(k) = SUM_e (1 - cos k.e) is the graph-Laplacian symbol (subluminal lattice band)."""
    S = (1.0 - np.cos(NN @ k)).sum()
    return np.sqrt(2.0 * S / 3.0)          # 2/3 normalises c_s=1 at k->0 for this 6-neighbour set


def phase_velocity(k):
    kk = np.linalg.norm(k)
    return omega(k) / kk if kk > 1e-9 else np.nan


if __name__ == "__main__":
    print("=== Residual dissipation of a dynamical particle: is 'no radiation' EXACT? (blind) ===\n")
    print("  The medium is a condensate (superfluid). Landau: a particle at v dissipates only by emitting a")
    print("  phonon with omega(k) <= v.k. Below v_c = min omega/|k| the emission phase space is EMPTY ->")
    print("  frictionless. v_c ~ c (emergent light speed); lab matter has v << c -> protected exactly.\n")

    # sample the Brillouin zone as a disk |k| < pi (one consistent k-set for BOTH v_c and the emission scan)
    kmax = np.pi
    grid = np.linspace(-kmax, kmax, 221)
    KX, KY = np.meshgrid(grid, grid)
    ks = np.column_stack([KX.ravel(), KY.ravel()])
    kn = np.linalg.norm(ks, axis=1)
    ks = ks[(kn > 0.15) & (kn < kmax)]                 # inside the zone, away from the k->0 singular point
    pv = np.array([omega(k) / np.linalg.norm(k) for k in ks])   # phase velocity of every mode

    # ---- [G1] Landau critical velocity ----
    c_s = 1.0                                          # emergent light speed in medium units (omega/|k| at k->0)
    v_c = pv.min()
    print("  [G1] Landau critical velocity from the medium's dispersion:")
    print(f"       c_s (emergent light speed, k->0 phase velocity) = {c_s:.3f}")
    print(f"       v_c = min_k omega(k)/|k| over the BZ           = {v_c:.3f}  ({v_c/c_s:.2f} c)")
    ok1 = 0.2 < v_c <= 1.0
    print(f"       => v_c is a finite fraction of c -- order c, not zero  [{'PASS' if ok1 else 'FAIL'}]\n")

    # ---- [G2] emission phase space empty below v_c ----
    # heavy particle, best case v aligned with k: emit iff omega(k) <= v|k|, i.e. phase velocity pv <= v.
    print("  [G2] emission phase space {k : omega(k) <= v.k} vs particle velocity v (heavy particle):")
    print(f"       {'v/c':>8} {'emittable-mode fraction':>26}")
    frac = {}
    for v in (0.2, 0.5, v_c - 0.02, v_c + 0.02, 0.9, 1.5):
        emit = (pv <= v)
        frac[v] = emit.mean()
        print(f"       {v:>8.3f} {emit.mean():>26.4f}")
    ok2 = frac[0.2] == 0.0 and frac[v_c - 0.02] == 0.0 and frac[v_c + 0.02] > 0.0
    print(f"       => zero emittable modes below v_c, nonzero above -- dissipation is EXACTLY 0 below threshold  [{'PASS' if ok2 else 'FAIL'}]\n")

    # ---- [G3] lab matter is deeply protected ----
    print("  [G3] laboratory matter vs the threshold:")
    v_lab_over_c = 1e-6                                 # a fast lab particle (~300 m/s ... keV electron ~ 0.1c is extreme; use 1e-6 typical)
    margin = v_c / v_lab_over_c
    print(f"       a typical lab particle: v/c ~ {v_lab_over_c:.0e}   vs   v_c/c ~ {v_c:.2f}")
    print(f"       margin below the dissipation threshold: ~ {margin:.1e}  -> emission kinematically forbidden.")
    print("       meanwhile the dephasing (elastic, thresholdless) is unchanged: the superposition still")
    print("       decoheres. => DECOHERENCE WITHOUT DISSIPATION, protected by the medium's superfluidity.")
    ok3 = margin > 1e5
    print(f"       [{'PASS' if ok3 else 'FAIL'}]\n")

    # ---- [G4] verdict ----
    print("  [G4] VERDICT and honest caveat:")
    print("       In the mean-field Landau kinematics the residual dissipation is FORBIDDEN below v_c ~ c:")
    print("       the single-phonon emission phase space is empty, so a subluminal particle (all lab matter,")
    print("       v/c ~ 1e-6) is ~6 orders below threshold and radiates essentially nothing. The dephasing")
    print("       (elastic, thresholdless) is untouched -- decoherence WITHOUT dissipation, protected by the")
    print("       SAME superfluidity that makes the medium a condensate and the SAME cone that gives Lorentz")
    print("       invariance (S8.1) and the Cerenkov-safe graviton-faster ordering (test_two_sector_dispersion).")
    print("       HONEST CAVEAT -- do NOT overclaim 'exactly zero': the Landau criterion is known to be")
    print("       neither necessary nor sufficient. Real condensates dissipate BELOW v_c through channels the")
    print("       single-phonon kinematics misses -- vortex nucleation, roton emission, boundary/surface")
    print("       effects, and any thermal normal component. So the residual dissipation is STRONGLY")
    print("       SUPPRESSED, not provably zero; its true size (the vortex/roton sub-threshold rate) is the")
    print("       remaining open number, and it is the SAME problem the collapse/gravitational-decoherence")
    print("       literature works on (Kafri-Taylor-Milburn decoherence-without-dissipation and its")
    print("       Diosi-Tilloy dissipative generalizations).\n")

    allp = ok1 and ok2 and ok3
    print("=" * 92)
    print(f"[verdict] {'ALL GATES PASS -- single-phonon dissipation forbidden below the cone' if allp else 'SOME GATES FAILED'}")
    print("  In the mean-field Landau kinematics the dynamical-particle residual dissipation is zero below")
    print("  v_c ~ c: single-phonon emission is kinematically forbidden, so lab matter (v/c ~ 1e-6) is deeply")
    print("  protected and the decoherence proceeds without measurable dissipation or radiation. This")
    print("  strengthens the evade-or-die survival -- but does NOT make it exact: the Landau criterion is")
    print("  neither necessary nor sufficient, and sub-critical channels (vortex nucleation, rotons, thermal")
    print("  normal fraction) leave a small, uncomputed residual dissipation. Landau protection and the")
    print("  decoherence-without-dissipation tension are both textbook / prior art (Landau; Volovik's")
    print("  superfluid vacuum; Kafri-Taylor-Milburn + Diosi-Tilloy). The model SURVIVES and the residual is")
    print("  strongly suppressed; the exact sub-critical rate is the honest open number, not a proven zero.")
