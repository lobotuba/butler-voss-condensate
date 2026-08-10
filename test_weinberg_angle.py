"""An un-tuned constant of nature: sin^2 theta_W = 3/8 at the substrate scale, from induced couplings.

The epistemic audit's stated ceiling is that no measured dimensionless CONSTANT OF NATURE falls out of the
model unforced -- its reality-contact is thin. This supplies one, by combining two ingredients the model
already has and had not put together.

  (1) The gauge couplings are INDUCED, not fundamental (Sakharov / Volovik, S8.14, S8.5): integrating out
      the fermions gives 1/g_a^2 = (1/12 pi^2) Tr(T_a^2) ln(Lambda^2/mu^2) + ..., so the loop factor and the
      log are COMMON to all three gauge groups and the ratio of couplings is fixed purely by the fermion
      content:  1/g_a^2 : 1/g_b^2 = Tr(T_a^2) : Tr(T_b^2).
  (2) The hypercharges are DERIVED, not input (anomaly cancellation on the exact lattice gauge symmetry,
      S8.42): Y = 1/6, -2/3, 1/3, -1/2, 1 for Q, u^c, d^c, L, e^c.

For one Standard-Model generation with those hypercharges, Tr(T_a^2) is EQUAL for all three groups in GUT
normalisation (SU(3): 2, SU(2): 2, U(1): (3/5)Tr(Y^2) = 2). So the induced couplings UNIFY -- g1 = g2 = g3
at the substrate scale -- and

        sin^2 theta_W = Tr(T2^2) / (Tr(T2^2) + Tr(Y^2)) = 2 / (2 + 10/3) = 3/8 = 0.375,

with NO tuning, and WITHOUT assuming any grand-unified gauge group. This is the SU(5) value, reached by a
different mechanism (induced couplings, not a unifying group). Confronted with data by one-loop running, it
lands where the measured couplings actually go -- with the same virtues and the same non-SUSY blemishes as
the textbook minimal unification.

  [G1] Tr(T_a^2) is equal for the three groups (GUT norm) -- the induced couplings unify. Un-tuned.
  [G2] hence sin^2 theta_W(substrate) = 3/8, an un-tuned dimensionless output.
  [G3] the measured couplings run UP to sin^2 theta_W = 3/8 at ~10^13 GeV, where the three approximately
       meet (spread ~13%, the known non-SUSY near-miss); this validates the M_Z inputs used.
  [G4] confronted the other way -- predict sin^2 theta_W(M_Z) from unification + the measured alpha_em,
       alpha_s -- it comes out ~0.20 vs the measured 0.231: right to ~10-12%, the standard non-SUSY gap,
       and the unification scale ~10^13 GeV sits a few orders below the Planck substrate.

Honest scope: 3/8 is the standard GUT-normalisation value; the model's contribution is that it is UN-TUNED
here -- it follows from the model's own derived hypercharges and its induced-coupling mechanism, with no
GUT group assumed -- and it inherits the standard non-SUSY imperfections (the ~12% M_Z miss, and a
unification scale below Planck). One-loop; a genuine constant of nature that the model does not tune.
"""
from __future__ import annotations
import numpy as np

SM = [("Q", 3, 2, 1/6), ("uc", 3, 1, -2/3), ("dc", 3, 1, 1/3), ("L", 1, 2, -1/2), ("ec", 1, 1, 1)]
T_FUND = 0.5
MZ = 91.1876
b = np.array([41/10, -19/6, -7.0])                       # one-loop SM (b1 GUT-normalised)


def traces():
    """Tr(T_a^2) over one SM generation: (color, weak, hypercharge-physical)."""
    T3 = sum(w * (T_FUND if c == 3 else 0.0) for _, c, w, Y in SM)
    T2 = sum(c * (T_FUND if w == 2 else 0.0) for _, c, w, Y in SM)
    TY = sum(c * w * Y * Y for _, c, w, Y in SM)
    return T3, T2, TY


def run(ainv_MZ, mu):
    return ainv_MZ - b / (2 * np.pi) * np.log(mu / MZ)


def main():
    print("=== An un-tuned constant of nature: sin^2 theta_W = 3/8 (induced couplings + anomaly Y) ===\n")
    ok = True
    T3, T2, TY = traces()
    TYg = (3/5) * TY

    # [G1] induced couplings unify: Tr(T_a^2) equal in GUT normalisation
    print("  [G1] Tr(T_a^2) over one SM generation (induced-coupling ratio 1/g_a^2 ~ Tr(T_a^2)):")
    print(f"       SU(3) = {T3:.3f}   SU(2) = {T2:.3f}   U(1)_Y = {TY:.4f} -> GUT-norm (3/5)Tr(Y^2) = {TYg:.3f}")
    g1 = abs(T3 - 2) < 1e-9 and abs(T2 - 2) < 1e-9 and abs(TYg - 2) < 1e-9
    ok &= g1
    print(f"       => all three equal (2,2,2): the induced couplings UNIFY, g1=g2=g3  -> {'PASS' if g1 else 'FAIL'}\n")

    # [G2] sin^2 theta_W = 3/8
    s2w = T2 / (T2 + TY)
    g2 = abs(s2w - 3/8) < 1e-9
    ok &= g2
    print(f"  [G2] sin^2 theta_W(substrate) = Tr(T2^2)/(Tr(T2^2)+Tr(Y^2)) = {T2:.3f}/{T2+TY:.3f} = {s2w:.4f}")
    print(f"       = 3/8, exactly, with no tuning  -> {'PASS' if g2 else 'FAIL'}\n")

    # measured couplings at M_Z (PDG), GUT-normalised alpha_1
    aem_inv, s2w_mz, as_mz = 127.951, 0.23122, 0.1179
    a2 = s2w_mz * aem_inv
    a1 = (3/5) * (1 - s2w_mz) * aem_inv
    a3 = 1.0 / as_mz
    ainv0 = np.array([a1, a2, a3])
    s2w_check = a2 / (a2 + (5/3) * a1)
    print(f"  [G2 check] measured M_Z couplings reproduce sin^2 theta_W(M_Z) = {s2w_check:.4f} (PDG 0.2312)\n")

    # [G3] run up: sin^2 theta_W -> 3/8 near 10^13 GeV; three couplings approximately meet
    print("  [G3] running the measured couplings UP (one-loop): sin^2 theta_W climbs to 3/8:")
    print(f"       {'mu (GeV)':>9} {'a1^-1':>7} {'a2^-1':>7} {'a3^-1':>7} {'sin^2 thW':>10}")
    for mu in (1e2, 1e8, 1e13, 1e16, 1.2e19):
        r = run(ainv0, mu)
        print(f"       {mu:>9.1e} {r[0]:>7.2f} {r[1]:>7.2f} {r[2]:>7.2f} {r[1]/(r[1]+5/3*r[0]):>10.4f}")
    t12 = (a1 - a2) / ((b[0] - b[1]) / (2 * np.pi))
    mu12 = MZ * np.exp(t12)
    rat = run(ainv0, mu12)
    spread = abs(rat[2] - rat[0]) / rat[0]
    g3 = 1e12 < mu12 < 1e14 and spread < 0.2
    ok &= g3
    print(f"       => sin^2 theta_W = 3/8 exactly at mu = {mu12:.1e} GeV (where a1=a2); the three couplings")
    print(f"          meet there to {spread:.0%} (a3^-1={rat[2]:.1f} vs a1,2^-1={rat[0]:.1f}) -- the non-SUSY")
    print(f"          near-miss  -> {'PASS' if g3 else 'FAIL'}\n")

    # [G4] the back-prediction of sin^2 theta_W(M_Z) from minimal unification
    s2w_pred = 1/6 + (5/9) * (1/aem_inv) / as_mz          # standard one-loop minimal-unification formula
    g4 = abs(s2w_pred - s2w_mz) / s2w_mz < 0.15
    ok &= g4
    print(f"  [G4] the reverse test -- predict sin^2 theta_W(M_Z) from unification + measured alpha_em, alpha_s:")
    print(f"       sin^2 theta_W(M_Z) = 1/6 + (5/9) alpha_em/alpha_s = {s2w_pred:.4f}  vs measured {s2w_mz:.4f}")
    print(f"       => right to {abs(s2w_pred-s2w_mz)/s2w_mz:.0%} -- the textbook non-SUSY gap; the model shares it"
          f"  -> {'PASS' if g4 else 'FAIL'}\n")

    print("=" * 88)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  sin^2 theta_W = 3/8 is a genuine UN-TUNED output of the model: the gauge couplings are induced")
    print("  (1/g_a^2 ~ Tr(T_a^2), a common loop), and for one SM generation with the model's anomaly-derived")
    print("  hypercharges the three Tr(T_a^2) are equal in GUT normalisation, so the couplings unify and")
    print("  sin^2 theta_W = 3/8 -- with no tuning and no GUT group assumed. It is the SU(5) value reached by")
    print("  a different mechanism. Run to data it lands where the measured couplings go: 3/8 at a unification")
    print("  scale ~10^13 GeV, the three couplings meeting to ~13%, and a back-predicted sin^2 theta_W(M_Z) ~")
    print("  0.20 vs the measured 0.231 -- right to ~10-12%, with the standard non-SUSY blemishes (the M_Z")
    print("  miss, a scale below Planck). This is the first constant of nature the model produces unforced --")
    print("  a real, if imperfect, dent in the audit's reality-contact ceiling.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
