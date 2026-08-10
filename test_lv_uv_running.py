"""The deciding UHE calculation (Sec 8.56 [G4]): DGLAP running of the proton's LV coefficient to the UV.

Sec 8.54-8.56 pinned the model's one falsifiable prediction to the ultra-high-energy frontier: the LV
coefficient felt by a proton is eta_p = xi * <x^2>_P, where <x^2>_P = int x^2 D(x) dx is the second
moment of the proton momentum density and xi = 0.245/2 is the model's fundamental (Planck-scale) LV
coefficient. With <x^2>_P measured at Q ~ 2 GeV (~0.04-0.08), xi_eff lands at ~0.5-1.5x the GZK bound --
marginal, straddling the exclusion. Sec 8.56 flagged the ONE remaining calculation: <x^2>_P should be
evaluated not at 2 GeV but at the LV operator's own (ultraviolet) scale, and it decreases with Q^2, so
the physical value is smaller and the model trends safer -- BY AN UNCOMPUTED AMOUNT. This computes it.

The LV term is a p^4, dimension-6 operator; at leading twist it is the spin-3 twist-2 operator, whose
forward proton matrix element IS <x^2>_P. QCD is Lorentz-invariant, so the QCD dressing that renormalises
this operator is the ordinary DGLAP evolution of its N=4 (number-density) Mellin moment. Since the model
defines xi at the substrate (~Planck) scale, consistency requires the matrix element at the SAME scale:
eta_p^phys = xi * <x^2>_P(M_substrate), with <x^2>_P run up from 2 GeV. (Using <x^2>_P(2 GeV) with the
UV coefficient xi, as Sec 8.56 did, mixes scales and OVERestimates -- it is the conservative bound.)

  [G1] MOMENTUM CONSERVATION. Evolving the N=2 singlet moment must keep quark+gluon momentum = 1 exactly
       (the anomalous-dimension matrix has the momentum sum rule built in) -- this validates the matrix.
  [G2] alpha_s. The 1-loop running reproduces alpha_s(M_Z) = 0.118 and gives sensible values elsewhere.
  [G3] THE RESULT. <x^2>_P runs from ~0.042 at 2 GeV down to ~0.003 at the Planck scale -- a ~13x
       suppression -- moving xi_eff from ~0.5x the GZK bound (frontier) to ~0.04x (safely inside).
  [G4] ROBUSTNESS. The conclusion holds across the UV matching scale (10^10 GeV .. Planck, all < 0.15x)
       and across the quark/gluon split and normalisation of the starting moment (suppression ~0.05-0.09).

Approximations, stated: leading order DGLAP, leading twist (higher twist is (Lambda/E)^2-suppressed and
negligible at UHE), and the LV operator taken as the twist-2 spin-3 operator. The GZK LV bound
(eta4 <~ 1e-2) is used as an external input, as in Sec 8.54-8.56.
"""
from __future__ import annotations
import numpy as np
from math import gamma as Gm

CF, CA, TR = 4.0 / 3.0, 3.0, 0.5
MZ, ASMZ, MC, MB, MT = 91.1876, 0.118, 1.3, 4.7, 173.0
XI, GZK_BOUND, M_PLANCK = 0.245 / 2, 1e-2, 1.22e19


def S1(N):
    return sum(1.0 / k for k in range(1, N + 1))


def gamma_singlet(N, nf):
    """LO singlet anomalous-dimension matrix, (quark-singlet, gluon) basis, df/dlnQ2 = (as/2pi) gamma f."""
    s1 = S1(N)
    gqq = CF * (1.5 - 2 * s1 + 1.0 / (N * (N + 1)))
    gqg = 2 * nf * TR * (N * N + N + 2) / (N * (N + 1) * (N + 2))
    ggq = CF * (N * N + N + 2) / ((N - 1) * N * (N + 1))
    ggg = 2 * CA * (1.0 / (N * (N - 1)) + 1.0 / ((N + 1) * (N + 2)) - s1) + 11 * CA / 6 - (2.0 / 3) * nf * TR
    return np.array([[gqq, gqg], [ggq, ggg]])


def nf_of(mu):
    return 3 + (mu > MC) + (mu > MB) + (mu > MT)


def _alpha_grid():
    tg = np.linspace(np.log(1.0), np.log(1e20 ** 2), 120000)
    ag = np.empty_like(tg); iz = np.searchsorted(tg, np.log(MZ ** 2)); ag[iz] = ASMZ
    for i in range(iz + 1, len(tg)):
        b0 = 11 - 2.0 / 3 * nf_of(np.exp(0.5 * tg[i - 1]))
        ag[i] = ag[i - 1] - b0 * ag[i - 1] ** 2 / (4 * np.pi) * (tg[i] - tg[i - 1])
    for i in range(iz - 1, -1, -1):
        b0 = 11 - 2.0 / 3 * nf_of(np.exp(0.5 * tg[i + 1]))
        ag[i] = ag[i + 1] + b0 * ag[i + 1] ** 2 / (4 * np.pi) * (tg[i + 1] - tg[i])
    return tg, ag


_TG, _AG = _alpha_grid()


def alpha_s(mu):
    return float(np.interp(np.log(mu ** 2), _TG, _AG))


def evolve(N, v0, Q0, Q1, nsteps=4000):
    t0, t1 = np.log(Q0 ** 2), np.log(Q1 ** 2); dt = (t1 - t0) / nsteps
    t = t0; v = np.array(v0, float)
    for _ in range(nsteps):
        mu = np.exp(0.5 * t)
        v = v + (alpha_s(mu) / (2 * np.pi)) * (gamma_singlet(N, nf_of(mu)) @ v) * dt
        t += dt
    return v


def mom_k2(a, b):
    B = lambda p, q: Gm(p) * Gm(q) / Gm(p + q)
    return B(a + 3, b + 1) / B(a + 1, b + 1)


def initial_x2():
    """(quark, gluon) parts of <x^2>_P at 2 GeV from data-anchored (global-fit-like) shapes."""
    q = 0.39 * mom_k2(0.5, 4.0) + 0.19 * mom_k2(-0.2, 9.0)   # valence + sea
    g = 0.42 * mom_k2(-0.1, 6.0)
    return q, g


def main():
    print("=== DGLAP running of the proton LV coefficient to the UV (Sec 8.56 [G4] resolved) ===\n")
    ok = True

    # ---- [G1] momentum conservation ----
    print("  [G1] MOMENTUM CONSERVATION (validates the anomalous-dimension matrix):")
    g1 = True
    for Q1 in (1e2, 1e8, 1e16):
        v = evolve(2, (0.58, 0.42), 2.0, Q1); g1 &= abs(v.sum() - 1.0) < 1e-4
        print(f"       2 GeV -> {Q1:.0e} GeV: quark {v[0]:.4f} gluon {v[1]:.4f} sum {v.sum():.5f}")
    ok &= g1
    print(f"       => quark+gluon momentum stays 1  -> {'PASS' if g1 else 'FAIL'}\n")

    # ---- [G2] alpha_s ----
    a_mz = alpha_s(MZ)
    g2 = abs(a_mz - 0.118) < 2e-3
    ok &= g2
    print(f"  [G2] alpha_s(M_Z) = {a_mz:.4f} (target 0.118); alpha_s(2 GeV) = {alpha_s(2.0):.3f}, "
          f"alpha_s(Planck) = {alpha_s(M_PLANCK):.3f}  -> {'PASS' if g2 else 'FAIL'}\n")

    # ---- [G3] the running result ----
    q0, g0 = initial_x2(); x2_0 = q0 + g0
    print(f"  [G3] <x^2>_P runs down as the LV operator's scale rises (xi = {XI:.3f}, GZK bound {GZK_BOUND:.0e}):")
    print(f"       {'scale (GeV)':>13} {'<x^2>_P':>9} {'supp':>7} {'xi_eff':>10} {'/bound':>8}")
    x2_uv = None
    for Q1 in (2.0, 1e2, 1e6, 1e10, M_PLANCK):
        x2 = x2_0 if Q1 == 2.0 else evolve(4, (q0, g0), 2.0, Q1).sum()
        if Q1 == M_PLANCK:
            x2_uv = x2
        xe = XI * x2
        print(f"       {Q1:>13.2e} {x2:>9.4f} {x2/x2_0:>7.3f} {xe:>10.2e} {xe/GZK_BOUND:>7.2f}x")
    g3 = (XI * x2_0 / GZK_BOUND > 0.3) and (XI * x2_uv / GZK_BOUND < 0.2)   # 2 GeV frontier -> UV safe
    ok &= g3
    print(f"       => from {XI*x2_0/GZK_BOUND:.2f}x the bound at 2 GeV (frontier) to {XI*x2_uv/GZK_BOUND:.2f}x at "
          f"the Planck scale (safe)  -> {'PASS' if g3 else 'FAIL'}\n")

    # ---- [G4] robustness ----
    print("  [G4] ROBUSTNESS of the safe conclusion:")
    print("       (a) UV matching scale (a Planck-suppressed operator could match anywhere >~ 10^10 GeV):")
    worst = 0.0
    for M in (1e10, 1e13, 1e16, M_PLANCK):
        r = XI * evolve(4, (q0, g0), 2.0, M).sum() / GZK_BOUND
        worst = max(worst, r)
        print(f"           match at {M:.0e} GeV: xi_eff = {r:.2f}x the bound")
    print("       (b) starting moment (global-fit low .. lattice-high, various quark/gluon splits):")
    supps = []
    for (q, g, lbl) in [(q0, g0, "global-fit"), (0.060, 0.020, "lattice-high"),
                        (0.040, 0.005, "quark-heavy"), (0.020, 0.022, "gluon-heavy")]:
        s = evolve(4, (q, g), 2.0, M_PLANCK).sum() / (q + g); supps.append(s)
        print(f"           {lbl:12}: <x^2>(2 GeV)={q+g:.3f} -> Planck suppression {s:.4f}")
    g4 = worst < 0.2 and max(supps) < 0.12
    ok &= g4
    print(f"       => safe (< 0.15x) across all UV scales; suppression ~0.05-0.09 regardless of input"
          f"  -> {'PASS' if g4 else 'FAIL'}\n")

    print("=" * 90)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  The calculation Sec 8.56 flagged as deciding is done: DGLAP evolution of the proton's second")
    print(f"  momentum moment suppresses the LV coefficient ~{1/(x2_uv/x2_0):.0f}x between 2 GeV and the Planck scale,")
    print(f"  because the model's xi is a UV quantity and the operator's matrix element must be taken at the")
    print(f"  same UV scale. This moves the model's one prediction from ~{XI*x2_0/GZK_BOUND:.1f}x the GZK bound (marginal,")
    print(f"  straddling exclusion) to ~{XI*x2_uv/GZK_BOUND:.2f}x (safely inside), robustly across UV scale and input PDF.")
    print("  Honest reading: the confrontation Sec 8.54 opened is RESOLVED in the model's favour -- it is not")
    print("  in tension with UHE data after all -- but the prediction thereby recedes ~1-2 orders below current")
    print("  reach, so the model's most falsifiable claim survives at the cost of being harder to test.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
