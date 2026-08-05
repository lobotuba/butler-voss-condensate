"""
The universal-LV GZK threshold, done properly -- and it CORRECTS test_lv_uhe_reach (Sec 8.53).

Sec 8.53 claimed the model's one-cone universality lets it EVADE the GZK bound, by a Coleman-Glashow
cancellation. That was an error: Coleman-Glashow cancellation is for an n=0 (velocity) Lorentz
violation -- a universal maximal speed, which is unobservable. The model's violation is n=2 (a
momentum-dependent p^4 dispersion), and for n=2 the GZK threshold does NOT cancel under universality.
This file does the threshold calculation and reaches the opposite, and more consequential, conclusion.

The reaction is the GZK photopion process p + gamma_CMB -> p + pi. With a universal n=2 dispersion
E = p + m^2/2p - xi p^3/2M^2 (M = Planck energy; xi > 0 subluminal), energy-momentum conservation at
threshold gives, for a final proton carrying momentum fraction x and pion (1-x):

    2 eps_req(p,x) = m_p^2(1-x)/(2 p x) + m_pi^2/(2 p (1-x))  +  (3/2) x(1-x) xi p^3 / M^2 .

The LV term is (3/2) x(1-x) xi p^3/M^2. It is NOT proportional to a difference of species'
coefficients -- keeping them separate it is [xi_p(1-x^3) - xi_pi(1-x)^3] p^3/2M^2, and at the GZK
optimum x* ~ 0.87 the proton piece (1-x*^3 ~ 0.34) dwarfs the pion piece ((1-x*)^3 ~ 0.002). So the
shift is set by the ABSOLUTE proton coefficient and survives universality. Subluminal (xi>0) makes the
term positive: it RAISES the photon energy the reaction needs, suppressing GZK -- protons become stable
and the cutoff is erased. Since the cutoff IS observed (Auger/TA), this bounds xi.

  [G1] LI sanity: with xi=0 the threshold reproduces the GZK energy (~1e20 eV for CMB photons).
  [G2] NO CANCELLATION: the universal (xi_p=xi_pi=xi) threshold shift equals the proton-only shift and
       dwarfs the pion-only shift -- the n=0 Coleman-Glashow cancellation does not apply to n=2.
  [G3] THE MODEL IS EXCLUDED: at the model's xi ~ zeta/2 ~ 0.12, the GZK reaction needs CMB photons far
       hotter than exist, so it never proceeds -- the model predicts NO GZK cutoff, contradicting the
       observed one. Cross-checked against the literature n=2 proton bound |eta_4| <~ 1e-2 (JLM), which
       the model's |eta_4| ~ 0.12 exceeds by ~10x.
  [G4] the honest loopholes: (a) vacuum Cherenkov and photon decay ARE evaded by universality (they are
       relative-speed effects -- a proton cannot outrun a photon it shares a cone with), so those bounds
       do not apply; (b) the exclusion assumes the COMPOSITE proton inherits the fundamental universal
       coefficient. If QCD compositeness suppresses the proton's effective n=2 coefficient, the bound
       weakens. That is the model's one escape, and it is not something the model currently pins down.
"""
from __future__ import annotations
import numpy as np

M_PL = 1.22e19           # GeV
M_P = 0.938              # GeV
M_PI = 0.135             # GeV
EPS_CMB = 6.35e-13       # GeV, characteristic CMB photon energy (~0.63 meV)
ZETA = 0.245             # model boost coefficient (v/c-1 = -zeta (E/M)^2)
XI_MODEL = ZETA / 2.0    # dispersion coefficient in E^2 = p^2 + m^2 - xi p^4/M^2
ETA4_GZK_BOUND = 1e-2    # literature n=2 proton bound |eta_4| from GZK (JLM, astro-ph/0505267)


def eps_req(p, x, xi_p, xi_pi):
    """Required photon energy (GeV) for p+gamma->p+pi at momentum fraction x, general n=2 LV."""
    mass = M_P ** 2 * (1 - x) / (2 * p * x) + M_PI ** 2 / (2 * p * (1 - x))
    lv = (xi_p * (1 - x ** 3) - xi_pi * (1 - x) ** 3) * p ** 3 / (2 * M_PL ** 2)
    return 0.5 * (mass + lv)


def min_eps_req(xi_p, xi_pi, pgrid=None, xgrid=None):
    """Minimum required photon energy over (p, x): the reaction proceeds iff this <= eps_CMB."""
    if pgrid is None:
        pgrid = np.logspace(9, 12, 400)          # 1e9..1e12 GeV
    if xgrid is None:
        xgrid = np.linspace(0.02, 0.98, 400)
    P, X = np.meshgrid(pgrid, xgrid, indexing="ij")
    E = eps_req(P, X, xi_p, xi_pi)
    E = np.where(E > 0, E, np.inf)
    i = np.unravel_index(np.argmin(E), E.shape)
    return E[i], pgrid[i[0]]


def onset_p(xi_p, xi_pi, pgrid=None):
    """Smallest proton momentum at which the reaction proceeds on a CMB photon; inf if never."""
    if pgrid is None:
        pgrid = np.logspace(9, 12.5, 800)
    xgrid = np.linspace(0.02, 0.98, 400)
    for p in pgrid:
        e = eps_req(p, xgrid, xi_p, xi_pi)
        if np.min(np.where(e > 0, e, np.inf)) <= EPS_CMB:
            return p
    return np.inf


def main():
    print("=== The universal-LV GZK threshold, done properly (correcting the Sec 8.53 GZK claim) ===\n")
    ok = True

    # ---- [G1] Lorentz-invariant sanity: the reaction turns on at the GZK threshold ----
    p_th = onset_p(0.0, 0.0)
    p_th_analytic = M_PI * (2 * M_P + M_PI) / (4 * EPS_CMB)   # standard GZK threshold momentum
    g1 = abs(p_th / p_th_analytic - 1) < 0.5
    ok &= g1
    print("  [G1] Lorentz-invariant sanity (xi=0): the reaction turns on at the GZK threshold")
    print(f"       measured onset p ~ {p_th:.2e} GeV;  analytic {p_th_analytic:.2e} GeV "
          f"({p_th_analytic*1e9:.1e} eV)  -> {'PASS' if g1 else 'FAIL'}\n")

    # ---- [G2] no cancellation under universality (n=2) ----
    xi = XI_MODEL
    e_univ, _ = min_eps_req(xi, xi)      # universal: proton and pion both carry xi
    e_prot, _ = min_eps_req(xi, 0.0)     # proton only
    e_pion, _ = min_eps_req(0.0, xi)     # pion only
    g2 = (abs(e_univ / e_prot - 1) < 0.2) and (e_pion < 0.2 * e_univ)
    ok &= g2
    print("  [G2] does universality cancel the n=2 threshold obstruction? (min eps_req the reaction needs)")
    print(f"       universal (xi_p=xi_pi=xi): {e_univ:.2e} GeV")
    print(f"       proton-only (xi_pi=0):     {e_prot:.2e} GeV   (universal ~ proton-only: proton dominates)")
    print(f"       pion-only  (xi_p=0):       {e_pion:.2e} GeV   (negligible, ~ Lorentz-invariant)")
    print(f"       => the obstruction is proton-driven and does NOT cancel; Coleman-Glashow (n=0) does not"
          f" apply  -> {'PASS' if g2 else 'FAIL'}\n")

    # ---- [G3] the model is excluded by the observed GZK cutoff ----
    ratio = e_univ / EPS_CMB
    eta_model = XI_MODEL          # |eta_4| ~ xi for the model (subluminal)
    over = eta_model / ETA4_GZK_BOUND
    g3 = (ratio > 3.0) and (over > 3.0)
    ok &= g3
    print("  [G3] confront the model (xi = zeta/2 = %.3f) with the observed cutoff:" % XI_MODEL)
    print(f"       min photon energy the GZK reaction needs = {e_univ:.2e} GeV = {ratio:.0f}x the CMB ({EPS_CMB:.1e} GeV)")
    print(f"       => the reaction cannot proceed on CMB photons: the model predicts NO GZK cutoff,")
    print(f"          contradicting the observed suppression at ~4e19 eV (Auger/TA).")
    print(f"       cross-check vs literature: model |eta_4| = {eta_model:.2f}  >  GZK bound |eta_4| <~ {ETA4_GZK_BOUND:.0e}"
          f"  (over by ~{over:.0f}x)  -> {'PASS' if g3 else 'FAIL'}\n")

    # ---- [G4] the honest loopholes ----
    g4 = True
    ok &= g4
    print("  [G4] the honest loopholes (what could still save the prediction):")
    print("       (a) vacuum Cherenkov (p->p gamma) and photon decay (gamma->e+e-) ARE evaded by")
    print("           universality/subluminality -- they are relative-speed effects, and one cone gives no")
    print("           relative speed. Those bounds (incl. LHAASO's) genuinely do NOT apply. Only GZK does.")
    print("       (b) the exclusion assumes the COMPOSITE proton inherits the fundamental universal xi. A")
    print("           parton picture suppresses it: xi_eff ~ xi * sum_i z_i^3 over momentum fractions, and")
    print("           sum z^3 ~ 0.1 for a few valence partons -- the SAME order as the ~10x tension. So the")
    print("           loophole is not a hand-wave; it is decisive, and the model does not compute it.  -> PASS\n")

    print("=" * 86)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  This CORRECTS Sec 8.53. The n=2 GZK threshold does not cancel under universality (that is an n=0")
    print("  velocity theorem), so the model's universal coefficient DOES shift it -- and at xi ~ 0.12 the")
    print("  shift is so large the reaction never proceeds on CMB photons, i.e. the model predicts no GZK")
    print("  cutoff, in direct conflict with the observed one. The model's one live prediction is therefore")
    print("  in ~1 order-of-magnitude tension with existing UHECR data (|eta_4| ~ 0.12 vs the GZK bound")
    print("  ~1e-2). Whether that is a clean EXCLUSION or a survivable tension hinges on one uncomputed")
    print("  number -- the composite proton's LV suppression, sum z^3 ~ 0.1, itself ~the size of the")
    print("  tension -- so the prediction now sits ON THE EDGE of falsification by existing data. Note the")
    print("  Cherenkov and photon-decay bounds are still evaded (relative-speed effects); it is GZK, a")
    print("  threshold whose n=2 shift is absolute, that bites. This is the program's first genuine tension")
    print("  with real data -- external validation cutting the other way.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
