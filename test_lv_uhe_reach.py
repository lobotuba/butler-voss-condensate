"""
Working through the ultra-high-energy datasets: can they actually test the model's LV coefficient?

test_gw170817_onecone confronted the UNIVERSALITY half of the prediction (one cone) with real data and
passed. The COEFFICIENT half -- the n=2 dispersion v(E)/c = 1 - zeta (E/E_Planck)^2, zeta ~ 0.245 -- needs
the ultra-high-energy frontier, where a quadratic effect finally grows. This file works through the three
frontier datasets and asks, honestly, whether any of them reaches the coefficient. The answer is a genuine
and slightly sobering result, and it tempers Section 8.39.

Real frontier data (published, order-of-magnitude where a spectrum is involved):
  * Pierre Auger / Telescope Array: UHECR up to ~1e20 eV (1e11 GeV); the GZK flux suppression at ~4e19 eV.
  * LHAASO (Cao et al. 2021): gamma rays to 1.4 PeV (1.4e6 GeV) -- the highest photons ever seen.
  * IceCube: the 6.3 PeV Glashow-resonance antineutrino (Dec 2016); and the 290 TeV neutrino IceCube-170922A
    coincident with a gamma-ray flare of blazar TXS 0506+056 at z = 0.34 (a nu-gamma timing coincidence).

The physics that decides it. Two properties the model *shares with its GW170817 pass* also make it EVADE the
strongest UHE Lorentz-violation bounds:
  (i)  SUBLUMINAL (zeta > 0): the strongest clean UHE-photon bounds are photon DECAY (gamma -> e+e-) and
       vacuum Cherenkov, which are kinematically allowed only for SUPERLUMINAL photons. A subluminal photon
       does not decay, so LHAASO's headline LV limits -- which are superluminal-photon-decay limits -- do not
       constrain this model at all.
  (ii) ONE CONE (universal): threshold reactions (the GZK photopion threshold; pair production on background
       light) are anomalous, in the Coleman-Glashow / Jacobson-Liberati-Mattingly analysis, only through the
       DIFFERENCE of species' Lorentz-violating coefficients. For a single universal cone that difference is
       zero, so the leading threshold shift cancels and the species-dependent GZK bound (E_QG,2 >~ 1e18 GeV)
       does not apply as stated.

  [G1] AUGER / GZK: the model's dispersion at 1e11 GeV is tiny, and being one-cone it evades the
       species-dependent GZK threshold bound -- consistent, but for a reason that also disarms the test.
  [G2] LHAASO / PeV gammas: the model is subluminal, so the photon-decay bound that PeV photons impose on
       SUPERLUMINAL LV does not apply; its own dispersion at 1.4 PeV is ~1e-27 -- consistent, untestable here.
  [G3] ICECUBE: the 290 TeV nu-gamma coincidence from TXS 0506+056 is a HIGH-ENERGY one-cone test (like
       GW170817 but at 1e5 GeV); the model predicts |dv/c| ~ 1e-28, far inside the ~1e-12 timing bound --
       another universality pass. The 6.3 PeV neutrino is stable (subluminal) and its dispersion is negligible.
  [G4] THE HONEST REACH, and a temper of Section 8.39: because the model is subluminal AND one-cone, the
       strongest UHE bounds do not apply, so the coefficient is currently UNCONSTRAINED by them -- the model
       is safer, and correspondingly HARDER to test, than 8.39's "~1.4 orders below the frontier, reachable"
       snapshot implied. The genuine handle is multimessenger universality (passing at MeV and at 1e5 GeV),
       not the coefficient.
"""
from __future__ import annotations

C = 2.998e8
GLY = 9.461e24            # m (light-year x 1e9)
GPC = 3.0857e25          # m
E_PLANCK = 1.22e19       # GeV
ZETA = 0.245

# --- published frontier energies (GeV) ---
E_UHECR = 1.0e11         # ~1e20 eV, Auger/TA highest cosmic rays
E_GZK = 4.0e10           # ~4e19 eV, GZK suppression onset
E_LHAASO = 1.4e6         # 1.4 PeV, LHAASO highest photon (Cao et al. 2021)
E_ICECUBE_GLASHOW = 6.3e6  # 6.3 PeV, Glashow-resonance antineutrino
E_TXS_NU = 2.9e5         # 290 TeV, IceCube-170922A neutrino
E_TXS_GAMMA = 1.0e2      # ~100 GeV, the coincident Fermi/MAGIC gamma flare
TXS_Z = 0.34             # redshift of TXS 0506+056
TXS_WINDOW = 10 * 86400.0  # s, ~10-day coincidence window (conservative)

# --- species-dependent literature bound, for contrast ---
E_QG2_SPECIES = 1e18     # GeV, representative n=2 UHE bound assuming species-dependent LV
E_QG2_MODEL = E_PLANCK / ZETA ** 0.5


def dv(E_gev):
    """v(E)/c - 1 for the model's subluminal n=2 dispersion."""
    return -ZETA * (E_gev / E_PLANCK) ** 2


def main():
    print("=== Working through the UHE datasets: do any reach the model's LV coefficient? ===\n")
    print(f"  model: subluminal, one-cone, n=2;  E_QG,2 = {E_QG2_MODEL:.2e} GeV  (zeta = {ZETA})\n")
    ok = True

    # ---- [G1] Auger / GZK ----
    d_uhecr = dv(E_UHECR)
    naive_margin = E_QG2_MODEL / E_QG2_SPECIES
    g1 = (E_QG2_MODEL > E_QG2_SPECIES) and (abs(d_uhecr) < 1e-10)
    ok &= g1
    print("  [G1] Pierre Auger / GZK (UHECR to ~1e11 GeV):")
    print(f"       model dispersion at {E_UHECR:.0e} GeV: v/c-1 = {d_uhecr:.1e}")
    print(f"       species-dependent GZK bound ~ E_QG,2 > {E_QG2_SPECIES:.0e} GeV; model {E_QG2_MODEL:.1e} "
          f"(nominally safe x{naive_margin:.0f})")
    print(f"       BUT one cone => Coleman-Glashow threshold difference = 0 => the GZK bound does not apply as")
    print(f"       stated; the model evades it -> consistent, test disarmed  -> {'PASS' if g1 else 'FAIL'}\n")

    # ---- [G2] LHAASO / PeV photons ----
    d_lhaaso = dv(E_LHAASO)
    subluminal = ZETA > 0
    g2 = subluminal and (abs(d_lhaaso) < 1e-20)
    ok &= g2
    print("  [G2] LHAASO / PeV gamma rays (to 1.4e6 GeV):")
    print(f"       model dispersion at {E_LHAASO:.1e} GeV: v/c-1 = {d_lhaaso:.1e}")
    print(f"       LHAASO's strong LV limits are photon-DECAY (gamma->e+e-) limits = SUPERLUMINAL only;")
    print(f"       the model is SUBLUMINAL (zeta>0={subluminal}) => photon decay forbidden => bound N/A")
    print(f"       -> consistent, and the headline UHE-photon bound does not apply  -> {'PASS' if g2 else 'FAIL'}\n")

    # ---- [G3] IceCube: universality at high energy + stable UHE neutrino ----
    D_txs = TXS_Z * (C * (1.0 / 2.27e-18))   # crude Hubble distance z*c/H0, H0~70 km/s/Mpc ~2.27e-18 /s
    bound_txs = TXS_WINDOW / (D_txs / C)
    dv_txs = abs(dv(E_TXS_NU) - dv(E_TXS_GAMMA))
    d_glashow = dv(E_ICECUBE_GLASHOW)
    g3 = (dv_txs < bound_txs) and (abs(d_glashow) < 1e-20)
    ok &= g3
    print("  [G3] IceCube (TXS 0506+056 nu-gamma coincidence, and the 6.3 PeV Glashow neutrino):")
    print(f"       TXS: 290 TeV nu vs ~100 GeV gamma over z={TXS_Z} (~{D_txs/GLY:.1f} Gly), ~10-day window")
    print(f"            timing bound |dv/c| < {bound_txs:.1e};  model predicts {dv_txs:.1e} -> one-cone passes at 1e5 GeV")
    print(f"       Glashow: 6.3 PeV antineutrino, dispersion {d_glashow:.1e}, stable (subluminal)"
          f"  -> {'PASS' if g3 else 'FAIL'}\n")

    # ---- [G4] the honest reach + temper of 8.39 ----
    g4 = True
    ok &= g4
    print("  [G4] honest reach -- and a temper of Section 8.39:")
    print("       The two properties that PASSED GW170817 -- subluminal, and one universal cone -- are exactly")
    print("       what makes the model EVADE the strongest UHE bounds: photon-decay limits need superluminal,")
    print("       threshold (GZK) limits need species-dependent LV. Neither applies. So the coefficient zeta is")
    print("       currently UNCONSTRAINED by UHE data -- the model is safer, and harder to test, than 8.39's")
    print("       '~1.4 orders below the frontier, reachable by next-gen' snapshot implied (that compared E_QG")
    print("       to species-dependent bounds). The real, working handle is multimessenger UNIVERSALITY, which")
    print(f"       the model passes at MeV (GW170817) and at 1e5 GeV (TXS 0506+056).  -> PASS\n")

    print("=" * 84)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  Worked through Auger (GZK), LHAASO (PeV gammas) and IceCube (PeV neutrinos + the TXS coincidence).")
    print("  The model is CONSISTENT with every one -- but the honest finding is why: its subluminal sign kills")
    print("  the photon-decay bounds and its one-cone universality kills the species-dependent threshold bounds,")
    print("  so the coefficient zeta is not actually reached by any of them. The frontier tests the UNIVERSALITY")
    print("  (which passes, now at both MeV and 1e5 GeV), not the coefficient. Net: Section 8.39's reachability")
    print("  was too optimistic; the model's one live prediction is confirmed in its structural half and, in its")
    print("  quantitative half, currently beyond reach for a reason internal to the prediction itself.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
