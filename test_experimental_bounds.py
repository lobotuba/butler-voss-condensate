"""
The model against the data: its two live predictions vs current experimental bounds.

A model earns its keep by being falsifiable. This project now carries two genuine predictions and one
inherited fine-tuning claim, and none of them has been confronted with real numbers. This file does
that, and the confrontation forces an honest correction to an earlier result.

  [A] n=2 LORENTZ VIOLATION (test_lv_prediction): v(E)/c = 1 - zeta (E/E_Planck)^2, subluminal,
      cross-species universal, quadratic. Compared against photon time-of-flight limits on quadratic
      dispersion (GRB/AGN; the literature quotes E_QG,2 of order 1e10-1e12 GeV).

  [B] SCALE-DEPENDENT gamma (test_two_gravities): with a massless graviton and a GAPPED amplitude
      mode, gamma(r) = g_t / (g_t + g_s exp(-r/lambda)) climbs from below 1 up to 1 across the
      amplitude Compton wavelength lambda = 1/m_A. Equivalently the amplitude mode is a Yukawa
      ADDITION to gravity of range lambda. Compared against Cassini's solar-system gamma and against
      short-range (Eot-Wash torsion-balance) tests of the inverse-square law.

  [C] THE 1e122 CRITICALITY TUNING (test_scale_fixing) — which this section RETRACTS. That result
      read gravity's RANGE as 1/m_A, so any bound on gravity being 1/r^2 out to solar-system or
      larger scales forced 1/m_A astronomically large, hence m_A^2 = 2a tiny, hence a medium tuned to
      about one part in 1e122. That chain assumed the SCALAR (amplitude-mode) picture of gravity.
      The tensor arc (test_deconfinement, test_induced_sign, test_spin2_dynamical) replaced it: the
      long-range force is carried by the MASSLESS deconfined graviton, and the amplitude mode is only
      a short-range correction. The constraint on m_A is therefore a LOWER bound, not an upper one —
      and an untuned m_A satisfies it with enormous margin. The fine-tuning is gone.

Scope and honesty. The experimental numbers used here are standard literature ORDERS OF MAGNITUDE
(Cassini's gamma is the one precise input); the comparison is an order-of-magnitude confrontation,
not a likelihood analysis. Two model-side assumptions are made explicit: the scalar and tensor
couplings are taken comparable (g_s ~ g_t), and the Yukawa strength relative to Newtonian gravity is
taken of order unity (alpha ~ 1), both natural since each couples to energy. Loosening either
weakens the bound in [B] proportionally.
"""
from __future__ import annotations
import numpy as np

HBARC = 1.9733e-7        # eV*m
E_PLANCK = 1.22e19       # GeV
AU = 1.496e11            # m
L_PLANCK = 1.616e-35     # m
M_PLANCK_EV = 1.22e28    # eV
DARK_ENERGY_MEV = 2.4    # (rho_Lambda)^(1/4), meV
CASSINI = 2.3e-5         # 1-sigma on |gamma - 1|


def gamma_of_r(r, lam, gs_over_gt=1.0):
    """gamma(r) = g_t / (g_t + g_s exp(-r/lambda)) -- 1 at long range, suppressed inside lambda."""
    return 1.0 / (1.0 + gs_over_gt * np.exp(-r / lam))


if __name__ == "__main__":
    print("=== The model against the data: two predictions and one retraction ===\n")

    # ---------------- [A] Lorentz violation ----------------
    print("  [A] n=2 LORENTZ VIOLATION:  v(E)/c = 1 - zeta (E/E_Planck)^2,  zeta ~ O(1)")
    print(f"      {'photon bound E_QG,2':>20} {'case':>12} {'effect/sensitivity':>19} "
          f"{'orders (energy)':>16} {'orders (effect)':>16}")
    for EQG, tag in ((1e10, "conservative"), (1e11, "typical"), (1e12, "optimistic")):
        supp = (EQG / E_PLANCK) ** 2
        print(f"      {EQG:>16.0e} GeV {tag:>12} {supp:>19.1e} "
              f"{np.log10(E_PLANCK / EQG):>16.1f} {np.log10(1 / supp):>16.0f}")
    print("      => SAFE against every current bound, by ~16 orders in the effect. But that cuts")
    print("         both ways: the prediction is NOT currently falsifiable. Reaching it needs ~8")
    print("         orders of improvement in quadratic-dispersion sensitivity. An honest tempering")
    print("         of the claim that this is the model's 'testable' prediction: it is safe, and out")
    print("         of reach. (A quadratic, Planck-suppressed effect is intrinsically hard; a LINEAR")
    print("         n=1 signal, which the model does NOT predict, is what current experiments probe.)\n")

    # ---------------- [B] short-range gravity ----------------
    print("  [B] SCALE-DEPENDENT gamma: the amplitude mode as a Yukawa addition of range 1/m_A")
    lam_cassini = AU / np.log(1.0 / CASSINI)
    print(f"      Cassini: |gamma-1| < {CASSINI:.1e} at 1 AU  ->  lambda < {lam_cassini:.2e} m "
          f"(a weak bound)")
    print(f"      {'Eot-Wash reach':>16} {'-> m_A bound':>14} {'vs dark-energy scale':>22}")
    for lam_um in (56, 40, 30):
        mA_meV = HBARC / (lam_um * 1e-6) * 1e3
        print(f"      {lam_um:>13} um {mA_meV:>11.2f} meV {mA_meV / DARK_ENERGY_MEV:>19.2f} x")
    mA_bound = HBARC / 50e-6 * 1e3
    print(f"      => the BINDING constraint is short-range gravity, beating Cassini by "
          f"{np.log10(lam_cassini / 50e-6):.0f} orders.")
    print(f"         Taking alpha ~ 1, the amplitude gap must satisfy m_A >~ {mA_bound:.1f} meV.")
    print(f"         That sits within a factor {mA_bound/DARK_ENERGY_MEV:.1f} of the dark-energy scale "
          f"({DARK_ENERGY_MEV} meV) -- the")
    print("         well-known coincidence that makes ~100 um the frontier of these experiments.")
    print("         THIS is the model's genuinely testable prediction: a Yukawa of gravitational")
    print("         strength just below the current reach would show up as gamma < 1 at short range.")
    print("      gamma(r) profile for lambda = 50 um (g_s = g_t):")
    lam = 50e-6
    print(f"        {'r':>12} {'gamma(r)':>10}")
    for r, lab in ((5e-6, "0.1 lam"), (5e-5, "1 lam"), (2.5e-4, "5 lam"), (5e-4, "10 lam"), (AU, "1 AU")):
        print(f"        {lab:>12} {gamma_of_r(r, lam):>10.6f}")
    print("      => gamma is GR (=1) to far beyond experimental precision at solar-system scales, and")
    print("         departs only inside the sub-millimetre window. Consistent, and falsifiable there.\n")

    # ---------------- [C] the retraction ----------------
    print("  [C] RETRACTION: the 1e122 criticality tuning of test_scale_fixing does NOT survive.")
    a_if_at_bound = (mA_bound * 1e-3 / M_PLANCK_EV) ** 2
    print(f"      OLD (scalar gravity): gravity's range = 1/m_A, so 1/m_A had to be astronomical")
    print(f"        -> m_A^2 = 2a tiny -> medium tuned to ~1 part in 1e122. A severe fine-tuning.")
    print(f"      NEW (tensor gravity): the long-range force is the MASSLESS deconfined graviton;")
    print(f"        the amplitude mode is only a short-range correction. The data give a LOWER bound,")
    print(f"        m_A >~ {mA_bound:.1f} meV, with NO upper bound -- a larger gap is only safer.")
    print(f"      => an UNTUNED medium (m_A ~ M_Planck) clears the bound by "
          f"{np.log10(M_PLANCK_EV / (mA_bound*1e-3)):.0f} orders.")
    print(f"        (Only if m_A sat exactly at the bound would a/a0 ~ {a_if_at_bound:.0e} be needed;")
    print("         nothing requires that.) The fine-tuning is DISSOLVED, not merely reduced.\n")

    print("[verdict] the model meets the data, and sheds a fine-tuning doing it:")
    print("  * LORENTZ VIOLATION: safe by ~16 orders in the effect -- and, honestly, out of reach.")
    print("    A quadratic Planck-suppressed signal is not falsifiable with present sensitivity; the")
    print("    earlier framing of this as the model's testable prediction was too generous.")
    print("  * SHORT-RANGE GRAVITY: this is the real testable one. A massless graviton plus a gapped")
    print("    amplitude mode predicts gamma < 1 inside 1/m_A, and current torsion-balance limits")
    print(f"    already require m_A >~ {mA_bound:.1f} meV -- coincidentally within a factor "
          f"{mA_bound/DARK_ENERGY_MEV:.1f} of the dark-")
    print("    energy scale, exactly the sub-millimetre window these experiments are now pushing into.")
    print("  * AND the confrontation RETRACTS test_scale_fixing's 1e122 criticality tuning: that")
    print("    number came from reading gravity's range as 1/m_A, which the tensor arc replaced. The")
    print(f"    surviving constraint is a lower bound an untuned medium clears by "
          f"~{np.log10(M_PLANCK_EV / (mA_bound*1e-3)):.0f} orders. Together")
    print("    with the cosmological-constant result (which dissolved the OTHER 1e122 tuning by the")
    print("    condensate's equilibrium thermodynamics), the model has now shed BOTH of its 10^122")
    print("    fine-tunings -- each one removed by a structural result rather than a fitted parameter.")
