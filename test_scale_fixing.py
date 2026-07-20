"""
Fixing the medium's ABSOLUTE SCALE from measured constants.

*** STATUS UPDATE -- the 1e122 criticality tuning below is RETRACTED (test_experimental_bounds).
    What STANDS: a0 = l_Planck. Matching the measured G with an order-unity coupling still forces
    the node spacing to the Planck length, and that consistency result is unaffected.
    What FAILS: the inference that the medium must be tuned to ~1 part in 1e122 of criticality. That
    chain assumed gravity's RANGE is 1/m_A -- i.e. that the amplitude mode IS gravity's mediator --
    so any bound on gravity being 1/r^2 out to large distances forced 1/m_A astronomically large and
    hence m_A^2 = 2a absurdly small. The tensor-gravity arc replaced that premise: the long-range
    force is carried by the MASSLESS deconfined graviton (test_deconfinement, test_induced_sign,
    test_spin2_dynamical), and the amplitude mode is only a short-range correction. The surviving
    experimental constraint is therefore a LOWER bound, m_A >~ 4 meV from short-range gravity, with
    NO upper bound -- a larger gap is only safer -- and an UNTUNED medium (m_A ~ M_Planck) clears it
    by ~30 orders. The fine-tuning is dissolved, not merely reduced.
    Also superseded: "the model predicts a Yukawa, not a pure 1/r^2". Long-range gravity is now the
    massless graviton's exact 1/r^2; the Yukawa survives only as the sub-millimetre correction that
    test_experimental_bounds turns into the model's one genuinely testable prediction. ***

Until now every number in this project has been in lattice units. Robert's instinct was the
right one: find a real observable that pins the node spacing a0 to metres. His route -- use
black-hole masses -- does not work (see the companion test_collapse.py, and the note below),
but the gravity result of test_critical_gravity hands us a better anchor, because it makes
gravity's RANGE a measurable prediction of the medium.

The model now supplies:
    c     : the medium's wave speed              (the emergent light speed)
    hbar  : via test_emergent_qm, hbar/2m = c^2/(2 Omega)   -- a material property
    G     : via test_critical_gravity, the amplitude-mode exchange
              E_int(R) = -(2 g phi0)^2 Q1 Q2 exp(-m_A R) / (4 pi R)
            which, compared with -G M1 M2 / R and with gravitational charge = mass, gives
              G = (2 g phi0)^2 / (4 pi mu^2)
    range : lambda_g = 1 / m_A                  -- gravity is a YUKAWA, not a pure 1/r^2

Two things then follow, and they are the point of this file.

(1) NATURALNESS FIXES a0. In lattice units the only length is a0, so the gravitational
    coupling per unit ENERGY is 1/M_* where M_* is the lattice energy scale. Matching the
    measured G forces M_* = M_Planck, i.e. a0 = l_Planck. Note what this is: test_lv_prediction
    ASSUMED a0 = the Planck length. It is now a CONSISTENCY RESULT -- the model reproduces G
    with an order-unity coupling if and only if the nodes sit at the Planck length. (Honest
    status: this is naturalness + dimensional analysis, not a hard derivation, because hbar is
    not independently derived -- the model gives relations among hbar, c and Omega, and
    quantisation is still imposed, cf. test_quantization.)

(2) THE OBSERVED RANGE OF GRAVITY MEASURES THE MEDIUM'S DISTANCE FROM CRITICALITY. This is the
    real calibration, and it is data, not assumption. Gravity is long-range, so the amplitude
    gap is tiny: m_A = 1/lambda_g. In lattice units m_A^2 = 2a, so
              a = (a0 / lambda_g)^2 / 2,
    and a must be tuned to that fraction of its natural (order-unity) value. Every graviton-mass
    bound therefore converts directly into a number telling us how finely the condensate is
    tuned. That number is the model's hierarchy problem, quantified.

FALSIFIABLE. The model does NOT predict a pure inverse-square law: it predicts a YUKAWA, with a
finite range. A measured nonzero graviton mass would not refute the model -- it would MEASURE
the medium's distance from its critical point.
"""
from __future__ import annotations

# --- measured constants (SI) ---
C = 2.99792458e8            # m/s
HBAR = 1.054571817e-34      # J s
G = 6.67430e-11             # m^3 kg^-1 s^-2
EV = 1.602176634e-19        # J

L_PLANCK = (HBAR * G / C ** 3) ** 0.5
M_PLANCK = (HBAR * C / G) ** 0.5                 # kg
E_PLANCK = M_PLANCK * C ** 2 / EV                # eV

# --- graviton-mass / gravity-range bounds (Compton wavelength lambda_g, metres) ---
BOUNDS = [
    ("LIGO-Virgo GW dispersion",        1.6e16),
    ("Solar-system (Yukawa) tests",     2.8e15),
    ("Galaxy-cluster dynamics",         6.2e22),
    ("Cosmological / large-scale",      1.0e26),
]


def m_A_eV(lam):
    """Graviton (amplitude-mode) mass from its Compton wavelength."""
    return HBAR * C / lam / EV


if __name__ == "__main__":
    print("=== Fixing the medium's absolute scale from measured constants ===\n")

    print("  [1] NATURALNESS FIXES THE NODE SPACING")
    print(f"      Planck length  l_P = sqrt(hbar G / c^3) = {L_PLANCK:.3e} m")
    print(f"      Planck energy  E_P                      = {E_PLANCK:.3e} eV")
    print("      The model's gravity is amplitude-mode exchange with coupling (2 g phi0); matching")
    print("      the measured G with an ORDER-UNITY lattice coupling requires the lattice energy")
    print("      scale to be the Planck scale, i.e.")
    print(f"          a0 = l_Planck = {L_PLANCK:.3e} m")
    print("      This was an ASSUMPTION in test_lv_prediction. It is now a CONSISTENCY RESULT:")
    print("      the model reproduces Newton's constant without an unnatural coupling exactly")
    print("      when the nodes sit at the Planck length -- and that is the same a0 the")
    print("      Lorentz-violation prediction was built on, so the two hang together.\n")

    print("  [2] THE RANGE OF GRAVITY MEASURES THE DISTANCE FROM CRITICALITY")
    print("      m_A = 1/lambda_g, and in lattice units m_A^2 = 2a, so a = (a0/lambda_g)^2 / 2.")
    print("      'tuning' = a / a_natural, with a_natural = O(1) in lattice units.\n")
    print(f"      {'bound on gravity range':<30} {'lambda_g (m)':>11} {'m_A (eV)':>11} "
          f"{'a0/lambda_g':>12} {'a (lattice)':>13} {'tuning':>10}")
    for name, lam in BOUNDS:
        r = L_PLANCK / lam
        a = r ** 2 / 2
        print(f"      {name:<30} {lam:>11.1e} {m_A_eV(lam):>11.1e} {r:>12.1e} {a:>13.1e} "
              f"{'1 in 1e%d' % round(-__import__('math').log10(a)):>10}")

    strongest = max(BOUNDS, key=lambda t: t[1])
    r = L_PLANCK / strongest[1]
    a = r ** 2 / 2
    print(f"\n      => taking the strongest bound ({strongest[0]}, lambda_g > {strongest[1]:.0e} m):")
    print(f"         the condensate must sit within ~1 part in 1e{round(-__import__('math').log10(a))}"
          " of its critical point.")
    print("         For scale: the cosmological-constant problem is the famous ~1e120. This is the")
    print("         SAME KIND of fine-tuning, and very nearly the same size -- which is a strong")
    print("         hint that the two are the same problem wearing different clothes. (The model")
    print("         already grew its own cosmological term, cf. the induced Pi_+(0,0) of")
    print("         test_induced_gravity.)\n")

    print("[verdict] the medium now has an absolute scale, and one honest new prediction.")
    print("  * a0 = l_Planck = 1.6e-35 m, not by assumption but by requiring an order-unity")
    print("    gravitational coupling to reproduce the measured G. The Lorentz-violation")
    print("    prediction of test_lv_prediction, which assumed this, is therefore self-consistent.")
    tune = round(-__import__('math').log10(a))
    print(f"  * The price is explicit and enormous: the condensate must be critical to ~1 part in")
    print(f"    1e{tune}. That IS gravity's weakness, restated as a property of the medium -- and it")
    print("    is the model's hierarchy problem, now quantified rather than hidden. That it lands")
    print("    on ~1e120, the cosmological-constant number, is either a coincidence or the point.")
    print("  * FALSIFIABLE, and this is the useful part: the model does NOT predict a pure 1/r^2.")
    print("    It predicts a YUKAWA of finite range lambda_g = 1/m_A. A measured graviton mass")
    print("    would not kill the model -- it would MEASURE how far the medium sits from its")
    print("    critical point, via a = (a0/lambda_g)^2/2. Gravity's range is a thermometer for")
    print("    the condensate.")
    print("\n  Note on the black-hole route. Using black-hole masses to fix a0 cannot work, for a")
    print("  reason worth recording: a 'node spacing -> 0' criterion is a DENSITY condition, so it")
    print("  predicts ONE universal density for black holes. Real horizon density goes as 1/M^2 --")
    print("  ~1.8e17 kg/m^3 for a 10-solar-mass hole, but ~4e-3 kg/m^3 for TON 618, three hundred")
    print("  times THINNER THAN AIR. Twenty orders of magnitude. A horizon is about COMPACTNESS")
    print("  (2GM/rc^2 = 1), not density. See test_collapse.py for what the medium's collapse")
    print("  threshold actually is, and why it is not a black hole.")
