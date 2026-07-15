"""
Two gravities, and which one wins at long range. The scalar-vs-tensor showdown.

The model carries TWO universal attractions, and they disagree on gamma:
  * SCALAR  -- the amplitude (Higgs) mode. Mass couples to it at tree level (g rho |chi|^2).
    It gives the working attractive force of test_critical_gravity. It is spin-0, so gamma = 0
    (Nordstrom): no light bending. Its range is 1/m_A, and m_A^2 = 2a is UNPROTECTED (a radial
    mode is never symmetry-protected), so it is generically GAPPED -- short range -- and only
    became long-range in test_critical_gravity by TUNING a -> 0 (criticality).
  * TENSOR  -- the incompatible/graviton sector. Spin-2, so gamma = 1 (GR, light bending x2). Its
    mass is PROTECTED to zero IF linearised diffeomorphism invariance emerges (test_graviton_ward:
    argued from the conserved IR stress tensor, not yet cleanly measured).

Which dominates the long-range force? The decisive fact is not the coupling strength -- it is the
MASS. A massless mediator's 1/r ALWAYS beats a massive mediator's e^{-m r}/r at large r
(power law beats exponential), for ANY couplings. So:

    long-range gravity  =  the LIGHTER mediator.

If the graviton is massless (protected) and the amplitude mode is gapped (unprotected), then the
GRAVITON wins at long range -> gamma -> 1 (real GR asymptotically), while the scalar survives only
inside its Compton wavelength 1/m_A as a gamma=0 Yukawa correction. The observable is then a
SCALE-DEPENDENT gamma: GR at long range, scalar-contaminated below 1/m_A.

The sting, and the honest crux: this REINTERPRETS test_critical_gravity. Its "working gravity"
was the scalar made long-range by tuning m_A -> 0 -- which is exactly the wrong move for GR,
because a long-range scalar DOMINATES and forces gamma = 0. The correct picture leaves the
amplitude mode GAPPED (its natural state) and lets the massless graviton be gravity. But the
graviton only exists if the incompatible sector is DECONFINED into a massless mode -- and in the
pure medium that sector is CONFINING (test_disclination_force: energy grows with separation) and
must be deconfined by the Sakharov-induced term. So:

    the two-gravities hurdle  ==  deconfine the spin-2 sector into a massless graviton.

Until that is done, the model's only demonstrated long-range attraction is the gamma=0 scalar.
This file makes the showdown quantitative: the potentials, the long-range winner as a function of
the two masses, and the gamma(r) profile that results.
"""
from __future__ import annotations
import numpy as np


def phi_scalar(r, gs, mA):
    """Amplitude-mode (Yukawa) time-potential. Contributes only to g_00 (gamma = 0)."""
    return gs * np.exp(-mA * r) / r


def phi_tensor(r, gt, mg):
    """Graviton time-potential. Being spin-2 it contributes EQUALLY to g_00 and g_ij (gamma = 1)."""
    return gt * np.exp(-mg * r) / r


if __name__ == "__main__":
    print("=== Two gravities: which wins at long range? ===\n")
    print("  scalar (amplitude mode): spin-0, gamma=0, mass m_A UNPROTECTED (gapped unless tuned)")
    print("  tensor (graviton):       spin-2, gamma=1, mass m_g PROTECTED to 0 IF diff-invariant")
    print("  Long-range winner = the LIGHTER mediator (1/r beats e^{-mr}/r for any coupling).\n")

    r = np.linspace(2, 200, 400)
    gs, gt = 1.0, 1.0          # equal couplings, to isolate the MASS effect

    # ---- A. graviton massless, scalar gapped (the natural case) ----
    mA = 0.05                  # amplitude gap: range 1/m_A = 20
    mg = 0.0                   # protected massless graviton
    Ps, Pt = phi_scalar(r, gs, mA), phi_tensor(r, gt, mg)
    gamma = Pt / (Ps + Pt)     # gamma(r) = Psi/Phi = tensor / (scalar + tensor)
    print("  [A] graviton massless (m_g=0), amplitude mode gapped (m_A=0.05, range 20):")
    print(f"      {'r':>6} {'Phi_scalar':>11} {'Phi_tensor':>11} {'gamma(r)':>9}")
    for rr in (5, 10, 20, 40, 80, 160):
        i = np.argmin(np.abs(r - rr))
        print(f"      {rr:>6} {Ps[i]:>11.3e} {Pt[i]:>11.3e} {gamma[i]:>9.3f}")
    print("      => gamma climbs 0 -> 1 across the amplitude Compton wavelength 1/m_A = 20:")
    print("         GR (gamma=1) at long range, scalar-contaminated (gamma<1) at short range.")
    print("         The MASSLESS graviton wins asymptotically -- real gravity is spin-2 at range.\n")

    # ---- B. the test_critical_gravity tuning: scalar made massless too ----
    Ps2 = phi_scalar(r, gs, 0.0)
    Pt2 = phi_tensor(r, gt, 0.0)
    gamma2 = Pt2 / (Ps2 + Pt2)
    print("  [B] test_critical_gravity's tuning -- amplitude mode ALSO driven massless (m_A -> 0):")
    print(f"      both now 1/r; gamma = {gamma2[0]:.2f} at ALL r (scalar and tensor tie).")
    print("      Tuning the scalar long-range forces gamma = 1/2, NOT 1. A long-range scalar is a")
    print("      permanent gamma-spoiler. The right move is to leave the amplitude mode GAPPED.\n")

    # ---- C. the actual pure-medium situation ----
    print("  [C] but in the PURE medium the tensor sector is CONFINING, not massless:")
    print("      test_disclination_force: two curvature charges have energy GROWING with R")
    print("      (repulsive, ~R^2). A confining sector contributes NOTHING to the long-range")
    print("      attraction -- so the pure-medium long-range force is 100% the scalar => gamma = 0.")
    print("      The massless graviton of [A] exists ONLY if the Sakharov-induced Einstein term")
    print("      DECONFINES the incompatible sector (turns +R into -1/r). That is the hurdle.\n")

    print("[verdict] the showdown has a clean logic and a hard crux:")
    print("  * WHICH WINS AT LONG RANGE IS SET BY MASS, not coupling: the lighter mediator's 1/r")
    print("    always beats the heavier's Yukawa. So IF the graviton is massless and the amplitude")
    print("    mode is gapped, gravity is spin-2 (gamma=1) at long range and scalar (gamma=0) only")
    print("    within 1/m_A -- a SCALE-DEPENDENT gamma, a genuine falsifiable prediction (a Yukawa")
    print("    deviation from GR below the amplitude-mode Compton wavelength, the regime short-")
    print("    range gravity experiments probe).")
    print("  * REINTERPRETATION: test_critical_gravity's 'working gravity' was the scalar tuned")
    print("    long-range (m_A -> 0) -- which is exactly what makes gamma = 0 and is the WRONG move")
    print("    for GR. The amplitude mode should stay gapped (a short-range correction); gravity")
    print("    proper is the graviton.")
    print("  * THE HURDLE, sharpened: the spin-2 sector is CONFINING in the pure medium, so it is")
    print("    not yet a massless graviton at all. Clearing the two-gravities hurdle = DECONFINING")
    print("    the incompatible sector into a massless spin-2 via the induced (Sakharov) term. That")
    print("    is the same fermion-loop question as the graviton Ward identity -- now with a sharp")
    print("    physical target: turn the confining +R potential of test_disclination_force into a")
    print("    Newtonian -1/r. If it deconfines, the mass hierarchy does the rest and gamma -> 1.")
