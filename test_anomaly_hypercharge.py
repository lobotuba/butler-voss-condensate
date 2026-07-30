"""
Anomaly cancellation fixes the hypercharges: the model constrains the Standard Model's charges.

test_yang_mills and test_electroweak realized the SM gauge dynamics and its breaking, but took the
fermion hypercharges (and hence the electric charges) as inputs. This file asks whether the model can
instead PREDICT them. The lever is anomaly cancellation: a chiral gauge symmetry is quantum-mechanically
consistent only if its anomalies cancel, and the model's emergent gauge symmetries are EXACT -- realized
as Wilson-link lattice gauge invariances (test_yang_mills) with an exact photon Ward identity
(test_lattice_ward). An exact gauge symmetry cannot be anomalous: you cannot put an anomalous chiral
gauge theory on a lattice while keeping gauge invariance exact. So in this model anomaly cancellation is
not optional -- it is forced -- and it constrains the fermion hypercharges.

Given the observed representations of one generation -- the quark doublet Q = (3, 2), the quark singlets
u^c = (3bar, 1) and d^c = (3bar, 1), the lepton doublet L = (1, 2) and the lepton singlet e^c = (1, 1),
all left-handed Weyl, with unknown hypercharges (Y1..Y5) -- the four anomaly conditions are

   [SU(3)^2 U(1)]  : 2 Y1 + Y2 + Y3 = 0                                (color triplets)
   [SU(2)^2 U(1)]  : 3 Y1 + Y4 = 0                                    (weak doublets)
   [grav^2 U(1)]   : 6 Y1 + 3 Y2 + 3 Y3 + 2 Y4 + Y5 = 0               (all, by multiplicity)
   [U(1)^3]        : 6 Y1^3 + 3 Y2^3 + 3 Y3^3 + 2 Y4^3 + Y5^3 = 0     (cubic)

The three linear conditions leave a two-parameter family; the cubic collapses it to a unique ratio (up to
the trivial u<->d relabelling), so the hypercharges are fixed up to one overall normalisation. Fixing that
scale by the observed electric charge Q = T3 + Y reproduces the Standard Model exactly: Y = 1/6, -2/3,
1/3, -1/2, 1, hence quark charges +2/3 and -1/3 and lepton charges 0 and -1 -- and charge quantisation,
the exact cancellation of the proton and electron charge, as the unique consistent assignment.

  [A] the four anomaly functionals vanish on the SM hypercharges (each machine zero).
  [B] solving the conditions for unknown hypercharges: the linear system plus the cubic fix them uniquely
      (up to scale and the u<->d swap) -- the SM assignment, derived rather than input.
  [C] the consequence: the electric charges Q = T3 + Y come out quantised and SM-valued; the proton and
      electron charges cancel exactly because that is the anomaly-free condition.
"""
from __future__ import annotations
import numpy as np

# one generation, left-handed Weyl: (name, SU(3) dim, SU(2) dim, multiplicity = dSU3*dSU2, physical Y)
FIELDS = [("Q  (3,2)", 3, 2, 1.0 / 6), ("u^c (3bar,1)", 3, 1, -2.0 / 3), ("d^c (3bar,1)", 3, 1, 1.0 / 3),
          ("L  (1,2)", 1, 2, -1.0 / 2), ("e^c (1,1)", 1, 1, 1.0)]
Y_SM = np.array([f[3] for f in FIELDS])


def anomalies(Y):
    """The four gauge/gravitational anomaly coefficients for hypercharges Y = (Y1..Y5)."""
    dSU3 = np.array([3, 3, 3, 1, 1]); dSU2 = np.array([2, 1, 1, 2, 1]); mult = dSU3 * dSU2
    T3col = 0.5 * (dSU3 == 3)                       # SU(3) index T=1/2 for (anti)triplets, 0 for singlets
    T2wk = 0.5 * (dSU2 == 2)                        # SU(2) index T=1/2 for doublets
    A_33 = np.sum(dSU2 * T3col * Y)                 # [SU(3)^2 U(1)]
    A_22 = np.sum(dSU3 * T2wk * Y)                  # [SU(2)^2 U(1)]
    A_grav = np.sum(mult * Y)                       # [grav^2 U(1)]
    A_cubic = np.sum(mult * Y ** 3)                 # [U(1)^3]
    return np.array([A_33, A_22, A_grav, A_cubic])


if __name__ == "__main__":
    print("=== Anomaly cancellation fixes the hypercharges: the model constrains the SM charges ===\n")
    print("  The model's emergent gauge symmetries are EXACT (Wilson-link gauge invariance, exact photon")
    print("  Ward identity). An exact chiral gauge symmetry cannot be anomalous, so anomaly cancellation")
    print("  is mandatory here -- and it constrains the fermion hypercharges.\n")

    # ---------- [A] the SM hypercharges are anomaly-free ----------
    print("  [A] THE FOUR ANOMALY CONDITIONS on the observed one-generation content:")
    names = ["[SU(3)^2 U(1)]", "[SU(2)^2 U(1)]", "[grav^2 U(1)]", "[U(1)^3 cubic]"]
    A = anomalies(Y_SM)
    for n, a in zip(names, A):
        print(f"      {n:>16} = {a:+.2e}")
    print(f"      => all four vanish on the SM hypercharges (max |A| = {np.abs(A).max():.1e}, machine zero).\n")

    # ---------- [B] solve for unknown hypercharges: uniqueness ----------
    print("  [B] SOLVING for unknown hypercharges (Y1..Y5) from the conditions, given only the reps:")
    print("      linear: [SU2^2U1] -> Y4 = -3 Y1;  [grav] with [SU3^2U1] -> Y5 = 6 Y1;  Y2 + Y3 = -2 Y1.")
    print("      cubic collapses the remaining freedom: with r = Y2/Y1, the [U(1)^3] condition becomes")
    print("      r^2 + 2 r - 8 = 0, so r = -4 or r = +2.")
    roots = np.roots([1, 2, -8])
    for r in sorted(roots):
        Y1 = 1.0 / 6                                 # fix the overall scale by the observed Q = T3 + Y
        Y = np.array([Y1, r * Y1, (-2 - r) * Y1, -3 * Y1, 6 * Y1])
        tag = "  <- Standard Model" if abs(r + 4) < 1e-9 else "  (= the u<->d relabelling of the SM)"
        print(f"      r = {r:+.0f}:  Y = [{', '.join(f'{y:+.3f}' for y in Y)}]   max|anomaly| = "
              f"{np.abs(anomalies(Y)).max():.0e}{tag}")
    print("      => the two branches are the same theory with u^c and d^c swapped. Up to that relabelling")
    print("         and the overall scale, the hypercharges are UNIQUE: the Standard-Model assignment.\n")

    # ---------- [C] the consequence: quantised electric charges ----------
    print("  [C] THE ELECTRIC CHARGES Q = T3 + Y that follow (T3 = +/-1/2 in a doublet, 0 in a singlet):")
    print(f"      {'field':>14} {'Y':>8} {'Q (upper)':>11} {'Q (lower)':>11}")
    charges = {}
    for (name, d3, d2, Yv) in FIELDS:
        if d2 == 2:
            qU, qL = Yv + 0.5, Yv - 0.5
            print(f"      {name:>14} {Yv:>+8.3f} {qU:>+11.3f} {qL:>+11.3f}")
            charges[name] = (qU, qL)
        else:
            print(f"      {name:>14} {Yv:>+8.3f} {Yv:>+11.3f} {'--':>11}")
            charges[name] = (Yv,)
    q_up, q_down = charges["Q  (3,2)"]
    q_e = FIELDS[4][3]                                # e^c charge = +1, so electron (e_L) charge = -1
    q_proton = 2 * q_up + q_down                      # proton = uud
    print(f"      => quarks at +2/3 and -1/3, leptons at 0 and -1. Proton (uud) charge = "
          f"{q_proton:+.3f}, electron charge = {-q_e:+.3f}:")
    print(f"         proton + electron charge = {q_proton + (-q_e):+.1e} -- exact cancellation, i.e. CHARGE")
    print("         QUANTISATION, as the unique anomaly-free assignment. The neutron (udd) is neutral.\n")

    print("[verdict] anomaly cancellation, forced by the model's exact gauge invariance, fixes the SM charges:")
    print("  * The model's emergent gauge symmetries are EXACT lattice gauge invariances (test_yang_mills,")
    print("    exact photon Ward identity in test_lattice_ward), and an exact chiral gauge symmetry cannot")
    print("    be anomalous. So anomaly cancellation is mandatory, not an extra assumption.")
    print("  * Given only the observed representations of one generation, the four anomaly conditions fix")
    print("    the five hypercharges uniquely -- up to an overall scale and the trivial u<->d relabelling --")
    print("    to the Standard-Model values 1/6, -2/3, 1/3, -1/2, 1 [A][B]. The hypercharges are DERIVED,")
    print("    not input.")
    print("  * Hence the electric charges Q = T3 + Y are quantised and SM-valued, and the proton and")
    print("    electron charges cancel exactly [C] -- charge quantisation as a consistency requirement,")
    print("    not a coincidence. This is the sharpest Standard-Model prediction the model has made.")
    print("  * HONEST scope: anomaly cancellation fixes the hypercharges GIVEN the representations (which")
    print("    SU(3)xSU(2) rep each fermion sits in) and the one-generation content -- those remain inputs,")
    print("    as does the number of generations and the gauge group itself. What is newly derived is the")
    print("    U(1) hypercharge assignment, the one continuous freedom, forced to the SM values. And it")
    print("    assumes the emergent gauge symmetry is a genuine exact 4D symmetry; the domain-wall")
    print("    construction (test_anomaly_inflow) can park an anomaly in the bulk, so a truly 4D gauge")
    print("    symmetry -- no accessible bulk -- is the setting in which the constraint bites.")
