"""
Cone universality across STATISTICS: is the fermion cone the boson cone?

test_lorentz_unified.py locked the medium and the (bosonic) field to a single
universal light cone by governing both with ONE operator. test_dirac.py then added
relativistic fermions -- with their own Fermi velocity v_F. That is a THIRD cone,
and nothing so far forces v_F = c_boson. If they differ, Lorentz invariance is
broken BETWEEN the fermion and boson sectors, and the earlier "emergent Lorentz"
result was a within-sector statement only.

Both cones live on the SAME honeycomb structure factor f(k) = sum_j exp(i k.delta_j):
  * FERMION (tight-binding, hopping t): bands E = +/- t|f|. Near a Dirac point |f|
    vanishes linearly -> v_F = (3/2) t a.
  * BOSON (vector-Hooke springs, stiffness K, mass m): acoustic branch
    omega^2 = (K/m)(z - |f|), z = 3. Near Gamma, z-|f| ~ (3/4)|k|^2 a^2
    -> c_B = (sqrt(3)/2) sqrt(K/m) a.
So v_F / c_B = sqrt(3) * t / sqrt(K/m): a FREE ratio, set by independent couplings.
Equality is a fine-tuning, not a symmetry.
"""
from __future__ import annotations
import numpy as np

D = np.array([[0.0, 1.0], [np.sqrt(3) / 2, -0.5], [-np.sqrt(3) / 2, -0.5]])   # nn vectors, a=1
KPT = np.array([4 * np.pi / (3 * np.sqrt(3)), 0.0])                            # a Dirac point


def f(k):
    return np.sum(np.exp(1j * (D @ np.asarray(k, float))))


def fermi_velocity(t=1.0, rho=1e-4, ndir=16):
    """v_F = slope of E = t|f| just off the Dirac point (isotropic)."""
    v = [abs(f(KPT + rho * np.array([np.cos(a), np.sin(a)]))) / rho
         for a in np.linspace(0, 2 * np.pi, ndir, endpoint=False)]
    return t * float(np.mean(v))


def boson_speed(K=1.0, m=1.0, rho=1e-4, ndir=16):
    """Acoustic branch of the vector-Hooke honeycomb: omega^2 = (K/m)(3 - |f|)."""
    c = []
    for a in np.linspace(0, 2 * np.pi, ndir, endpoint=False):
        k = rho * np.array([np.cos(a), np.sin(a)])
        c.append(np.sqrt((K / m) * max(3.0 - abs(f(k)), 0.0)) / rho)
    return float(np.mean(c))


if __name__ == "__main__":
    print("=== Cone universality across statistics: fermion v_F vs boson c_B ===\n")
    print(f"  {'t':>5} {'K':>5} {'m':>5} | {'v_F (fermion)':>14} {'c_B (boson)':>12} {'v_F / c_B':>10}")
    rows = [(1.0, 1.0, 1.0), (1.0, 3.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 3.0)]
    for t, K, m in rows:
        vF, cB = fermi_velocity(t), boson_speed(K, m)
        print(f"  {t:>5.2f} {K:>5.2f} {m:>5.2f} | {vF:>14.4f} {cB:>12.4f} {vF/cB:>10.4f}")

    vF, cB = fermi_velocity(1.0), boson_speed(1.0, 1.0)
    print(f"\n  at t=K=m=1:  v_F/c_B = {vF/cB:.4f}  (= sqrt(3) = {np.sqrt(3):.4f})")
    print("  scaling: v_F/c_B = sqrt(3) * t / sqrt(K/m) -- an arbitrary, tunable ratio.")
    Ktune = 3.0
    print(f"  to force one cone you must TUNE, e.g. K/m = 3 t^2 gives v_F/c_B = "
          f"{fermi_velocity(1.0)/boson_speed(Ktune,1.0):.4f}")

    print("\n  => the fermion and boson cones are set by INDEPENDENT couplings (t vs K/m).")
    print("     Equality is a fine-tuning, not a symmetry: generically TWO cones, i.e. Lorentz")
    print("     violation BETWEEN statistics. (Real graphene is the cautionary case: its Dirac")
    print("     fermions have v_F ~ c/300 while its phonons run ~100x slower still -- only the")
    print("     fermion sector is even approximately relativistic.)")
    print("\n  Honest correction: the emergent-Lorentz result (test_lorentz*) is a WITHIN-SECTOR")
    print("  statement -- one round universal cone for the medium + bosonic field. Adding fermions")
    print("  opens a new cone that is not automatically locked. A genuinely Lorentz-invariant world")
    print("  needs ALL excitations to descend from ONE structure (bosons as collective modes of the")
    print("  fermions, as near a Fermi point), or a symmetry that relates them -- supersymmetry is")
    print("  precisely the boson-fermion symmetry that would lock the two cones together.")
