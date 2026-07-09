"""
Locking the cones: a boson made OF the fermions rides the fermion light cone.

test_cone_universality.py found the gap: an INDEPENDENT lattice boson (springs,
stiffness K) has cone c_B = (sqrt3/2) sqrt(K/m), while the Dirac fermions have
v_F = (3/2) t. The ratio v_F/c_B = sqrt(3) t / sqrt(K/m) is a free, tunable number
-- equality is a fine-tuning, not a symmetry. Two cones = Lorentz violation between
statistics.

The structural fix (Volovik's "everything from one structure"; Sakharov's induced
dynamics): do not put the boson in by hand. Let it be a COLLECTIVE MODE of the
fermions -- a fermion bilinear (particle-hole) excitation, or a gauge field induced
by integrating the fermions out. Then its dispersion is inherited from the fermion
sector, and the only velocity available is v_F.

The decisive, computable statement: every composite (particle-hole) boson of
momentum q costs at least the LOWER EDGE of the interband continuum,
    omega_min(q) = min_k [ E_+(k+q) - E_-(k) ] = min_k v_F(|k+q| + |k|) = v_F |q|,
so the composite-boson light cone IS the fermion light cone -- exactly, with no
tuning. A massless collective mode rides that edge.
"""
from __future__ import annotations
import numpy as np

D = np.array([[0.0, 1.0], [np.sqrt(3) / 2, -0.5], [-np.sqrt(3) / 2, -0.5]])   # nn, a=1
KPT = np.array([4 * np.pi / (3 * np.sqrt(3)), 0.0])                            # Dirac point
T = 1.0                                                                        # hopping


def absf(k):
    """|f(k)| on a grid of k (shape (...,2)) -> band magnitude E_+ = t|f|."""
    ph = k @ D.T                                   # (...,3)
    return np.abs(np.exp(1j * ph).sum(-1))


def ph_edge(q, n=501, span=3.0):
    """Lower edge of the interband particle-hole continuum at momentum transfer q:
    min_k [ E_+(k+q) + E_+(k) ]  (since E_- = -E_+)."""
    R = max(span * np.linalg.norm(q), 1e-3)
    g = np.linspace(-R, R, n)
    X, Y = np.meshgrid(g, g, indexing="ij")
    k = KPT + np.stack([X, Y], -1)
    return T * float((absf(k + q) + absf(k)).min())


def fermi_velocity(rho=1e-4, ndir=16):
    v = [absf(KPT + rho * np.array([np.cos(a), np.sin(a)])) / rho
         for a in np.linspace(0, 2 * np.pi, ndir, endpoint=False)]
    return T * float(np.mean(v))


def boson_speed(K=1.0, m=1.0, rho=1e-4):
    k = rho * np.array([1.0, 0.0])
    return float(np.sqrt((K / m) * max(3.0 - absf(k), 0.0)) / rho)


if __name__ == "__main__":
    print("=== Locking the cones: composite bosons ride the fermion light cone ===\n")
    vF = fermi_velocity()
    print(f"  fermion cone  v_F = {vF:.4f}   (= 3/2 * t, t=1)\n")

    print("[A] lower edge of the interband particle-hole continuum (the composite-boson cone)")
    print(f"  {'|q|':>7} {'direction':>12} {'omega_min(q)':>13} {'omega_min/|q|':>14} {'/ v_F':>8}")
    for qm in (0.02, 0.05, 0.10, 0.20):
        for lbl, u in (("[1,0]", (1.0, 0.0)), ("[0.6,0.8]", (0.6, 0.8))):
            u = np.array(u); u = u / np.linalg.norm(u)
            w = ph_edge(qm * u)
            print(f"  {qm:>7.2f} {lbl:>12} {w:>13.5f} {w/qm:>14.4f} {w/qm/vF:>8.4f}")
    print("\n  => omega_min(q)/|q| -> v_F as q -> 0 (0.995 at |q|=0.02), the same in every")
    print("     direction; the residual at larger q is band curvature beyond the linear Dirac")
    print("     cone. So the composite (particle-hole) bosons live on the FERMION light cone: a")
    print("     massless collective mode rides that edge, and its speed IS v_F -- inherited, not")
    print("     chosen, with no free parameter to tune.\n")

    print("[B] contrast: an INDEPENDENT lattice boson has its own, tunable cone")
    for K, m in ((1.0, 1.0), (3.0, 1.0), (1.0, 3.0)):
        cB = boson_speed(K, m)
        print(f"   springs K={K}, m={m}:  c_B = {cB:.4f},  v_F/c_B = {vF/cB:.4f}  (free parameter)")

    print("\n[C] verdict")
    print("  * An independent boson field brings its OWN cone -> v_F/c_B is a tunable ratio, and")
    print("    equality is a fine-tuning. That was the Lorentz gap found in test_cone_universality.")
    print("  * A COMPOSITE boson -- a fermion bilinear, or a gauge field induced by integrating the")
    print("    fermions out (Sakharov) -- has no cone of its own: it inherits v_F exactly. ONE cone,")
    print("    no tuning, protected by the structure rather than by a coincidence of couplings.")
    print("  => the cure is the project's recurring lesson at its deepest level: ALL excitations")
    print("     must descend from ONE structure. Emergent Lorentz invariance across statistics is")
    print("     then automatic -- and, in Volovik's picture, the emergent gauge field and graviton")
    print("     (tetrad fluctuations of the fermion cone) ride that same cone too.")
