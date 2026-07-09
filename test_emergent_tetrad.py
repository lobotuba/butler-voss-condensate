"""
The capstone: one structure -> fermions + an emergent photon + an emergent graviton,
all on the same light cone.

Near a Dirac node the fermion Hamiltonian is, to linear order,
    H = e^a_i  sigma_a  (k_i - A_i),
where
  * A_i  = the POSITION of the Dirac node  -> an emergent U(1) GAUGE FIELD (photon),
  * e^a_i = the SHAPE of the Dirac cone     -> an emergent TETRAD / metric (graviton).
Both are features of the fermion dispersion itself, so neither has a light cone of
its own: they ride the fermion cone by construction. This is Volovik's mechanism,
and it is exactly what test_cone_universality showed was needed.

Here we show the model's OWN medium supplies them. Perturb the three nearest-
neighbour bonds of the honeycomb, t_j -> t(1 + u_j) -- the medium's elastic degrees
of freedom -- and read off, from the fermion bands alone:
  A   = shift of the Dirac node            (the emergent gauge field)
  G_ij = Re(M_i conj(M_j)),  M = grad f    (the emergent metric; E(q)^2 = q.G.q)
        trace part  -> conformal factor (v_F^2)
        traceless part h_ij -> the spin-2 graviton polarisations
"""
from __future__ import annotations
import numpy as np

D = np.array([[0.0, 1.0], [np.sqrt(3) / 2, -0.5], [-np.sqrt(3) / 2, -0.5]])   # nn, a=1
KPT = np.array([4 * np.pi / (3 * np.sqrt(3)), 0.0])


def f_and_grad(k, t):
    ph = np.exp(1j * (D @ k))                       # (3,)
    f = float_complex = (t * ph).sum()
    M = 1j * (t * ph) @ D                           # dF/dk  (complex 2-vector)
    return f, M


def find_node(t, k0=KPT, iters=60):
    """Newton-solve f(k)=0 (two real equations) for the Dirac-node position."""
    k = np.array(k0, float)
    for _ in range(iters):
        f, M = f_and_grad(k, t)
        F = np.array([f.real, f.imag])
        J = np.array([[M[0].real, M[1].real], [M[0].imag, M[1].imag]])
        if abs(np.linalg.det(J)) < 1e-12:
            break
        k = k - np.linalg.solve(J, F)
    return k


def emergent_fields(u):
    """Return (A, v2, h) for bond perturbation u=(u1,u2,u3)."""
    t = 1.0 + np.asarray(u, float)
    K = find_node(t)
    A = K - KPT                                     # emergent gauge field
    _, M = f_and_grad(K, t)
    G = np.real(np.outer(M, np.conj(M)))            # metric: E^2 = q.G.q
    G = 0.5 * (G + G.T)
    v2 = 0.5 * np.trace(G)                          # conformal factor = v_F^2
    h = G - v2 * np.eye(2)                          # traceless: the graviton
    return A, v2, h


def show(name, u):
    A, v2, h = emergent_fields(u)
    hmag = np.sqrt(h[0, 0] ** 2 + h[0, 1] ** 2)
    print(f"  {name:>26}  A = ({A[0]:+.5f},{A[1]:+.5f})  |A|={np.linalg.norm(A):.5f}"
          f"   v_F^2 = {v2:7.4f}   |h| = {hmag:.5f}")
    return A, v2, h


if __name__ == "__main__":
    print("=== One structure: bonds -> fermions + emergent photon + emergent graviton ===\n")
    print("  Perturb the medium's three nn bonds, t_j -> t(1+u_j), and read the fermion bands:")
    print("    A  = Dirac-node shift      (emergent GAUGE FIELD, spin-1)")
    print("    h  = traceless part of G   (emergent METRIC / GRAVITON, spin-2)\n")

    show("unperturbed (0,0,0)", (0, 0, 0))
    print("    -> isotropic cone, v_F^2 = 2.25 = (3/2)^2; no photon, no graviton.\n")

    show("uniform  (u,u,u)", (0.05, 0.05, 0.05))
    print("    -> node does NOT move (A=0): a uniform bond stretch is a pure CONFORMAL rescaling")
    print("       of the cone (v_F changes, h=0). No photon, no graviton -- only a scale change.\n")

    show("doublet  (2u,-u,-u)", (0.06, -0.03, -0.03))
    show("doublet  (0,u,-u)", (0.0, 0.03, -0.03))
    show("single bond (u,0,0)", (0.05, 0.0, 0.0))
    print("    -> the doublet (E-representation) bond fluctuations move the node (a PHOTON, A != 0)")
    print("       AND deform the cone anisotropically (a GRAVITON, h != 0). One microscopic")
    print("       degree of freedom -- the bonds -- sources both.\n")

    # linearity: both A and h are linear in u at small amplitude
    print("  linearity check (both fields are linear responses to the bond fluctuation):")
    for eps in (0.02, 0.01, 0.005):
        A, _, h = emergent_fields((2 * eps, -eps, -eps))
        hm = np.sqrt(h[0, 0] ** 2 + h[0, 1] ** 2)
        print(f"    u={eps:.3f}:  |A|/u = {np.linalg.norm(A)/eps:.4f}   |h|/u = {hm/eps:.4f}")

    print("\n  => H_eff = e^a_i sigma_a (k_i - A_i): the photon (A) and the graviton (e, via G)")
    print("     are BOTH read off the fermion dispersion. Neither is an independent field, so")
    print("     neither carries a light cone of its own -- both ride the fermion cone, exactly")
    print("     as test_cone_lock and test_induced_action require.")
    print("  => ONE structure (the medium's bonds) yields the fermions, the gauge field and the")
    print("     metric, all Lorentz-invariant on a single cone. This is the construction the")
    print("     whole program pointed to: the emergent photon and the emergent graviton are not")
    print("     added to the medium -- they ARE the medium's fluctuations, seen by its fermions.")
