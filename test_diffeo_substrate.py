"""
Which substrate protects emergent diffeomorphism invariance? The rank-4 tensor decides, and a
crystal fails where an isotropic (icosahedral or amorphous) medium would not.

test_graviton_transversality left the gravity arc on a sharp, unexplained asymmetry. Emergent LORENTZ
invariance is IRRELEVANT-protected: the fermion cone's anisotropy is a dimension-6 operator that flows
away, which is why the model's one live prediction is a tiny (E/E_Planck)^2 effect and why Lorentz
symmetry emerges cleanly. Emergent DIFFEOMORPHISM invariance is only MARGINAL: the induced graviton's
rotational anisotropy holds a flat ratio as q -> 0 (a converged 12.4% split of the two transverse-
traceless polarisations) and does NOT flow away, which is why the realized theory is Nordstrom,
gamma = 0. Section 8.37 recorded the asymmetry; it did not explain it. This file does, and in doing so
identifies exactly what a substrate would need to change to make gamma = 1 reachable.

The explanation is one object: the substrate's rank-4 neighbour tensor A_ijkl = <n_i n_j n_k n_l>. Every
anisotropy the emergent fields can feel is built from the substrate's invariant tensors, and the
lowest anisotropic one is rank 4. The point of leverage is that the two fields couple to it at
DIFFERENT orders:

  * the FERMION (a scalar/spinor dispersion) feels A_ijkl only at O(k^4) in its energy -- a
    dimension-6, IRRELEVANT operator. So even on a crystal the fermion's Lorentz breaking flows away.
  * the GRAVITON (spin 2) contracts A_ijkl with its polarisation indices, so the anisotropy enters its
    TWO-derivative kinetic term directly -- a dimension-4, MARGINAL operator. On a crystal it does not
    flow away. That is the whole of the asymmetry: spin 2 promotes the same rank-4 anisotropy from
    irrelevant (for the fermion) to marginal (for the graviton).

So the obstruction to gamma = 1 is not the field theory (in the continuum it is Einstein, by Weinberg
-- see test_induced_transversality) but the substrate's rank-4 anisotropy. A substrate with NO rank-4
anisotropy removes it. Two exist:

  [A] THE RANK-4 OBSTRUCTION, and the graviton split it produces. Cubic is rank-4 anisotropic; the
      ICOSAHEDRAL shell is rank-4 EXACTLY isotropic (to machine precision), with its first anisotropy
      deferred to rank 6. The spin-2 h+/hx split, read straight off A_ijkl, is order-one for cubic and
      machine-zero for icosahedral -- the geometric origin of the 12.4%.
  [B] THE RG ORDER, measured on a tight-binding cone. Cubic dispersion anisotropy ~ k^2 (so the
      graviton's two-derivative term is anisotropic: MARGINAL); icosahedral ~ k^4 (deferred two orders:
      the graviton term is isotropic and anisotropy is dimension-6: IRRELEVANT, the same footing as
      emergent Lorentz).
  [C] AMORPHOUS works too. A statistically isotropic (random) medium has rank-4 anisotropy ~ 1/sqrt(N)
      -> 0 in the continuum. Space as a condensate is more naturally amorphous than crystalline, and an
      amorphous substrate protects diffeomorphism invariance the same way an icosahedral one does.

SCOPE, stated plainly. This identifies the substrate PROPERTY that controls gamma and shows the model's
crystalline substrate (triangular / cubic) is exactly the wrong one -- its rank-4 anisotropy is the
marginal diffeomorphism breaking measured in test_graviton_transversality. It does NOT build a working
emergent-Dirac-plus-graviton on an icosahedral or amorphous substrate; those have no Bloch theorem, and
that construction is a major open problem. What is established is the target: gamma = 0 is an artifact
of assuming a crystal, and an isotropic substrate defers the breaking to irrelevant order -- the one
concrete route by which the model's gravity could become Einstein.
"""
from __future__ import annotations
import numpy as np

PHI = (1 + np.sqrt(5)) / 2


def icosahedral_shell():
    """12 vertices of a regular icosahedron: cyclic perms of (0, +-1, +-phi), unit-normalised."""
    v = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            v += [(0, s1, s2 * PHI), (s1, s2 * PHI, 0), (s2 * PHI, 0, s1)]
    v = np.array(v, float)
    return v / np.linalg.norm(v, axis=1, keepdims=True)


CUBIC = np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)], float)
ICO = icosahedral_shell()


def rank_tensor(V, r):
    """A_{i...} = (1/N) sum_n n_i n_j ... (r indices)."""
    T = np.zeros((3,) * r)
    for v in V:
        term = v
        for _ in range(r - 1):
            term = np.multiply.outer(term, v)
        T += term
    return T / len(V)


def rank4_anisotropy(V):
    """||A - (isotropic projection)|| / ||A||, zero iff the rank-4 tensor is isotropic."""
    A = rank_tensor(V, 4)
    I = np.eye(3)
    iso = (np.einsum('ij,kl->ijkl', I, I) + np.einsum('ik,jl->ijkl', I, I) + np.einsum('il,jk->ijkl', I, I))
    c = np.sum(A * iso) / np.sum(iso * iso)
    return np.linalg.norm(A - c * iso) / np.linalg.norm(A)


def rank6_anisotropy(V):
    A = rank_tensor(V, 6)
    return abs(A[0, 0, 0, 0, 0, 0] - 5 * A[0, 0, 0, 0, 1, 1]) / abs(A[0, 0, 0, 0, 0, 0])


def graviton_split(V):
    """The spin-2 h+/hx anisotropy read off A_ijkl: the graviton's two-derivative rotational split."""
    A = rank_tensor(V, 4)
    hp = A[0, 0, 0, 0] + A[1, 1, 1, 1] - 2 * A[0, 0, 1, 1]     # h+ = (xx - yy)
    hx = 4 * A[0, 1, 0, 1]                                     # hx = (xy)
    return abs(hp - hx) / ((abs(hp) + abs(hx)) / 2 + 1e-30)


def dispersion_anisotropy(V, kmag, dirs):
    vals = np.array([np.sum(1 - np.cos(kmag * (d @ V.T))) for d in dirs])
    return (vals.max() - vals.min()) / vals.mean()


if __name__ == "__main__":
    print("=== Which substrate protects emergent diffeomorphism invariance? The rank-4 tensor decides ===\n")
    print("  Emergent Lorentz is IRRELEVANT-protected (fermion feels rank-4 anisotropy at O(k^4)); the")
    print("  graviton (spin 2) feels the SAME anisotropy in its two-derivative term -> MARGINAL. So the")
    print("  substrate's rank-4 tensor is the whole obstruction to gamma = 1. Remove it and it flows away.\n")

    # ---------- [A] the rank-4 obstruction and the graviton split ----------
    print("  [A] THE RANK-4 TENSOR, and the graviton anisotropy it produces:")
    print(f"      {'substrate':>14} {'rank-4 aniso':>13} {'rank-6 aniso':>13} {'graviton h+/hx split':>21}")
    for name, V in (("cubic", CUBIC), ("icosahedral", ICO)):
        print(f"      {name:>14} {rank4_anisotropy(V):>13.2e} {rank6_anisotropy(V):>13.2e} "
              f"{graviton_split(V):>21.2e}")
    print("      => cubic: rank-4 anisotropic, order-one graviton split -- this IS the marginal breaking")
    print("         test_graviton_transversality measured as 12.4%. Icosahedral: rank-4 isotropic to")
    print("         machine precision, graviton split machine-zero, first anisotropy only at rank 6.\n")

    # ---------- [B] the RG order, on a tight-binding cone ----------
    rng = np.random.default_rng(0)
    dirs = rng.normal(size=(400, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    ks = np.array([0.15, 0.25, 0.40, 0.60])
    print("  [B] THE RG ORDER, from a tight-binding dispersion eps(k) = sum_n [1 - cos(k.n)]:")
    print(f"      {'|k|':>6} {'cubic aniso':>13} {'slope':>7} {'ico aniso':>13} {'slope':>7}")
    ac = [dispersion_anisotropy(CUBIC, k, dirs) for k in ks]
    ai = [dispersion_anisotropy(ICO, k, dirs) for k in ks]
    for i, k in enumerate(ks):
        sc = np.log(ac[i] / ac[i - 1]) / np.log(ks[i] / ks[i - 1]) if i else np.nan
        si = np.log(ai[i] / ai[i - 1]) / np.log(ks[i] / ks[i - 1]) if i else np.nan
        print(f"      {k:>6.2f} {ac[i]:>13.2e} {sc:>7.2f} {ai[i]:>13.2e} {si:>7.2f}")
    print("      => cubic anisotropy ~ k^2: the graviton's two-derivative term is anisotropic = MARGINAL.")
    print("         icosahedral ~ k^4: deferred two orders, the two-derivative term is isotropic and the")
    print("         first anisotropy is dimension-6 = IRRELEVANT -- the same footing as emergent Lorentz.\n")

    # ---------- [C] amorphous substrate ----------
    print("  [C] AMORPHOUS (statistically isotropic) substrate: rank-4 anisotropy -> 0 in the continuum.")
    print(f"      {'N points':>10} {'rank-4 anisotropy':>18}")
    rng2 = np.random.default_rng(1)
    for N in (50, 200, 800, 3200):
        W = rng2.normal(size=(N, 3))
        W /= np.linalg.norm(W, axis=1, keepdims=True)
        print(f"      {N:>10} {rank4_anisotropy(W):>18.2e}")
    print("      => falls as ~1/sqrt(N). A random isotropic medium has no rank-4 anisotropy in the")
    print("         continuum limit, so it protects diffeomorphism invariance like the icosahedron does.")
    print("         Space as a condensate is more naturally amorphous than crystalline.\n")

    print("[verdict] gamma = 0 is an artifact of the CRYSTALLINE substrate, not of the field theory:")
    print("  * Every anisotropy the emergent fields feel is built from the substrate's rank-4 tensor.")
    print("    The fermion feels it at O(k^4) (dimension-6, IRRELEVANT), so emergent Lorentz works even")
    print("    on a crystal. The graviton, being spin 2, feels the SAME tensor in its two-derivative")
    print("    term (dimension-4, MARGINAL), so emergent diffeomorphism fails on a crystal. That single")
    print("    fact -- spin 2 promoting rank-4 anisotropy from irrelevant to marginal -- is the entire")
    print("    asymmetry Section 8.37 recorded but did not explain.")
    print("  * The model's substrate is a crystal (triangular / cubic), with an order-one rank-4")
    print("    anisotropy: exactly the marginal breaking test_graviton_transversality measured. In the")
    print("    continuum the induced action is Einstein (Weinberg, test_induced_transversality); it is")
    print("    the crystal that spoils it.")
    print("  * A substrate with NO rank-4 anisotropy defers the breaking to rank 6 = dimension-6 =")
    print("    IRRELEVANT, putting emergent diffeomorphism on the same footing as emergent Lorentz.")
    print("    Two exist: the ICOSAHEDRAL shell (exactly isotropic rank-4, measured to 1e-16) and any")
    print("    AMORPHOUS medium (rank-4 anisotropy ~ 1/sqrt(N) -> 0). Either is a concrete substrate on")
    print("    which the model's gravity could become Einstein.")
    print("  * SCOPE. This identifies the property and the target; it does not build the emergent Dirac")
    print("    and graviton on such a substrate -- neither has a Bloch theorem, and that is the open")
    print("    construction. What is settled: gamma = 0 follows from assuming a crystal, and an isotropic")
    print("    substrate is the one concrete route by which it could become gamma = 1.")
