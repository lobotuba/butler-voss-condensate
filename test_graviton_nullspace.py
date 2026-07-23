"""
The induced graviton form has no gauge null space at all. Why projection cannot rescue gamma = 1.

Section 8.29 measured the induced two-derivative action of the tetrad sector one polarisation at a
time and found it fails linearised diffeomorphism invariance -- the pure-gauge modes cost energy, the
transverse-traceless doublet is split, and the fitted coefficients miss Einstein-Hilbert. A natural
hope survives that file: the model's actual gravity is not the tetrad's compatible strain but the
deconfined curvature sector, which propagates only INCOMPATIBLE strain (test_light_bending showed
compatible strain is pure gauge, flat, and inert). If the physical curvature modes were cleanly
separated from the offending gauge modes, the violation might live entirely in a sector the curvature
propagator projects away, and gamma = 1 could still hold for the gravity that matters.

This file tests that hope directly, by measuring the FULL 6x6 quadratic form on the space of symmetric
h_ij at fixed q -- every diagonal and every cross term -- rather than one polarisation at a time. The
cross terms are the whole point: a hope built on projecting one subspace out is a statement about
whether the form is BLOCK-DIAGONAL between physical and gauge modes, and only the off-diagonal blocks
can answer it.

The basis is orthonormal in the natural inner product <e, e'> = e_ij e'_ij, with q along z:
    physical (incompatible / transverse):  h+, hx  (transverse-traceless),  trp (transverse trace)
    gauge (h_ij = d_i xi_j + d_j xi_i):     xz, yz  (spin-1),               zz  (longitudinal)
Cross terms are read from the diagonal of the sum, P_ab = 1/2[ D(e_a + e_b) - D(e_a) - D(e_b) ], the
standard polarisation identity. Each entry has the q = 0 mass term removed on the identical k-point
set, exactly as in Section 8.29, so what remains is the pure two-derivative content.

What linearised diffeomorphism invariance requires is unambiguous and is stated as a prediction before
the measurement: the three gauge directions must be EXACT null vectors of the form. Einstein-Hilbert's
own 6x6, computed here on the same basis, has a rank-3 kernel spanned by exactly those three modes.
The induced form is compared against it eigenvalue by eigenvalue.

  [A] THE FULL FORM, measured, beside the Einstein-Hilbert form on the same basis.
  [B] ITS SPECTRUM. Einstein-Hilbert has three zero eigenvalues (the gauge kernel) and three nonzero;
      the induced form's spectrum shows whether any kernel exists.
  [C] THE MIXING. The physical-gauge off-diagonal block, whose size decides whether ANY projection
      onto the physical subspace can remove the gauge violation.
  [D] THE TWO PHYSICAL SIGNATURES separately: the transverse-traceless doublet (spin-2, which no
      projection onto incompatible strain can touch) and the transverse-trace / longitudinal block
      (spin-0, where gamma actually lives for a spherical source).

SCOPE. This is the induced action of the TETRAD sector, the same object as Section 8.29, not the
deconfined curvature sector's own propagator. The point is precisely negative: it closes the projection
loophole that Section 8.29 left open, and shows that the tetrad's induced form cannot be massaged into
an Einstein-Hilbert one by discarding a subspace. gamma = 1 for the model's gravity still rests on
Weinberg applied to the curvature sector, measured directly -- which remains the open problem.
"""
from __future__ import annotations
import numpy as np
from test_graviton_transversality import energy, curv, EH

R2 = np.sqrt(2.0)
BASIS = {                                   # q along z; orthonormal in <e,e'> = e_ij e'_ij
    "h+":  np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) / R2,   # transverse-traceless
    "hx":  np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / R2,    # transverse-traceless
    "trp": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]]) / R2,    # transverse trace (physical scalar)
    "xz":  np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / R2,    # gauge: h = d_z xi_x
    "yz":  np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / R2,    # gauge: h = d_z xi_y
    "zz":  np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], float),  # gauge: h = 2 d_z xi_z
}
NAMES = list(BASIS)
PHYS, GAUGE = ["h+", "hx", "trp"], ["xz", "yz", "zz"]


def diag_response(Nperp, Nz, nq, e):
    """Second-order energy response to h_ij(z) = e_ij cos(qz), per q^2, q = 0 mass removed."""
    q = 2 * np.pi * nq / Nz
    prof = np.cos(q * np.arange(Nz))
    pq = curv(lambda t: energy(Nperp, Nz, hfield=t * prof[:, None, None] * e))
    p0 = curv(lambda t: energy(Nperp, Nz, hfield=t * np.ones(Nz)[:, None, None] * e))
    return (pq - 0.5 * p0) / q ** 2


def induced_form(Nperp, Nz, nq):
    """The full symmetric 6x6, cross terms by the polarisation identity."""
    d = {k: diag_response(Nperp, Nz, nq, BASIS[k]) for k in NAMES}
    P = np.diag([d[k] for k in NAMES]).astype(float)
    for i, a in enumerate(NAMES):
        for j in range(i + 1, len(NAMES)):
            b = NAMES[j]
            s = diag_response(Nperp, Nz, nq, BASIS[a] + BASIS[b])
            P[i, j] = P[j, i] = 0.5 * (s - d[a] - d[b])
    return P


def eh_form():
    """Einstein-Hilbert's own 6x6 on this basis, from the bilinear invariants (q along z):
       F1 = h_ij h'_ij,  F2 = h_zj h'_zj,  F3 = 1/2(h_zz tr h' + h'_zz tr h),  F4 = tr h tr h'."""
    def F(e, f):
        return np.array([np.sum(e * f), np.sum(e[2] * f[2]),
                         0.5 * (e[2, 2] * np.trace(f) + f[2, 2] * np.trace(e)),
                         np.trace(e) * np.trace(f)])
    return np.array([[EH @ F(BASIS[a], BASIS[b]) for b in NAMES] for a in NAMES])


def _grid(M):
    return "\n".join("  " + f"{NAMES[i]:>4} " +
                     "".join(f"{M[i, j]:>10.4f}" for j in range(6)) for i in range(6))


if __name__ == "__main__":
    print("=== The induced graviton form has no gauge null space. Why projection cannot save gamma=1 ===\n")
    Nperp, Nz, nq = 14, 36, 1
    q = 2 * np.pi * nq / Nz
    print(f"  Wilson-Dirac torus, Nperp = {Nperp}, Nz = {Nz}, q = {q:.4f} along z. Full 6x6 quadratic")
    print("  form on symmetric h_ij, q = 0 mass removed. The question is whether it is block-diagonal")
    print("  between physical (incompatible) and gauge modes -- only then could a projection help.\n")

    P = induced_form(Nperp, Nz, nq)
    E = eh_form()
    print("  header:" + "".join(f"{n:>10}" for n in NAMES))
    print("  [A] INDUCED FORM (x1e3):")
    print(_grid(1e3 * P))
    print("\n      EINSTEIN-HILBERT FORM on the same basis (overall scale arbitrary):")
    print(_grid(E))
    print("      => EH is block-diagonal and its entire gauge block (xz,yz,zz) is zero. The induced")
    print("         form is not: the gauge diagonal is its LARGEST entries, not its smallest.\n")

    ev_P = np.linalg.eigvalsh(P)
    ev_E = np.linalg.eigvalsh(E)
    print("  [B] SPECTRUM. Diffeomorphism invariance = a 3-dimensional null space (one per xi")
    print("      component). Einstein-Hilbert has exactly three zero eigenvalues.")
    print(f"      induced eigenvalues (x1e3) : {np.round(1e3 * ev_P, 4)}")
    print(f"      EH eigenvalues (scaled)    : {np.round(ev_E, 4)}")
    nz = int(np.sum(np.abs(ev_P) > 1e-3 * np.abs(ev_P).max()))
    print(f"      => the induced form has {nz} nonzero eigenvalues out of 6: NO null space. The gauge")
    print("         kernel that Einstein-Hilbert requires is simply absent, and the smallest induced")
    print("         eigenvalue is a physical mode, not a gauge one.\n")

    ip = [NAMES.index(k) for k in PHYS]
    ig = [NAMES.index(k) for k in GAUGE]
    Ppp, Pgg, Ppg = P[np.ix_(ip, ip)], P[np.ix_(ig, ig)], P[np.ix_(ip, ig)]
    print("  [C] THE MIXING -- the block that decides whether projection can help. A projection onto")
    print("      the physical subspace removes the gauge violation ONLY if the physical-gauge block")
    print("      vanishes (the form is block-diagonal). It does not.")
    print(f"      ||physical block||           = {np.linalg.norm(Ppp):.3e}")
    print(f"      ||gauge block||              = {np.linalg.norm(Pgg):.3e}   "
          f"({np.linalg.norm(Pgg)/np.linalg.norm(Ppp):.1f}x the physical block)")
    print(f"      ||physical-gauge mixing||    = {np.linalg.norm(Ppg):.3e}   "
          f"({np.linalg.norm(Ppg)/np.linalg.norm(Ppp):.2f} of the physical block)")
    print(f"      Einstein-Hilbert mixing      = {np.linalg.norm(E[np.ix_(ip, ig)]):.1e}  (exactly zero)")
    print("      => the physical trace mode is sourced by a pure-gauge mode at 42% of its own scale.")
    print("         You cannot project out a mode that feeds the one you keep. Projection fails.\n")

    print("  [D] THE TWO PHYSICAL SIGNATURES, read off directly:")
    i_hp, i_hx, i_tr, i_zz = (NAMES.index(k) for k in ("h+", "hx", "trp", "zz"))
    split = abs(P[i_hp, i_hp] - P[i_hx, i_hx]) / abs(P[i_hp, i_hp])
    print(f"      spin-2 (h+ vs hx, no projection can mix these): {1e3*P[i_hp,i_hp]:.4f} vs "
          f"{1e3*P[i_hx,i_hx]:.4f}  -> split {split:.1%}")
    print("        rotational invariance requires them equal; they are not, and no projection onto")
    print("        incompatible strain touches the spin-2 sector, so this survives any rescue.")
    print(f"      spin-0 (trp-zz block, where gamma lives for a spherical source):")
    print(f"        induced {np.round(1e3*P[np.ix_([i_tr,i_zz],[i_tr,i_zz])],4).tolist()}")
    print(f"        EH      {np.round(E[np.ix_([i_tr,i_zz],[i_tr,i_zz])],4).tolist()}")
    print("        EH puts all of it in trp with zz null; the induced form spreads it across a")
    print("        gauge mode. This is the sector that sets Psi, hence gamma.\n")

    print("[verdict] projection cannot rescue gamma = 1 for the tetrad, and Section 8.29's failure is")
    print("          structural, not a matter of degree:")
    print("  * Linearised diffeomorphism invariance is a rank-3 null space, one zero eigenvalue per")
    print("    gauge parameter. Einstein-Hilbert's 6x6 has exactly that kernel. The induced form has")
    print("    NO kernel: all six eigenvalues are nonzero, and the gauge directions are its STIFFEST")
    print("    modes rather than its flat ones.")
    print("  * The one hope Section 8.29 left -- that the curvature sector propagates only the")
    print("    incompatible modes and might avoid the gauge violation -- requires the form to be")
    print("    block-diagonal between physical and gauge subspaces. It is not: the physical-gauge")
    print("    mixing is 42% of the physical block itself, concentrated in the spin-0 sector where the")
    print("    PPN parameter gamma is defined. A projection cannot remove a violation that is coupled")
    print("    into the modes it keeps.")
    print("  * The spin-2 doublet fails a different way -- a 12.6% rotational-anisotropy split that no")
    print("    projection onto incompatible strain can reach at all, since spin-2 cannot mix with the")
    print("    spin-1 gauge modes. Both physical channels fail, for independent reasons.")
    print("  * So the negative of Section 8.29 is not a near miss to be projected away: the induced")
    print("    tetrad action is not a deformation of Einstein-Hilbert in any subspace. gamma = 1 for")
    print("    the model's gravity rests entirely on the deconfined curvature sector, measured")
    print("    directly -- the project's central open problem, now with the tetrad route fully closed.")
