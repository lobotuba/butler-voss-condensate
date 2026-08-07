"""Is the induced graviton anisotropy the substrate's rank-4 tensor? Measured directly -- and it is not.

Section 8.58 (test_diffeo_substrate) argued that the marginal diffeomorphism breaking of Section 8.37 --
the converged 12.4% split of the induced tetrad-graviton's two transverse-traceless polarisations -- is
the substrate's rank-4 neighbour tensor A_ijkl = <n_i n_j n_k n_l>, and concluded that a substrate with
NO rank-4 anisotropy (icosahedral, or amorphous) would remove it. That section measured only the
substrate GEOMETRY (the rank-4 tensor of the bare neighbour vectors and a small-k tight-binding
dispersion). It never measured the induced GRAVITON on a substrate whose rank-4 tensor had been tuned
away. This file does, on the exact instrument that produced the 12.4%, and the specific mechanism fails.

The instrument is the sea-energy method of test_graviton_transversality (z in real space, xy in k-space,
the metric entering through the tetrad), generalised so the Wilson-Dirac hopping runs over an arbitrary
neighbour shell. The knob is the weight w of a (110)-type diagonal shell added to the nearest-neighbour
cubic shell. This keeps the lattice periodic (Bloch holds, the instrument is unchanged) while sweeping
the effective shell's rank-4 anisotropy through EXACT zero: the cubic shell has A_1111 = 2, A_1122 = 0;
the (110) shell has A_1111 = 8, A_1122 = 4; the combination is rank-4 isotropic (A_1111 = 3 A_1122) at
w* = 1/2. Section 8.58's claim, tested directly, is that the induced h+/hx split should collapse at w*.

  [A] CALIBRATION. Pure nearest-neighbour cubic reproduces Section 8.37's converged split (~12.4%),
      confirming the generalised instrument is the same measurement.

  [B] THE SWEEP -- the refutation. As w runs 0 -> large, the substrate's rank-4 anisotropy DIPS to
      exact zero at w* = 1/2 and climbs back. The measured induced split does NOTHING of the kind: it
      falls MONOTONICALLY through w*, with no feature there, and its smallest values occur where the
      rank-4 anisotropy is LARGE. At w* (rank-4 identically zero) the split is still ~8%; at w = 2
      (rank-4 anisotropy larger than the cubic value) it is under 1%. A single point settles it: rank-4
      exactly zero, split manifestly nonzero. No weighting of the neighbour rank-4 tensor can track a
      monotone curve with a non-monotone one.

  [C] CONVERGENCE. The ~8% residual at w* is a property of the theory, not the grid: it converges in the
      transverse momentum count and is flat in q, the controls Section 8.37 insisted on.

WHAT CONTROLS IT, and why this UNIFIES 8.57 and 8.58. The induced two-derivative graviton term is a
whole-Brillouin-zone integral, and test_induced_transversality (Section 8.57) already showed its
MAGNITUDE is a Sakharov term fixed at the substrate (cutoff) scale, with no long-wavelength finite part.
Its ANISOTROPY is the same kind of object: a whole-zone, substrate-scale quantity, not the small-k
rank-4 neighbour moment. That is exactly why tuning a single low-order moment to zero does nothing --
the small-k expansion is the wrong scale -- and why adding more hopping directions (large w), which makes
the WHOLE dispersion sample the sphere more evenly, is what actually shrinks the split.

SCOPE, and what survives. This corrects Section 8.58's MECHANISM: the induced graviton anisotropy is not
the substrate's rank-4 neighbour tensor, and a crystal tuned to rank-4 isotropy is still Nordstrom. Its
CONCLUSION survives and is sharpened. Removing a whole-zone anisotropy needs a substrate isotropic at
ALL orders, not at rank 4 alone -- which is precisely a statistically isotropic AMORPHOUS medium (no
Brillouin zone, no preferred direction at any scale), the naturally-condensate substrate. The icosahedral
route is now expected to be imperfect, not exact: its rank-4 tensor is isotropic but its whole-zone
integral carries its rank-6 anisotropy, so it would leave a small residual split rather than the machine
zero Section 8.58 read off the geometry. The one clean route to gamma = 1 is an amorphous substrate, and
the reason is stronger than 8.58 stated: the obstruction lives at the substrate scale, where 8.57 put it.
"""
from __future__ import annotations
import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
AL = [np.block([[Z2, s], [s, Z2]]) for s in (sx, sy, sz)]     # alpha_x, alpha_y, alpha_z
BE = np.block([[I2, Z2], [Z2, -I2]])                          # beta
M0, RW = 0.6, 1.0


def tetrad(h):
    """Symmetric square root of the inverse metric g = 1 + h (removes the O(h^2) remapping)."""
    w, U = np.linalg.eigh(np.eye(3) + h)
    return (U * (1.0 / np.sqrt(w))) @ U.T


def energy(Nperp, Nz, shells, M=M0, hfield=None):
    """Filled-sea energy per site for a Wilson-Dirac fermion whose hopping runs over `shells`.

    shells : list of (n, t_cone, t_wilson), n an integer 3-vector (one representative per +/- pair).
    hfield : h_ij(z) background; each hop direction n rotates as n -> V(z) n through the tetrad, so a
             diagonal-shell hop feels the metric exactly as the nearest-neighbour cone does.
    """
    gp = (np.arange(Nperp) + 0.5) / Nperp * 2 * np.pi
    KX, KY = np.meshgrid(gp, gp, indexing="ij")
    kx, ky = KX.ravel(), KY.ravel()
    P = kx.size
    V = (np.broadcast_to(np.eye(3), (Nz, 3, 3)) if hfield is None
         else np.array([tetrad(hfield[z]) for z in range(Nz)]))

    H = np.zeros((P, 4 * Nz, 4 * Nz), dtype=complex)
    onsite = M + sum(tw for (_, _, tw) in shells)         # each (1-cos) shell rep adds tw on-site
    for z in range(Nz):
        H[:, 4 * z:4 * z + 4, 4 * z:4 * z + 4] += onsite * BE

    for (n, tc, tw) in shells:
        nx, ny, nz = n
        phase = kx * nx + ky * ny                          # k.n restricted to the periodic (x,y) part
        for z in range(Nz):
            b = 4 * z
            Vmid = V[z] if nz == 0 else 0.5 * (V[z] + V[(z + nz) % Nz])
            Vn = Vmid @ np.array(n, float)
            A = Vn[0] * AL[0] + Vn[1] * AL[1] + Vn[2] * AL[2]   # sum_i alpha_i (V n)_i
            if nz == 0:
                H[:, b:b + 4, b:b + 4] += tc * np.sin(phase)[:, None, None] * A[None]
                H[:, b:b + 4, b:b + 4] += (-tw) * np.cos(phase)[:, None, None] * BE[None]
            else:
                bp = 4 * ((z + nz) % Nz)
                ph = np.exp(1j * phase)[:, None, None]
                Mhop = (tc * A / (2j))[None] * ph + (-0.5 * tw * BE)[None] * ph
                H[:, b:b + 4, bp:bp + 4] += Mhop
                H[:, bp:bp + 4, b:b + 4] += np.conj(np.transpose(Mhop, (0, 2, 1)))
    w = np.linalg.eigvalsh(H)
    return float(np.sum(w[w < 0])) / (P * Nz)


def curv(f, st=1e-2):
    a, b, c, d, g = f(-2 * st), f(-st), f(0.0), f(st), f(2 * st)
    return (-a + 16 * b - 30 * c + 16 * d - g) / (12 * st ** 2)


TThp = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], float)     # transverse-traceless, on the axes
TThx = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], float)      # transverse-traceless, at 45 deg


def response(Nperp, Nz, nq, shells, e):
    """Second-order sea-energy response to h_ij(z) = e cos(qz), per q^2, q = 0 mass term removed."""
    q = 2 * np.pi * nq / Nz
    prof = np.cos(q * np.arange(Nz))
    pq = curv(lambda t: energy(Nperp, Nz, shells, hfield=t * prof[:, None, None] * e))
    p0 = curv(lambda t: energy(Nperp, Nz, shells, hfield=t * np.ones(Nz)[:, None, None] * e))
    return (pq - 0.5 * p0) / q ** 2


CUBIC_REPS = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
DIAG_REPS = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]


def shells_for(w):
    s = [(n, 1.0, RW) for n in CUBIC_REPS]
    if w != 0.0:
        s += [(n, w, w) for n in DIAG_REPS]
    return s


def rank4_aniso(w):
    """Cubic anisotropy ||A - iso||/||A|| of the cone-weighted rank-4 hopping tensor. Zero at w = 1/2."""
    A = np.zeros((3, 3, 3, 3))
    for reps, t in ((CUBIC_REPS, 1.0), (DIAG_REPS, w)):
        for n in reps:
            for s1 in (1, -1):
                v = s1 * np.array(n, float)
                A += t * np.einsum('i,j,k,l->ijkl', v, v, v, v)
    I = np.eye(3)
    iso = (np.einsum('ij,kl->ijkl', I, I) + np.einsum('ik,jl->ijkl', I, I) + np.einsum('il,jk->ijkl', I, I))
    c = np.sum(A * iso) / np.sum(iso * iso)
    return np.linalg.norm(A - c * iso) / np.linalg.norm(A)


def split(Nperp, Nz, nq, w):
    sh = shells_for(w)
    p = response(Nperp, Nz, nq, sh, TThp)
    x = response(Nperp, Nz, nq, sh, TThx)
    return p, x, abs(p - x) / abs(p)


if __name__ == "__main__":
    print("=== Is the induced graviton anisotropy the substrate rank-4 tensor? Measured -- it is not ===\n")
    print("  Instrument: the sea-energy graviton of Section 8.37, with the Wilson-Dirac hopping run over")
    print("  a tunable shell. Adding a (110) shell of weight w sweeps the effective rank-4 anisotropy")
    print("  through EXACT zero at w* = 1/2, while the lattice stays periodic. Section 8.58 predicts the")
    print("  induced h+/hx split collapses at w*.\n")
    Nperp = 14

    # ---------- [A] calibration against Section 8.37 ----------
    print("  [A] CALIBRATION: pure nearest-neighbour cubic must reproduce Section 8.37's ~12.4% split.")
    print(f"      {'q':>7} {'TT h+':>11} {'TT hx':>11} {'anisotropy':>11}")
    for (nz, nq) in ((24, 2), (36, 1), (48, 1)):
        p, x, s = split(Nperp, nz, nq, 0.0)
        print(f"      {2*np.pi*nq/nz:>7.4f} {p:>11.6f} {x:>11.6f} {s:>10.2%}")
    print("      => converged ~12.4%, flat in q: the generalised instrument is the same measurement.\n")

    # ---------- [B] the sweep: the refutation ----------
    print("  [B] THE SWEEP. Rank-4 anisotropy dips to EXACT zero at w* = 1/2 and climbs back; the")
    print("      measured split falls MONOTONICALLY through w*, with no feature, smallest where rank-4")
    print("      is LARGE. (Section 8.58 predicted a collapse AT w*.)")
    print(f"      {'w':>6} {'rank-4 aniso':>13} {'TT h+':>11} {'TT hx':>11} {'split':>8}")
    for w in (0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0):
        p, x, s = split(Nperp, 36, 1, w)
        tag = "  <- rank-4 = 0" if w == 0.5 else ""
        print(f"      {w:>6.2f} {rank4_aniso(w):>13.4f} {p:>11.6f} {x:>11.6f} {s:>7.2%}{tag}")
    print("      => at w* (rank-4 identically zero) the split is still large; at w = 2, where the shell")
    print("         is MORE rank-4-anisotropic than the cubic lattice, the split is far smaller. The")
    print("         neighbour rank-4 tensor does not control the induced graviton anisotropy.\n")

    # ---------- [C] convergence control at w* ----------
    print("  [C] CONVERGENCE at w* = 1/2: the residual split is the theory's, not the grid's.")
    print(f"      {'Nperp':>6} {'Nz':>4} {'nq':>3} {'q':>7} {'split':>8}")
    for Np in (10, 20, 40):
        _, _, s = split(Np, 24, 2, 0.5)
        print(f"      {Np:>6} {24:>4} {2:>3} {2*np.pi*2/24:>7.4f} {s:>7.2%}")
    for (Nz, nq) in ((24, 2), (36, 1), (48, 1)):
        _, _, s = split(20, Nz, nq, 0.5)
        print(f"      {20:>6} {Nz:>4} {nq:>3} {2*np.pi*nq/Nz:>7.4f} {s:>7.2%}")
    print("      => converges to ~8% and is flat in q. Rank-4 isotropy leaves a large split.\n")

    print("[verdict] The induced graviton anisotropy is a whole-Brillouin-zone (substrate-scale)")
    print("  quantity, not the substrate's rank-4 neighbour tensor. This corrects Section 8.58's")
    print("  mechanism and unifies it with Section 8.57:")
    print("  * DIRECT REFUTATION. Tuning the neighbour rank-4 tensor to EXACT zero (w* = 1/2) leaves the")
    print("    induced h+/hx split near 8%; a shell that is MORE rank-4-anisotropic (w = 2) has a split")
    print("    under 1%. The split is monotone in w while the rank-4 anisotropy is not -- no weighting of")
    print("    the neighbour rank-4 tensor can track it. Section 8.58 measured the substrate geometry and")
    print("    assumed the graviton inherited it; measured on the induced action, it does not.")
    print("  * WHY. The induced two-derivative term is a whole-zone integral. Section 8.57 showed its")
    print("    MAGNITUDE is a Sakharov term fixed at the substrate scale with no long-wavelength part;")
    print("    its ANISOTROPY is the same kind of object. A small-k moment is the wrong scale, which is")
    print("    why nulling it does nothing, and why adding hopping directions -- making the WHOLE")
    print("    dispersion sample the sphere evenly -- is what actually shrinks the split.")
    print("  * WHAT SURVIVES. Section 8.58's CONCLUSION holds and is sharpened: removing a whole-zone")
    print("    anisotropy needs a substrate isotropic at ALL orders, i.e. an AMORPHOUS (statistically")
    print("    isotropic, no Brillouin zone) medium -- the naturally-condensate substrate. Even an")
    print("    icosahedral quasicrystal is now expected to leave a small rank-6 residual, not the machine")
    print("    zero 8.58 read off the geometry. The route to gamma = 1 is amorphous, for the stronger")
    print("    reason that the obstruction lives at the substrate scale, exactly where 8.57 put it.")
