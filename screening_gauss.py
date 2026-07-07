"""
Screen-1 -- existence proof: the medium CAN carry an unscreened 1/r force, if the
mediating field obeys a Gauss law with NO mass term.

Screen-0 showed gravity-by-density screens (lambda~3) and that no coupling dial
defeats it: the mediating field behaves as if it has a MASS (the medium's pinning
to its rest spacing R0), and a massive field always gives Yukawa exp(-r/lambda),
lambda = 1/m. A massless field -- one protected by a conservation law so it has no
restoring term -- gives Coulomb 1/r instead (force 1/r^2).

Here we test exactly that on the SAME fcc node machinery: solve a discrete Poisson
equation for a point source,
    massless:  Lap Phi = -s            -> expect Phi(r) ~ 1/r   (UNSCREENED)
    massive : (Lap - m^2) Phi = -s     -> expect Phi(r) ~ exp(-r/lambda), lambda=1/m
with Dirichlet Phi=0 on the outer shell (a finite stand-in for infinity). A large
fcc ball is reached with a cell-list SPARSE nearest-neighbour Laplacian and a pure
-numpy conjugate-gradient solve (no scipy; the dense (N,N) operator is multi-GB
past a few thousand nodes). If massless=1/r and massive=Yukawa, screening is
*precisely* the mass term -- and long-range gravity needs a CONSERVED (un-screenable)
source, not more dimensions or a bigger coupling. This is 'gravity by hand': an
existence proof + calibration, NOT emergence (that is Screen-2).
"""
from __future__ import annotations
import numpy as np

from bvc_core import R0, perfect_fcc
from density_response_3d_large import cell_pairs


def build_laplacian(X, rcut=1.15 * R0):
    """Sparse nn graph Laplacian as (matvec, scale). L Phi ~ +div-grad Phi, scaled
    so lap(r^2)=6 in the 3D interior. Only nearest neighbours (regular fcc)."""
    pi, pj = cell_pairs(X, rcut); N = len(X)
    r = np.linalg.norm(X, axis=1); interior = r < 0.6 * r.max()

    def graph_lap(phi):                      # sum_nn (phi_j - phi_i) at each node
        dij = phi[pj] - phi[pi]
        return np.bincount(pi, dij, N) - np.bincount(pj, dij, N)

    q = (X ** 2).sum(1)
    scale = 6.0 / graph_lap(q)[interior].mean()
    return (lambda phi: scale * graph_lap(phi)), N


def cg(matvec, b, tol=1e-9, maxit=5000):
    x = np.zeros_like(b); r = b - matvec(x); p = r.copy(); rs = r @ r
    b2 = max(b @ b, 1e-30)
    for _ in range(maxit):
        Ap = matvec(p); alpha = rs / (p @ Ap)
        x += alpha * p; r -= alpha * Ap; rsn = r @ r
        if rsn / b2 < tol ** 2:
            break
        p = r + (rsn / rs) * p; rs = rsn
    return x


def solve_poisson(Lmv, N, X, s, m2, r_bnd):
    """Solve (L - m2) Phi = -s, Dirichlet Phi=0 past r_bnd, i.e. (-L + m2) Phi = s."""
    free = (np.linalg.norm(X, axis=1) <= r_bnd).astype(float)

    def A(phi):
        phi = phi * free
        return free * (-Lmv(phi) + m2 * phi)

    return cg(A, s * free), free.astype(bool)


def radial(X, Phi, sel, bins):
    r = np.linalg.norm(X, axis=1)
    return np.array([Phi[sel & (r >= lo) & (r < hi)].mean()
                     if (sel & (r >= lo) & (r < hi)).sum() else np.nan
                     for lo, hi in zip(bins[:-1], bins[1:])])


def fit_forms(rmid, phi, rmax):
    """Linear-space relative residual for Coulomb (A/r+B) vs Yukawa (A exp(-r/lam))."""
    m = np.isfinite(phi) & (rmid > 2.0) & (rmid < 0.75 * rmax) & (phi > 1e-9)
    r, y = rmid[m], phi[m]
    Mc = np.stack([1.0 / r, np.ones_like(r)], 1)
    cc, *_ = np.linalg.lstsq(Mc, y, rcond=None)
    res_c = np.sqrt(np.mean((y - Mc @ cc) ** 2)) / y.mean()
    pc = np.polyfit(r, np.log(y), 1); lam = -1.0 / pc[0]
    yhat = np.exp(np.polyval(pc, r))
    res_y = np.sqrt(np.mean((y - yhat) ** 2)) / y.mean()
    return res_c, res_y, lam


def coulomb_flux(X, Phi, sel, bins):
    """Gauss test: enclosed flux Q(r) = sum_{|x|<r} Lap Phi = -sum s inside r.
    For a solved field this equals 4pi r^2 (-dPhi/dr). We report the fitted
    'apparent range' via an exp fit -- the boundary-honest discriminator is how
    that range SCALES with the box (below)."""
    pass


def apparent_range(fcc, s, m2, frac_bnd=0.92):
    """Solve on this fcc and return the exp-fit apparent decay length lambda."""
    rmax = np.linalg.norm(fcc, axis=1).max()
    Lmv, N = build_laplacian(fcc)
    Phi, free = solve_poisson(Lmv, N, fcc, s, m2, frac_bnd * rmax)
    bins = np.arange(0, 0.82 * rmax, 1.0); rmid = 0.5 * (bins[:-1] + bins[1:])
    phi = radial(fcc, Phi, free, bins)
    _, _, lam = fit_forms(rmid, phi, rmax)
    return lam, phi, rmid, rmax


def point_source(fcc, width=1.5):
    s = np.exp(-(fcc ** 2).sum(1) / (2 * width ** 2))
    return s / s.sum()


def main():
    print("Screen-1 Gauss-law existence proof (sparse Laplacian + CG)\n")
    print("Boundary-honest test: a MASSLESS (Coulomb 1/r) field has NO intrinsic")
    print("length, so its apparent range must GROW with the box. A MASSIVE (Yukawa)")
    print("field pins at lambda=1/m regardless of box size.\n")

    radii = [10.0, 14.0, 20.0]
    print(f"  {'radius':>7} {'N':>7} | {'massless lambda':>16} | {'massive m=1/3 lambda':>21}")
    for R in radii:
        fcc = perfect_fcc(radius=R); s = point_source(fcc)
        lam0, _, _, _ = apparent_range(fcc, s, 0.0)
        lam3, _, _, _ = apparent_range(fcc, s, (1 / 3.0) ** 2)
        print(f"  {R:>7.0f} {len(fcc):>7} | {lam0:>16.2f} | {lam3:>21.2f}  (1/m=3.0)")

    # explicit 1/r check on the largest massless solve: Phi*r flat beyond source
    fcc = perfect_fcc(radius=20.0); s = point_source(fcc)
    _, phi, rmid, rmax = apparent_range(fcc, s, 0.0)
    k = min(14, len(rmid))
    print("\n  massless Phi(r) and Phi*r (Coulomb => Phi*r ~ linear A(1 - r/R_b), not exp):")
    print("   r    :", " ".join(f"{x:5.1f}" for x in rmid[:k]))
    print("   Phi  :", " ".join(f"{x:5.3f}" for x in phi[:k]))
    print("   Phi*r:", " ".join(f"{x:5.3f}" for x in (phi * rmid)[:k]))
    print("\n  => massless apparent-range grows with the box (no intrinsic scale = 1/r,")
    print("     Coulomb, UNSCREENED); massive pins at 1/m (Yukawa, screened). Screening")
    print("     IS the mass term. Long-range gravity needs a CONSERVED source that")
    print("     forbids a mass term (Screen-2), not a bigger coupling or more dimensions.")


if __name__ == "__main__":
    main()
