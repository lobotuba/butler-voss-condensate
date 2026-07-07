"""
Phase 3D-2 -- sparse scale-up: is 3D gravity long-range (power law) or screened?

3D-1 showed the ordered-fcc density response barely decays out to r~7 (lambda~14),
but radius-9 was too small to tell long-lambda / power-law / finite-size apart.
Here we scale to a large fcc medium with a CELL-LIST sparse force (O(N*nbrs)), so
we can measure Drho(r) far from both the source and the free surface.

A gently strained lattice barely changes neighbours, so the pair list is built
ONCE from the perfect fcc and reused through the (cooled) relaxation.
Sparse forces are verified against the dense engine at small N first.
"""
from __future__ import annotations
import math
from collections import defaultdict
import numpy as np

from bvc_core import R0, EPS, SIGMA, RCUT, perfect_fcc, pairwise
import integration_phase3_variational as V

H, RCUT_FIELD, GMIN, GMAX, BETA = R0, 2.0 * R0, 0.02, 1.95, 60.0


# ------------------------------------------------------ cell-list pair finder --
def cell_pairs(X, cutoff):
    lo = X.min(0)
    cell = np.floor((X - lo) / cutoff).astype(np.int64)
    buckets = defaultdict(list)
    for idx in range(len(X)):
        buckets[tuple(cell[idx])].append(idx)
    for k in buckets:
        buckets[k] = np.array(buckets[k])
    offs = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    PI, PJ = [], []
    for k, mem in buckets.items():
        blk = np.concatenate([buckets[(k[0]+o[0], k[1]+o[1], k[2]+o[2])]
                              for o in offs if (k[0]+o[0], k[1]+o[1], k[2]+o[2]) in buckets])
        for i in mem:
            d2 = ((X[blk] - X[i]) ** 2).sum(1)
            j = blk[(d2 < cutoff ** 2) & (blk > i)]
            if len(j):
                PI.append(np.full(len(j), i)); PJ.append(j)
    return np.concatenate(PI), np.concatenate(PJ)


def _scatter(pi, pj, vals, N):                    # +vals on i, -vals on j (antisym force)
    return np.bincount(pi, vals, N) - np.bincount(pj, vals, N)


def lj_force_sparse(X, pi, pj, N):
    d = X[pi] - X[pj]; r2 = (d ** 2).sum(1); w = r2 < RCUT ** 2
    inv2 = np.where(w, SIGMA ** 2 / r2, 0.0); inv6 = inv2 ** 3; inv12 = inv6 ** 2
    coef = np.where(w, 24 * EPS * (2 * inv12 - inv6) / np.where(w, r2, 1.0), 0.0)
    fp = coef[:, None] * d
    return np.stack([_scatter(pi, pj, fp[:, c], N) for c in range(3)], 1)


def coupling_sparse(X, u, pi, pj, N, rho0):
    d = X[pi] - X[pj]; r2 = (d ** 2).sum(1)
    W = np.where(r2 < RCUT_FIELD ** 2, np.exp(-r2 / (2 * H ** 2)), 0.0)
    rho = np.bincount(pi, W, N) + np.bincount(pj, W, N)   # sum_{j!=i} W_ij (matches dense _geom)
    Q = (u[pi] - u[pj]) ** 2
    A = np.bincount(pi, W * Q, N) + np.bincount(pj, W * Q, N)
    s = 1.0 / (1.0 + np.exp(BETA * (rho / rho0 - 1.0)))
    g = GMIN + (GMAX - GMIN) * s
    gp = (GMAX - GMIN) * (-s * (1 - s)) * (BETA / rho0)
    bracket = (g[pi] + g[pj]) * Q + gp[pi] * A[pi] + gp[pj] * A[pj]
    coef = (1.0 / (2 * H ** 2)) * bracket * W
    fp = coef[:, None] * d
    F = np.stack([_scatter(pi, pj, fp[:, c], N) for c in range(3)], 1)
    return F, rho


def density_sparse(X, pi, pj, N):
    r2 = ((X[pi] - X[pj]) ** 2).sum(1)
    W = np.where(r2 < RCUT_FIELD ** 2, np.exp(-r2 / (2 * H ** 2)), 0.0)
    return np.bincount(pi, W, N) + np.bincount(pj, W, N)   # sum_{j!=i} W_ij (matches dense _geom)


def relax(X0, u, pi, pj, rho0, steps=2500, dt=0.002, cool=0.97):
    X = X0.copy(); Vn = np.zeros_like(X); N = len(X)
    F = lj_force_sparse(X, pi, pj, N) + coupling_sparse(X, u, pi, pj, N, rho0)[0]
    for _ in range(steps):
        X += Vn * dt + 0.5 * F * dt ** 2
        Fn = lj_force_sparse(X, pi, pj, N) + coupling_sparse(X, u, pi, pj, N, rho0)[0]
        Vn += 0.5 * (F + Fn) * dt; Vn *= cool; F = Fn
    return X


def verify_small():
    """Sparse coupling force must match the dense variational engine at small N."""
    fcc = perfect_fcc(radius=5.0)
    N = len(fcc); pi, pj = cell_pairs(fcc, 2.9)
    u = 0.8 * np.exp(-(fcc ** 2).sum(1) / (2 * 3.0 ** 2))
    rho0 = float(np.median(density_sparse(fcc, pi, pj, N)))
    Fs, _ = coupling_sparse(fcc, u, pi, pj, N, rho0)
    f = V.VariationalCoupled(fcc, beta=BETA, m2=1.0, g_min=GMIN, damping=1.0)
    f.u = u.copy(); f.pi = np.zeros(N)
    _, Fd, _ = f.forces()
    Fd = Fd - __import__("bvc_core").lj_forces_energy(fcc)[0]   # dense coupling part only
    err = np.linalg.norm(Fs - Fd) / (np.linalg.norm(Fd) + 1e-9)
    print(f"sparse-vs-dense coupling force rel-err @N={N}: {err:.2e}  "
          f"({'OK' if err < 1e-3 else 'MISMATCH'})")


def main():
    verify_small()
    fcc = perfect_fcc(radius=13.0)
    N = len(fcc)
    print(f"\nlarge fcc: N={N}  radius~{np.linalg.norm(fcc, axis=1).max():.1f}")
    pi, pj = cell_pairs(fcc, 2.9)
    print(f"pairs: {len(pi)}  (~{2*len(pi)/N:.0f} neighbours/node)")
    rho0 = float(np.median(density_sparse(fcc, pi, pj, N)))

    Xref = relax(fcc, np.zeros(N), pi, pj, rho0)
    src = 0.8 * np.exp(-(fcc ** 2).sum(1) / (2 * 3.5 ** 2))
    Xsrc = relax(fcc, src, pi, pj, rho0)

    r = np.linalg.norm(fcc, axis=1)
    drho_node = density_sparse(Xsrc, pi, pj, N) - density_sparse(Xref, pi, pj, N)
    bins = np.arange(0, 12.1, 1.0); rmid = 0.5 * (bins[:-1] + bins[1:])
    drho = np.array([drho_node[(r >= lo) & (r < hi)].mean() if ((r >= lo) & (r < hi)).sum() else np.nan
                     for lo, hi in zip(bins[:-1], bins[1:])])
    print("\n3D-2 large-fcc Drho(r):")
    print("  r   :", " ".join(f"{x:5.1f}" for x in rmid))
    print("  Drho:", " ".join(f"{x:+6.4f}" for x in drho))
    m = np.isfinite(drho) & (np.abs(drho) > 1e-5) & (rmid > 2.5) & (rmid < 10.5)
    if m.sum() >= 3:
        rf, y = rmid[m], np.abs(drho[m])
        lam = -1 / np.polyfit(rf, np.log(y), 1)[0]
        pn = np.polyfit(np.log(rf), np.log(y), 1); n = -pn[0]
        se = np.sum((np.log(y) - np.polyval(np.polyfit(rf, np.log(y), 1), rf)) ** 2)
        sp = np.sum((np.log(y) - np.polyval(pn, np.log(rf))) ** 2)
        print(f"\n  fit: exp lambda={lam:.2f} (SS={se:.3f})   power n={n:.2f} (SS={sp:.3f})")
        print("  => better:", "EXPONENTIAL/screened" if se < sp else f"POWER-LAW  ~1/r^{n:.1f}")
        print("  (Newtonian 3D force would be ~1/r^2; the *potential/well* ~1/r)")


if __name__ == "__main__":
    main()
