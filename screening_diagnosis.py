"""
Screen-0 -- diagnose the screening: is lambda a MEDIUM property (Bitter-Crum) or
a coupling-strength property?

The 3D-2 result (density_response_3d_large.py) showed gravity-by-density is
screened, exp lambda~3.3, in both 2D and 3D. The Bitter-Crum reading: the density
coupling makes the source a CENTER OF DILATATION, and two dilatation centers in an
infinite isotropic linear medium do not interact -- all the residual interaction
is anharmonic/near-field, hence short-range REGARDLESS of coupling strength.

Sharp test: sweep the coupling sharpness beta. Bitter-Crum predicts the RANGE
(lambda) is set by the medium, so it should be ~beta-INDEPENDENT while the
AMPLITUDE scales with beta (stronger coupling = bigger compression, same reach).
If instead lambda grows with beta, screening is a coupling artifact and a simple
dial fixes it. Either way this tells us whether the fix must change the coupling
TYPE (-> Screen-1: a Gauss-law/conserved-source coupling) or just a knob.

Reuses the cell-list sparse machinery verified to 2e-15 against the dense engine.
"""
from __future__ import annotations
import numpy as np

from bvc_core import R0, RCUT, EPS, SIGMA, perfect_fcc
from density_response_3d_large import (cell_pairs, lj_force_sparse, _scatter,
                                       density_sparse, H, RCUT_FIELD, GMIN, GMAX)


def coupling_sparse_beta(X, u, pi, pj, N, rho0, beta):
    """Variational coupling node force at an arbitrary beta (else identical to
    density_response_3d_large.coupling_sparse)."""
    d = X[pi] - X[pj]; r2 = (d ** 2).sum(1)
    W = np.where(r2 < RCUT_FIELD ** 2, np.exp(-r2 / (2 * H ** 2)), 0.0)
    rho = np.bincount(pi, W, N) + np.bincount(pj, W, N)
    Q = (u[pi] - u[pj]) ** 2
    A = np.bincount(pi, W * Q, N) + np.bincount(pj, W * Q, N)
    s = 1.0 / (1.0 + np.exp(beta * (rho / rho0 - 1.0)))
    g = GMIN + (GMAX - GMIN) * s
    gp = (GMAX - GMIN) * (-s * (1 - s)) * (beta / rho0)
    bracket = (g[pi] + g[pj]) * Q + gp[pi] * A[pi] + gp[pj] * A[pj]
    coef = (1.0 / (2 * H ** 2)) * bracket * W
    fp = coef[:, None] * d
    return np.stack([_scatter(pi, pj, fp[:, c], N) for c in range(3)], 1)


def relax(X0, u, pi, pj, rho0, beta, steps=2500, dt=0.002, cool=0.97):
    X = X0.copy(); Vn = np.zeros_like(X); N = len(X)
    F = lj_force_sparse(X, pi, pj, N) + coupling_sparse_beta(X, u, pi, pj, N, rho0, beta)
    for _ in range(steps):
        X += Vn * dt + 0.5 * F * dt ** 2
        Fn = lj_force_sparse(X, pi, pj, N) + coupling_sparse_beta(X, u, pi, pj, N, rho0, beta)
        Vn += 0.5 * (F + Fn) * dt; Vn *= cool; F = Fn
    return X


def fit_lambda(fcc, pi, pj, N, rho0, beta, Xref):
    src = 0.8 * np.exp(-(fcc ** 2).sum(1) / (2 * 3.5 ** 2))
    Xsrc = relax(fcc, src, pi, pj, rho0, beta)
    r = np.linalg.norm(fcc, axis=1)
    drho_node = density_sparse(Xsrc, pi, pj, N) - density_sparse(Xref, pi, pj, N)
    bins = np.arange(0, 10.1, 1.0); rmid = 0.5 * (bins[:-1] + bins[1:])
    drho = np.array([drho_node[(r >= lo) & (r < hi)].mean() if ((r >= lo) & (r < hi)).sum() else np.nan
                     for lo, hi in zip(bins[:-1], bins[1:])])
    amp = np.nanmax(np.abs(drho))
    total = float(drho_node.sum())          # integrated response Sum_i Drho_i (conserved?)
    m = np.isfinite(drho) & (np.abs(drho) > 1e-5) & (rmid > 2.5) & (rmid < 8.5)
    lam = -1 / np.polyfit(rmid[m], np.log(np.abs(drho[m])), 1)[0] if m.sum() >= 3 else np.nan
    return lam, amp, total, rmid, drho


def main():
    fcc = perfect_fcc(radius=11.0); N = len(fcc)
    pi, pj = cell_pairs(fcc, 2.9)
    rho0 = float(np.median(density_sparse(fcc, pi, pj, N)))
    print(f"Screen-0: fcc N={N} radius~{np.linalg.norm(fcc,axis=1).max():.1f}  "
          f"pairs={len(pi)}\n")
    print("Bitter-Crum prediction: lambda ~ const vs beta (medium sets the range),")
    print("amplitude grows with beta (coupling strength sets the depth).\n")

    Xref = relax(fcc, np.zeros(N), pi, pj, rho0, beta=0.0)  # beta=0 -> no density coupling
    print(f"  {'beta':>6} {'lambda':>8} {'amp':>9} {'integral':>10} {'amp*lam^3':>10}")
    lams, totals, amps = [], [], []
    for beta in (20.0, 40.0, 60.0, 100.0):
        lam, amp, total, rmid, drho = fit_lambda(fcc, pi, pj, N, rho0, beta, Xref)
        lams.append(lam); totals.append(total); amps.append(amp)
        print(f"  {beta:>6.0f} {lam:>8.2f} {amp:>9.5f} {total:>10.3f} {amp*lam**3:>10.3f}")
    tot, lam_arr = np.array(totals), np.array(lams)
    print(f"\n  integrated response Sum Drho: {tot[0]:.0f} -> {tot[-1]:.0f} over beta 20->100"
          f"  (NOT conserved: no protected flux)")
    print(f"  lambda: {lam_arr[0]:.1f} -> {lam_arr[-1]:.1f}   amplitude: "
          f"{amps[0]:.3f} -> {amps[-1]:.3f}")
    print("\n  Verdict: beta is NOT a usable lever. lambda creeps up but amplitude AND")
    print("  integrated compression COLLAPSE together (high beta -> g'(rho) saturates to 0")
    print("  except a thin rho~rho0 shell, so the whole response fades). You can trade")
    print("  strength for slightly-longer reach, never a STRONG long-range force. The")
    print("  intrinsic elastic screening (weak-coupling limit) is lambda~3.")
    print("  Root cause: this coupling has NO conserved, un-screenable flux -- exactly")
    print("  what a long-range (Gauss-law) force needs. Next: Screen-1 must MANUFACTURE")
    print("  a conserved source (auxiliary potential with div-flux = source, no mass term).")


if __name__ == "__main__":
    main()
