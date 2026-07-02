"""
3D density response (Phases 3D-0 / 3D-1) -- does 3D lengthen the gravity force?

2D screens the medium's strain hard (perfect hex: screening length lambda~5.7).
In 3D, elasticity is longer-range, so an ordered fcc medium might give a
much longer-range (or power-law) force -- Newtonian territory.

The variational engine (integration_phase3_variational.VariationalCoupled) is
dimension-agnostic (pairwise ops + an inline symmetric operator), so it runs in
3D unchanged. Here:
  3D-0  seed a lump on an fcc medium; check the engine runs + conserves energy.
  3D-1  static density response Drho(r) to a gentle frozen mass on ordered fcc,
        3D radial shells, fit exp vs power-law; compare to 2D hex (lambda~5.7).

Result (fcc radius 9, N~3055): Drho(r) barely decays over the measurable range
(factor ~1.6 to r=6.8, vs 2D's factor ~5) -- 3D is MUCH less screened
(lambda~14, ~2.4x 2D; power-law not excluded). Caveat: radius 9 is < half a
screening length and near the surface, so long-lambda vs power-law vs finite-size
can't be distinguished here -> motivates 3D-2 (sparse scale-up to radius ~15-20).
"""
from __future__ import annotations
import numpy as np

from bvc_core import R0, pairwise, perfect_fcc
import integration_phase3_variational as V


def sph_density(X, h=R0):
    _, r2 = pairwise(X)
    return np.exp(-r2 / (2 * h ** 2)).sum(1)


def radial(X, rho, bins):
    r = np.linalg.norm(X - X.mean(0), axis=1)
    return np.array([rho[(r >= lo) & (r < hi)].mean() if ((r >= lo) & (r < hi)).sum() else np.nan
                     for lo, hi in zip(bins[:-1], bins[1:])])


def relax_with(fcc, u_field, beta=60.0, steps=3000, dt=0.002):
    g = V.VariationalCoupled(fcc, beta=beta, m2=1.0, g_min=0.02, damping=1.0, dt=dt)
    g.u = u_field.copy(); g.pi = np.zeros(g.N)
    for _ in range(steps):
        _, FX, _ = g.forces(); g.Wn += 0.5*FX*dt; g.X += g.Wn*dt
        _, FX, _ = g.forces(); g.Wn += 0.5*FX*dt; g.Wn *= 0.97
    return g.X


def main():
    fcc = perfect_fcc(radius=9.0)
    print(f"3D fcc: N={len(fcc)}  radius~{np.linalg.norm(fcc, axis=1).max():.1f}  dim={fcc.shape[1]}")

    # 3D-0: engine runs + conserves in 3D?
    f = V.VariationalCoupled(fcc, beta=40, m2=1.0, g_min=0.02, damping=1.0, dt=0.002)
    f.seed_lump(amp=1.0, width=3.0); E0 = f.energy()
    for _ in range(1500): f.step()
    print(f"3D-0: t={f.time:.1f}  energy drift={100*(f.energy()-E0)/abs(E0):+.3f}%  "
          f"max|u|={abs(f.u).max():.2f}  => {'OK' if abs(f.u).max()<20 else 'BLEW UP'}")

    # 3D-1: static density response on ordered fcc
    Xref = relax_with(fcc, np.zeros(len(fcc)))
    src = 0.8 * np.exp(-(fcc ** 2).sum(1) / (2 * 3.5 ** 2))
    Xsrc = relax_with(fcc, src)
    bins = np.arange(0, 7.6, 0.8); rmid = 0.5 * (bins[:-1] + bins[1:])
    drho = radial(Xsrc, sph_density(Xsrc), bins) - radial(Xref, sph_density(Xref), bins)
    print("\n3D-1 ordered-fcc Drho(r):")
    print("  r   :", " ".join(f"{r:5.1f}" for r in rmid))
    print("  Drho:", " ".join(f"{d:+6.4f}" for d in drho))
    m = np.isfinite(drho) & (np.abs(drho) > 1e-5) & (rmid > 1.5) & (rmid < 6.5)
    if m.sum() >= 3:
        rf, y = rmid[m], np.abs(drho[m])
        lam = -1 / np.polyfit(rf, np.log(y), 1)[0]
        n = -np.polyfit(np.log(rf), np.log(y), 1)[0]
        print(f"  fit: exp lambda={lam:.2f}   power n={n:.2f}   (2D hex: lambda~5.7)")
        print("  => 3D is much less screened than 2D (radius 9 too small to fix the form)")


if __name__ == "__main__":
    main()
