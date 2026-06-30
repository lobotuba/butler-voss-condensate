"""
Integration Phases 0-1 : a field on a self-assembled medium
===========================================================

PROTOTYPE.  Brings the H10 mobile-node medium together with the H6-H9 complex
field, in two de-risking steps:

  PHASE 0  Infrastructure + the central go/no-go test.  MEASURE how badly the
           meshfree Laplacian degrades off a perfect lattice -- weighted-graph vs
           a least-squares (LSQ) reference, on a perfect hex lattice vs a relaxed
           (self-assembled) cloud.  If the cheap graph Laplacian holds up, the
           whole integration is viable on it; if not, the field needs LSQ.

  PHASE 1  Freeze a self-assembled medium and run the complex field on it.
           Seed a vortex and check it PERSISTS with a quantized, conserved
           winding -- i.e. a topological particle can live on a medium the nodes
           built themselves (vs the same test on a perfect hex lattice).

Pure numpy (no scipy).  Shared medium + operators come from bvc_core.
"""
from __future__ import annotations
import math
import numpy as np

from bvc_core import (R0, relax_medium, perfect_hex, pairwise,
                      laplacian_matrix, lsq_laplacian_values, lsq_laplacian_matrix,
                      interior_mask)


# ============================================================ PHASE 0 =========
def phase0():
    print("PHASE 0 -- meshfree Laplacian accuracy (go/no-go)")
    print("  test f = cos(k x), wavelength ~ 8 spacings; RMS rel. error on interior\n")
    hexX = perfect_hex()
    cloud = relax_medium()
    k = 2 * math.pi / (8 * R0)

    def errors(X, label):
        f = np.cos(k * X[:, 0])
        ana = -k ** 2 * f
        Ag = laplacian_matrix(X, rcut=1.45 * R0)
        Lg = Ag @ f
        Ll = lsq_laplacian_values(X, f, rcut=1.9 * R0)
        m = interior_mask(X) & np.isfinite(Ll)
        rms = lambda L: math.sqrt(np.mean((L[m] - ana[m]) ** 2) / np.mean(ana[m] ** 2))
        nn = np.sqrt(np.sort(pairwise(X)[1], axis=1)[:, 1]).mean()
        print(f"  {label:16} N={len(X):4d}  mean-spacing={nn:5.3f}  "
              f"graph-RMS={rms(Lg)*100:6.2f}%   LSQ-RMS={rms(Ll)*100:6.2f}%")
        return rms(Lg)

    e_hex = errors(hexX, "perfect hex")
    e_cloud = errors(cloud, "relaxed cloud")
    print(f"\n  -> weighted-graph Laplacian degradation, cloud vs hex: "
          f"{e_cloud*100:.2f}% vs {e_hex*100:.2f}%  "
          f"(factor {e_cloud/max(e_hex,1e-9):.1f}x)")
    verdict = ("GO: graph Laplacian usable for Phase 1" if e_cloud < 0.15
               else "MARGINAL: prefer LSQ operator for the field")
    print(f"  -> {verdict}\n")
    return cloud, hexX


# ============================================================ PHASE 1 =========
class FrozenField:
    """Complex field on a fixed (irregular) point set; Laplacian precomputed."""
    def __init__(self, X, laplacian="graph", rcut=None, v0=1.0, lam=1.0, c=1.0,
                 dt=0.02, damping=1.0):
        self.X, self.center = X, X.mean(0)
        if laplacian == "lsq":
            self.A = lsq_laplacian_matrix(X, rcut or 1.9 * R0)
        else:
            self.A = laplacian_matrix(X, rcut or 1.45 * R0)
        self.v0, self.lam, self.c, self.dt, self.damping = v0, lam, c, dt, damping
        self.psi = np.full(len(X), v0, dtype=np.complex128)
        self.psid = np.zeros(len(X), dtype=np.complex128)
        self.time = 0.0

    def seed_vortex(self, n=1, core=2.5):
        d = self.X - self.center
        r = np.linalg.norm(d, axis=1)
        phi = np.arctan2(d[:, 1], d[:, 0])
        self.psi = self.v0 * np.tanh(r / core) * np.exp(1j * n * phi)
        self.psid = np.zeros_like(self.psi)

    def step(self):
        rho = np.abs(self.psi) ** 2
        force = -self.lam * (rho - self.v0 ** 2) * self.psi
        accel = self.c ** 2 * (self.A @ self.psi) + force
        self.psid += accel * self.dt
        self.psid *= self.damping
        self.psi += self.psid * self.dt
        self.time += self.dt

    def _winding_at(self, frac, band=1.2):
        d = self.X - self.center
        r = np.linalg.norm(d, axis=1)
        R = frac * r.max()
        sel = np.where(np.abs(r - R) < band * R0)[0]
        if len(sel) < 8:
            return None
        ang = np.arctan2(d[sel, 1], d[sel, 0])
        loop = sel[np.argsort(ang)]
        p = self.psi[loop]
        dphi = np.angle(np.roll(p, -1) * np.conj(p))
        return int(round(dphi.sum() / (2 * math.pi)))

    def enclosed_winding(self, fracs=(0.40, 0.50, 0.60)):
        """Charge enclosed, voted over several loop radii (a single loop on an
        irregular mesh is noisy; the median is robust)."""
        vals = [w for f in fracs if (w := self._winding_at(f)) is not None]
        return int(np.median(vals)) if vals else None

    def core_depth(self, frac=0.12):
        d = np.linalg.norm(self.X - self.center, axis=1)
        inner = d < frac * d.max()
        return float(np.abs(self.psi[inner]).min() / self.v0) if inner.any() else 1.0


def phase1(cloud, hexX, steps=2500):
    print("PHASE 1 -- does a vortex persist on a self-assembled medium?")
    print("  seed n=+1 vortex; track enclosed winding & core depth over time\n")
    # (operator, dt, damping, steps) -- LSQ needs the smaller dt; mild damping
    # sheds the radiation from seeding a non-equilibrium vortex.  All reach t~50.
    cases = [
        ("perfect hex   (graph op, control)",   hexX,  "graph", 0.02, 1.000, 2500),
        ("relaxed cloud (graph op -- inadequate)", cloud, "graph", 0.02, 1.000, 2500),
        ("relaxed cloud (LSQ op + mild damping)",  cloud, "lsq",   0.01, 0.998, 5000),
    ]
    for label, X, op, dt, damp, ns in cases:
        f = FrozenField(X, laplacian=op, dt=dt, damping=damp)
        f.seed_vortex(n=1)
        w0 = f.enclosed_winding()
        traj = []
        for kk in range(ns):
            f.step()
            if kk % (ns // 6) == 0 or kk == ns - 1:
                traj.append((round(f.time, 1), f.enclosed_winding(), round(f.core_depth(), 2)))
        wF = f.enclosed_winding()
        maj = sum(1 for _, w, _ in traj if w == 1) / len(traj)
        end_core = traj[-1][2]
        # topological charge can only change if the zero exits or pair-annihilates;
        # matching start==end==1 (+ majority) means it was conserved, transient
        # mis-reads are the wandering core crossing a measurement loop.
        charge_ok = (w0 == 1 and wF == 1 and maj >= 0.7)
        core = "sharp" if end_core < 0.3 else "under-resolved (shallow core)"
        print(f"  {label}:  start winding={w0}")
        print("    " + "  ".join(f"t{t}:w={w},core={c}" for t, w, c in traj))
        print(f"    end winding={wF}  =>  charge "
              f"{'CONSERVED' if charge_ok else 'LOST'};  core {core}\n")


if __name__ == "__main__":
    print("=== Integration Phases 0-1 : field on a self-assembled medium ===\n")
    cloud, hexX = phase0()
    phase1(cloud, hexX)
