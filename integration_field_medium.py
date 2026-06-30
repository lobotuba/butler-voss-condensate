"""
Integration Phases 0-1 : a field on a self-assembled medium
===========================================================

PROTOTYPE.  Brings the H10 mobile-node medium together with the H6-H9 complex
field, in two de-risking steps:

  PHASE 0  Infrastructure + the central go/no-go test.  Build the meshfree
           operators on an irregular point cloud and MEASURE how badly the
           Laplacian degrades off a perfect lattice -- weighted-graph vs a
           least-squares (LSQ) reference, on a perfect hex lattice vs a relaxed
           (self-assembled) cloud.  If the weighted-graph Laplacian holds up,
           the whole integration is viable on the cheap operator.

  PHASE 1  Freeze a self-assembled medium and run the complex field on it.
           Seed a vortex and check it PERSISTS with a quantized, conserved
           winding -- i.e. a topological particle can live on a medium the nodes
           built themselves (vs the same test on a perfect hex lattice).

Pure numpy (no scipy).  Reuses the Lennard-Jones medium from
prototype_mobile_nodes and the field ideas from prototype_complex.
"""
from __future__ import annotations
import math
import numpy as np
from prototype_mobile_nodes import lj_forces_energy, R0, RCUT


# ============================================================ medium / cloud ==
def relax_medium(N=300, steps=6000, dt=0.005, cool=0.999, seed=1):
    """Self-assemble a cooled LJ droplet; return its node positions (quiet)."""
    rng = np.random.default_rng(seed)
    area_per, = (1.7,)
    R = math.sqrt(N * area_per / math.pi)
    X = []
    while len(X) < N:
        p = (rng.random(2) * 2 - 1) * R
        if math.hypot(*p) <= R and all(math.hypot(*(p - q)) > 0.95 for q in X):
            X.append(p)
    X = np.array(X)
    V = np.zeros_like(X)
    F, _ = lj_forces_energy(X)
    for _ in range(steps):
        X += V * dt + 0.5 * F * dt ** 2
        Fn, _ = lj_forces_energy(X)
        V += 0.5 * (F + Fn) * dt
        V *= cool
        F = Fn
    return X - X.mean(0)


def perfect_hex(radius_cells=9, a=R0):
    """Triangular (hex) lattice clipped to a disk; the isotropic control."""
    pts, Rmax = [], radius_cells * a
    for j in range(-2 * radius_cells, 2 * radius_cells + 1):
        for i in range(-2 * radius_cells, 2 * radius_cells + 1):
            x = a * (i + 0.5 * (j & 1))
            y = a * math.sqrt(3) / 2 * j
            if x * x + y * y <= Rmax * Rmax:
                pts.append((x, y))
    X = np.array(pts)
    return X - X.mean(0)


# ===================================================== meshfree operators =====
def _pairwise(X):
    d = X[:, None, :] - X[None, :, :]            # d[i,j] = X_i - X_j
    return d, (d ** 2).sum(-1)


def laplacian_matrix(X, rcut):
    """Consistent weighted-graph Laplacian as a matrix A (so L f = A @ f).
    Normalised by 2D / sum_j w_ij |r_ij|^2 so it reproduces div-grad for smooth
    f on an isotropic neighbourhood (D = 2)."""
    d, r2 = _pairwise(X)
    N = len(X)
    W = ((r2 > 1e-12) & (r2 < rcut ** 2)).astype(float)
    swr2 = (W * r2).sum(1)
    norm = (2 * 2) / np.where(swr2 > 0, swr2, 1.0)
    A = norm[:, None] * W
    A[np.arange(N), np.arange(N)] = -A.sum(1)
    return A


def lsq_laplacian_values(X, f, rcut):
    """Least-squares quadratic-fit Laplacian (the accurate reference)."""
    d, r2 = _pairwise(X)
    out = np.full(len(X), np.nan)
    for i in range(len(X)):
        sel = np.where((r2[i] > 1e-12) & (r2[i] < rcut ** 2))[0]
        if len(sel) < 6:
            continue
        dl = X[sel] - X[i]                        # delta = X_j - X_i
        B = np.stack([dl[:, 0], dl[:, 1],
                      0.5 * dl[:, 0] ** 2, dl[:, 0] * dl[:, 1], 0.5 * dl[:, 1] ** 2], 1)
        c, *_ = np.linalg.lstsq(B, f[sel] - f[i], rcond=None)
        out[i] = c[2] + c[4]                      # Hxx + Hyy = trace(Hessian)
    return out


def lsq_laplacian_matrix(X, rcut):
    """The LSQ Laplacian as a matrix A (L f = A @ f).  Each row is the linear
    functional that extracts trace(Hessian) from the local quadratic fit, so it
    stays accurate on irregular meshes (it removes the spurious gradient term
    that wrecks the plain graph Laplacian)."""
    d, r2 = _pairwise(X)
    N = len(X)
    A = np.zeros((N, N))
    for i in range(N):
        sel = np.where((r2[i] > 1e-12) & (r2[i] < rcut ** 2))[0]
        if len(sel) < 6:
            continue
        dl = X[sel] - X[i]
        B = np.stack([dl[:, 0], dl[:, 1],
                      0.5 * dl[:, 0] ** 2, dl[:, 0] * dl[:, 1], 0.5 * dl[:, 1] ** 2], 1)
        M = np.linalg.pinv(B)                     # (5, m) = (B^T B)^-1 B^T
        row = M[2] + M[4]                         # trace functional over neighbours
        A[i, sel] += row
        A[i, i] -= row.sum()
    return A


def interior_mask(X, frac=0.62):
    r = np.linalg.norm(X - X.mean(0), axis=1)
    return r < frac * r.max()


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
        nn = np.sqrt(np.sort(_pairwise(X)[1], axis=1)[:, 1]).mean()
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
