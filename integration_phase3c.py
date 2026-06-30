"""
Integration Phase 3c : close the loop on ONE lump (self-focusing)
=================================================================

PROTOTYPE.  3a/3b validated the two halves of gravity-by-density in isolation and
fixed the sign (denser => slower; field energy compresses the medium).  3c turns
BOTH couplings on at once, with a REAL evolving field, on a single lump:

    field energy  --(A)-->  compress nodes  -->  denser  --(B)-->  slower waves
         ^                                                              |
         +-------------------  lump disperses less, stays put  <--------+

A lump that would otherwise disperse should SELF-FOCUS.  The danger is runaway
(energy -> compression -> slower -> more focusing -> ...), the same feedback that
blew up cubic3d.  The natural brake is the Lennard-Jones short-range repulsion:
nodes cannot compress past ~r0, so the density well -- and the focusing -- should
SATURATE at a finite width rather than collapse.  First milestone: a STABLE run.

Couplings (both off in the control):
  A  node force  F += alpha * grad(e_field)        [3b sign: pulls nodes into energy]
  B  field speed c^2 * g(rho),  g = 1/(1+beta*(rho/rho0 - 1))   [3a sign: denser=slower]

Field: a massive real scalar (no self-trapping of its own), so any focusing comes
from the coupling, not the potential.
"""
from __future__ import annotations
import numpy as np

from bvc_core import lj_forces_energy, relax_medium, pairwise, R0, brookshaw_laplacian


# ----------------------------------------------- meshfree Laplacian + gradient -
def lsq_operators(X, rcut):
    """LSQ Laplacian A and gradient Gx, Gy as matrices (all from one local fit)."""
    _, r2 = pairwise(X)
    N = len(X)
    A = np.zeros((N, N)); Gx = np.zeros((N, N)); Gy = np.zeros((N, N))
    for i in range(N):
        sel = np.where((r2[i] > 1e-12) & (r2[i] < rcut ** 2))[0]
        if len(sel) < 6:
            continue
        dl = X[sel] - X[i]
        B = np.stack([dl[:, 0], dl[:, 1],
                      0.5 * dl[:, 0] ** 2, dl[:, 0] * dl[:, 1], 0.5 * dl[:, 1] ** 2], 1)
        M = np.linalg.pinv(B)
        for G, row in ((Gx, M[0]), (Gy, M[1]), (A, M[2] + M[4])):
            G[i, sel] += row
            G[i, i] -= row.sum()
    return A, Gx, Gy


def sph_density(X, h):
    _, r2 = pairwise(X)
    return np.exp(-r2 / (2 * h ** 2)).sum(1)


class CoupledLump:
    def __init__(self, X, alpha=0.0, beta=0.0, c=1.0, m2=1.0, h=None,
                 dt=0.005, field_damping=0.999, cool=0.999, g_min=0.08,
                 rcut=1.9 * R0, rebuild_every=25):
        self.X = X.copy(); self.V = np.zeros_like(X); self.N = len(X)
        self.alpha, self.beta, self.c, self.m2 = alpha, beta, c, m2
        self.h = h or R0
        self.dt, self.field_damping, self.cool = dt, field_damping, cool
        self.g_min = g_min
        self.rcut, self.rebuild_every = rcut, rebuild_every
        # Brookshaw (symmetric, stable) for the FIELD Laplacian; LSQ gradient
        # (Gx,Gy) for the node force (applied to smoothed energy, so it's fine).
        _, self.Gx, self.Gy = lsq_operators(self.X, rcut)
        self.A = brookshaw_laplacian(self.X, rcut=rcut)
        self.rho0 = float(np.median(sph_density(self.X, self.h)))
        self.u = np.zeros(self.N); self.ud = np.zeros(self.N)
        self.Fmed, _ = lj_forces_energy(self.X)
        self.time = 0.0

    def seed_lump(self, amp=1.2, width=4.0):
        r2 = (self.X ** 2).sum(1)
        self.u = amp * np.exp(-r2 / (2 * width ** 2))
        self.ud = np.zeros(self.N)

    def energy_density(self):
        return 0.5 * self.ud ** 2 + 0.5 * self.m2 * self.u ** 2

    def g_factor(self):
        rho = sph_density(self.X, self.h)
        g = 1.0 / (1.0 + self.beta * (rho / self.rho0 - 1.0))
        return np.clip(g, self.g_min, 2.0)

    def smooth_energy(self, e):
        """Kernel-average the field energy so its gradient (the node force) is
        smooth -- the raw mesh field is too rough and ejects nodes."""
        _, r2 = pairwise(self.X)
        W = np.exp(-r2 / (2 * self.h ** 2))
        return (W @ e) / W.sum(1)

    def step(self, k):
        e = self.energy_density()
        es = self.smooth_energy(e)
        # Coupling A: nodes pulled up the (smoothed) field-energy gradient
        Ffield = self.alpha * np.stack([self.Gx @ es, self.Gy @ es], 1)
        # medium velocity-Verlet (LJ + field force), gentle cooling
        Ftot = self.Fmed + Ffield
        self.X += self.V * self.dt + 0.5 * Ftot * self.dt ** 2
        Fmed_new, _ = lj_forces_energy(self.X)
        self.Fmed = Fmed_new
        self.V += 0.5 * (Ftot + Fmed_new + Ffield) * self.dt
        self.V *= self.cool
        if k % self.rebuild_every == 0:
            _, self.Gx, self.Gy = lsq_operators(self.X, self.rcut)
            self.A = brookshaw_laplacian(self.X, rcut=self.rcut)
        # Coupling B: field on the current geometry, wave speed slowed where dense
        g = self.g_factor()
        udd = self.c ** 2 * g * (self.A @ self.u) - self.m2 * self.u
        self.ud += udd * self.dt
        self.ud *= self.field_damping
        self.u += self.ud * self.dt
        self.time += self.dt

    # -- diagnostics --
    def lump_width(self):
        e = self.energy_density(); s = e.sum()
        if s < 1e-9:
            return float("nan")
        c = (self.X * e[:, None]).sum(0) / s
        return float(np.sqrt(((self.X - c) ** 2).sum(1) @ e / s))

    def central_density_ratio(self):
        e = self.energy_density(); s = e.sum()
        c = (self.X * e[:, None]).sum(0) / s if s > 1e-9 else np.zeros(2)
        near = np.linalg.norm(self.X - c, axis=1) < 1.5 * self.h
        return float(sph_density(self.X, self.h)[near].mean() / self.rho0)

    def field_energy(self):
        return float(self.energy_density().sum())


def run(alpha, beta, label, steps=4000, m2=1.0, cool=0.999, fd=0.999):
    cloud = relax_medium(N=260, seed=3)
    f = CoupledLump(cloud, alpha=alpha, beta=beta, m2=m2, cool=cool, field_damping=fd)
    f.seed_lump()
    w0 = f.lump_width()
    traj = []
    for k in range(steps):
        f.step(k)
        if k % (steps // 6) == 0 or k == steps - 1:
            traj.append((round(f.time, 1), round(f.lump_width(), 2),
                         round(f.central_density_ratio(), 2),
                         round(float(np.abs(f.u).max()), 2)))
    print(f"  {label}:")
    print("    " + "  ".join(f"t{t}:width={w},dens={d},|u|max={m}" for t, w, d, m in traj))
    wF = traj[-1][1]
    blew = traj[-1][3] > 50 or not np.isfinite(wF)
    trend = "DISPERSES" if wF > w0 + 0.5 else ("FOCUSES" if wF < w0 - 0.5 else "holds")
    print(f"    start width={w0:.2f} -> end {wF}  =>  "
          f"{'BLEW UP' if blew else trend}\n")


if __name__ == "__main__":
    print("=== Phase 3c : close the loop on one lump (self-focusing?) ===")
    print("  width = energy-weighted RMS radius; dens = central node density / rho0\n")
    # Brookshaw field operator is stable, so beta (well depth) can be pushed hard
    # AND the field run near-conservatively (fd~1) -- both impossible with LSQ.
    # Near-conservative is the decisive self-trapping test: width AND |u| holding.
    run(0.0, 0.0, "control (fd~1, no coupling)", fd=0.9998)
    run(10.0, 60.0, "Brookshaw coupled beta=60 (fd~1)", cool=0.99, fd=0.9998)
