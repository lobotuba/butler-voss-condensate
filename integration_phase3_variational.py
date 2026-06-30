"""
Integration Phase 3 : VARIATIONAL coupling (one energy, conserved by construction)
==================================================================================

PROTOTYPE.  In 3c the two gravity-by-density couplings were assembled by hand, so
the system was not Hamiltonian and energy LEAKED (the self-trapped lump slowly
faded).  Here both couplings come from a SINGLE energy functional, so the field
force and the node force are exact gradients of the same E -> energy is conserved
and the field operator is automatically symmetric (stable).

Energy functional (W_ij = exp(-r_ij^2/2h^2), rho_i = sum_j W_ij, gamma=g(rho)):

    E =  1/2 sum pi_i^2                       (field kinetic)
       + 1/4 sum_ij (gamma_i+gamma_j) W_ij (u_i-u_j)^2     (field "spring" energy)
       + 1/2 m^2 sum u_i^2                    (field mass)
       + 1/2 sum |Wn_i|^2  +  U_LJ(X)         (medium kinetic + Lennard-Jones)

  d/du  -> field force = symmetric Laplacian with weights (gamma_i+gamma_j)W_ij
           (stable; wave speed ~ gamma(rho), so denser=slower = Coupling B)
  d/dX  -> node force (Coupling A), derived analytically below; the same
           gamma'(rho)<0 that slows waves also compresses the medium.

gamma(rho)=1/(1+beta*(rho/rho0-1)).  beta=0 decouples the density response (the
field still tugs the nodes through W(r); that residual is the control).
"""
from __future__ import annotations
import numpy as np

from bvc_core import lj_forces_energy, relax_medium, R0


class VariationalCoupled:
    def __init__(self, X, beta=0.0, m2=1.0, h=None, rcut=None, dt=0.004, damping=1.0):
        self.X = X.copy(); self.Wn = np.zeros_like(X); self.N = len(X)
        self.beta, self.m2 = beta, m2
        self.h = h or R0
        self.rcut = rcut or 2.0 * R0
        self.dt, self.damping = dt, damping
        self.u = np.zeros(self.N); self.pi = np.zeros(self.N)
        self.rho0 = float(np.median(self._geom(self.X)[2]))
        self.time = 0.0

    def _geom(self, X):
        d = X[:, None, :] - X[None, :, :]
        r2 = (d ** 2).sum(-1)
        within = (r2 > 1e-12) & (r2 < self.rcut ** 2)
        W = np.exp(-r2 / (2 * self.h ** 2)) * within
        rho = W.sum(1)
        return d, W, rho

    def _gamma(self, rho):
        # exp form: bounded & positive for all rho (denser=>slower), no singularity
        # the 1/(1+b(rho/rho0-1)) form blows up for rarefied nodes at large beta.
        g = np.exp(-self.beta * (rho / self.rho0 - 1.0))
        gp = -(self.beta / self.rho0) * g               # dg/drho
        return g, gp

    def forces(self):
        d, W, rho = self._geom(self.X)
        g, gp = self._gamma(rho)
        Q = (self.u[:, None] - self.u[None, :]) ** 2      # (u_i-u_j)^2
        A = (W * Q).sum(1)                                # A_i = sum_j W_ij Q_ij
        Gamma = g[:, None] + g[None, :]                   # gamma_i + gamma_j
        Wg = Gamma * W                                    # symmetric field weights
        Fu = (Wg @ self.u) - Wg.sum(1) * self.u - self.m2 * self.u   # field force
        # node force = -dU/dX = (1/2h^2) sum_j bracket_ij W_ij (X_i - X_j)
        bracket = Gamma * Q + gp[:, None] * A[:, None] + gp[None, :] * A[None, :]
        M = bracket * W
        FX_var = (1.0 / (2 * self.h ** 2)) * (M.sum(1)[:, None] * self.X - M @ self.X)
        Flj, Elj = lj_forces_energy(self.X)
        return Fu, FX_var + Flj, (d, W, rho, g, Q, Elj)

    def energy(self, cache=None):
        if cache is None:
            _, _, c = self.forces(); _, W, rho, g, Q, Elj = c
        else:
            _, W, rho, g, Q, Elj = cache
        Gamma = g[:, None] + g[None, :]
        U_spring = 0.25 * (Gamma * W * Q).sum()
        return (0.5 * (self.pi ** 2).sum() + U_spring + 0.5 * self.m2 * (self.u ** 2).sum()
                + 0.5 * (self.Wn ** 2).sum() + Elj)

    def step(self):
        Fu, FX, _ = self.forces()
        self.pi += 0.5 * Fu * self.dt;  self.Wn += 0.5 * FX * self.dt
        self.u += self.pi * self.dt;    self.X += self.Wn * self.dt
        Fu, FX, _ = self.forces()
        self.pi += 0.5 * Fu * self.dt;  self.Wn += 0.5 * FX * self.dt
        self.pi *= self.damping;        self.Wn *= self.damping
        self.time += self.dt

    def seed_lump(self, amp=1.0, width=4.0):
        self.u = amp * np.exp(-(self.X ** 2).sum(1) / (2 * width ** 2))
        self.pi = np.zeros(self.N)

    def lump_width(self):
        e = 0.5 * self.pi ** 2 + 0.5 * self.m2 * self.u ** 2
        s = e.sum()
        if s < 1e-9:
            return float("nan")
        c = (self.X * e[:, None]).sum(0) / s
        return float(np.sqrt(((self.X - c) ** 2).sum(1) @ e / s))


def run(beta, label, steps=4000, damping=1.0, dt=0.004):
    cloud = relax_medium(N=170, seed=3)
    f = VariationalCoupled(cloud, beta=beta, damping=damping, dt=dt)
    f.seed_lump()
    E0 = f.energy(); w0 = f.lump_width()
    traj = []
    for k in range(steps):
        f.step()
        if k % (steps // 6) == 0 or k == steps - 1:
            traj.append((round(f.time, 1), round(f.lump_width(), 2),
                         round(float(np.abs(f.u).max()), 2),
                         round(100 * (f.energy() - E0) / abs(E0), 2)))
    print(f"  {label}:")
    print("    " + "  ".join(f"t{t}:width={w},|u|={m},dE={e}%" for t, w, m, e in traj))
    wF = traj[-1][1]
    print(f"    width {w0:.2f}->{wF}  energy drift {traj[-1][3]:+.2f}%  =>  "
          f"{'DISPERSES' if wF > w0 + 0.5 else ('FOCUSES' if wF < w0 - 0.5 else 'holds')}\n")


if __name__ == "__main__":
    print("=== Phase 3 variational : one energy functional, conserved by construction ===")
    print("  key check: energy drift dE near 0 (vs 3c's leakage); then width behaviour\n")
    run(0.0, "control (beta=0)", dt=0.003, steps=6000)
    run(10.0, "coupled beta=10 (damping=1: energy check)", dt=0.003, steps=6000)
    run(10.0, "coupled beta=10 (mild damping: settle)", dt=0.003, steps=6000, damping=0.999)
