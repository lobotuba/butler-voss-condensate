"""
Integration Phase 2 : the medium MOVES while the field rides on it
==================================================================

PROTOTYPE.  Phase 1 put a field on a FROZEN self-assembled medium.  Phase 2 lets
the nodes keep moving (one-way coupling: the medium pushes the field, the field
does NOT yet push the medium) and asks the decisive question:

    does the topological winding stay quantized when the lattice actively
    REARRANGES -- i.e. when nodes swap neighbors (reconnections)?

Design:
  * Field is LAGRANGIAN: psi_i is attached to node i and carried as it moves.
  * The LSQ meshfree Laplacian (Phase 0 showed it is required) is rebuilt from
    the current positions every `rebuild_every` steps.
  * A thermostat heats the medium to temperature T so it genuinely rearranges;
    we COUNT reconnections to prove the mesh topology is changing, then check the
    winding and whether a vortex core still exists anywhere.

Topology says winding can only change if the field zero leaves the droplet or
pair-annihilates -- thermal jostling and reconnections should NOT be able to
unwind it.  Phase 2 tests that numerically.
"""
from __future__ import annotations
import math
import numpy as np
from bvc_core import lj_forces_energy, relax_medium, lsq_laplacian_matrix, pairwise, R0


class MovingMediumField:
    def __init__(self, X, T=0.0, dt=0.005, field_damping=0.999,
                 rcut_field=1.9 * R0, rcut_bond=1.45 * R0, rebuild_every=20,
                 v0=1.0, lam=1.0, c=1.0, seed=0):
        self.X = X.copy()
        self.V = np.zeros_like(X)
        self.N = len(X)
        self.T, self.dt, self.field_damping = T, dt, field_damping
        self.rcut_field, self.rcut_bond, self.rebuild_every = rcut_field, rcut_bond, rebuild_every
        self.v0, self.lam, self.c = v0, lam, c
        self.rng = np.random.default_rng(seed)
        self.A = lsq_laplacian_matrix(self.X, rcut_field)
        self.psi = np.full(self.N, v0, dtype=np.complex128)
        self.psid = np.zeros(self.N, dtype=np.complex128)
        self.Fmed, _ = lj_forces_energy(self.X)
        self.bonds = self._bond_matrix()
        self.reconnections = 0
        self.time = 0.0
        # seed thermal velocities at temperature T
        if T > 0:
            self.V = self.rng.normal(0, math.sqrt(T), size=X.shape)
            self.V -= self.V.mean(0)

    # -- medium bookkeeping --------------------------------------------------
    def _bond_matrix(self):
        _, r2 = pairwise(self.X)
        return (r2 > 1e-12) & (r2 < self.rcut_bond ** 2)

    def _thermostat(self):
        self.V -= self.V.mean(0)                       # remove centre-of-mass drift
        if self.T > 0:
            ke = 0.5 * (self.V ** 2).sum() / self.N
            if ke > 1e-9:
                self.V *= math.sqrt(self.T / ke)       # isokinetic rescale

    # -- seeding -------------------------------------------------------------
    def seed_vortex(self, n=1, core=2.5):
        c = self.X.mean(0)
        d = self.X - c
        r = np.linalg.norm(d, axis=1)
        phi = np.arctan2(d[:, 1], d[:, 0])
        self.psi = self.v0 * np.tanh(r / core) * np.exp(1j * n * phi)
        self.psid = np.zeros_like(self.psi)

    # -- combined step -------------------------------------------------------
    def step(self, k):
        # medium velocity-Verlet under Lennard-Jones
        self.X += self.V * self.dt + 0.5 * self.Fmed * self.dt ** 2
        Fnew, _ = lj_forces_energy(self.X)
        self.V += 0.5 * (self.Fmed + Fnew) * self.dt
        self.Fmed = Fnew
        self._thermostat()
        # periodically rebuild operator + count reconnections
        if k % self.rebuild_every == 0:
            nb = self._bond_matrix()
            self.reconnections += int(np.triu(nb != self.bonds, 1).sum())
            self.bonds = nb
            self.A = lsq_laplacian_matrix(self.X, self.rcut_field)
        # field step on the current (moving) geometry; psi rides with the nodes
        rho = np.abs(self.psi) ** 2
        force = -self.lam * (rho - self.v0 ** 2) * self.psi
        accel = self.c ** 2 * (self.A @ self.psi) + force
        self.psid += accel * self.dt
        self.psid *= self.field_damping
        self.psi += self.psid * self.dt
        self.time += self.dt

    # -- diagnostics ---------------------------------------------------------
    def _winding_at(self, frac, band=1.2):
        c = self.X.mean(0)
        d = self.X - c
        r = np.linalg.norm(d, axis=1)
        R = frac * r.max()
        sel = np.where(np.abs(r - R) < band * R0)[0]
        if len(sel) < 8:
            return None
        loop = sel[np.argsort(np.arctan2(d[sel, 1], d[sel, 0]))]
        p = self.psi[loop]
        return int(round(np.angle(np.roll(p, -1) * np.conj(p)).sum() / (2 * math.pi)))

    def winding(self):
        vals = [w for f in (0.40, 0.50, 0.60) if (w := self._winding_at(f)) is not None]
        return int(np.median(vals)) if vals else None

    def core_exists(self):
        """global min |psi|/v0: ~0 if a vortex core exists somewhere, ~1 if the
        field has healed (vortex gone)."""
        return float(np.abs(self.psi).min() / self.v0)

    def gyration(self):
        d = self.X - self.X.mean(0)
        return float(np.sqrt((d ** 2).sum(1).mean()))


def run(T, label, steps=4000, seed=0):
    cloud = relax_medium(seed=1)
    f = MovingMediumField(cloud, T=T, seed=seed)
    f.seed_vortex(n=1)
    Rg0 = f.gyration()
    traj = []
    for k in range(steps):
        f.step(k)
        if k % (steps // 6) == 0 or k == steps - 1:
            traj.append((round(f.time, 1), f.winding(), round(f.core_exists(), 2),
                         f.reconnections))
    w0, wF = traj[0][1], traj[-1][1]
    maj = sum(1 for _, w, _, _ in traj if w == 1) / len(traj)
    cohesive = abs(f.gyration() - Rg0) / Rg0 < 0.25
    print(f"  {label} (T={T}):")
    print("    " + "  ".join(f"t{t}:w={w},core={c},reconn={n}" for t, w, c, n in traj))
    charge_ok = (w0 == 1 and wF == 1 and maj >= 0.7)
    print(f"    reconnections={f.reconnections}  cohesive={cohesive}  =>  charge "
          f"{'CONSERVED through a rearranging medium' if charge_ok else 'LOST'}\n")


if __name__ == "__main__":
    print("=== Phase 2 : winding survival on a MOVING, rearranging medium ===\n")
    print("  metrics: w=winding, core=global min|psi|/v0 (0=core exists), reconn=cumulative bond swaps\n")
    run(0.00, "frozen control", steps=4000)
    run(0.12, "warm (rearranging)", steps=4000)
    run(0.30, "hot (near melting)", steps=4000)
