"""
Butler-Voss Condensate -- COMPLEX-FIELD PROTOTYPE  (charge & spin from texture)
===============================================================================

PROTOTYPE / not yet merged into simulation.py.  Purpose: promote the node field
from a real scalar u to a COMPLEX order parameter psi, so the medium can carry
the two kinds of "charge" Robert wants to distinguish (the H6/H7 program):

  * NOETHER charge  Q = Im( sum conj(psi) * psi_dot )   -- charge from the U(1)
        symmetry of the potential (a Q-ball; conserved, continuous). This is
        "charge from the dynamics."
  * TOPOLOGICAL charge  n = (1/2pi) * sum of phase winding around each plaquette
        -- an INTEGER vortex number; conserved because you cannot unwind it.
        This is "charge from the SHAPE of the between-space" (H6).

Why a Mexican-hat vacuum: a protected winding needs the phase defined in the
vacuum, i.e. |psi| = v0 > 0.  The original "oscillon" potential (vacuum at 0)
supports non-topological lumps (Noether charge only); the "mexicanhat" potential
supports vortices (quantized winding).  Both are provided so the difference is
visible -- and so the particle census (H8) can sort structures by sector.

This prototype runs on a square 2D lattice (clean plaquettes for the winding
sum); the same field can later be carried on the hex/fcc lattices of the main
engine.  Dynamics mirror simulation.py: graph Laplacian + symplectic step.
"""
from __future__ import annotations
import math
import numpy as np


# ---------------------------------------------------------------- lattice ----
def square_lattice(rows, cols):
    """4-neighbor square grid; returns pos, neighbor_idx, valid, plaquettes."""
    N = rows * cols
    idx = lambda r, c: r * cols + c
    pos = np.array([(c, r) for r in range(rows) for c in range(cols)], float)
    nidx = np.zeros((N, 4), np.int64)
    valid = np.zeros((N, 4), float)
    offs = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for r in range(rows):
        for c in range(cols):
            i = idx(r, c)
            for k, (dr, dc) in enumerate(offs):
                rr, cc = r + dr, c + dc
                if 0 <= rr < rows and 0 <= cc < cols:
                    nidx[i, k] = idx(rr, cc); valid[i, k] = 1.0
                else:
                    nidx[i, k] = i
    # plaquettes: the 4 corners of each unit cell, in loop order (CCW)
    plaq = [(idx(r, c), idx(r, c + 1), idx(r + 1, c + 1), idx(r + 1, c))
            for r in range(rows - 1) for c in range(cols - 1)]
    return pos, nidx, valid, np.array(plaq, np.int64)


# ------------------------------------------------------------------ field ----
class ComplexFabric:
    def __init__(self, rows=72, cols=72, potential="mexicanhat",
                 v0=1.0, lam=1.0, mass2=1.0, focus=1.0, saturation=0.5,
                 dt=0.02, c=1.0, damping=1.0):
        self.pos, self.nidx, self.valid, self.plaq = square_lattice(rows, cols)
        self.N = self.pos.shape[0]
        self.rows, self.cols = rows, cols
        self.deg_norm = float(self.valid.sum(axis=1).max())
        self.ncount = np.clip(self.valid.sum(axis=1), 1.0, None)
        self.potential, self.v0, self.lam = potential, v0, lam
        self.mass2, self.focus, self.sat = mass2, focus, saturation
        self.dt, self.c, self.damping = dt, c, damping
        # vacuum: v0 for the broken potential, 0 for the oscillon
        base = v0 if potential == "mexicanhat" else 0.0
        self.psi = np.full(self.N, base, dtype=np.complex128)
        self.psid = np.zeros(self.N, dtype=np.complex128)
        self.time = 0.0

    # -- operators ----------------------------------------------------------
    def laplacian(self, f):
        neigh = (f[self.nidx] * self.valid).sum(axis=1)
        return (neigh - self.ncount * f) / self.deg_norm

    def potential_force(self):
        rho = np.abs(self.psi) ** 2
        if self.potential == "mexicanhat":
            return -self.lam * (rho - self.v0 ** 2) * self.psi
        # oscillon (vacuum at 0): -m^2 + focus*rho - sat*rho^2
        return (-self.mass2 + self.focus * rho - self.sat * rho ** 2) * self.psi

    def potential_energy(self):
        rho = np.abs(self.psi) ** 2
        if self.potential == "mexicanhat":
            return 0.25 * self.lam * (rho - self.v0 ** 2) ** 2
        return 0.5 * self.mass2 * rho - 0.25 * self.focus * rho ** 2 \
            + (self.sat / 6.0) * rho ** 3

    # -- dynamics (symplectic, identical scheme to simulation.py) -----------
    def step(self):
        accel = self.c ** 2 * self.laplacian(self.psi) + self.potential_force()
        self.psid += accel * self.dt
        self.psid *= self.damping
        self.psi += self.psid * self.dt
        self.time += self.dt

    # -- conserved quantities ----------------------------------------------
    def noether_charge(self):
        """U(1) charge Q = Im( sum conj(psi) psi_dot ).  Continuous, conserved."""
        return float(np.imag(np.sum(np.conj(self.psi) * self.psid)))

    def topological_charge(self, per_plaquette=False):
        """Integer winding summed over plaquettes (the vortex number)."""
        p = self.psi[self.plaq]                       # (P,4)
        nxt = np.roll(p, -1, axis=1)
        dphi = np.angle(nxt * np.conj(p))             # wrapped to (-pi,pi]
        wind = np.rint(dphi.sum(axis=1) / (2 * math.pi)).astype(int)
        if per_plaquette:
            return wind
        return int(wind.sum())

    def energy(self):
        kin = 0.5 * np.abs(self.psid) ** 2
        diff = (self.psi[self.nidx] - self.psi[:, None]) * self.valid
        grad = 0.5 * np.sum(np.abs(diff) ** 2, axis=1) / self.deg_norm
        return float(np.sum(kin + grad + self.potential_energy()))

    def localized_energy_fraction(self, radius_frac=0.25):
        """Fraction of (energy above vacuum) sitting in the central region --
        a crude 'did it stay a localized particle?' measure."""
        e = (0.5 * np.abs(self.psid) ** 2
             + self.potential_energy()
             - self.potential_energy_vacuum())
        e = np.clip(e, 0, None)
        cen = self.pos.mean(axis=0)
        r = np.linalg.norm(self.pos - cen, axis=1)
        R = radius_frac * r.max()
        tot = e.sum()
        return float(e[r < R].sum() / tot) if tot > 1e-9 else 0.0

    def potential_energy_vacuum(self):
        if self.potential == "mexicanhat":
            return 0.0
        return 0.0

    # -- seeding ------------------------------------------------------------
    def _center(self, fx, fy):
        lo, hi = self.pos.min(0), self.pos.max(0)
        return lo + np.array([fx, fy]) * (hi - lo)

    def set_vortices(self, vortices, core=3.0):
        """vortices = [(fx, fy, n), ...]; build psi with those windings on a
        |psi|->v0 background.  Total topological charge = sum of n."""
        psi = np.full(self.N, self.v0, dtype=np.complex128)
        for fx, fy, n in vortices:
            cx, cy = self._center(fx, fy)
            dx, dy = self.pos[:, 0] - cx, self.pos[:, 1] - cy
            r = np.hypot(dx, dy)
            phi = np.arctan2(dy, dx)
            psi = psi * np.tanh(r / core) * np.exp(1j * n * phi)
        self.psi, self.psid = psi, np.zeros(self.N, np.complex128)

    def set_qball(self, fx=0.5, fy=0.5, amp=1.4, sigma=5.0, omega=1.0):
        """Non-topological lump with internal rotation -> Noether charge,
        zero winding (only meaningful for the oscillon potential)."""
        cx, cy = self._center(fx, fy)
        r2 = (self.pos[:, 0] - cx) ** 2 + (self.pos[:, 1] - cy) ** 2
        env = amp * np.exp(-r2 / (2 * sigma ** 2))
        self.psi = env.astype(np.complex128)
        self.psid = (1j * omega) * self.psi          # rotating phase -> charge

    def set_breather(self, fx=0.5, fy=0.5, amp=0.6, sigma=5.0):
        """Neutral bump on the v0 vacuum: charge 0, no winding (a 'meson')."""
        cx, cy = self._center(fx, fy)
        r2 = (self.pos[:, 0] - cx) ** 2 + (self.pos[:, 1] - cy) ** 2
        self.psi = (self.v0 + amp * np.exp(-r2 / (2 * sigma ** 2))).astype(np.complex128)
        self.psid = np.zeros(self.N, np.complex128)


# ----------------------------------------------------- census demonstration --
def census(steps=800):
    """Seed several candidate 'species' and report their conserved quantities
    and persistence -- a first pass at H8 (how many distinct stable types?)."""
    specs = [
        ("breather  n=0",  "mexicanhat", lambda f: f.set_breather()),
        ("vortex    n=+1", "mexicanhat", lambda f: f.set_vortices([(0.5, 0.5, +1)])),
        ("antivortex n=-1", "mexicanhat", lambda f: f.set_vortices([(0.5, 0.5, -1)])),
        ("v-pair   n=+1,-1", "mexicanhat",
         lambda f: f.set_vortices([(0.35, 0.5, +1), (0.65, 0.5, -1)])),
        ("vortex    n=+2", "mexicanhat", lambda f: f.set_vortices([(0.5, 0.5, +2)])),
        ("Q-ball    (osc)", "oscillon",  lambda f: f.set_qball()),
    ]
    print(f"{'species':18}{'pot':11}{'topo0':>6}{'topoF':>6}{'Noether0':>10}"
          f"{'NoetherF':>10}{'E0':>9}{'EF':>9}{'loc%F':>7}")
    for name, pot, seed in specs:
        f = ComplexFabric(potential=pot, damping=1.0)
        seed(f)
        t0, q0, e0 = f.topological_charge(), f.noether_charge(), f.energy()
        for _ in range(steps):
            f.step()
        tF, qF, eF = f.topological_charge(), f.noether_charge(), f.energy()
        loc = f.localized_energy_fraction()
        print(f"{name:18}{pot:11}{t0:>6}{tF:>6}{q0:>10.3f}{qF:>10.3f}"
              f"{e0:>9.2f}{eF:>9.2f}{loc:>7.2f}")


if __name__ == "__main__":
    print("=== Butler-Voss Condensate :: complex-field prototype ===")
    print("Census of candidate species (H8): conserved charges + persistence\n")
    census()
    print("\nLegend: topo0/topoF = topological winding (start/end), conserved & integer;")
    print("        Noether0/F = U(1) charge; E0/EF = total energy;")
    print("        loc%F = fraction of excess energy still central at the end.")
