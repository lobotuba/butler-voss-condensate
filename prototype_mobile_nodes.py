"""
H10 -- Self-cohering medium: mobile, mutually-attracting nodes
=============================================================

PROTOTYPE.  In simulation.py and prototype_complex.py the nodes are FIXED and
only the field moves.  Robert's H10 conjecture: the nodes themselves should
attract each other (with short-range repulsion), so the medium self-organizes.
If true, three things should follow -- and this script checks each:

  1. The inter-node SPACING is no longer a free knob: it is set by the balance
     of attraction vs. short-range repulsion (the pair-potential minimum).
     -> answers "how much space exists between nodes?"
  2. The medium SELF-SELECTS an isotropic, close-packed arrangement (in 2D a
     triangular/hexagonal lattice; in 3D, fcc/hcp).  That is exactly the
     isotropic geometry we found numerically stable -- so the cubic/square
     anisotropy artifact would dissolve on its own.
  3. Mutual attraction makes the medium COHERE into a stable droplet rather
     than dispersing -> long-run stability.

The Lennard-Jones medium and its self-assembly live in bvc_core (`relax_medium`,
with report=True to log the ordering).  LJ minimum is at r0 = 2^(1/6) sigma; the
2D ground state is the triangular lattice -> 6 neighbors, hexagonal bond order.
"""
from __future__ import annotations
import numpy as np

from bvc_core import EPS, SIGMA, RCUT, R0, relax_medium


# ---------------------------------- why hexagonal: static energy comparison ---
def lattice_energy_comparison():
    """Energy per node of a perfect SQUARE vs TRIANGULAR(hex) patch at the
    spacing that minimizes each -- shows which arrangement the medium prefers."""
    def patch_energy(kind, a, m=7):
        pts = []
        for i in range(-m, m + 1):
            for j in range(-m, m + 1):
                if kind == "square":
                    pts.append((i * a, j * a))
                else:                                # triangular
                    pts.append((a * (i + 0.5 * (j & 1)), a * np.sqrt(3) / 2 * j))
        X = np.array(pts, float)
        c = np.array([0.0, 0.0])
        # energy of the CENTER node only (bulk), to avoid edge effects
        d = X - c
        r2 = (d ** 2).sum(1); r2 = r2[r2 > 1e-9]
        r2 = r2[r2 < RCUT ** 2]
        inv6 = (SIGMA ** 2 / r2) ** 3
        return float((4 * EPS * (inv6 ** 2 - inv6)).sum())

    print("   arrangement  best-spacing  energy/node (lower = preferred)")
    for kind in ("square", "triangular(hex)"):
        k = "square" if kind == "square" else "tri"
        aa = np.linspace(0.9 * R0, 1.3 * R0, 81)
        es = [patch_energy(k, a) for a in aa]
        i = int(np.argmin(es))
        print(f"   {kind:16} {aa[i]:>7.3f}     {es[i]:>9.4f}")
    print("   (2D ground state is triangular/hex -> 6 neighbors, isotropic;")
    print("    in 3D the same logic selects fcc/hcp ~ the stable fcc3d.)")


if __name__ == "__main__":
    print("=== H10 :: self-cohering medium of mutually-attracting nodes ===\n")
    print("1+2+3. Self-assembly from a disordered cloud:")
    relax_medium(N=200, seed=0, report=True)
    print("\nWhy it picks hexagonal (not square/cubic):")
    lattice_energy_comparison()
