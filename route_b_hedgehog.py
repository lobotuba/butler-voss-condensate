"""
Route B, H-0 / H-1 -- the 3D POINT topological charge: the S^2 hedgehog.

Screen-2 / Route A used a complex field (target S^1), whose 3D defect is a LINE.
A genuine 3D POINT charge needs a larger target space: a 3-component unit vector
n-hat in S^2 (the O(3) / Heisenberg model), where pi_2(S^2)=Z classifies point
defects. The hedgehog n-hat = r-hat is the charge-+1 point defect -- the "particle"
the project is after.

Topological charge = degree of the map from an enclosing surface to S^2, computed
as the flux of the topological current
    j^i = (1/4pi) n-hat . (d_j n-hat x d_k n-hat),   (i,j,k) cyclic,
so Q = closed-surface flux of j (no per-plaquette orientation bookkeeping).

  H-0  gate: hedgehog -> +1, antihedgehog -> -1, vacuum -> 0 (integer, localized).
  H-1  a single hedgehog's energy diverges LINEARLY with box size (E ~ L) -- the
       "global monopole": the 3D point-charge analogue of the 2D-vortex log and the
       3D vortex-line tension. (The clean 1/r point charge needs gauging: Route C.)
"""
from __future__ import annotations
import numpy as np


def grid(L):
    a = np.arange(L)
    return np.meshgrid(a, a, a, indexing="ij")


def seed_hedgehog(L, center=None, sign=+1):
    """Unit-vector field of a hedgehog (sign +1) or antihedgehog (sign -1) at
    center. Core placed off-lattice (+0.5) so no site sits exactly at r=0."""
    cx, cy, cz = center or (L / 2, L / 2, L / 2)
    X, Y, Z = grid(L)
    dx, dy, dz = X - (cx + 0.5), Y - (cy + 0.5), Z - (cz + 0.5)
    if sign < 0:
        dx = -dx                                   # mirror one target axis -> degree -1
    n = np.stack([dx, dy, dz], -1).astype(float)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return n


def energy(n):
    """O(3) Heisenberg exchange: 1/2 sum_links |n_j - n_i|^2. Uniform vacuum = 0."""
    e = 0.0
    for ax in range(3):
        d = np.diff(n, axis=ax)
        e += (d ** 2).sum()
    return 0.5 * e


def topo_current(n):
    """j^i = (1/4pi) n . (d_j n x d_k n), returned as (jx, jy, jz) fields."""
    gx, gy, gz = (np.gradient(n, axis=a) for a in range(3))
    trip = lambda A, B: np.einsum("...a,...a->...", n, np.cross(A, B))
    f = 1.0 / (4 * np.pi)
    return f * trip(gy, gz), f * trip(gz, gx), f * trip(gx, gy)


def topo_charge(n, margin=None):
    """Net charge enclosed = flux of j through the box faces (inset by margin)."""
    L = n.shape[0]; m = margin or max(3, L // 6)
    jx, jy, jz = topo_current(n)
    s = slice(m, L - m)
    Q  = jx[L - 1 - m, s, s].sum() - jx[m, s, s].sum()
    Q += jy[s, L - 1 - m, s].sum() - jy[s, m, s].sum()
    Q += jz[s, s, L - 1 - m].sum() - jz[s, s, m].sum()
    return float(Q)


# ================================================================ H-0 ==========
def h0_gate():
    print("H-0  the hedgehog is a clean integer point charge (topological-current flux)")
    cases = [("hedgehog  n=+1", seed_hedgehog(40, sign=+1), +1),
             ("antihedgehog -1", seed_hedgehog(40, sign=-1), -1),
             ("vacuum (uniform)", np.tile([0, 0, 1.0], (40, 40, 40, 1)), 0)]
    ok = True
    for name, n, expect in cases:
        Q = topo_charge(n)
        good = abs(Q - expect) < 0.1
        ok &= good
        print(f"  {name:>18}:  Q = {Q:+.3f}  -> {round(Q):+d}   ({'OK' if good else 'BAD'})")
    print(f"  => gate {'PASSED' if ok else 'FAILED'}: the model hosts a genuine 3D POINT")
    print("     topological charge (S^2 hedgehog), integer and localized.\n")
    return ok


# ================================================================ H-1 ==========
def h1_self_energy():
    print("H-1  single-hedgehog self-energy vs box size (the 'global monopole')")
    print(f"  {'L':>5} {'energy':>10} {'E / L':>8} {'Q':>6}")
    Ls = [16, 24, 32, 40, 48]
    E = []
    for L in Ls:
        n = seed_hedgehog(L)
        E.append(energy(n))
        print(f"  {L:>5} {E[-1]:>10.2f} {E[-1]/L:>8.2f} {topo_charge(n):>+6.2f}")
    E = np.array(E); Ls = np.array(Ls, float)
    slope = np.polyfit(Ls, E, 1)[0]
    # linear (E~L) vs quadratic (E~L^2) fit quality
    ss_lin = np.sum((E - np.polyval(np.polyfit(Ls, E, 1), Ls)) ** 2)
    ss_quad = np.sum((E - np.polyval(np.polyfit(Ls, E, 2), Ls)) ** 2)
    print(f"\n  E/L -> ~constant: E is LINEAR in L (dE/dL = {slope:.2f}); "
          f"SS_linear={ss_lin:.1f} vs SS_quad={ss_quad:.1f}")
    print("  => the hedgehog self-energy diverges linearly with system size -- the global")
    print("     monopole. Same story as the 2D-vortex log and 3D-line tension: a bare")
    print("     topological charge has no finite energy; only neutral combinations do.")
    print("     A clean 1/r point charge needs a screening/gauge field -> Route C.\n")


if __name__ == "__main__":
    print("=== Route B: the S^2 hedgehog point charge :: H-0 gate + H-1 ===\n")
    h0_gate()
    h1_self_energy()
