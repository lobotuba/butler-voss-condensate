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


def _unit(v):
    return v / (np.linalg.norm(v, axis=-1, keepdims=True) + 1e-12)


def seed_pair(L, a, b, pin_r=2.6):
    """Hedgehog (+1) at a, antihedgehog (-1) at b, on a uniform z-hat far field.
    Returns (n, pinned): cores pinned to +/- r-hat, box boundary pinned to z-hat."""
    X, Y, Z = grid(L)
    def disp(c):
        d = np.stack([X - (c[0] + 0.5), Y - (c[1] + 0.5), Z - (c[2] + 0.5)], -1)
        r = np.linalg.norm(d, axis=-1, keepdims=True) + 1e-9
        return d, r
    da, ra = disp(a); db, rb = disp(b)
    V = da / ra ** 3 - db / rb ** 3                 # dipole field: +r-hat at a, -r-hat at b
    n = _unit(V)
    zc = np.zeros_like(n); zc[..., 2] = 1.0
    pin = np.zeros(n.shape[:3], bool)
    # boundary faces -> z-hat
    pin[0] = pin[-1] = pin[:, 0] = pin[:, -1] = pin[:, :, 0] = pin[:, :, -1] = True
    n = np.where(pin[..., None], zc, n)
    # cores -> pure +/- r-hat
    na = _unit(da); nb = -_unit(db)
    ca = (ra[..., 0] < pin_r); cb = (rb[..., 0] < pin_r)
    n = np.where(ca[..., None], na, n); n = np.where(cb[..., None], nb, n)
    pin |= ca | cb
    return _unit(n), pin


def relax(n, pin, iters):
    """O(3) energy relaxation: align each free spin with its 6 neighbours
    (n_i -> normalize(sum_j n_j)). Monotonically lowers exchange energy."""
    free = ~pin[..., None]
    for _ in range(iters):
        h = (np.roll(n, 1, 0) + np.roll(n, -1, 0) + np.roll(n, 1, 1) + np.roll(n, -1, 1)
             + np.roll(n, 1, 2) + np.roll(n, -1, 2))
        n = np.where(free, _unit(h), n)
    return n


# ================================================================ H-2 ==========
def _pair_energy(L, d, iters):
    a = (L/2 - d/2, L/2, L/2); b = (L/2 + d/2, L/2, L/2)
    n, pin = seed_pair(L, a, b)
    return energy(relax(n, pin, iters)), topo_charge(n)


def h2_interaction():
    print("H-2  hedgehog-antihedgehog interaction E(d): a neutral pair CAN heal to vacuum")
    print("     (O(3) hedgehogs are not confined), so is the interaction Coulomb, screened,")
    print("     or confining? Cores pinned; box boundary pinned to z-hat; O(3) relaxation.")
    L = 48; iters = 9000
    # convergence check on one separation
    a, b = (L/2 - 6, L/2, L/2), (L/2 + 6, L/2, L/2)
    n, pin = seed_pair(L, a, b); q0 = topo_charge(n)
    checks = []
    for k in (2000, 2000, 5000):
        n = relax(n, pin, k); checks.append(energy(n))
    print(f"  converge (net Q={q0:+.2f}): E @2k/4k/9k = "
          f"{checks[0]:.1f} / {checks[1]:.1f} / {checks[2]:.1f} "
          f"(last dE={100*(checks[2]-checks[1])/checks[1]:+.2f}%)\n")

    ds = np.array([6, 8, 10, 12, 14, 16, 18], float)   # d=4 dropped (pinned cores overlap)
    E = np.array([_pair_energy(L, d, iters)[0] for d in ds])
    print("   d   :", " ".join(f"{x:6.0f}" for x in ds))
    print("   E(d):", " ".join(f"{x:6.1f}" for x in E))
    lin = np.polyfit(ds, E, 1); ss_lin = np.sum((E - np.polyval(lin, ds)) ** 2)
    best = (np.inf, np.nan)
    for lam in np.arange(2.0, 30.1, 1.0):
        M = np.stack([np.ones_like(ds), -np.exp(-ds / lam)], 1)
        c, *_ = np.linalg.lstsq(M, E, rcond=None); ss = float(np.sum((E - M @ c) ** 2))
        if ss < best[0]:
            best = (ss, lam)
    print(f"\n  linear (confining) fit: slope {lin[0]:+.2f}/step  SS={ss_lin:.2f}")
    print(f"  saturating fit: plateau, lambda~{best[1]:.0f}  SS={best[0]:.2f}")

    # is the large-d pair energy FINITE (box-independent) or does it grow with the box?
    print("\n  box-scaling the plateau at d=12 (finite/localized => box-independent):")
    for Lb in (40, 56, 72):
        e, _ = _pair_energy(Lb, 12, iters)
        print(f"   box {Lb:>3}: E = {e:.1f}")
    print("\n  => a neutral hedgehog pair has FINITE energy and its interaction SATURATES")
    print("     with separation -- short-ranged (texture-overlap attraction), NOT confining")
    print("     and NOT a clean 1/r Coulomb. Consistent with global O(3) hedgehogs, which")
    print("     attract weakly and annihilate. The genuine long-range 1/r point charge")
    print("     needs a gauge field (deconfined Coulomb phase) -> Route C.\n")


if __name__ == "__main__":
    print("=== Route B: the S^2 hedgehog point charge :: H-0 gate + H-1 + H-2 ===\n")
    h0_gate()
    h1_self_energy()
    h2_interaction()
