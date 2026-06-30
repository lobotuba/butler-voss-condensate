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

Mechanics: classic Lennard-Jones molecular dynamics (velocity-Verlet) with
cooling, on a free 2D cluster (LJ attraction binds it, so no box is needed).
LJ minimum is at r0 = 2^(1/6) sigma; the 2D ground state is the triangular
lattice -> 6 neighbors, hexagonal bond order.
"""
from __future__ import annotations
import numpy as np

EPS, SIGMA, RCUT = 1.0, 1.0, 2.5
R0 = 2.0 ** (1.0 / 6.0) * SIGMA          # pair-potential minimum (target spacing)


# ---------------------------------------------------------- LJ interactions --
def lj_forces_energy(X, fcap=200.0):
    """Vectorized pairwise Lennard-Jones forces (N,2) and total energy."""
    d = X[:, None, :] - X[None, :, :]               # (N,N,2)  r_i - r_j
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)                    # no self-interaction
    within = r2 < RCUT ** 2
    inv2 = np.where(within, SIGMA ** 2 / r2, 0.0)
    inv6 = inv2 ** 3
    inv12 = inv6 ** 2
    # F(r)/r = 24 eps (2 inv12 - inv6) / r^2  -> force vector = (F/r) * d
    coef = 24 * EPS * (2 * inv12 - inv6) / np.where(within, r2, 1.0)
    coef = np.clip(coef, -fcap, fcap)
    F = (coef[:, :, None] * d).sum(axis=1)          # (N,2)
    E = 0.5 * (4 * EPS * (inv12 - inv6))[within].sum()
    return F, float(E)


# --------------------------------------------------------------- order params -
def neighbors_within(X, rmax):
    d = X[:, None, :] - X[None, :, :]
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    return r2, (r2 < rmax ** 2)


def mean_nn_spacing(X):
    r2, _ = neighbors_within(X, np.inf)
    return float(np.sqrt(r2.min(axis=1)).mean())


def coordination(X, rmax=1.35 * R0):
    _, near = neighbors_within(X, rmax)
    return near.sum(axis=1)


def hex_order(X, rmax=1.35 * R0):
    """Global 2D bond-orientational order psi6 (|<exp(6 i theta)>|).
    ~1 for a triangular/hexagonal lattice, ~0 for square or disordered."""
    d = X[:, None, :] - X[None, :, :]
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    near = r2 < rmax ** 2
    theta = np.arctan2(d[:, :, 1], d[:, :, 0])
    psi = np.zeros(len(X), dtype=complex)
    coord = near.sum(axis=1)
    for i in range(len(X)):
        if coord[i] >= 3:                            # interior nodes only
            psi[i] = np.mean(np.exp(6j * theta[i, near[i]]))
    interior = coord >= 3
    return float(np.abs(psi[interior].mean())) if interior.any() else 0.0


# ------------------------------------------------------------- the simulation -
def self_assemble(N=200, steps=6000, dt=0.005, cool=0.999, seed=0):
    rng = np.random.default_rng(seed)
    # disordered start: random in a disk, sized looser than close packing,
    # with a minimum separation so we don't start inside the LJ wall.
    area_per = 1.7                                   # > hex optimum (~0.97) => room to order
    R = np.sqrt(N * area_per / np.pi)
    X = []
    while len(X) < N:
        p = (rng.random(2) * 2 - 1) * R
        if np.hypot(*p) <= R and all(np.hypot(*(p - q)) > 0.95 for q in X):
            X.append(p)
    X = np.array(X)
    V = np.zeros_like(X)

    print(f"   N={N}  target spacing r0={R0:.3f}")
    print(f"   {'t':>6} {'KE/N':>9} {'spacing':>8} {'coord':>6} {'psi6':>6}")
    F, _ = lj_forces_energy(X)
    snap = lambda k: k % (steps // 8) == 0 or k == steps - 1
    for k in range(steps):
        X += V * dt + 0.5 * F * dt ** 2
        Fnew, E = lj_forces_energy(X)
        V += 0.5 * (F + Fnew) * dt
        V *= cool
        F = Fnew
        if snap(k):
            ke = 0.5 * (V ** 2).sum() / N
            print(f"   {k*dt:>6.1f} {ke:>9.4f} {mean_nn_spacing(X):>8.3f} "
                  f"{coordination(X).mean():>6.2f} {hex_order(X):>6.3f}")
    return X


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
    self_assemble()
    print("\nWhy it picks hexagonal (not square/cubic):")
    lattice_energy_comparison()
