"""
bvc_core -- shared primitives for the Butler-Voss Condensate prototypes
=======================================================================

Single home for the pieces that the H10 / integration prototypes all need, so
they stop being copied (or imported prototype-from-prototype):

  * the Lennard-Jones medium and its self-assembly (`relax_medium`)
  * medium order parameters (spacing, coordination, hexagonal bond order)
  * lattice builders (`perfect_hex`, `square_lattice`)
  * meshfree differential operators on an irregular point cloud
    (weighted-graph and least-squares Laplacians)

Pure numpy.  The field engines themselves (ComplexFabric / FrozenField /
MovingMediumField) stay with their own scripts -- they are genuinely different
engines, not duplication.
"""
from __future__ import annotations
import math
import numpy as np


# ============================================================ Lennard-Jones ===
EPS, SIGMA, RCUT = 1.0, 1.0, 2.5
R0 = 2.0 ** (1.0 / 6.0) * SIGMA          # pair-potential minimum (target spacing)


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


# ========================================================= medium order params =
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


# =========================================================== medium assembly ===
def relax_medium(N=300, steps=6000, dt=0.005, cool=0.999, seed=1, report=False):
    """Self-assemble a cooled LJ droplet from a disordered cloud; return its
    (centered) node positions.  With report=True, log the ordering as it freezes
    (the H10 demo); otherwise run quiet (the integration medium)."""
    rng = np.random.default_rng(seed)
    area_per = 1.7                                   # > hex optimum (~0.97) => room to order
    R = math.sqrt(N * area_per / math.pi)
    X = []
    while len(X) < N:
        p = (rng.random(2) * 2 - 1) * R
        if math.hypot(*p) <= R and all(math.hypot(*(p - q)) > 0.95 for q in X):
            X.append(p)
    X = np.array(X)
    V = np.zeros_like(X)
    if report:
        print(f"   N={N}  target spacing r0={R0:.3f}")
        print(f"   {'t':>6} {'KE/N':>9} {'spacing':>8} {'coord':>6} {'psi6':>6}")
    F, _ = lj_forces_energy(X)
    for k in range(steps):
        X += V * dt + 0.5 * F * dt ** 2
        Fn, _ = lj_forces_energy(X)
        V += 0.5 * (F + Fn) * dt
        V *= cool
        F = Fn
        if report and (k % (steps // 8) == 0 or k == steps - 1):
            ke = 0.5 * (V ** 2).sum() / N
            print(f"   {k*dt:>6.1f} {ke:>9.4f} {mean_nn_spacing(X):>8.3f} "
                  f"{coordination(X).mean():>6.2f} {hex_order(X):>6.3f}")
    return X - X.mean(0)


# ================================================================= lattices ====
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


# ======================================================= meshfree operators ====
def pairwise(X):
    """(d, r2) where d[i,j] = X_i - X_j and r2 = |d|^2."""
    d = X[:, None, :] - X[None, :, :]
    return d, (d ** 2).sum(-1)


def laplacian_matrix(X, rcut):
    """Consistent weighted-graph Laplacian as a matrix A (so L f = A @ f).
    Normalised by 2D / sum_j w_ij |r_ij|^2 so it reproduces div-grad for smooth
    f on an isotropic neighbourhood (D = 2).  Cheap, but degrades badly off a
    regular lattice -- see lsq_laplacian_matrix."""
    d, r2 = pairwise(X)
    N = len(X)
    W = ((r2 > 1e-12) & (r2 < rcut ** 2)).astype(float)
    swr2 = (W * r2).sum(1)
    norm = (2 * 2) / np.where(swr2 > 0, swr2, 1.0)
    A = norm[:, None] * W
    A[np.arange(N), np.arange(N)] = -A.sum(1)
    return A


def lsq_laplacian_values(X, f, rcut):
    """Least-squares quadratic-fit Laplacian of f (the accurate reference)."""
    d, r2 = pairwise(X)
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
    d, r2 = pairwise(X)
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


def brookshaw_laplacian(X, h=None, rcut=None):
    """Symmetric SPH (Brookshaw) Laplacian: L f|_i = sum_j w_ij (f_j - f_i),
    w_ij = (2/h^2) * Vij * W(r_ij), symmetric node volumes Vij=(V_i+V_j)/2,
    V_i = 1/rho_i.  Symmetric (real spectrum) and negative-semidefinite (stable)
    by construction -- unlike the LSQ Laplacian, which is accurate but indefinite.
    A single global scale enforces consistency (reproduces lap(x^2+y^2)=4)."""
    h = h or R0
    rcut = rcut or 1.9 * R0
    _, r2 = pairwise(X)
    within = (r2 > 1e-12) & (r2 < rcut ** 2)
    W = np.exp(-r2 / (2 * h ** 2)) * within
    rho = W.sum(1) + 1.0
    V = 1.0 / rho
    Vij = 0.5 * (V[:, None] + V[None, :])
    w = (2.0 / h ** 2) * Vij * W
    L = w.copy()
    np.fill_diagonal(L, -w.sum(1))
    q = (X ** 2).sum(1)
    r = np.linalg.norm(X - X.mean(0), axis=1)
    m = r < 0.6 * r.max()
    return L * (4.0 / np.mean((L @ q)[m]))
