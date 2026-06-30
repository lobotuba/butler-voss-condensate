"""
Symmetric field-operator scorecard
==================================

The Phase 3c blocker: the LSQ meshfree Laplacian is accurate but INDEFINITE
(max Re eigenvalue > 0) and asymmetric, so it blows up the field before the
gravity coupling can be pushed.  Here we implement symmetric candidates and score
each on the metrics that decide a field operator:

  symmetry  ||A - A^T|| / ||A||        (precondition for a real spectrum)
  STABILITY max Re(eigenvalue)         (<= 0 required; this IS the growth rate)
  complex   max |Im(eigenvalue)|       (0 for symmetric)
  ACCURACY  RMS error vs analytic Laplacian on the relaxed cloud
  CFL       dt < 2 / sqrt(|lambda_min|)
  + dynamic: energy drift & max|u| in a damping=1 field run (the end-to-end check)

Candidates:
  LSQ            -- current operator (accurate, indefinite)
  Brookshaw      -- symmetric SPH Laplacian, density/volume weighted
  stabilized-LSQ -- symmetrize LSQ then clip positive eigenvalues to 0 (NSD)
"""
from __future__ import annotations
import math
import numpy as np

from bvc_core import relax_medium, lsq_laplacian_matrix, pairwise, R0


# ----------------------------------------------------- candidate operators ----
def brookshaw_laplacian(X, h=R0, rcut=1.9 * R0):
    """Symmetric SPH (Brookshaw) Laplacian: L f|_i = sum_j w_ij (f_j - f_i),
    w_ij = (2/h^2) * Vij * W(r_ij), with symmetric node volumes Vij=(V_i+V_j)/2
    and V_i = 1/rho_i.  Symmetric (=> real spectrum) and negative-semidefinite
    (=> stable) by construction; a single global scale enforces consistency."""
    _, r2 = pairwise(X)
    within = (r2 > 1e-12) & (r2 < rcut ** 2)
    W = np.exp(-r2 / (2 * h ** 2)) * within
    rho = W.sum(1) + 1.0                       # SPH density (W(0)=1 self term)
    V = 1.0 / rho
    Vij = 0.5 * (V[:, None] + V[None, :])
    w = (2.0 / h ** 2) * Vij * W               # symmetric positive weights, w_ii=0
    L = w.copy()
    np.fill_diagonal(L, -w.sum(1))             # row sums to 0; NSD
    # consistency: scale so L reproduces lap(x^2+y^2)=4 on the interior
    q = (X ** 2).sum(1)
    r = np.linalg.norm(X - X.mean(0), axis=1)
    m = r < 0.6 * r.max()
    return L * (4.0 / np.mean((L @ q)[m]))


def stabilized_lsq(X, rcut=1.9 * R0):
    """Symmetrize the LSQ Laplacian, then project its spectrum onto (-inf, 0]
    (clip spurious positive eigenvalues).  Symmetric + NSD by construction."""
    A = lsq_laplacian_matrix(X, rcut)
    As = 0.5 * (A + A.T)
    w, Vv = np.linalg.eigh(As)
    w = np.minimum(w, 0.0)
    return (Vv * w) @ Vv.T


# ------------------------------------------------------------- scorecard ------
def score(X, A, name):
    asym = np.linalg.norm(A - A.T) / np.linalg.norm(A)
    ev = np.linalg.eigvals(A)
    lam_max, imax, lam_min = ev.real.max(), np.abs(ev.imag).max(), ev.real.min()
    k = 2 * math.pi / (8 * R0)
    f = np.cos(k * X[:, 0]); ana = -k ** 2 * f
    r = np.linalg.norm(X - X.mean(0), axis=1); m = r < 0.6 * r.max()
    Af, a = (A @ f)[m], ana[m]
    s = float(Af @ a / (Af @ Af))                 # best global scale: judge SHAPE, not calibration
    rms = math.sqrt(np.mean((s * Af - a) ** 2) / np.mean(a ** 2))
    dt = 2 / math.sqrt(abs(lam_min)) if lam_min < 0 else float("inf")
    stab = "STABLE" if lam_max <= 1e-6 else "unstable"
    print(f"  {name:16} asym={asym:7.1e}  maxRe={lam_max:+7.3f} [{stab:8}]  "
          f"max|Im|={imax:5.2f}  acc-RMS={rms*100:5.1f}%  CFL dt<{dt:.2f}")


def dynamic(X, A, name, c=1.0, m2=0.2, dt=0.1, steps=500, seed=0):
    """Field run at damping=1 (conservative): does it blow up, and is the
    discrete energy conserved?  A good operator -> bounded |u|, ~0 energy drift."""
    rng = np.random.default_rng(seed)
    u = np.exp(-(X ** 2).sum(1) / (2 * (3 * R0) ** 2)) * (1 + 0.05 * rng.standard_normal(len(X)))
    ud = np.zeros(len(X))

    def energy():
        return 0.5 * (ud ** 2).sum() + 0.5 * c ** 2 * (-(u @ (A @ u))) + 0.5 * m2 * (u ** 2).sum()

    E0 = energy(); umax = abs(u).max()
    for _ in range(steps):
        ud += (c ** 2 * (A @ u) - m2 * u) * dt
        u += ud * dt
        umax = max(umax, float(np.abs(u).max()))
        if umax > 1e3:
            break
    drift = 100 * (energy() - E0) / abs(E0)
    print(f"  {name:16} max|u| {abs(np.exp(0)):.0f}->{umax:8.2f}   energy drift {drift:+8.1f}%   "
          f"{'BLOWS UP' if umax > 50 else 'bounded'}")


if __name__ == "__main__":
    print("=== Symmetric field-operator scorecard ===\n")
    X = relax_medium(N=220, seed=3)
    ops = {
        "LSQ (current)": lsq_laplacian_matrix(X, 1.9 * R0),
        "Brookshaw": brookshaw_laplacian(X),
        "stabilized-LSQ": stabilized_lsq(X),
    }
    print("STATIC scorecard (relaxed cloud):")
    for name, A in ops.items():
        score(X, A, name)
    print("\nDYNAMIC test (damping=1, m2=0.2, to t=50):")
    for name, A in ops.items():
        dynamic(X, A, name)
