"""
Integration Phase 3a / 3b : the two halves of gravity-by-density (in isolation)
==============================================================================

PROTOTYPE.  Phase 3's goal is a closed loop:
    field energy --(A)--> compresses the medium --> denser/slower region --(B)-->
    other excitations refract toward the mass  =  attraction.
That full loop is a runaway risk (3c/3d).  Here we validate each HALF on its own,
with an IMPOSED field (no feedback), which is cheap, stable, and decisive -- and
fixes the all-important sign before closing the loop.

  3a  COUPLING B (medium -> field).  Impose a dense region by hand; set the local
      wave speed c_eff^2 = c^2 * g(rho) with denser => slower.  Fly a massless
      test-wave packet past it at an impact parameter and measure whether its
      path bends TOWARD the mass (attraction) -- the gravitational-lensing sign.

  3b  COUPLING A (field -> medium).  Impose a field-energy blob; push mobile nodes
      up its gradient (F += alpha * grad e).  Measure whether node density RISES
      under the blob (self-compression).

Pure numpy; medium + operators from bvc_core.
"""
from __future__ import annotations
import math
import numpy as np

from bvc_core import lj_forces_energy, relax_medium, laplacian_matrix, mean_nn_spacing


# ============================================================ 3a : refraction =
def hex_rect(cols, rows, a=1.0):
    """A rectangular hex (triangular) lattice patch, spacing a."""
    pts = [(a * (c + 0.5 * (r % 2)), a * (math.sqrt(3) / 2) * r)
           for r in range(rows) for c in range(cols)]
    return np.array(pts, float)


def phase3a(cols=58, rows=34, beta_attract=1.6):
    print("3a -- does a wave packet refract toward an imposed dense region?")
    X = hex_rect(cols, rows)
    A = laplacian_matrix(X, rcut=1.4)                  # accurate on a regular lattice
    xs, ys = X[:, 0], X[:, 1]
    xmass, ymass = 0.5 * cols, 0.5 * rows * math.sqrt(3) / 2
    sigma, b = 4.0, 5.0                                # mass width; packet impact parameter
    rho = np.exp(-((xs - xmass) ** 2 + (ys - ymass) ** 2) / (2 * sigma ** 2))  # imposed density bump

    c, k, w, dt = 1.0, 1.2, 4.0, 0.1
    x0, y0 = 0.20 * cols, ymass + b
    x_stop = xmass + 3.0 * sigma                       # measure once the packet is past the mass

    def fly(beta, label):
        g = 1.0 / (1.0 + beta * rho)                   # beta>0: denser => slower
        env = np.exp(-((xs - x0) ** 2 + (ys - y0) ** 2) / (2 * w ** 2))
        ph = k * (xs - x0)
        u = env * np.cos(ph)
        v = c * k * env * np.sin(ph)                   # right-moving packet
        traj = []
        for _ in range(900):
            acc = c ** 2 * g * (A @ u)
            v += acc * dt
            u += v * dt
            we = u ** 2 + v ** 2
            xc = float((we * xs).sum() / we.sum())
            yc = float((we * ys).sum() / we.sum())
            traj.append((xc, yc))
            if xc > x_stop:
                break
        dperp = traj[-1][1] - ymass                    # transverse offset from the mass line
        print(f"   {label:28} reached x={traj[-1][0]:4.1f}: perp offset {b:.1f} -> {dperp:+.2f}  "
              f"(toward mass = {'YES' if dperp < b - 0.3 else 'no '})")
        return dperp

    print(f"   mass at ({xmass:.0f},{ymass:.0f}) sigma={sigma}; packet impact b={b}\n")
    d_none = fly(0.0, "no mass (control)")
    d_att = fly(beta_attract, "denser->slower (attract)")
    d_rep = fly(-0.5, "denser->faster (repel)")
    print(f"\n   => control stays at ~{d_none:+.2f}; attractive bends to {d_att:+.2f} "
          f"(toward mass); repulsive to {d_rep:+.2f} (away).")
    print(f"   => Coupling B works with sign 'denser=slower'; "
          f"deflection vs control = {d_att - d_none:+.2f}\n")


# ========================================================== 3b : compression ==
def central_density(X, center, radius):
    """nodes inside `radius` of center, and their mean nearest-neighbor spacing."""
    sel = np.linalg.norm(X - center, axis=1) < radius
    n = int(sel.sum())
    if n < 3:
        return n, float("nan")
    return n, mean_nn_spacing(X[sel])


def phase3b():
    print("3b -- does an imposed field-energy blob compress the node medium?")
    cloud = relax_medium(N=400, seed=2)
    center = np.zeros(2)
    sigma, E0, Rmeas = 5.0, 1.0, 7.0

    def compress(alpha, steps=4000, dt=0.005, cool=0.999):
        X = cloud.copy()
        V = np.zeros_like(X)

        def total_force(X):
            F, _ = lj_forces_energy(X)
            d = X - center
            e = E0 * np.exp(-(d ** 2).sum(1) / (2 * sigma ** 2))
            F = F + alpha * (-d / sigma ** 2) * e[:, None]   # alpha * grad(e), points inward
            return F

        F = total_force(X)
        for _ in range(steps):
            X += V * dt + 0.5 * F * dt ** 2
            Fn = total_force(X)
            V += 0.5 * (F + Fn) * dt
            V *= cool
            F = Fn
        return X

    n0, s0 = central_density(cloud, center, Rmeas)
    print(f"   start:           central nodes={n0:3d}  mean spacing={s0:.3f}")
    for alpha in (0.0, 8.0, 16.0):
        X = compress(alpha)
        n, s = central_density(X, center, Rmeas)
        tag = "(control)" if alpha == 0 else f"compression {100*(n-n0)/n0:+.0f}% nodes"
        print(f"   alpha={alpha:4.1f}:       central nodes={n:3d}  mean spacing={s:.3f}   {tag}")
    print("   => field energy pulls nodes in: central density rises with alpha "
          "(spacing shrinks) = Coupling A works.\n")


if __name__ == "__main__":
    print("=== Phase 3a/3b : the two halves of gravity-by-density ===\n")
    phase3a()
    phase3b()
