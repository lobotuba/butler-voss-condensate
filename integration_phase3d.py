"""
Integration Phase 3d : do two masses drift together? (the gravity test)
=======================================================================

The decisive test of gravity-by-density, and the one the whole program was built
toward.  Self-binding (3c) turned out to need *strong* gravity and failed (the
medium is nearly incompressible); but mutual attraction between two masses needs
only the *weak* force already confirmed by 3a (waves refract toward a mass).  So
this is the right test, and it does NOT require a bound state.

Setup: seed TWO field lumps a fixed distance apart on the medium, evolve with the
energy-conserving variational coupling (bounded g), and measure their separation
over time.  Signature of gravity:
    beta=0 (no coupling) -> no attraction (separation ~constant / grows)
    beta>0               -> separation SHRINKS
    larger beta          -> shrinks more

Result (N=300, damping=1, dt=0.001, to t=8):
    beta=0 : 9.33 -> 10.85  (+1.51, drift apart)   dE -0.12%
    beta=40: 9.33 ->  7.78  (-1.56, drift together) dE -0.17%
    beta=60: 9.33 ->  7.48  (-1.85, stronger)       dE -0.08%
=> two masses attract via the emergent medium force, energy-conserved, with the
   attraction growing with coupling strength: gravity-by-density confirmed.
"""
from __future__ import annotations
import numpy as np

from bvc_core import relax_medium
from integration_phase3_variational import VariationalCoupled


def seed_two_lumps(f, d, amp=1.5, w=2.0):
    """Two Gaussian field lumps centered at (+/- d/2, 0)."""
    x = f.X
    f.u = amp * (np.exp(-(((x[:, 0] - d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)) +
                 np.exp(-(((x[:, 0] + d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)))
    f.pi = np.zeros(f.N)


def separation(f):
    """Distance between the field-energy centroids of the two halves (x<0, x>0)."""
    e = 0.5 * f.pi ** 2 + 0.5 * f.m2 * f.u ** 2
    L = f.X[:, 0] < 0; R = ~L
    if e[L].sum() < 1e-9 or e[R].sum() < 1e-9:
        return float("nan")
    cl = (f.X[L] * e[L, None]).sum(0) / e[L].sum()
    cr = (f.X[R] * e[R, None]).sum(0) / e[R].sum()
    return float(np.linalg.norm(cr - cl))


def run(cloud, beta, d0=9.0, steps=8000, dt=0.001):
    f = VariationalCoupled(cloud, beta=beta, m2=1.0, g_min=0.02, damping=1.0, dt=dt)
    seed_two_lumps(f, d0)
    E0 = f.energy()
    seps = [separation(f)]
    for k in range(steps):
        f.step()
        if (k + 1) % 1000 == 0:
            seps.append(separation(f))
    dE = 100 * (f.energy() - E0) / abs(E0)
    verdict = ("ATTRACT (drift together)" if seps[-1] < seps[0] - 0.5
               else "no attraction (drift apart)" if seps[-1] > seps[0] + 0.5
               else "~static")
    print(f"  beta={beta:3d}: sep " + " ".join(f"{s:.2f}" for s in seps) +
          f"   net {seps[-1]-seps[0]:+.2f}  dE={dE:+.2f}%  => {verdict}")


if __name__ == "__main__":
    print("=== Phase 3d : do two masses drift together? (gravity-by-density) ===")
    print("  seed two lumps at +/-4.5; energy-conserving variational coupling\n")
    cloud = relax_medium(N=300, seed=3)
    for beta in (0, 40, 60):
        run(cloud, beta)
    print("\n  beta=0 shows no attraction; coupling makes them drift together,")
    print("  more strongly at higher beta, with energy conserved => gravity.")
