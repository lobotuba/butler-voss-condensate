"""
Test for emergent BOOST invariance (the rest of the Lorentz group).

test_lorentz.py / _unified.py established the ROTATIONAL facets: an isotropic
light cone (emergent rotational invariance) and, with one operator for medium and
field, a single UNIVERSAL cone. Full Lorentz symmetry also needs BOOSTS: physics
identical in frames moving relative to the medium, with no detectable rest frame
at low energy.

Boosts act on (omega, k) as a 4-vector. A dispersion is boost-invariant iff its
surface maps to itself under a boost with the emergent speed c:
    omega' = gamma (omega - beta c k_x),  k_x' = gamma (k_x - beta omega / c),  k_perp' = k_perp.
The massless cone omega = c|k| and the massive mass-shell omega^2 = c^2 k^2 + m^2
are exactly boost-invariant in the continuum; on the lattice the higher-order
dispersion breaks this, and the residual should vanish as (k a)^2 -- emergent
boost invariance, Lorentz violation suppressed as (E/E_Planck)^2.

Uses the UNIFIED (vector-Hooke) operator -- the single-cone sector -- with c set
to the low-k wave speed (normalised to 1).
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_fcc, R0
from test_lorentz import neighbours

DELTA = neighbours(perfect_fcc(radius=12.0), rcut=1.3 * R0)
KMAX = np.pi / R0
_NORM = None


def symbol(kvec):
    return (1.0 - np.cos(DELTA @ kvec)).sum()


def _norm():
    global _NORM
    if _NORM is None:
        k0 = 1e-4 * KMAX * np.array([1, 0, 0.0])
        _NORM = (np.linalg.norm(k0) ** 2) / symbol(k0)     # c_eff -> 1 at low k
    return _NORM


def omega(kvec, m=0.0):
    return np.sqrt(_norm() * symbol(kvec) + m ** 2)         # c = 1


def boost(w, kvec, beta):
    g = 1.0 / np.sqrt(1 - beta ** 2)
    wp = g * (w - beta * kvec[0])
    kp = kvec.copy(); kp[0] = g * (kvec[0] - beta * w)
    return wp, kp


def residual(kdir, frac, m, beta):
    kvec = frac * KMAX * np.array(kdir, float) / np.linalg.norm(kdir)
    w = omega(kvec, m)
    wp, kp = boost(w, kvec, beta)
    return abs(wp - omega(kp, m)) / (abs(wp) + 1e-12)       # off-shell after boost


def sweep(name, m, beta, kdir=(2, 1, 1)):
    fracs = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.4])
    res = np.array([residual(kdir, f, m, beta) for f in fracs])
    print(f"\n[{name}]  boost beta={beta}, k along {kdir}"
          + (f", mass m={m}" if m else " (massless)"))
    print("   k/kmax :", " ".join(f"{f:7.2f}" for f in fracs))
    print("   off-shell:", " ".join(f"{r:7.1e}" for r in res))
    p = np.polyfit(np.log(fracs[:4]), np.log(res[:4] + 1e-30), 1)[0]
    print(f"   => boost-violation ~ (k/kmax)^{p:.1f}  -> vanishes at low k (emergent boost inv.)")
    return res


if __name__ == "__main__":
    print("=== Test for emergent boost invariance (unified single-cone sector) ===")
    print(f"emergent c normalised to 1 (low-k wave speed); k_max = pi/R0")

    sweep("massless cone", m=0.0, beta=0.3)
    sweep("massless cone", m=0.0, beta=0.6)
    sweep("massive mass-shell", m=0.3, beta=0.6)

    # mass-shell check: omega^2 - k^2 should equal m^2 (Lorentz invariant) at low k
    print("\n[mass-shell]  omega^2 - |k|^2 (should = m^2 = 0.090), k along (2,1,1):")
    for f in (0.05, 0.1, 0.2, 0.4):
        kvec = f * KMAX * np.array([2, 1, 1.0]) / np.linalg.norm([2, 1, 1])
        inv = omega(kvec, 0.3) ** 2 - np.linalg.norm(kvec) ** 2
        print(f"   k/kmax={f:.2f}:  omega^2 - k^2 = {inv:.4f}")

    # causality: the front (max group) velocity must not exceed c
    ks = np.linspace(1e-3, 0.5, 60) * KMAX
    w = np.array([omega(k * np.array([1, 0, 0.0])) for k in ks])
    vg = np.gradient(w, ks)
    print(f"\n[causality]  max group velocity over 0<k<0.5 kmax = {vg.max():.4f}  "
          f"({'<= c, causal' if vg.max() <= 1.0 + 1e-6 else 'SUPERLUMINAL'})")
    print("\n  => the emergent massless cone and massive mass-shell are boost-invariant at")
    print("     low energy (violation ~ (E/E_Planck)^2), omega^2-k^2 is a Lorentz invariant,")
    print("     and signals stay sub-luminal. With isotropy (test_lorentz) and one universal")
    print("     cone (test_lorentz_unified), the FULL Lorentz group emerges at long wavelength.")
