"""
The sharpest barrier: long-range spin-2 gravity. What it needs vs what we have.

The forces program proved the model's gravity-by-density SCREENS (Yukawa, lambda~3):
it couples to energy density, which the medium relaxes (Bitter-Crum). Real gravity
is (i) long-range 1/r^2, (ii) universally ATTRACTIVE (couples to mass/energy, one
sign), and (iii) SPIN-2 (a massless symmetric-tensor graviton -- needed for light
bending and gravitational waves, not just Newton). This test checks each requirement
against the model.

  A. LONG-RANGE + UNIVERSAL: a MASSLESS mediator coupled to mass density gives a
     1/r^2 force, and because mass density is intrinsically POSITIVE it is
     universally attractive (unlike EM's +/- charges). The machinery can carry this
     Newtonian (scalar) piece -- IF the mediator is massless.
  B. SPIN-2: the graviton needs helicity +/-2. The medium's displacement field is a
     VECTOR: its phonons are 1 longitudinal (helicity 0) + 2 transverse (helicity
     +/-1) -- spin-0 and spin-1, like a photon. NO spin-2 mode exists among them.
  C. Verdict: gravity is doubly out of reach -- the mass-coupling is MASSIVE
     (screened, so not even scalar-Newtonian survives long-range), and there is no
     spin-2 degree of freedom (and Weinberg-Witten forbids an emergent graviton
     coupling to the full T_mu_nu). Precisely the sharpest open problem.
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_fcc, R0
from test_lorentz import neighbours, lj_derivs

# ------------------------------------------- A. long-range, universal attraction
L = 48
_k1 = 2 * np.pi * np.fft.fftfreq(L)
KX, KY, KZ = np.meshgrid(_k1, _k1, _k1, indexing="ij")
SYM = 2 * ((1 - np.cos(KX)) + (1 - np.cos(KY)) + (1 - np.cos(KZ)))     # ~ k^2 (massless)


def blob(c):
    x = np.arange(L)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r2 = (X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2
    return np.exp(-r2 / 2.0)


def potential(rho, m2=0.0):
    rk = np.fft.fftn(rho); Gk = np.zeros_like(rk); msk = (SYM + m2) > 1e-9
    Gk[msk] = rk[msk] / (SYM[msk] + m2)
    return np.fft.ifftn(Gk).real


def part_A():
    print("[A] gravity mediator: massless -> long-range 1/r; massive -> screened (the model's case)")
    c = L // 2
    Phi0 = potential(blob((c, c, c)))                 # massless mediator
    Phim = potential(blob((c, c, c)), m2=0.25)        # massive mediator (m=0.5), = gravity-by-density
    r = np.arange(2, 13); g0 = Phi0[c + r, c, c]; gm = Phim[c + r, c, c]
    fit = (r >= 2) & (r <= 8)
    p0 = -np.polyfit(np.log(r[fit]), np.log(g0[fit]), 1)[0]
    xim = -1.0 / np.polyfit(r[fit], np.log(np.abs(gm[fit])), 1)[0]
    print(f"   massless potential Phi(r) ~ 1/r^{p0:.2f}   (Newtonian 1/r; long-range, unscreened)")
    print("     r   :", " ".join(f"{x:6d}" for x in r[:8]))
    print("     Phi :", " ".join(f"{x:6.3f}" for x in g0[:8]))
    print(f"   massive potential Phi(r) ~ e^(-r/xi)/r, xi={xim:.1f}  (SCREENED -- the model's own")
    print("     gravity-by-density mediator is massive by Bitter-Crum: this is why gravity is short-range)")
    print("   => a MASSLESS mass-mediator would give Newtonian 1/r, and since mass density is")
    print("      intrinsically POSITIVE the force is UNIVERSALLY ATTRACTIVE (unlike EM's +/- charges).")
    print("      But the model's mediator is the massive one -> screened. Masslessness is the missing piece.\n")


# ------------------------------------------------------- B. spin content (no spin-2)
def part_B():
    print("[B] the medium's massless modes: what spin? (graviton needs helicity +/-2)")
    fcc = perfect_fcc(radius=10.0); delta = neighbours(fcc, rcut=2.5)
    khat = np.array([2, 1, 1.0]); khat /= np.linalg.norm(khat)
    kvec = 0.1 * (np.pi / R0) * khat
    D = np.zeros((3, 3))
    for dj in delta:
        r = np.linalg.norm(dj); u = dj / r; p1, p2 = lj_derivs(r)
        D += (p2 * np.outer(u, u) + (p1 / r) * (np.eye(3) - np.outer(u, u))) * (1 - np.cos(kvec @ dj))
    w2, evec = np.linalg.eigh(D)
    print(f"   {'branch':>7} {'omega':>8} {'|e.k^|':>8} {'character':>14}")
    for i in range(3):
        proj = abs(evec[:, i] @ khat)
        kind = "longitudinal" if proj > 0.9 else "transverse"
        print(f"   {i:>7} {np.sqrt(max(w2[i],0))/np.linalg.norm(kvec):>8.3f} {proj:>8.3f} {kind:>14}")
    print("   => 1 longitudinal (helicity 0 = spin-0) + 2 transverse (helicity +/-1 = spin-1).")
    print("      A vector displacement field cannot carry helicity +/-2: there is NO spin-2")
    print("      (graviton) mode among the phonons. The medium offers at most a photon-like")
    print("      massless sector, not a graviton.\n")


if __name__ == "__main__":
    print("=== Long-range spin-2 gravity: what it needs vs what the model has ===\n")
    part_A()
    part_B()
    print("[C] verdict -- gravity is the sharpest, still-unsolved barrier:")
    print("  * long-range + universal attraction: ACHIEVABLE in principle (massless mediator,")
    print("    positive mass density) -- but the model's mass-coupling (gravity-by-density) is")
    print("    MASSIVE/screened (Bitter-Crum, lambda~3), so even scalar-Newtonian gravity does")
    print("    NOT survive at long range. That screening is the core contradiction.")
    print("  * spin-2: MISSING -- the medium's phonons are spin-0 + spin-1 only; a graviton")
    print("    needs a fundamental symmetric-tensor d.o.f. the displacement medium lacks, and")
    print("    Weinberg-Witten forbids an emergent massless spin-2 coupling to the full T_mu_nu.")
    print("  => Long-range spin-2 gravity is NOT achieved. To get it the model must add a")
    print("     tensor d.o.f. whose mass-coupling is protected massless (gauged diffeomorphisms)")
    print("     -- the analog of how gauging gave EM its 1/r^2 (Route C), but for spin-2. This")
    print("     is the frontier the whole forces + fundamental-physics program points to.")
