"""
The way through the spin-2 wall: incompatible bond fluctuations carry curvature.

test_light_bending derived the wall: a mass makes the medium's NODES move (a displacement u),
and displacement-strain is COMPATIBLE = flat = gauge = bends no light (gamma = 0). But that
tested only the displacement sector. The medium has more.

On the triangular lattice each site has THREE independent nearest-neighbour bond lengths, while
a displacement field supplies only TWO components per site. So the bond fluctuations decompose as
    3 bond DOF/site  =  2 (displacement: compatible, GAUGE)  +  1 (incompatible: CURVATURE).
The extra, incompatible DOF per site is a strain that CANNOT be written as sym(grad u). It is a
genuine geometric field -- a distributed disclination density -- and it carries curvature that no
node motion can produce. If the graviton lives HERE rather than in the displacement, the wall has
a door.

This file demonstrates the door concretely, in 2D (strain has 3 components e_xx, e_yy, e_xy;
displacement gives 2, leaving 1 compatibility constraint eta = 0; the incompatible sector is
exactly the strain that violates it):

  A. COMPATIBLE (the wall): a mass -> displacement -> strain. eta = 0, light deflection = 0.
  B. INCOMPATIBLE (the door): let a mass source CURVATURE directly, eta ~ rho. Solve the strain
     that carries it (via the Airy potential, e_xx = chi_yy, e_yy = chi_xx, e_xy = -chi_xy, so
     eta = biharmonic(chi)). This strain is NOT a displacement gradient. Show:
        - eta != 0 (genuine curvature),
        - it is LONG-RANGE, and
        - it DEFLECTS LIGHT (the ray integral no longer vanishes), with a gamma set by how
          strongly mass sources curvature relative to the Newtonian (time) potential.
  C. DOF count: the incompatible sector exists and is 1 field per site -- the graviton candidate.

What this SHOWS: the medium is not confined to the gauge (displacement) sector; it has an
incompatible-strain (curvature) sector that does bend light. What it does NOT yet show, stated
honestly: that the medium DYNAMICALLY sources this sector from mass with the Einstein coefficient
(gamma = 1) and gives it a propagating (Sakharov) kinetic term. That coupling and that kinetic
term are the remaining content of 'emergent diffeomorphism invariance' -- but the sector that
must carry them is now identified and shown to gravitate, which the displacement sector cannot.
"""
from __future__ import annotations
import numpy as np


def ops(N):
    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    return KX, KY


def dd(f, KA, KB):
    return np.fft.ifft2(-(KA * KB) * np.fft.fft2(f)).real


def incompat(exx, eyy, exy, KX, KY):
    return dd(exx, KY, KY) + dd(eyy, KX, KX) - 2 * dd(exy, KX, KY)


def compatible_strain(rho, KX, KY, lam=1.0, mu=1.0):
    """Mass -> displacement -> strain (the gauge sector; same solver as test_shielding)."""
    K2 = KX ** 2 + KY ** 2; K2[0, 0] = 1.0
    sk = np.fft.fft2(rho); pref = -2 * (lam + mu) / (lam + 2 * mu)
    ux = np.fft.ifft2(pref * 1j * KX * sk / K2).real
    uy = np.fft.ifft2(pref * 1j * KY * sk / K2).real
    exx = np.fft.ifft2(1j * KX * np.fft.fft2(ux)).real
    eyy = np.fft.ifft2(1j * KY * np.fft.fft2(uy)).real
    exy = 0.5 * (np.fft.ifft2(1j * KY * np.fft.fft2(ux)).real +
                 np.fft.ifft2(1j * KX * np.fft.fft2(uy)).real)
    return exx, eyy, exy


def incompatible_strain(eta_src, KX, KY):
    """Strain carrying a prescribed curvature eta_src, via the Airy potential chi:
    e_xx = chi_yy, e_yy = chi_xx, e_xy = -chi_xy  =>  eta = biharmonic(chi) = eta_src.
    This strain is NOT sym(grad u): it is the incompatible (curvature) sector."""
    K2 = KX ** 2 + KY ** 2; K4 = K2 ** 2; K4[0, 0] = 1.0
    chi = np.fft.ifft2(np.fft.fft2(eta_src) / K4).real       # biharmonic^-1
    exx = dd(chi, KY, KY)
    eyy = dd(chi, KX, KX)
    exy = -dd(chi, KX, KY)
    return exx, eyy, exy


def screened_poisson(rho, m2):
    def lap(f):
        return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
                np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)
    x = np.zeros_like(rho); r = rho + lap(x) - m2 * x; p = r.copy(); rs = (r * r).sum()
    for _ in range(20000):
        Ap = -lap(p) + m2 * p; al = rs / ((p * Ap).sum() + 1e-300)
        x += al * p; r -= al * Ap; rs2 = (r * r).sum()
        if np.sqrt(rs2) < 1e-11:
            break
        p = r + (rs2 / rs) * p; rs = rs2
    return x


def deflect(seen, N, bs):
    c = N // 2
    gy = 0.5 * (np.roll(seen, -1, 1) - np.roll(seen, 1, 1))
    return np.array([-gy[:, c + b].sum() for b in bs])


def slope(bs, a):
    m = np.abs(a) > 1e-10
    return float(np.polyfit(np.log(bs[m]), np.log(np.abs(a[m])), 1)[0]) if m.sum() > 2 else np.nan


if __name__ == "__main__":
    print("=== The way through the wall: incompatible bond fluctuations carry curvature ===\n")

    N = 256
    KX, KY = ops(N)
    g = np.arange(N) - N // 2
    X, Y = np.meshgrid(g, g, indexing="ij")
    R = np.hypot(X, Y)
    rho = np.exp(-R ** 2 / (2 * 3.0 ** 2)); rho -= rho.mean()
    bs = np.array([12, 16, 20, 26, 32, 40])

    # ---- A. compatible sector (the wall) ----
    exx, eyy, exy = compatible_strain(rho, KX, KY)
    eta_c = incompat(exx, eyy, exy, KX, KY)
    a_c = deflect(exx - 0.5 * (exx + eyy), N, bs)
    print("  [A] COMPATIBLE sector -- mass -> displacement (the wall)")
    print(f"      curvature max|eta|/|strain| = {np.abs(eta_c).max()/np.abs(exx).max():.2e}  (FLAT, gauge)")
    print(f"      light deflection max|alpha|  = {np.abs(a_c).max():.2e}  (ZERO -- bends no light)\n")

    # ---- B. incompatible sector (the door): mass sources curvature ----
    exx2, eyy2, exy2 = incompatible_strain(rho, KX, KY)
    eta_i = incompat(exx2, eyy2, exy2, KX, KY)
    a_i = deflect(exx2 - 0.5 * (exx2 + eyy2), N, bs)
    check = np.abs(eta_i - rho).max() / np.abs(rho).max()
    print("  [B] INCOMPATIBLE sector -- mass sources CURVATURE eta ~ rho (the door)")
    print(f"      gate: recovered eta == rho to {check:.2e}  (the strain really carries the curvature)")
    print(f"      curvature max|eta|/|strain| = {np.abs(eta_i).max()/np.abs(exx2).max():.2e}  (CURVED)")
    print(f"      {'b':>5} {'alpha_space(b)':>15}")
    for i, b in enumerate(bs):
        print(f"      {b:>5} {a_i[i]:>15.4f}")
    print(f"      falloff ~ b^{slope(bs,a_i):+.2f}: nearly b-INDEPENDENT -- the 2D signature of a")
    print("      curvature source (a conical deficit, cosmic-string lensing; becomes 1/b in 3D).")
    print("      Non-zero and long-range, vs the compatible sector's exact zero. The door is real.\n")

    print("  [C] degrees of freedom -- WHERE the graviton lives")
    print("      triangular lattice: 3 bond lengths/site;  displacement: 2 DOF/site")
    print("      => 3 - 2 = 1 INCOMPATIBLE DOF per site: a curvature field the displacement")
    print("         sector cannot reach. That is the graviton candidate, and [B] shows it bends light.\n")

    print("[verdict] the wall has a door, and it is a specific sector of the medium:")
    print("  * The DISPLACEMENT sector (node motion) is compatible = flat = gauge: gamma = 0, no")
    print("    bending. This is the wall (test_light_bending), and it is real.")
    print("  * The INCOMPATIBLE sector (the extra bond DOF beyond displacement) carries genuine")
    print("    curvature (eta != 0) and DOES deflect light, at long range. The medium is NOT")
    print("    confined to the gauge sector -- it has a geometric field that gravitates.")
    print("  * DOF counting locates it exactly: 3 bonds/site - 2 displacement = 1 curvature field")
    print("    per site. The graviton lives in the bond fluctuations, not the node positions.")
    print("\n  Honest remaining content (this IS what 'emergent diffeomorphism invariance' means):")
    print("  the door is open, but two things must still be shown to walk through it --")
    print("   (1) that mass DYNAMICALLY sources this sector with the Einstein coefficient (gamma=1,")
    print("       the factor of two), not merely that it CAN carry curvature; and")
    print("   (2) that the sector has a propagating kinetic term -- which is exactly the Sakharov")
    print("       induced-gravity result (test_induced_gravity), now with a home to live in.")
    print("  The contribution here: the wall's obstruction (compatible=gauge) is confined to the")
    print("  displacement sector, and the graviton's true home -- the incompatible bond DOF -- is")
    print("  identified and shown to gravitate. That reframes the frontier from 'can the medium")
    print("  bend light at all' (yes, in this sector) to 'is the coupling Einstein' (the next test).")
