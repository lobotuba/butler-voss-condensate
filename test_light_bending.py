"""
Does the model bend light? The PPN gamma, and WHY the wall is where it is.

*** STATUS UPDATE (conclusion superseded by the tensor-gravity arc). The measurement and its
    derivation below are correct AND still load-bearing: mass deforms the medium by a DISPLACEMENT,
    displacement-derived strain is COMPATIBLE, compatible strain is flat (pure gauge), and it bends
    no light. That is why the ELASTIC/scalar sector has gamma = 0.
    What is superseded is the headline "the model gives gamma = 0". It does not, at long range. The
    curvature (incompatible) sector -- which this file correctly identified as the only source of
    genuine bending, but dismissed because it REPELLED -- was later shown to be CONFINING rather
    than repulsive in the usual sense, and to DECONFINE into an attractive massless graviton once
    the fermion loop supplies a positive induced Einstein term (test_deconfinement, with mu > 0
    measured in test_induced_sign). That graviton is dynamical and healthy in 3+1D
    (test_spin2_dynamical), and gamma -> 1 follows from Weinberg on the conserved IR stress tensor
    (test_lattice_ward). So: gamma = 0 for the compatible/scalar sector (this file, still true),
    gamma -> 1 for the model's actual long-range gravity. ***

The spin-2 wall's one observable crux is light bending, set by gamma = Psi/Phi (space/time
curvature): gamma = 0 is scalar gravity (no Eddington bending), gamma = 1 is GR (factor of two).
A first attempt built Psi from a stand-in and got a meaningless number; a second, using the
medium's real strain, found the shear channel returned exactly ZERO and I did not yet see why.

Here is why, and it is the whole result. A mass can only make the medium respond ONE way: its
nodes move, i.e. a displacement field u(x). The spatial metric a ruler (or a light ray) then
sees is g_ij = delta_ij + 2 e_ij with e_ij = 1/2 (d_i u_j + d_j u_i). But a strain derived from
a single-valued displacement is COMPATIBLE, and a compatible strain is a FLAT metric -- it is
literally the coordinate relabeling x -> x + u, a diffeomorphism, pure gauge. Flat space has no
curvature and bends no light. That is why the shear integrated to zero: it was gauge all along.

In 2D the curvature of a strain field is its incompatibility (the St-Venant / Gaussian-curvature
operator)
        eta = d2(e_xx)/dy2 + d2(e_yy)/dx2 - 2 d2(e_xy)/dxdy,
which is IDENTICALLY zero for any e from a displacement. Genuine curvature (eta != 0) requires a
strain that CANNOT be written as a displacement -- a topological DISCLINATION. And that is
exactly the sector test_disclination_force measured: it carries real curvature, but its like
charges REPEL with a force that grows with distance. So the medium's only source of genuine
spatial curvature has the wrong sign for gravity, while its response to actual mass is pure
gauge. That is the spin-2 wall, derived from one principle rather than asserted.

This file demonstrates all three legs:
  A. A mass -> displacement -> strain. Its incompatibility eta ~ 0 (FLAT, gauge), and the light
     deflection from the spatial channel is ~ 0. The spatial sector cannot bend light.
  B. A disclination -> incompatible strain. eta != 0: genuine curvature exists, but ONLY for the
     topological defect (which repels -- test_disclination_force).
  C. The amplitude-mode potential Phi is a genuine scalar field (not a displacement gradient), so
     it DOES bend light -- but through g_00 only, giving gamma = 0. Scalar gravity, no factor of 2.
"""
from __future__ import annotations
import numpy as np


def d(f, ax, k):
    fk = np.fft.fft2(f)
    return np.fft.ifft2((1j * k[ax]) ** 1 * fk).real if False else None


def spectral_ops(N):
    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    return KX, KY


def dd(f, KA, KB):
    """second derivative d^2 f / dA dB, spectral."""
    return np.fft.ifft2(-(KA * KB) * np.fft.fft2(f)).real


def elastic_strain(rho, KX, KY, lam=1.0, mu=1.0):
    """Displacement u sourced by a mass (centre of dilatation, eigenstrain ~ rho), then the
    strain e_ij = 1/2(d_i u_j + d_j u_i). Same solver as test_shielding."""
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0
    sk = np.fft.fft2(rho)
    pref = -2 * (lam + mu) / (lam + 2 * mu)
    ux = np.fft.ifft2(pref * 1j * KX * sk / K2).real
    uy = np.fft.ifft2(pref * 1j * KY * sk / K2).real
    exx = np.fft.ifft2(1j * KX * np.fft.fft2(ux)).real
    eyy = np.fft.ifft2(1j * KY * np.fft.fft2(uy)).real
    exy = 0.5 * (np.fft.ifft2(1j * KY * np.fft.fft2(ux)).real +
                 np.fft.ifft2(1j * KX * np.fft.fft2(uy)).real)
    return ux, uy, exx, eyy, exy


def incompat(exx, eyy, exy, KX, KY):
    """eta = e_xx,yy + e_yy,xx - 2 e_xy,xy  -- the 2D curvature (strain incompatibility)."""
    return dd(exx, KY, KY) + dd(eyy, KX, KX) - 2 * dd(exy, KX, KY)


def screened_poisson(rho, m2, tol=1e-11, itmax=20000):
    def lap(f):
        return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
                np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)
    x = np.zeros_like(rho); r = rho + lap(x) - m2 * x; p = r.copy(); rs = (r * r).sum()
    for _ in range(itmax):
        Ap = -lap(p) + m2 * p
        al = rs / ((p * Ap).sum() + 1e-300); x += al * p; r -= al * Ap
        rs2 = (r * r).sum()
        if np.sqrt(rs2) < tol:
            break
        p = r + (rs2 / rs) * p; rs = rs2
    return x


def deflect(seen, N, bs):
    """alpha(b) = -integral d(seen)/dy dx along +x rays at offset b (leading-order geodesic)."""
    c = N // 2
    gy = 0.5 * (np.roll(seen, -1, 1) - np.roll(seen, 1, 1))
    return np.array([-gy[:, c + b].sum() for b in bs])


if __name__ == "__main__":
    print("=== Does the model bend light? gamma, and WHY the wall is where it is ===\n")

    N = 256
    KX, KY = spectral_ops(N)
    g = np.arange(N) - N // 2
    X, Y = np.meshgrid(g, g, indexing="ij")
    R = np.hypot(X, Y)
    rho = np.exp(-R ** 2 / (2 * 3.0 ** 2)); rho -= rho.mean()

    # ---- A. the mass's spatial response is a displacement => FLAT (pure gauge) ----
    ux, uy, exx, eyy, exy = elastic_strain(rho, KX, KY)
    eta_mass = incompat(exx, eyy, exy, KX, KY)
    strain_norm = np.abs(exx).max()
    print("  [A] a MASS: the medium responds with a displacement u, so its strain is COMPATIBLE.")
    print(f"      curvature (incompatibility) eta:  max|eta|/|strain| = "
          f"{np.abs(eta_mass).max()/strain_norm:.2e}  -> FLAT (pure gauge, a diffeomorphism)")
    bs = np.array([12, 16, 20, 26, 32, 40])
    a_space = deflect(exx - 0.5 * (exx + eyy), N, bs)      # what a +x ray sees from the shear
    print(f"      light deflection from the SPATIAL (strain) channel: "
          f"max|alpha| = {np.abs(a_space).max():.2e}  -> ZERO. Gauge bends no light.\n")

    # ---- B. genuine curvature exists only for a DISCLINATION (incompatible strain) ----
    # a wedge disclination: impose an incompatible isotropic strain theta*delta_ij with
    # eta = lap(theta) != 0 (theta from a single-valued field, but entering as the PHYSICAL
    # strain, not as 1/2(du+du^T) -- that is precisely what a topological defect does).
    theta = np.exp(-R ** 2 / (2 * 4.0 ** 2))
    eta_disc = incompat(theta, theta, np.zeros_like(theta), KX, KY)   # = lap(theta)
    print("  [B] a DISCLINATION: incompatible strain, NOT derivable from any displacement.")
    print(f"      curvature eta:  max|eta|/|strain| = {np.abs(eta_disc).max()/np.abs(theta).max():.2e}"
          "  -> NONZERO. Genuine curvature -- but this is the topological sector, whose like")
    print("      charges REPEL with a growing force (test_disclination_force). Wrong sign for gravity.\n")

    # ---- C. what DOES bend light: the amplitude-mode scalar Phi (but gamma = 0) ----
    Phi = screened_poisson(4 * np.pi * rho, 0.05 ** 2)
    a_time = deflect(Phi, N, bs)
    print("  [C] the amplitude-mode potential Phi is a genuine SCALAR field (not a displacement).")
    print(f"      it bends light through g_00: max|alpha_time| = {np.abs(a_time).max():.3e}  -> NONZERO,")
    print("      but it is the TIME part only. gamma = (spatial bending)/(time bending) =")
    gamma = np.abs(a_space).max() / np.abs(a_time).max()
    print(f"          gamma = {gamma:.2e}  ~ 0   (scalar gravity; GR would be 1, the factor of two)\n")

    print("[verdict] gamma = 0, and now it is DERIVED, not just measured:")
    print("  * A mass can only move the medium's nodes -- a displacement field u. Any strain from a")
    print("    displacement is COMPATIBLE, i.e. a FLAT metric, i.e. the coordinate change x -> x+u:")
    print("    pure gauge. Flat space carries no curvature and bends no light. That is why the")
    print("    long-range tetrad shear -- real as a field -- does ZERO gravitational work (and why")
    print("    test_tetrad_force found no force, and this finds no bending): it is a diffeomorphism.")
    print("  * Genuine spatial curvature exists ONLY as an INCOMPATIBLE strain -- a topological")
    print("    disclination -- and that sector REPELS and confines (test_disclination_force). The")
    print("    medium's one source of real curvature has the wrong sign for gravity.")
    print("  * So the only thing that attracts at range is the SCALAR amplitude mode, which bends")
    print("    light through g_00 alone: gamma = 0, no Eddington factor of two.")
    print("\n  => The spin-2 wall, derived from one principle: a fixed-background medium responds to")
    print("     mass by a DISPLACEMENT (compatible strain = flat = gauge), so it cannot generate the")
    print("     genuine spatial curvature GR requires. GR needs a metric perturbation that is NOT a")
    print("     displacement gradient -- a dynamical, diffeomorphism-invariant field -- which is")
    print("     exactly what the lattice lacks. This is not a numerical limit to push past; it is")
    print("     structural. The route through it is to make the metric dynamical in its own right")
    print("     (incompatible strain as a propagating field), not to refine the medium.")
