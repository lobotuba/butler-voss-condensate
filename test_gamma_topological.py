"""
The topological channel behind gamma = 1: does mass carry net curvature charge? Gauss-Bonnet says no.

*** STATUS UPDATE -- the dimensional caveat below is now RESOLVED, and the 3+1D answer is the same.
    The verdict here notes that Gauss-Bonnet is a 2D theorem and defers whether the 3+1D medium evades
    it the way general relativity does. test_gamma_3d does that calculation, by direct ray tracing in
    3D, and the medium does NOT reach gamma = 1. The reason turns out to be deeper than Gauss-Bonnet
    and independent of dimension: general relativity sources the spatial metric by a Poisson equation,
    lap(Psi) = 4 pi G rho, so Psi is the long-range POTENTIAL and gamma = Psi/Phi = 1; the medium's
    compression sets Psi = theta* = rho ALGEBRAICALLY, so Psi is LOCAL, and gamma = Psi/Phi falls to
    zero at range (a general-relativity control with Psi = the potential holds gamma = 1 at every
    impact parameter, confirming the ray tracer). The 2D zero-charge result was a symptom of this. So
    the caveat does not rescue gamma = 1; the 3+1D answer matches the 2D one, for a cleaner reason. ***

test_gamma_source closed the SMOOTH route to gamma = 1: a static mass is a scalar energy density, and
its induced coupling to the spin-2 spatial stress vanishes by a selection rule, so the smooth loop
sources the Newtonian potential Phi and no curvature. But that selection rule is specifically about
spin: it forbids a scalar from sourcing the spin-2 part of h_ij. Curvature itself -- the 2D Ricci
scalar, the incompatibility eta -- is a SCALAR, and test_incompatible_gravity located a genuine scalar
curvature degree of freedom in the medium (the third bond fluctuation per site, beyond the two
displacement modes). A scalar mass sourcing a scalar curvature is not forbidden by any spin rule. So
the topological/curvature channel is the one route left that could still give gamma = 1 directly, and
this file measures it.

The mechanism is the model's own: gravity here is COMPRESSION -- field energy changes the medium's
preferred local density (Phase 3). That is an EIGENSTRAIN theta*(x) ~ rho(x): a local, isotropic,
stress-free change of preferred bond length, exactly the non-uniform thermal expansion that buckles a
heated plate. Unlike a body force, which relaxes to a compatible displacement (test_einstein_source
[A], eta ~ 1e-15), an inhomogeneous eigenstrain leaves real Gaussian curvature that light, riding the
bond metric (test_incompatible_gravity), must see:

        eta = incompat(theta* delta_ij) = lap(theta*) ~ lap(rho).

This is nonzero, so mass DOES curve the medium locally. The question is whether that curvature is the
Einstein shape. gamma = 1 needs the light-bending curvature charge to be the mass itself, eta ~ rho,
so that the spatial potential Psi tracks the Newtonian Phi. Compression gives eta ~ lap(rho) instead,
and the difference is not cosmetic -- it is Gauss-Bonnet.

  [A] Mass as an eigenstrain sources curvature eta = lap(theta*), nonzero and local.
  [B] But its NET curvature charge vanishes identically: integral(eta) = integral(lap theta*) = 0, a
      total derivative, to machine precision and for every profile. A smooth localized compression
      makes a curvature DIPOLE -- a dome in the compressed core, a saddle in the surrounding ring --
      with no net deficit angle. This is the Gauss-Bonnet theorem: smooth deformations carry no
      topological curvature charge.
  [C] Zero net charge means no long-range bending. The deflection from the compression curvature dies
      with distance, while a genuine curvature charge (eta ~ rho, a disclination density proportional
      to the mass) bends light at range. gamma from compression is 0 at long range.
  [D] The complete picture across all channels, and what gamma = 1 actually requires.

This is the last of the accessible direct channels. Every one the model exposes at low dimension --
the smooth induced loop (test_gamma_source, 2+1D), the elastic body force (test_einstein_source, 2D),
and now the compression eigenstrain (2D) -- sources zero net curvature charge, each for its own clean
reason: a spin selection rule, strain compatibility, and Gauss-Bonnet. gamma = 1 requires mass to
carry net curvature charge, i.e. to nucleate a topological DISCLINATION density proportional to
itself, which is exactly the Einstein source coupling test_einstein_source isolated and which no
mechanism in the model supplies. The caveat that keeps this from being a proof is dimensional and is
stated in full in the verdict: Gauss-Bonnet is a 2D theorem, and in 3+1D a smooth mass sources smooth
curvature without any topological defect, so whether the 3+1D medium evades the obstruction the way
general relativity does is a further calculation. What that leaves for gamma = 1 is the same status as the graviton Ward identity in
test_lattice_ward: an INFRARED-EMERGENT identity, argued from Weinberg's theorem on the conserved IR
stress tensor, with no direct in-model realization -- and now, unlike the model's directly confirmed
emergent Lorentz invariance, a series of direct results against it at the accessible scale.
"""
from __future__ import annotations
import numpy as np


def ops(N):
    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    return KX, KY


def dd(f, KA, KB):
    return np.fft.ifft2(-(KA * KB) * np.fft.fft2(f)).real


def rlap(f):
    """Real-space Laplacian, so the total-derivative cancellation in [B] is exact, not spectral."""
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
            np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)


def poisson(rho, KX, KY):
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0
    out = np.fft.ifft2(np.fft.fft2(rho) / K2).real
    return out - out.mean()


def incompatible_strain(eta, KX, KY):
    """Strain carrying curvature eta via the Airy potential: biharm(chi) = eta,
    e_xx = chi_yy, e_yy = chi_xx, e_xy = -chi_xy. The strain that is NOT a displacement gradient."""
    K2 = KX ** 2 + KY ** 2
    K4 = K2 ** 2
    K4[0, 0] = 1.0
    chi = np.fft.ifft2(np.fft.fft2(eta) / K4).real
    return dd(chi, KY, KY), dd(chi, KX, KX), -dd(chi, KX, KY)


def deflect(seen, N, bs):
    """Deflection of +x rays at impact parameters bs, from the spatial field they integrate."""
    c = N // 2
    gy = 0.5 * (np.roll(seen, -1, 1) - np.roll(seen, 1, 1))
    return np.array([-gy[:, c + b].sum() for b in bs])


if __name__ == "__main__":
    print("=== The topological channel behind gamma = 1: does mass carry net curvature charge? ===\n")
    print("  The smooth route died on a spin rule (test_gamma_source): a scalar mass can't source the")
    print("  spin-2 spatial stress. But curvature is itself a SCALAR, so a scalar mass sourcing scalar")
    print("  curvature is allowed. The model's gravity is COMPRESSION -- mass as an eigenstrain")
    print("  theta*(x) ~ rho -- and an inhomogeneous eigenstrain leaves curvature eta = lap(theta*).\n")

    N = 384
    KX, KY = ops(N)
    g = np.arange(N) - N // 2
    X, Y = np.meshgrid(g, g, indexing="ij")
    R = np.hypot(X, Y)
    rho = np.exp(-R ** 2 / (2 * 3.0 ** 2))
    bs = np.array([16, 24, 32, 44, 60, 80])

    # ---------- [A] mass as eigenstrain sources curvature ----------
    theta = rho.copy()                                     # preferred local compression ~ mass
    eta = rlap(theta)                                      # eta = lap(theta*), the buckling curvature
    print("  [A] MASS AS EIGENSTRAIN (compression) sources curvature eta = lap(theta*):")
    print(f"      max|eta| / max|theta*| = {np.abs(eta).max()/theta.max():.3f}  -- nonzero, so mass")
    print("      does curve the medium locally (the heated-plate buckle). Not the flat compatible")
    print("      displacement of test_einstein_source [A]: this is genuine local Gaussian curvature.\n")

    # ---------- [B] Gauss-Bonnet: net curvature charge is zero ----------
    print("  [B] BUT ITS NET CURVATURE CHARGE VANISHES -- Gauss-Bonnet. integral(eta) = integral(lap")
    print("      theta*) is a total derivative, zero for any localized compression, at every width:")
    print(f"      {'width':>7} {'integral(eta)/max|theta*|':>26}")
    for w in (2.0, 3.0, 5.0, 8.0):
        th = np.exp(-R ** 2 / (2 * w ** 2))
        print(f"      {w:>7.1f} {rlap(th).sum()/th.max():>26.2e}")
    print("      => zero to machine precision. A smooth localized compression makes a curvature")
    print("         DIPOLE -- a dome in the core, a saddle in the ring -- with no net deficit angle.")
    print("         Net curvature charge is topological, and a smooth deformation carries none.\n")

    # ---------- [C] no net charge -> no long-range bending ----------
    print("  [C] ZERO NET CHARGE MEANS NO LONG-RANGE BENDING. Deflection from the compression")
    print("      curvature, beside a genuine curvature charge (eta ~ rho, a disclination density")
    print("      proportional to the mass -- what Einstein needs):")
    exx_c, eyy_c, exy_c = incompatible_strain(eta, KX, KY)
    a_comp = deflect(exx_c - 0.5 * (exx_c + eyy_c), N, bs)
    eta_E = rho - rho.mean()
    exx_E, eyy_E, exy_E = incompatible_strain(eta_E, KX, KY)
    a_E = deflect(exx_E - 0.5 * (exx_E + eyy_E), N, bs)
    print(f"      {'b':>5} {'compression |alpha|':>20} {'Einstein charge |alpha|':>24}")
    for i, b in enumerate(bs):
        print(f"      {b:>5} {abs(a_comp[i]):>20.4f} {abs(a_E[i]):>24.4f}")
    print("      => the compression deflection is zero at every impact parameter; the genuine")
    print("         curvature charge bends light and stays long-range. Compression gives gamma = 0.\n")

    # ---------- [D] the complete picture ----------
    print("  [D] WHAT gamma = 1 REQUIRES, and why every channel misses it:")
    print("      gamma = 1 needs the light-bending curvature charge to BE the mass: eta ~ rho, a net")
    print("      curvature (disclination density) proportional to rho. The model's channels give:")
    print("        smooth induced loop   : coupling to spatial stress = 0   (spin selection rule)")
    print("        elastic body force    : compatible displacement, eta = 0 (test_einstein_source)")
    print("        compression eigenstrain: eta = lap(theta*), net charge = 0 (Gauss-Bonnet, here)")
    print("      Each sources ZERO net curvature charge, for its own reason. Net curvature charge is")
    print("      carried only by a topological DISCLINATION, and mass does not nucleate one.\n")

    print("[verdict] no direct mechanism gives gamma = 1; every channel sources zero net curvature.")
    print("  * Compression -- the model's own gravity -- does curve the medium (eta = lap theta* =/= 0),")
    print("    but by Gauss-Bonnet its net curvature charge is exactly zero: a smooth localized")
    print("    eigenstrain is a curvature dipole with no deficit angle, so it bends no light at range.")
    print("    gamma from compression is 0, and the reason is a theorem, not a small number.")
    print("  * The smooth loop gave 0 by a spin selection rule (test_gamma_source), the elastic")
    print("    force gave 0 by strain compatibility (test_einstein_source), and compression gives 0")
    print("    by Gauss-Bonnet. gamma = 1 requires mass to carry net curvature charge -- to nucleate")
    print("    a disclination density proportional to itself -- the Einstein source coupling no")
    print("    sector supplies.")
    print("  * HONEST DIMENSIONAL CAVEAT, because the result is a negative. Gauss-Bonnet is a 2D")
    print("    theorem and this measurement is 2D -- the setting where the project's incompatible /")
    print("    disclination picture of curvature is defined (test_incompatible_gravity,")
    print("    test_light_bending, test_disclination_force), and test_gamma_source is likewise 2+1D.")
    print("    In 3+1D general relativity a smooth mass sources smooth Ricci curvature with NO")
    print("    topological defect -- that IS gamma = 1 -- so the 2D obstruction is not automatically")
    print("    the final word in 3+1D. The 3D incompatibility of an isotropic eigenstrain is the")
    print("    Einstein-tensor-shaped delta_ij lap(theta*) - d_i d_j theta*, whose curvature content")
    print("    is a separate calculation not done here. What is settled is the 2D compression channel;")
    print("    what is not is whether the 3+1D medium evades Gauss-Bonnet the way GR does.")
    print("  * So the model's gravity, as directly realized, couples to mass as a SCALAR (the")
    print("    Newtonian T00 alone): Nordstrom, gamma = 0, not Einstein. gamma = 1 survives only as an")
    print("    INFRARED-EMERGENT claim -- Weinberg's theorem on a massless spin-2 coupled to the")
    print("    conserved IR stress tensor -- on the same footing as the emergent graviton Ward")
    print("    identity (test_lattice_ward). Unlike the model's emergent Lorentz invariance, which is")
    print("    directly confirmed on the cone, that emergence has no positive in-model evidence and now")
    print("    a consistent set of direct results against it at the accessible scale. That is the")
    print("    sharpest honest statement of where the gravitational arc stands: Newtonian gravity is")
    print("    real and healthy; the Einstein completion is argued, not demonstrated.")
