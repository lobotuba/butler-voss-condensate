"""
Does the photon's WIDTH bend light? Testing whether the zero-width ray assumption hid the factor of two.

Robert's question: Einstein noted a photon has no length, but it does have a transverse WIDTH -- could the
bending of light (the Eddington factor of two, gamma = 1) be a width effect that a zero-width ray misses?
Every gamma measurement in this project (test_light_bending, test_gamma_3d) ray-traced light as a
zero-width geodesic, so "light is a thin ray" is a genuine buried premise. This file tests it directly by
computing the deflection of a photon of FINITE width -- a Gaussian beam of transverse width w -- past a
mass, and comparing it to the zero-width ray. The centroid deflection is exact via Ehrenfest's theorem
(d<k_x>/dz = -<d_x V>_beam), so it isolates the width effect cleanly: the beam bends by its
profile-AVERAGED force, independent of diffraction (which only spreads the profile further).

The physics first, so the numbers have a frame. The factor of two is not a property of the photon; it is a
property of the GEOMETRY it moves through. Light deflection is (1 + gamma) x Newton: the "1" is the time
warp g00 (the potential slowing the wave, a scalar refractive index n = 1 - Phi), the "gamma" is the
SPATIAL curvature Psi (n picks up a second -Psi). Light samples BOTH, and samples them EQUALLY, for one
reason -- it moves at c: on a null path dx = c dt, so the space and time parts of the metric enter with
equal weight. That is the origin of the factor of two: null-ness, not width. The model already has null
light (emergent Lorentz), but its mass sources only compression -> n = 1 - Phi with NO -Psi term
(Psi = 0, measured in every channel, Sections 8.32-8.36). So there is no spatial curvature for a photon of
ANY width to feel. The centroid of any beam obeys Ehrenfest exactly -- it bends by the beam-averaged
grad(n) -- so a finite width only averages grad(n) over the profile, an O((w/b)^2) tidal correction, never
a factor of two. The factor of two lives in n (the metric), not in the probe's width.

Measured here, three ways:
  [A] The zero-width RAY deflection (the Born integral), for the model's index n = 1 - Phi (gamma = 0) and
      for a GR index n = 1 - 2 Phi (gamma = 1), to set the two reference deflections.
  [B] A finite-width PHOTON (centroid deflection, exact via Ehrenfest) through the MODEL's index,
      deflection vs width w: it tracks the ray (gamma = 0) with only an O((w/b)^2) correction. Width
      does not add the factor of two.
  [C] The SAME finite-width photon through the GR index -- the control: it deflects by ~2x (gamma = 1) at
      every width, proving the beam method SEES the factor of two when the geometry has it. So the null in
      [B] is physics (no Psi to feel), not a limitation of the ray picture.
"""
from __future__ import annotations
import numpy as np

# numpy >= 2.0 renamed np.trapz -> np.trapezoid (same signature, incl. axis=)
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ----------------------------------------------------------------------------- field
MASS = 6.0            # weak-lensing strength (linear regime, so Born ~ exact)
SOFT = 1.5            # core softening of the potential


def phi(x, z):
    """Newtonian potential of a point mass at the origin (2D lensing slice)."""
    return -MASS / np.sqrt(x * x + z * z + SOFT * SOFT)


def dphidx(x, z):
    r2 = x * x + z * z + SOFT * SOFT
    return MASS * x / r2 ** 1.5


def ray_deflection(b, zmax=400.0, nz=40001, gamma_index=0.0):
    """Zero-width ray: the Born deflection alpha = -integral d_x(n-1) dz along the straight path x=b.
    n - 1 = -(1 + gamma_index) Phi, so alpha = (1 + gamma_index) integral d_x Phi dz."""
    z = np.linspace(-zmax, zmax, nz)
    alpha_newton = _trapz(dphidx(b, z), z)         # the (1+gamma)=1 piece = "Newtonian" reference
    return (1.0 + gamma_index) * alpha_newton, alpha_newton


def beam_deflection(b, w, gamma_index, zmax=400.0, nz=6001, nx=2001):
    """Finite-width photon: the centroid deflection of a collimated Gaussian beam of width w.
    Ehrenfest's theorem makes the centroid exact and diffraction-independent --
        d<k_x>/dz = -<d_x V>_beam,   V = -(1 + gamma_index) Phi,
    so the beam bends by its transverse-PROFILE-AVERAGED force, integrated along the path. This is
    the clean finite-width generalisation of the ray integral (which is the w -> 0 delta limit); it
    isolates the width effect from diffraction (which only spreads the profile further, never adds
    a factor of two)."""
    x = np.linspace(b - 6 * w, b + 6 * w, nx)                  # transverse profile support
    G = np.exp(-(x - b) ** 2 / (2 * w * w)); G /= _trapz(G, x)
    z = np.linspace(-zmax, zmax, nz)
    XX, ZZ = np.meshgrid(x, z, indexing="ij")
    avg_force = _trapz(G[:, None] * dphidx(XX, ZZ), x, axis=0)   # <d_x Phi>_beam (z)
    return (1.0 + gamma_index) * _trapz(avg_force, z)


if __name__ == "__main__":
    print("=== Does the photon's WIDTH bend light? Testing the zero-width ray assumption ===\n")
    print("  The factor of two (gamma=1) is the SPATIAL curvature Psi, which light samples because it")
    print("  moves at c (null path: dx = c dt weights space and time equally) -- not because it is wide.")
    print("  The model's mass makes only compression: n = 1 - Phi (Psi = 0). Test if width changes that.\n")

    b = 24.0

    # ---------- [A] zero-width ray references ----------
    a_model, a_newton = ray_deflection(b, gamma_index=0.0)
    a_gr, _ = ray_deflection(b, gamma_index=1.0)
    print(f"  [A] ZERO-WIDTH RAY deflection at impact parameter b = {b:.0f} (Born integral):")
    print(f"      model index  n = 1 - Phi   : alpha = {a_model:.5f}   (gamma = 0, the '1' = time warp only)")
    print(f"      GR index     n = 1 - 2 Phi : alpha = {a_gr:.5f}   (gamma = 1, the factor of two)")
    print(f"      ratio GR/model = {a_gr / a_model:.4f}  <- the factor of two lives in the INDEX (in Psi).\n")

    # ---------- [B] finite-width photon through the MODEL index ----------
    print("  [B] FINITE-WIDTH PHOTON (centroid, exact via Ehrenfest) through the MODEL index n = 1 - Phi.")
    print("      If width supplied the factor of two, the implied gamma would climb toward 1 as w grows:")
    print(f"      {'width w':>9} {'w/b':>7} {'alpha_beam':>12} {'alpha/alpha_ray':>16} {'implied gamma':>14}")
    for w in (0.5, 2.0, 4.0, 6.0, 8.0):
        ab = beam_deflection(b, w, gamma_index=0.0)
        print(f"      {w:>9.1f} {w / b:>7.2f} {ab:>12.5f} {ab / a_model:>16.4f} {ab / a_newton - 1:>14.4f}")
    print("      => alpha stays at the ray value (implied gamma ~ 0) for every width; the only change is")
    print("         a small O((w/b)^2) tidal correction from averaging grad(n) over the profile. Giving")
    print("         the photon a width does NOT create the factor of two.\n")

    # ---------- [C] finite-width photon through a GR index -- the control ----------
    print("  [C] CONTROL -- the SAME finite-width photon through a GR index n = 1 - 2 Phi (Psi = Phi):")
    print(f"      {'width w':>9} {'w/b':>7} {'alpha_beam':>12} {'alpha/alpha_ray':>16} {'implied gamma':>14}")
    for w in (0.5, 4.0, 8.0):
        ab = beam_deflection(b, w, gamma_index=1.0)
        print(f"      {w:>9.1f} {w / b:>7.2f} {ab:>12.5f} {ab / a_gr:>16.4f} {ab / a_newton - 1:>14.4f}")
    print("      => the beam deflects by ~2x (implied gamma ~ 1) at every width -- the width average")
    print("         SEES the factor of two WHEN the geometry has it. So the null in [B] is physics, not")
    print("         the ray approximation hiding something. (Diffraction, omitted here, only spreads the")
    print("         profile further and washes out width effects more -- it cannot add the factor of two.)\n")

    print("[verdict] the photon's width is not the missing factor of two.")
    print("  * The factor of two is the SPATIAL curvature Psi, carried in the refractive index n = 1 - Phi")
    print("    - Psi. Light samples it because it moves at c (a null path weights space and time equally),")
    print("    not because of any transverse width. A wide photon and a thin ray both bend by (1+gamma)x.")
    print("  * A finite-width photon through the model's index bends by exactly the ray amount (gamma = 0)")
    print("    up to an O((w/b)^2) tidal correction [B]; through a GR index it bends by 2x (gamma = 1) at")
    print("    every width [C]. Width changes the correction, never the factor of two.")
    print("  * So the model's gamma = 0 is robust to giving the photon a width: the mass sources no Psi")
    print("    for a probe of ANY size to feel. The zero-width ray was not hiding the factor of two -- the")
    print("    spatial curvature that would carry it is simply absent (Sections 8.32-8.36). Einstein's")
    print("    factor of two is about light's SPEED, not its width; the model has the speed but not the")
    print("    spatial curvature.")
