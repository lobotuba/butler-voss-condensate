"""
Is a graviton mass induced? Testing the assumption the whole gravity arc rests on.

*** STATUS UPDATE -- the gap this file names at the end has since been CLOSED, and it turned out not
    to be a separate question. The verdict below leaves the DECONFINED CURVATURE sector's q = 0 term
    unmeasured and calls it the sharpest load-bearing assumption remaining in the arc.
    test_curvature_mass shows it is the COSMOLOGICAL CONSTANT in different words. A vacuum energy
    enters the action as sqrt(g) Lambda, and expanding sqrt(g) = 1 + tr h/2 + (tr h)^2/8 - tr(h^2)/4
    gives a quadratic-in-h term carrying NO derivatives -- which is what a mass term is:
    S ⊃ -(Lambda/8)[2 h_ij h_ij - (tr h)^2]. For transverse-traceless h that is -(Lambda/4) h_ij h_ij,
    so m^2 is proportional to Lambda for the propagating spin-2 modes (verified against the closed
    form to 8e-10, with the TT cases at exactly -0.5).
    Bare, it is order unity -- as fatal as the tetrad result measured below. What removes it is the
    self-sustained condensate vacuum of test_cosmological_constant: the gravitating quantity is the
    grand potential -P, which cancels to the available precision across 122 decades, against a
    rigid-vacuum control that retains 0.75. Residual mass ~2e-8 of the bare scale, some 1000x below
    the exclusion bound of test_ir_fixed_point.
    Nothing measured below is retracted: the tetrad is still unprotected, and masslessness still
    cannot be inherited from it. What has changed is that the curvature sector does not need to
    inherit it. HONEST: that protection is an EQUILIBRIUM CONDITION, not a symmetry, and is the same
    Volovik mechanism as the cosmological-constant result, so the two stand or fall together. ***

*** SECOND STATUS UPDATE -- the verdict below says the tetrad mode is unprotected but stops at the
    MASS. test_graviton_transversality measured the rest of its induced two-derivative action and the
    picture is worse, and more coherent: the tetrad action also violates linearised diffeomorphism
    invariance MARGINALLY (a flat ratio to the invariant term as q -> 0, so it never flows away),
    breaks rotational invariance by a converged ~12.4%, and fits Einstein-Hilbert coefficients
    (1, -2, 2, -1) so badly that one of them comes out with the wrong sign. The instrument used there
    reproduces the -0.2746 measured below to ten digits by an independent code path, which is the
    cross-calibration for both files.
    The door this file left open is therefore now shut in both directions: the tetrad can supply
    neither masslessness NOR Einstein structure. That strengthens rather than weakens the reading
    below, since the project's gravity was never the tetrad. ***

Section 8.26 established that a graviton mass is the one RELEVANT deformation: it grows with
distance and destroys the inverse-square tail, so the infrared attractor stands or falls on its
absence. But that section only bounded a mass it had inserted BY HAND into the propagator, and the
propagator used everywhere else in this project, κ∇⁴ - μ∇², simply has no mass term written in it.
Masslessness has therefore been an ASSUMPTION of the operator form throughout, never a measurement.
The induced-gravity calculation makes the omission explicit: test_induced_sign extracts μ as the q²
coefficient of Π(q) - Π(0), discarding Π(0) as a contact term. That discarded number is precisely the
mass candidate, and nobody had asked whether it vanishes.

This file asks. The question is posed so that no perturbative bookkeeping can hide the answer: rather
than assembling a bubble and a seagull and hoping the set is complete, the ground-state energy of the
filled fermion sea is computed as a function of a CONSTANT (q = 0) deformation. A constant deformation
is exactly the q → 0 limit, and the sea energy contains every order at once -- bubble, seagull, and
everything above. If a uniform deformation changes the energy at all, the corresponding boson has an
unprotected q = 0 term. There is nothing left to leave out.

Two deformations are compared, in the same regulator (a periodic Brillouin-zone torus with a gapped
Wilson-Dirac spectrum, the regulator whose photon Ward identity closes to 1e-16 in test_lattice_ward):

  [A] A CONSTANT GAUGE FIELD, k → k + A. This is the photon's mass term, and gauge invariance
      requires it to vanish exactly: on a torus a constant A merely re-samples a complete period of a
      periodic function. This is the control, and it shows what a PROTECTED mode looks like.
  [B] A CONSTANT TRACELESS CONE ANISOTROPY, v = (1+ε, 1-ε, 1), which is the h₊ polarisation of the
      tetrad metric -- the emergent graviton of Section 8.5, whose degrees of freedom are the shape of
      the Dirac cone. Its curvature at ε = 0 is the induced q = 0 term -- the mass candidate.

The contrast is the result, and the sign of the answer was not assumed in advance.
"""
from __future__ import annotations
import numpy as np


def sea_energy(N, M, v=(1.0, 1.0, 1.0), A=(0.0, 0.0, 0.0), r=1.0):
    """Filled Dirac-sea energy per site on the 3D Brillouin torus.

    E(k) = sqrt( Σ_i v_i² sin²(k_i+A_i) + M(k+A)² ),  M(k) = M + r Σ_i (1 - cos k_i).
    v_i deforms the cone (the tetrad / metric); A_i is a constant gauge potential.
    """
    g = (np.arange(N) + 0.5) / N * 2 * np.pi
    kx, ky, kz = np.meshgrid(g, g, g, indexing="ij")
    s2 = 0.0
    mass = M
    for ki, vi, ai in ((kx, v[0], A[0]), (ky, v[1], A[1]), (kz, v[2], A[2])):
        s2 = s2 + (vi * np.sin(ki + ai)) ** 2
        mass = mass + r * (1.0 - np.cos(ki + ai))
    return -float(np.mean(np.sqrt(s2 + mass ** 2)))          # filled negative-energy band


def curvature(f, h=1e-3):
    """d²f/dx² at 0 by a symmetric five-point stencil (cancels the linear term exactly)."""
    fm2, fm1, f0, fp1, fp2 = (f(-2 * h), f(-h), f(0.0), f(h), f(2 * h))
    return (-fm2 + 16 * fm1 - 30 * f0 + 16 * fp1 - fp2) / (12 * h ** 2)


if __name__ == "__main__":
    print("=== Is a graviton mass induced? The assumption under the whole gravity arc ===\n")
    print("  A constant (q=0) deformation is the q → 0 limit, and the sea energy contains every")
    print("  order at once -- bubble, seagull and above. If a uniform deformation changes the energy")
    print("  at all, that boson has an unprotected q = 0 term. Nothing can be left out.\n")

    N, M = 40, 0.6
    print(f"  Regulator: periodic BZ torus, N = {N}³, Wilson-Dirac gap M = {M}.\n")

    print("  [A] PHOTON CONTROL -- constant gauge field k → k + A. Gauge invariance requires the")
    print("      energy to be flat: on a torus a constant A only re-samples a complete period.")
    print(f"      {'direction':>12} {'d²E/dA²':>16}")
    ph = []
    for lab, vec in (("A along x", (1, 0, 0)), ("A along y", (0, 1, 0)),
                     ("A diagonal", (1, 1, 1))):
        c = curvature(lambda a, vv=vec: sea_energy(N, M, A=tuple(a * t for t in vv)))
        ph.append(abs(c))
        print(f"      {lab:>12} {c:>16.3e}")
    print("      => flat to ~1e-10 (and improving with refinement): the photon mass is ZERO, and")
    print("         zero because a symmetry forbids it, not because a number happened to cancel.")
    print("         This is what a protected mode looks like.\n")

    print("  [B] TETRAD GRAVITON -- constant traceless cone anisotropy v = (1+ε, 1-ε, 1), the h₊")
    print("      polarisation. Its curvature at ε = 0 is the induced q = 0 term.")
    print(f"      {'M (gap)':>9} {'d²E/dε² (mass²)':>18} {'vs photon':>14}")
    gr = []
    for MM in (0.4, 0.6, 1.0, 1.5):
        c = curvature(lambda e, m=MM: sea_energy(N, m, v=(1 + e, 1 - e, 1.0)))
        gr.append(c)
        print(f"      {MM:>9.2f} {c:>18.6f} {abs(c)/max(max(ph),1e-30):>14.1e}")
    print("      => nonzero, O(1), and eight orders above the photon control. Note the SIGN: it is")
    print("         NEGATIVE, so a uniform shear LOWERS the sea energy rather than costing it. The")
    print("         fermion loop does not merely give the cone a mass, it DESTABILISES the symmetric")
    print("         cone, and the medium's own shear rigidity has to overcome that to keep the")
    print("         undeformed cone stable at all. Either sign settles the question at issue: the")
    print("         q = 0 term is O(1), nothing forbids it, and the tetrad mode is unprotected.\n")

    print("  [C] IS IT JUST THE REGULATOR? Refinement check that the contrast is physical rather")
    print("      than a discretisation artifact.")
    print(f"      {'N':>5} {'photon d²E/dA²':>16} {'graviton d²E/dε²':>18}")
    for NN in (24, 32, 40, 48):
        pa = curvature(lambda a: sea_energy(NN, M, A=(a, 0.0, 0.0)))
        gg = curvature(lambda e: sea_energy(NN, M, v=(1 + e, 1 - e, 1.0)))
        print(f"      {NN:>5} {pa:>16.3e} {gg:>18.6f}")
    print("      => the photon's zero and the graviton's O(1) value are both stable under refinement:")
    print("         the contrast is a property of the theory, not of the grid.\n")

    print("[verdict] the tetrad mode is UNPROTECTED -- and that is consistent, not fatal:")
    print("  * The photon is massless because an exact lattice symmetry forbids the mass; the tetrad")
    print("    graviton has no such protection and duly acquires an O(1) term. That is the concrete")
    print("    content of Section 8.26's finding that the infrared attractor is EMPIRICAL rather than PROTECTED.")
    print("  * It independently corroborates a result the project already had by another route. A")
    print("    tetrad mode with an O(1) quadratic term at q=0 is not a massless long-range mediator,")
    print("    which is exactly why test_tetrad_force measured NO long-range force from the tetrad's")
    print("    1/r² field and why the elastic route was declared dead (Eshelby-Crum). Two unrelated")
    print("    calculations agree, and the negative sign further explains why that route could not be")
    print("    rescued by tuning: the fermion contribution pushes the wrong way.")
    print("  * WHAT THIS DOES NOT SHOW. The gravity the project actually claims is NOT the tetrad")
    print("    mode: it is the DECONFINED CURVATURE sector (Sections 8.11-8.12), whose masslessness")
    print("    comes from the biharmonic structure of the curvature field rather than from any")
    print("    symmetry protecting the cone shape. This file does not measure THAT sector's mass, and")
    print("    the assumption therefore survives here -- narrowed, not removed. What is now known is")
    print("    that masslessness cannot be inherited from the tetrad, so the curvature sector must")
    print("    supply it, and that remains the load-bearing and still-unverified assumption.")
