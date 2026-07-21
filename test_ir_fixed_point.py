"""
Does the infrared fixed-point claim survive scrutiny? Measuring what "fixed point" actually asserts.

The report's largest single claim is that general relativity is reached as an INFRARED FIXED POINT.
It is also the least supported, and the gap is worth naming precisely rather than defending. "Fixed
point" is a renormalisation-group term. Until now nothing in this project performed a
renormalisation-group analysis: no operator was classified as relevant or irrelevant, no statement was
made about what the infrared forgets, and the phrase was carrying authority from a calculation that
had not been done. This file does that calculation, in the only form the model can actually support.

The operational content of an infrared fixed point is two statements, both measurable:

  (i)  UNIVERSALITY. Long-distance observables must become insensitive to the ultraviolet couplings.
       If the far field still remembers the lattice-scale coefficient, there is no fixed point.
  (ii) OPERATOR CLASSIFICATION. Deformations must sort into irrelevant ones (which the infrared
       forgets) and relevant ones (which destroy it). A claim of a fixed point is empty unless that
       sorting is exhibited.

Both are computed here from the exact Green's function of the deconfined sector. Adding a graviton
mass term to the biharmonic-plus-Einstein propagator, κq⁴ + μq² + μm², and factorising
κ(q² + A)(q² + B), the three-dimensional transform is closed-form,

    G(r) = [ exp(-√B r) - exp(-√A r) ] / ( 4π κ (A - B) r ),   A + B = μ/κ,  AB = μm²/κ,

which reduces at m = 0 to the (1 - exp(-r/ℓ))/(4πμr) of Section 8.11. No box, no lattice, no
periodic images: the quantity being tested is not contaminated by the tool testing it.

What is measured:
  [A] THE INFRARED FORGETS THE ULTRAVIOLET. The higher-derivative coefficient κ is varied over four
      decades at fixed μ. The far-field Green's function is unchanged, and the residual dependence
      falls as exp(-r/ℓ) — the infrared forgets the ultraviolet EXPONENTIALLY, not merely as a power.
      This is universality, and it is the substance of the claim.
  [B] THE OPERATOR SORTING, MEASURED. The same observable is exposed to two deformations. Multiplying
      κ by 10000 leaves the far field untouched (irrelevant). A graviton mass smaller than the
      Einstein scale by many orders destroys it outright (relevant). That contrast is the
      renormalisation-group content the phrase was asserting without support.
  [C] HOW MASSLESS, QUANTITATIVELY. Since a mass is the one relevant deformation, the claim stands or
      falls on its absence. The measured inverse-square tail is converted into an exclusion bound on
      any induced graviton mass, expressed against the model's own crossover scale ℓ.

What this does NOT establish, and the report should not say otherwise. The fixed point demonstrated
here is EMPIRICAL, not PROTECTED. Section 8.12 established that diffeomorphism invariance is not an
exact lattice symmetry, so nothing forbids a graviton mass from being induced; masslessness is
measured (in the tail exponent here, and independently in the radiation speed and the two
polarisations of Section 8.20) but it is not guaranteed by a symmetry the way the photon's masslessness
is guaranteed by the exactly-closing U(1) Ward identity. A fixed point that is observed rather than
protected is a weaker object, and the distinction is the honest core of the claim. Nor does anything
here reach NONLINEAR general relativity: the sector analysed is the linearised propagator, so what is
shown is that the linearised Einstein term is the attractor of the long-distance theory — not that the
full Einstein equations are.
"""
from __future__ import annotations
import numpy as np


def green(r, kappa, mu, m=0.0):
    """Exact G(r) for (κ∇⁴ - μ∇² + μm²)G = δ, by partial fractions on κ(q²+A)(q²+B)."""
    S = mu / kappa                      # A + B
    P = mu * m ** 2 / kappa             # A * B
    disc = S ** 2 - 4 * P
    if disc <= 0:
        raise ValueError("mass too large: poles collide (m > sqrt(mu/4kappa))")
    rt = np.sqrt(disc)
    A, B = 0.5 * (S + rt), 0.5 * (S - rt)
    if B == 0.0:
        return (1.0 - np.exp(-np.sqrt(A) * r)) / (4 * np.pi * mu * r)
    return ((np.exp(-np.sqrt(B) * r) - np.exp(-np.sqrt(A) * r))
            / (4 * np.pi * kappa * (A - B) * r))


def force_slope(kappa, mu, m=0.0, r1=200.0, r2=400.0, n=400):
    """log-log slope of the force -dG/dr over a far window (Newtonian tail => exactly -2)."""
    r = np.linspace(r1, r2, n)
    dr = r[1] - r[0]
    G = green(r, kappa, mu, m)
    F = -(G[2:] - G[:-2]) / (2 * dr)
    return np.polyfit(np.log(r[1:-1]), np.log(np.abs(F)), 1)[0]


if __name__ == "__main__":
    print("=== Is the infrared fixed-point claim earned? Measuring what it asserts ===\n")
    print("  'Fixed point' is a renormalisation-group statement. Its operational content is that the")
    print("  infrared FORGETS the ultraviolet, and that deformations SORT into relevant and")
    print("  irrelevant. Neither had been measured in this project. Both are measured here.\n")

    mu = 1.0
    print("  [A] THE INFRARED FORGETS THE ULTRAVIOLET. κ (the higher-derivative, lattice-scale")
    print("      coefficient) varied over four decades at fixed μ = 1; the Newtonian prediction is")
    print("      G → 1/(4πμr), with NO dependence on κ at all.")
    print("      κ enters ONLY through the crossover ℓ = √(κ/μ); probing each at the same physical")
    print("      depth into the far field (r = 50ℓ) isolates whether the LAW itself remembers κ.")
    print(f"      {'κ':>10} {'ℓ=√(κ/μ)':>10} {'r = 50ℓ':>10} {'G·4πμr':>14} {'force slope':>12}")
    for kappa in (1e-2, 1e-1, 1e0, 1e1, 1e2):
        ell = np.sqrt(kappa / mu)
        rp = 50.0 * ell
        g = green(np.array([rp]), kappa, mu)[0] * 4 * np.pi * mu * rp
        print(f"      {kappa:>10.0e} {ell:>10.3f} {rp:>10.2f} {g:>14.10f} "
              f"{force_slope(kappa, mu, 0.0, 200*ell, 400*ell):>12.5f}")
    print("      => across 10000× in the ultraviolet coefficient the far-field law is Newton to ten")
    print("         figures with exponent -2 to five. κ sets only WHERE the crossover sits, never")
    print("         what lies beyond it: the residual is exp(-r/ℓ), so the infrared forgets the")
    print("         ultraviolet EXPONENTIALLY. That is universality, and it is the claim's substance.\n")

    print("  [B] THE OPERATOR SORTING, MEASURED. Relevance is a statement about the r → ∞ LIMIT, not")
    print("      about any one radius, so each deformation is tracked outward. Deviation is")
    print("      |G·4πμr - 1|; it must fall to 0 for an irrelevant operator and rise to 1 for a")
    print("      relevant one.")
    print(f"      {'deformation':>32} {'r=3ℓ*':>10} {'r=10ℓ*':>10} {'r=30ℓ*':>10} {'r=100ℓ*':>10}"
          f"   verdict")
    ell_ref = np.sqrt(1.0 / mu)
    for lab, kap, mm in (("κ × 1        (reference)", 1.0, 0.0),
                         ("κ × 100      (higher-deriv.)", 100.0, 0.0),
                         ("κ × 10000    (higher-deriv.)", 10000.0, 0.0),
                         ("graviton mass m = 1e-4", 1.0, 1e-4),
                         ("graviton mass m = 1e-3", 1.0, 1e-3)):
        # scale radii by that deformation's OWN crossover, so each is probed at equal depth
        lref = np.sqrt(kap / mu)
        devs = []
        for f in (3, 10, 30, 100):
            rp = f * lref
            devs.append(abs(green(np.array([rp]), kap, mu, mm)[0] * 4 * np.pi * mu * rp - 1.0))
        trend = "irrelevant" if devs[-1] < devs[0] * 1e-3 else "RELEVANT"
        print(f"      {lab:>32} " + " ".join(f"{d:>10.2e}" for d in devs) + f"   {trend}")
    print("      => the two behave oppositely in the limit. Higher-derivative structure dies as")
    print("         exp(-r/ℓ) and is gone by 100ℓ however large κ is: IRRELEVANT. A graviton mass")
    print("         ten thousand times smaller than the Einstein scale instead GROWS outward until")
    print("         it removes the tail entirely: RELEVANT. That sorting is the renormalisation-")
    print("         group content the phrase 'fixed point' was asserting without support.\n")

    print("  [C] HOW MASSLESS, QUANTITATIVELY. The claim rests entirely on the relevant operator")
    print("      being absent, so the inverse-square tail is turned into an exclusion bound on any")
    print("      induced graviton mass. Tolerance: the fitted force exponent departing from -2.")
    print(f"      {'tolerance |slope+2|':>21} {'excluded above m':>18} {'m·ℓ':>12} "
          f"{'Compton/ℓ':>12}")
    ell1 = np.sqrt(1.0 / mu)
    for tol in (1e-2, 1e-3, 1e-4):
        lo, hi = 1e-8, 1e-1
        for _ in range(60):                       # bisect on the tolerance crossing
            mid = np.sqrt(lo * hi)
            if abs(force_slope(1.0, mu, mid) + 2.0) > tol:
                hi = mid
            else:
                lo = mid
        print(f"      {tol:>21.0e} {hi:>18.3e} {hi*ell1:>12.3e} {1.0/(hi*ell1):>12.3e}")
    print("      => the measured tail excludes any induced graviton mass above ~3e-5 in units of the")
    print("         Einstein scale — a Compton wavelength exceeding the crossover scale ℓ by about")
    print("         four orders. The graviton is massless to the precision the model can resolve,")
    print("         which is a bound, not a proof of exact masslessness.\n")

    print("[verdict] the fixed-point claim is now measured — and it is narrower than it sounded:")
    print("  * WHAT IS EARNED. The infrared forgets the ultraviolet exponentially (universality), and")
    print("    the operators sort as a fixed point requires: higher-derivative structure irrelevant,")
    print("    a graviton mass relevant. The linearised Einstein term is the ATTRACTOR of the long-")
    print("    distance theory. That is a real renormalisation-group statement and it is measured.")
    print("  * WHAT IS NOT. The fixed point is EMPIRICAL, not PROTECTED. Diffeomorphism invariance is")
    print("    not an exact lattice symmetry (Section 8.12), so nothing FORBIDS an induced graviton")
    print("    mass — the one relevant deformation. Masslessness is measured, not guaranteed, and")
    print("    that is a weaker object than the photon's symmetry-protected masslessness.")
    print("  * SCOPE. This is the LINEARISED propagator. Nothing here reaches the nonlinear Einstein")
    print("    equations, and γ = 1 remains ARGUED (Weinberg on the conserved infrared stress tensor)")
    print("    rather than measured, because the direct check is regulator-limited (Section 8.12).")
    print("  * The report should therefore say that the LINEARISED Einstein term is reached as an")
    print("    infrared attractor, with the relevant deformation measured to be absent but not")
    print("    symmetry-protected — not that 'general relativity is reached as an infrared fixed")
    print("    point', which claims the nonlinear theory and a protection neither of which is shown.")
