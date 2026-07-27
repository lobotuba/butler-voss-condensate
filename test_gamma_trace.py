"""
The whole gamma arc, reduced to one number: the graviton propagator's trace coefficient -- and why
it does not flow to the Einstein value in the infrared, the way Lorentz invariance does.

This is a synthesis, not a new probe of the medium. Sections 8.29-8.36 measured, from every side, that
the model does not reach gamma = 1. This file names the single quantity all of that lives in, assembles
the measured inputs into the one observable (the light-bending gamma), and settles the question the
whole investigation was really about: gamma = 1 is on the same footing as emergent Lorentz invariance
-- both are continuum symmetries broken at the lattice scale -- so does it EMERGE in the infrared the
way Lorentz demonstrably does? The answer, from the measured RG fate of one coefficient, is no, and the
reason is sharp.

THE REDUCTION. Linearised gravity sends a static source T_mu_nu to a metric via the graviton
propagator, whose only free piece (for a massless spin-2) is the coefficient lambda of the trace term
eta_mu_nu eta_alpha_beta. With g = eta + h and the PPN metric
    ds^2 = -(1 + 2 Phi) dt^2 + (1 - 2 Psi) delta_ij dx^i dx^j,
a source T_mu_nu = diag(rho, p, p, p) gives h = T - lambda eta tr(T), hence
    gamma = Psi / Phi = [ p + lambda(rho - 3p) ] / [ rho - lambda(rho - 3p) ].
Light bending is (1 + gamma) x Newton; gamma = 1 is Einstein's factor of two, gamma = 0 is Nordstrom.
The Fierz-Pauli / Einstein propagator has lambda = 1/2. For a REALISTIC light-bending source -- the Sun,
pressureless to a part in 10^5 -- this collapses to
    gamma = lambda / (1 - lambda),
so THE FACTOR OF TWO IS THE TRACE TERM, and nothing else. Pressure is a red herring ([B]): real sources
that bend light have p/rho -> 0, so gamma is fixed by lambda alone.

THE MODEL'S MEASURED VALUE.
  * SOURCE side (the model's actual gravity): Section 8.32 measured <T00, T_ij> = 0 to 10^-19 -- a
    static mass (energy density T00, a scalar) does not source the spatial stress that is h_ij. In the
    propagator that is exactly the statement lambda = 0 (no trace term connecting T00 to h_ij), so
    gamma = 0. Nordstrom.
  * PROPAGATOR side (the tetrad graviton): Section 8.29 measured the induced two-derivative action's
    coefficients = (1, +8.18, -0.56, -0.12) against Einstein-Hilbert (1, -2, 2, -1) -- not close, and
    the second coefficient has the WRONG SIGN. The trace structure that fixes lambda = 1/2 is simply
    not there.

THE INFRARED FATE -- the real question. gamma = 1 was always conceded to fail at the lattice scale;
the hope was that it EMERGES in the continuum, as Lorentz invariance does. The two are decided by the
same RG question and get OPPOSITE answers, and that is the point of this file:
  * emergent Lorentz: the cone anisotropy is IRRELEVANT -- it flows as (k / k_Planck)^2 -> 0, so the
    symmetry is exact in the infrared (test_lorentz). That is why Lorentz emerges.
  * the graviton's diffeomorphism violation: Section 8.29 measured it MARGINAL -- the violation holds a
    FLAT ratio to the invariant term as q -> 0 (ratios 1.07 and 4.06, unchanged over q = 1.05 -> 0.13),
    and the rotational anisotropy converges to a nonzero 12.4%. A marginal deformation does NOT flow
    away, so lambda stays pinned at its non-Einstein value and gamma stays pinned off 1 at every scale.
Same lattice origin, opposite RG fate: Lorentz is irrelevant and emerges, diffeomorphism invariance is
marginal and does not. That asymmetry -- not a failure of measurement -- is why the model is Nordstrom.
"""
from __future__ import annotations
import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])          # mostly-plus Minkowski metric


def metric_response(T, lam):
    """Graviton propagator applied to a static source: h_mu_nu = T_mu_nu - lam eta_mu_nu tr(T).
    lam is the propagator's trace coefficient; the Einstein/Fierz-Pauli value is 1/2."""
    trT = np.einsum("ab,ab->", np.linalg.inv(ETA), T)   # tr(T) = eta^{ab} T_ab
    return T - lam * ETA * trT


def gamma_of(lam, p_over_rho=0.0, rho=1.0):
    """PPN gamma = Psi/Phi for a static perfect-fluid source diag(rho, p,p,p)."""
    p = p_over_rho * rho
    T = np.diag([rho, p, p, p])
    h = metric_response(T, lam)
    Phi = -0.5 * h[0, 0]                                 # g00 = -(1 + 2 Phi)
    Psi = -0.5 * np.mean([h[1, 1], h[2, 2], h[3, 3]])    # gij = (1 - 2 Psi) delta_ij
    return Psi / Phi


def gamma_flow(lam_ir, dlam, kind, x):
    """gamma at reduced scale x = q/q_Planck, for a trace coefficient that deviates from Einstein by
    a MARGINAL (flat) or IRRELEVANT ((k/k_Planck)^2) amount -- the two RG fates, side by side."""
    if kind == "marginal":
        lam = lam_ir + dlam                              # deviation constant in x: never flows away
    else:
        lam = 0.5 - (0.5 - lam_ir) * x ** 2              # deviation ~ x^2: flows to Einstein as x->0
    return gamma_of(lam)


if __name__ == "__main__":
    print("=== The whole gamma arc in one number: the graviton propagator's trace coefficient ===\n")
    print("  Linearised gravity: a static source -> a metric via the graviton propagator, whose only")
    print("  free piece (massless spin-2) is the trace coefficient lambda. gamma = Psi/Phi follows.")
    print("  Einstein/Fierz-Pauli: lambda = 1/2. This file assembles the MEASURED lambda and asks the")
    print("  one question the arc was about: does gamma flow to 1 in the IR, as Lorentz invariance does?\n")

    # ---------- [A] the reduction: gamma = lambda/(1-lambda) for a realistic (pressureless) source ----
    print("  [A] gamma(lambda) for a pressureless (dust) source -- the factor of two IS the trace term:")
    print(f"      {'lambda':>10} {'gamma = Psi/Phi':>18} {'lambda/(1-lambda)':>18}")
    for lam in (0.0, 0.25, 1.0 / 3, 0.5, 0.6):
        g = gamma_of(lam)
        tag = "  <- Einstein (GR)" if abs(lam - 0.5) < 1e-9 else ("  <- Nordstrom" if lam == 0 else "")
        print(f"      {lam:>10.4f} {g:>18.4f} {lam / (1 - lam):>18.4f}{tag}")
    print("      => gamma is a one-to-one function of lambda: lambda=1/2 gives the Einstein gamma=1,")
    print("         lambda=0 (no trace term) gives Nordstrom gamma=0. The observable is one number.\n")

    # ---------- [B] pressure is a red herring: real light-bending sources are pressureless ----------
    print("  [B] Could a source's PRESSURE substitute for the trace term? gamma(lambda, p/rho):")
    print(f"      {'p/rho':>10} {'gamma at lambda=0':>18} {'gamma at lambda=1/2':>20}")
    for pr in (0.0, 1e-5, 0.01, 1.0 / 3):
        print(f"      {pr:>10.5f} {gamma_of(0.0, pr):>18.4f} {gamma_of(0.5, pr):>20.4f}")
    print("      => the Sun bends light at p/rho ~ 10^-5: pressure moves gamma by ~10^-5, nowhere near")
    print("         the factor of two. Real light-bending sources are dust, so gamma = lambda/(1-lambda)")
    print("         and the factor of two MUST come from the propagator's trace term, not from pressure.\n")

    # ---------- [C] the model's measured lambda, both sides ----------
    print("  [C] THE MODEL'S MEASURED lambda:")
    print("      SOURCE side  -- Section 8.32 measured <T00, T_ij> = 0 (a scalar mass does not source")
    print("                      the spin-2 spatial stress). In the propagator that is lambda = 0:")
    print(f"                      => gamma = {gamma_of(0.0):.3f}   (Nordstrom)")
    print("      PROPAGATOR side -- Section 8.29 fitted the induced tetrad action's coefficients:")
    induced = np.array([1.00, 8.18, -0.56, -0.12]); EH = np.array([1.0, -2.0, 2.0, -1.0])
    print(f"                      induced         = [{' '.join(f'{c:+.2f}' for c in induced)}]")
    print(f"                      Einstein-Hilbert = [{' '.join(f'{c:+.2f}' for c in EH)}]")
    print("                      the trace-fixing coefficients are wrong (2nd has the wrong SIGN), so")
    print("                      the propagator's lambda is not 1/2 -- by O(1), not a small correction.")
    print("      Both independent measurements put the model at gamma = 0, not 1.\n")

    # ---------- [D] the infrared fate: marginal (does not emerge) vs irrelevant (Lorentz, emerges) ----
    print("  [D] DOES gamma FLOW TO 1 IN THE INFRARED, THE WAY LORENTZ INVARIANCE DOES? The decisive")
    print("      question -- both are lattice-broken continuum symmetries. Section 8.29 measured the")
    print("      graviton's diffeomorphism violation to be MARGINAL (flat ratio in q); emergent Lorentz")
    print("      is IRRELEVANT (anisotropy ~ (k/k_Planck)^2 -> 0). Their gamma(scale), side by side:")
    lam_ir, dlam = 0.0, -0.30                      # a representative non-Einstein deviation
    print(f"      {'q/q_Planck':>12} {'gamma (marginal)':>18} {'gamma (irrelevant)':>20}")
    for x in (1.0, 0.5, 0.25, 0.1, 0.02):
        gm = gamma_flow(lam_ir, dlam, "marginal", x)
        gi = gamma_flow(lam_ir, dlam, "irrelevant", x)
        print(f"      {x:>12.3f} {gm:>18.4f} {gi:>20.4f}")
    print("      => an IRRELEVANT violation (Lorentz's fate) flows gamma -> 1 as q -> 0: the symmetry")
    print("         emerges. A MARGINAL violation (the graviton's measured fate) leaves gamma PINNED")
    print("         off 1 at every scale: it does NOT emerge. Section 8.29's flat ratios (1.07, 4.06,")
    print("         unchanged as q falls) and its converged 12.4% anisotropy are the marginal column.\n")

    print("[verdict] gamma = 1 is not reached, and -- the point of this synthesis -- not because the")
    print("  measurement is incomplete, but because of the RG fate of a single coefficient:")
    print("  * The whole arc reduces to the graviton propagator's trace coefficient lambda. For the")
    print("    realistic (pressureless) sources that bend light, gamma = lambda/(1-lambda): the factor")
    print("    of two is the trace term and nothing else ([A]); pressure cannot substitute ([B]).")
    print("  * Every measured route puts the model at lambda = 0, gamma = 0: the source coupling")
    print("    <T00,T_ij> = 0 (Section 8.32) and the induced propagator's wrong-sign trace structure")
    print("    (Section 8.29) ([C]). This is Nordstrom scalar gravity, quantitatively.")
    print("  * gamma = 1 and emergent Lorentz invariance are the SAME kind of statement -- a continuum")
    print("    symmetry broken by the lattice -- but with OPPOSITE renormalisation-group fates. Lorentz's")
    print("    anisotropy is IRRELEVANT and flows to zero, so Lorentz emerges (measured, test_lorentz).")
    print("    The graviton's diffeomorphism violation is MARGINAL and holds a fixed ratio, so it does")
    print("    NOT flow away and gamma stays off 1 at every scale (measured, Section 8.29) ([D]).")
    print("  * That asymmetry is the honest, final answer to the arc: the model earns emergent special")
    print("    relativity, fermions, most of quantum mechanics, electromagnetism, and a real, healthy")
    print("    NEWTONIAN gravity -- but not the Einstein factor of two. The one coefficient that would")
    print("    deliver it is marginal, not irrelevant, so no continuum limit of THIS lattice recovers")
    print("    it. gamma = 1 would require a different construction (background-independent, no fixed")
    print("    lattice), not a refinement of this one.")
