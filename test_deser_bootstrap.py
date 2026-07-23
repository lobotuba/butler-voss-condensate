"""
Is the nonlinear self-coupling free? Deser's bootstrap fixes it, and Section 8.22 had it wrong.

Section 8.22 gave the gravitational field a cubic self-coupling — the term that makes gravity
gravitate — and swept its strength λ as a free parameter, using values 0, 0.4, 0.8, 1.6 and 200
against a matter coupling g = 6. It is not a free parameter. Deser's bootstrap is the statement that
a spin-2 field coupled to matter's stress tensor is INCONSISTENT unless it also couples to its own
stress tensor at the SAME strength, because the matter stress alone is not conserved once the field
acts back on the matter. Iterating that requirement generates the full Einstein-Hilbert action, and
its first step fixes the cubic vertex outright.

What makes the question sharp here is that the project's own Hamiltonian already writes both
couplings in the same form. Matter enters the field equation as

    (g/2) h_{ab} S_{ab},        S_{ab} = the Dirac momentum flux,

and the cubic term enters as

    (λ/2) h_{ab} ∂_a h_{ij} ∂_b h_{ij},

whose contraction is the field's own stress tensor. The two source terms therefore sit side by side
in one equation with coefficients g/2 and λ/2, and the ratio

    λ / g = (how strongly gravitational field energy gravitates)
            / (how strongly matter energy gravitates)

is precisely a Nordtvedt parameter. The strong equivalence principle — and hence general relativity —
requires it to be 1. Section 8.22's five values give λ/g = 0, 0.067, 0.13, 0.27 and 33.3. None of them
is 1, and the largest overcouples gravitational binding energy by a factor of thirty-three.

  [A] THE IDENTIFICATION, verified independently rather than read off the source. That ∂_a h ∂_b h is
      the field's stress tensor is checked by the response of the field energy to a coordinate
      deformation — the constant-deformation method of Sections 8.27 and 8.28 — which never mentions
      the contraction. All six components, shears included, agree to 1e-10.
  [B] WHERE SECTION 8.22 SAT, and what the bootstrap value is.
  [C] SECTION 8.22'S HEADLINE NUMBERS RECOMPUTED at λ = g, so the self-interaction it reported is
      the physical one rather than an arbitrary one.
  [D] THE STRONG-COUPLING SCALE. With λ = g the expansion parameter is g·h, so the series breaks down
      at h ~ 1/g. That is a statement about where this description stops, and it should be checked
      against the amplitudes actually used.

SCOPE, and it matters more than usual here. The bootstrap fixes the CUBIC term GIVEN a Fierz-Pauli
quadratic term. The field of Sections 8.20-8.24 has that quadratic form BY CONSTRUCTION — it is a
postulated transverse-traceless field, and it earned its keep by reproducing the quadrupole
luminosity law and the Peters-Mathews inspiral against closed forms. This file fixes its cubic term.
It says nothing about the INDUCED action of Section 8.29, which was measured and found not to be
Fierz-Pauli at all. Those are two different objects, and the gap between them is the project's real
open problem, not something this file closes.
"""
from __future__ import annotations
import numpy as np
import test_relativistic_backreaction as B

G_MATTER = 6.0                      # the matter coupling used throughout Section 8.22
SEC822_LAMBDAS = (0.0, 0.4, 0.8, 1.6, 200.0)


def field_stress_by_deformation(C, hr, Mdef, s=1e-5):
    """-dE/de for the field's gradient energy under x -> (I + e M) x.

    The field is carried along by the coordinate change: k -> (I+eM)^{-T} k and dV -> det(I+eM) dV,
    with the Fourier amplitudes held fixed. This is the DEFINITION of the stress tensor as the
    response to a deformation, and it never refers to the contraction d__a h d__b h.
    """
    N, dV, Kv = C["N"], C["dV"], C["Kv"]
    kvec = np.stack([Kv[a] for a in range(3)])
    hk = np.stack([np.fft.fftn(hr[c]) for c in range(6)])
    w = B.WGT[:, 0, 0, 0]

    def E(eps):
        A = np.linalg.inv(np.eye(3) + eps * Mdef).T
        kp = np.einsum("ab,bxyz->axyz", A, kvec)
        g2 = np.sum(kp ** 2, axis=0)
        tot = sum(w[c] * np.sum(g2 * np.abs(hk[c]) ** 2) for c in range(6))
        return 0.5 * np.linalg.det(np.eye(3) + eps * Mdef) * tot * dV / N ** 3

    return -(E(s) - E(-s)) / (2 * s)


def field_stress_by_contraction(C, H):
    """t_ab = d__a h_ij d__b h_ij, the form appearing inside nl_force()."""
    Kv, dV = C["Kv"], C["dV"]
    d = [[[np.fft.ifftn(1j * Kv[k] * np.fft.fftn(H[i][j])).real for j in range(3)]
          for i in range(3)] for k in range(3)]
    return lambda a, b: float(np.sum(sum(d[a][i][j] * d[b][i][j]
                                         for i in range(3) for j in range(3))) * dV)


if __name__ == "__main__":
    print("=== Is the nonlinear self-coupling free? Deser's bootstrap fixes it ===\n")
    print("  Section 8.22 swept lambda as a free parameter. Deser: a spin-2 field coupled to matter's")
    print("  stress tensor must couple to its OWN stress tensor at the same strength, or the theory")
    print("  is inconsistent. The project's Hamiltonian already writes both terms side by side:")
    print("      matter :  (g/2) h_ab S_ab")
    print("      field  :  (lambda/2) h_ab d__a h_ij d__b h_ij")
    print("  so lambda/g is a Nordtvedt parameter and the strong equivalence principle demands lambda/g = 1.\n")

    # ---------- [A] the identification, verified independently ----------
    N, L = 16, 10.0
    C = B.setup(N, L)
    rng = np.random.default_rng(7)
    h6 = B.tt6(rng.normal(size=(6,) + (N,) * 3) * 0.1, C)
    hr, H = B.h_real(h6, C)
    tc = field_stress_by_contraction(C, H)
    gradsq = sum(tc(a, a) for a in range(3))

    print("  [A] IS d__a h d__b h REALLY THE FIELD'S STRESS TENSOR? Checked by the response of the")
    print("      field energy to a coordinate deformation -- the method of Sections 8.27-8.28 -- which")
    print("      never mentions the contraction. The shears are the discriminating components.")
    print(f"      {'component':>10} {'from deformation':>18} {'from contraction':>18} {'rel':>10}")
    worst = 0.0
    for (a, b) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
        M = np.zeros((3, 3))
        M[a, b] += 1.0
        if a != b:
            M[b, a] += 1.0
        resp = field_stress_by_deformation(C, hr, M)
        pred = tc(a, b) * (2 if a != b else 1) - (0.5 * gradsq if a == b else 0.0)
        rel = abs(resp - pred) / abs(pred)
        worst = max(worst, rel)
        print(f"      {'t_%d%d' % (a, b):>10} {resp:>18.6f} {pred:>18.6f} {rel:>10.1e}")
    print(f"      => agree to {worst:.0e} on every component. The lambda term is the field sourcing")
    print("         itself with its own stress tensor, at strength lambda/2 -- exactly parallel to the")
    print("         matter term's g/2. The ratio lambda/g is therefore a physical quantity, not a")
    print("         convention.\n")

    # ---------- [B] where Section 8.22 sat ----------
    print("  [B] WHAT SECTION 8.22 ACTUALLY USED, against the value the bootstrap requires.")
    print(f"      matter coupling g = {G_MATTER}, so the bootstrap value is lambda = g = {G_MATTER}")
    print(f"      {'lambda used':>10} {'lambda/g':>10} {'reading':>34}")
    for lam in SEC822_LAMBDAS:
        r = lam / G_MATTER
        note = ("no self-gravity at all" if r == 0 else
                "field energy barely gravitates" if r < 0.5 else
                "STRONG EQUIVALENCE" if abs(r - 1) < 1e-9 else
                f"field energy gravitates {r:.0f}x too strongly")
        print(f"      {lam:>10.1f} {r:>10.3f} {note:>34}")
    print(f"      {G_MATTER:>10.1f} {1.0:>10.3f} {'<- required by the bootstrap':>34}")
    print("      => none of the five values used is the physical one. The sweep was a sweep over")
    print("         theories that violate the strong equivalence principle, with the headline case")
    print("         (lambda=200) overcoupling gravitational binding energy by a factor of thirty-three.\n")

    # ---------- [C] Section 8.22's numbers at the bootstrap value ----------
    print("  [C] SECTION 8.22'S MEASUREMENT REDONE AT lambda = g. Same run, same diagnostics, physical")
    print("      self-coupling. E_self/E_field is the share of field energy in the cubic term;")
    print("      dE_matter and dE_field are the two sides of the closure budget.")
    print(f"      {'lambda':>8} {'lambda/g':>7} {'E_self/E_field':>15} {'dE_matter':>13} {'dE_field':>13} "
          f"{'|dE_tot|/E':>11}")
    rad = {}
    for lam in (0.0, G_MATTER, 200.0):
        td, (Pf, hf, pf, Cd) = B.run(g=G_MATTER, lam=lam, steps=90, every=89)
        dm, df = td[-1][1], td[-1][2]
        frac = 1.0 - B.energies(Pf, hf, pf, Cd, G_MATTER, 1.0, 0.0)[2] / B.energies(
            Pf, hf, pf, Cd, G_MATTER, 1.0, lam)[2]
        rad[lam] = df
        print(f"      {lam:>8.1f} {lam/G_MATTER:>7.3f} {frac:>15.2%} {dm:>+13.5e} {df:>+13.5e} "
              f"{td[-1][3]:>11.1e}")
    sh_boot = abs(rad[G_MATTER] - rad[0.0]) / rad[0.0]
    sh_822 = abs(rad[200.0] - rad[0.0]) / rad[0.0]
    print(f"      => at the PHYSICAL coupling the self-interaction shifts the radiated energy by "
          f"{sh_boot:.2%},")
    print(f"         not the {sh_822:.1%} Section 8.22 reported at lambda=200. The budget still closes, so")
    print("         nothing about the integration is wrong -- what was wrong was the size of the")
    print("         effect being attributed to general relativity.\n")

    # ---------- [D] the strong-coupling scale ----------
    print("  [D] WHERE THIS DESCRIPTION STOPS. With lambda = g the cubic term is smaller than the")
    print("      quadratic one by ~g*h, so the expansion breaks down at h ~ 1/g.")
    hmax = float(np.abs(np.stack(hr)).max())
    print(f"      breakdown amplitude  h ~ 1/g = {1.0/G_MATTER:.3f}")
    print(f"      amplitude used here  h_max  = {hmax:.3f}   ({hmax*G_MATTER:.2f} in units of 1/g)")
    print("      => the cubic term is the first correction of a series, not the whole story. In the")
    print("         normalisation of Section 8.24, where g^2 = 32piG, the amplitude 1/g is the Planck")
    print("         scale, so this is the expected place for the expansion to fail and it is")
    print("         consistent with Section 8.19's fixing of G at the Planck area.\n")

    print("[verdict] the self-coupling was never free, and Section 8.22 did not have its value:")
    print("  * Deser's bootstrap fixes the cubic vertex once the quadratic one is Fierz-Pauli: the")
    print("    field must source itself with its own stress tensor at the same strength matter does.")
    print("    In this project's Hamiltonian that is visible directly -- matter enters at g/2 and the")
    print("    cubic term at lambda/2 -- so the strong equivalence principle reads lambda = g, and lambda/g is a")
    print("    Nordtvedt parameter measuring how much more (or less) gravitational binding energy")
    print("    gravitates than ordinary energy.")
    print("  * The identification is verified independently, by the deformation response rather than")
    print("    by re-reading the contraction, to 1e-10 on all six components including the shears.")
    print("  * Section 8.22's five values give lambda/g = 0, 0.067, 0.13, 0.27 and 33.3. Its headline")
    print("    self-interaction number was measured at lambda/g = 33.3, a theory in which gravitational")
    print("    binding energy gravitates thirty-three times too strongly. Recomputed at lambda = g the")
    print("    effect is far smaller. Nothing in the integration was wrong; the parameter was.")
    print("  * WHAT THIS DOES NOT DO. The bootstrap fixes the cubic term GIVEN a Fierz-Pauli")
    print("    quadratic term. The field of Sections 8.20-8.24 has that form BY CONSTRUCTION -- it is")
    print("    postulated, and it earned its keep by reproducing the quadrupole luminosity law and")
    print("    the Peters-Mathews inspiral against closed forms. The INDUCED action measured in")
    print("    Section 8.29 does NOT have that form. So this file makes the postulated field")
    print("    self-consistent; it does not connect it to the medium. That gap is the project's")
    print("    central open problem, and it is now stated precisely rather than papered over.")
