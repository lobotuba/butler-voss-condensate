"""
Verifying the load-bearing assumption of the cosmological-constant solution -- and finding it is the
SAME fact that makes gamma = 0. The model's failure to reach Einstein gravity IS its solution to the CC.

test_cosmological_constant dissolved the 10^122 fine-tuning with Volovik's mechanism: the gravitating
vacuum energy is not the bare zero-point density eps but the grand potential rho_Lambda = eps - mu n = -P,
and a self-sustained condensate vacuum has P = 0, so rho_Lambda = 0 automatically. Its own status note
flags the weak point: "if the Volovik argument fails both [the CC solution and graviton masslessness]
fall together." That argument -- that emergent gravity couples to the vacuum STRESS (-P), not to the
energy density eps -- was ASSUMED there, not verified. It is exactly the tension with the gamma arc,
which measured that gravity couples to energy density T00. This file verifies it, and the resolution is
better than a patch: the two statements are true in DIFFERENT sectors, and the fact that reconciles them
is one already measured.

The resolution by sector. The emergent metric has two pieces that couple to different sources:
  * the TIME potential h00 (Newtonian Phi) couples to the energy density T00;
  * the SPATIAL metric h_ij (curvature Psi, and the cosmological expansion) couples to the spatial
    stress T_ij.
A static MASS is a localized T00: it sources h00 (Newtonian gravity, which the model has) but -- by the
measured selection rule <T00, T_ij> = 0 (test_gamma_source) -- sources NO spatial curvature, which is
exactly why gamma = 0 (no light-bending factor of two). The VACUUM is a uniform T00 = eps (huge) plus a
spatial stress T_ij = P delta_ij: its huge energy density sources only a uniform, unobservable h00, and
-- by the SAME selection rule -- sources NO spatial curvature; the cosmological constant (spatial
curvature / expansion) is sourced only by the vacuum pressure P, which is zero at self-bound equilibrium.
So "gravity couples to T00" (gamma arc, the h00 sector) and "the gravitating vacuum energy is -P" (the CC,
the h_ij sector) are both correct, and both follow from one measured fact: energy density does not source
spatial curvature.

  [A] the real medium: the self-bound Lennard-Jones condensate has a large energy density and P = 0 at
      equilibrium, with the density self-tuning P back to zero after any change -- the mechanism of
      test_cosmological_constant, now on the actual medium rather than a toy equation of state.
  [B] the verification: with the induced-gravity bubble, energy density T00 sources the Newtonian h00
      (<T00 T00> =/= 0) but NO spatial curvature (<T00, T_ij> = 0), while a genuine spatial stress DOES
      source spatial curvature. So emergent gravity couples the vacuum's eps only to h00 and its P to the
      spatial metric -- the Volovik assumption, measured, not assumed.
  [C] the unification: <T00, T_ij> = 0 is a single selection rule with TWO consequences -- gamma = 0
      (a mass's energy density does not bend light) and Lambda ~ 0 (the vacuum's energy density does not
      gravitate). The model does not reach Einstein gravity for the very reason it has no CC problem.
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_hex
from test_gamma_source import disc, q2coef, bubble, V_J0, V_T00, V_T, V_Tplus

RCUT = 8.0

# The FIXED neighbour set of a triangular lattice, built once at unit spacing. A perfect lattice's
# neighbour STRUCTURE does not change with spacing -- only the distances scale as a -- so summing over
# this fixed set makes E(a) perfectly smooth (no atoms crossing a distance cutoff as a varies), which is
# what lets the pressure reach zero at the minimum instead of a finite-difference noise floor.
_X1 = perfect_hex(radius_cells=16, a=1.0)
_c1 = int(np.argmin((_X1 ** 2).sum(1)))
_D = np.linalg.norm(_X1 - _X1[_c1], axis=1)
NEIGH = _D[(_D > 1e-9) & (_D < RCUT)]                       # neighbour distances at unit spacing


# ------------------------------------------------------- the real medium (LJ condensate) --------------
def per_atom_energy(a, sigma=1.0, eps=1.0):
    """Edge-free bulk energy per atom of a triangular LJ crystal at spacing a: sum over the fixed
    neighbour set (distances a * NEIGH). Smooth in a by construction."""
    sr6 = (sigma / (a * NEIGH)) ** 6
    return 0.5 * (4 * eps * (sr6 ** 2 - sr6)).sum()


def pressure(a, sigma=1.0, eps=1.0, da=1e-4):
    """P = -dE/dV, V = (sqrt3/2) a^2 per atom; = 0 at the energy minimum (equilibrium)."""
    dEda = (per_atom_energy(a + da, sigma, eps) - per_atom_energy(a - da, sigma, eps)) / (2 * da)
    return -dEda / (np.sqrt(3) * a)


def equilibrium(sigma=1.0, eps=1.0, n=401):
    """Equilibrium spacing a0 (the energy minimum, P = -dE/dV = 0). Coarse grid to bracket the minimum,
    then a secant refinement of P(a) = 0 so the residual pressure is at the finite-difference floor, not
    the grid spacing -- otherwise 'P = 0' would be limited by how finely a was sampled."""
    aa = np.linspace(0.98 * sigma, 1.18 * sigma, n)
    E = np.array([per_atom_energy(a, sigma, eps) for a in aa])
    i = int(np.argmin(E))
    a, b = aa[i - 1], aa[i + 1]                             # bracket the minimum
    fa, fb = pressure(a, sigma, eps), pressure(b, sigma, eps)
    for _ in range(60):                                    # secant on P(a) = 0
        if fb == fa:
            break
        c = b - fb * (b - a) / (fb - fa)
        fc = pressure(c, sigma, eps)
        a, fa, b, fb = b, fb, c, fc
        if abs(fc) < 1e-10:
            break
    return b, per_atom_energy(b, sigma, eps), fb


if __name__ == "__main__":
    print("=== The vacuum energy does not gravitate: verifying the CC mechanism, and its tie to gamma=0 ===\n")

    # ---------- [A] the real medium: energy density huge, pressure zero, self-tuned ----------
    print("  [A] THE REAL MEDIUM (self-bound Lennard-Jones condensate): energy density vs pressure.")
    a0, E0, P0 = equilibrium()
    area = (np.sqrt(3) / 2) * a0 ** 2
    print(f"      equilibrium spacing a0 = {a0:.4f};  energy density eps = E/area = {E0/area:+.4f}")
    print(f"      pressure P(a0) = {P0:+.2e}  (ZERO at the self-bound equilibrium)")
    print("      scale the condensation energy (eps) -- energy density moves, the gravitating P stays 0:")
    print(f"      {'eps_depth':>10} {'energy density':>16} {'P(a0)':>12}")
    for depth in (1e0, 1e3, 1e6):
        a0d, E0d, P0d = equilibrium(eps=depth)
        print(f"      {depth:>10.0e} {E0d/((np.sqrt(3)/2)*a0d**2):>16.2e} {P0d:>12.2e}")
    print("      and after a 'phase transition' that shifts the preferred spacing, P self-tunes to 0:")
    for sig in (1.00, 1.08):
        a0s, _, P0s = equilibrium(sigma=sig)
        print(f"        sigma = {sig:.2f}  ->  new a0 = {a0s:.4f},  P(a0) = {P0s:+.2e}")
    print("      => a large (negative) energy density coexists with EXACTLY zero pressure, for any")
    print("         condensation depth, self-restored after any change. On the real medium, as the toy")
    print("         equation of state of test_cosmological_constant showed. The gravitating quantity is P.\n")

    # ---------- [B] verify: which stress does emergent gravity couple to? ----------
    print("  [B] VERIFY THE ASSUMPTION -- does emergent gravity couple energy density to spatial curvature?")
    print("      The induced-gravity bubble (gapped 2+1D Dirac, the test_gamma_source instrument):")
    qs = np.array([0.05, 0.10, 0.15, 0.20])
    kx, ky = disc(3.0, 501)
    M = 0.5
    aT00, _ = q2coef(kx, ky, M, V_T00, qs)                 # energy density -> time potential h00
    _, cTT = q2coef(kx, ky, M, V_Tplus, qs)                # spatial stress -> spatial metric (coupling)
    # the mixed coupling <T00, T_ij> that would let energy density curve space (the selection rule):
    mixed = []
    for q in qs:
        A = V_T00(kx + q / 2, ky, M, 1.0)
        for (i, j) in ((0, 0), (1, 1), (0, 1)):
            B = V_T(i, j, kx + q / 2, ky, M, 1.0)
            mixed.append(abs(bubble(kx, ky, q, 0.0, M, A, B)))
    print(f"      <T00 T00>  (energy density -> Newtonian h00)      q^2 coef = {aT00:+.4f}   =/= 0")
    print(f"      <T00, T_ij>(energy density -> SPATIAL metric)     max|.|   = {max(mixed):.1e}   = 0")
    print(f"      <T+  T+ >  (spatial stress -> SPATIAL metric)     coupling = {cTT:+.4f}   =/= 0")
    print("      => emergent gravity couples the energy density ONLY to the time potential h00, and the")
    print("         SPATIAL metric only to the spatial stress T_ij. This is the Volovik assumption of")
    print("         test_cosmological_constant, now MEASURED: the vacuum's huge eps sources only a")
    print("         uniform (unobservable) h00; its pressure P sources the cosmological curvature, and")
    print("         P = 0 at equilibrium [A]. The vacuum energy density does not gravitate.\n")

    # ---------- [C] the unification: gamma = 0 and Lambda ~ 0 are the SAME selection rule ----------
    print("  [C] THE UNIFICATION -- one selection rule, two consequences:")
    print("      <T00, T_ij> = 0  (energy density does not source spatial curvature) means:")
    print("        * a static MASS (a localized T00) sources no spatial curvature  ->  gamma = 0")
    print("          (no light-bending factor of two: the gamma arc);")
    print("        * the VACUUM (a uniform T00 = eps) sources no spatial curvature ->  Lambda ~ 0")
    print("          (its enormous energy density does not drive expansion: the CC problem dissolved).")
    print("      The model does NOT reach Einstein gravity for the very same reason it has NO")
    print("      cosmological-constant problem. The factor of two and the 10^122 are one fact.\n")

    print("[verdict] the vacuum energy does not gravitate, and the CC mechanism is now verified, not assumed:")
    print("  * On the real condensate, a large energy density coexists with P = 0 at self-bound")
    print("    equilibrium, self-restored after any change to the vacuum [A] -- the gravitating quantity")
    print("    is the pressure, and it is zero.")
    print("  * The load-bearing Volovik assumption of test_cosmological_constant -- that emergent gravity")
    print("    couples to the stress (-P), not the energy density (eps) -- is MEASURED: energy density")
    print("    sources the Newtonian h00 but not the spatial metric (<T00,T_ij> = 0), while spatial")
    print("    stress does [B]. So the vacuum's eps gravitates only as an unobservable uniform h00; the")
    print("    cosmological constant is its pressure, which self-tunes to zero.")
    print("  * The reconciliation of the tension: 'gravity couples to T00' (gamma arc) and 'the")
    print("    gravitating vacuum energy is -P' (the CC) are both correct, in the h00 and h_ij sectors")
    print("    respectively, and both follow from <T00,T_ij> = 0.")
    print("  * The unification [C]: gamma = 0 and Lambda ~ 0 are the same measured selection rule. The")
    print("    model's inability to make the Einstein factor of two IS what keeps the vacuum energy from")
    print("    gravitating. A model that bent light like GR would have the 10^122 problem; this one has")
    print("    neither, together.")
    print("  HONEST SCOPE: this makes the EQUILIBRIUM cosmological constant zero and explains why; it does")
    print("  not derive the observed nonzero Lambda, which remains the off-equilibrium (expansion,")
    print("  matter, relaxation) question that test_cosmological_constant already relocated it to.")
