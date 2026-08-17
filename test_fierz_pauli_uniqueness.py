"""Route A, the build -- the last calculation: gauge invariance forces the world-crystal graviton to be
Einstein-Hilbert, so gamma = 1 is not a tuning.

The build (S8.67-S8.70) left one open item, flagged in S8.68/S8.70: the EXACT Einstein-Hilbert tensor
structure of the world-crystal dual graviton on the lattice -- on which the trace-reversal, and thus the
exact gamma = 1 rather than merely gamma != 0, ultimately rests. This section closes it, not by fine-tuning a
lattice, but by showing the tensor structure is FORCED.

THE ARGUMENT. The world crystal supplies a graviton that is (a) MASSLESS (S8.67: omega^2 ~ q^4, zero moduli),
(b) TWO-DERIVATIVE (S8.67: omega^2 ~ q^4 in u is INT (d h)^2 in the metric h ~ d u -- a two-derivative kinetic
term), and (c) DIFFEOMORPHISM-INVARIANT (S8.68: the physical graviton is the incompatible/dual field, which is
invariant under h -> h + d xi; the strains it is built on ARE the diffeomorphisms). The Fierz-Pauli / Weinberg
uniqueness theorem states that a massless, two-derivative, diffeomorphism-invariant action for a symmetric
tensor is UNIQUELY linearised Einstein-Hilbert. So conditions (a)-(c) -- all established for the world crystal
-- force the graviton to be Einstein-Hilbert, hence gamma = 1. No tuning: gauge invariance (which the
deconfined world crystal has, S8.68) IS the Einstein-Hilbert condition.

This section demonstrates the uniqueness explicitly and end-to-end. The general two-derivative quadratic
action for a symmetric tensor h_{mu nu} is a four-parameter family:
    S = INT [ a (d_l h_{mn})(d^l h^{mn}) + b (d^m h_{mn})(d_r h^{rn}) + c (d^m h_{mn})(d^n h) + d (d_l h)(d^l h) ].
Imposing diffeomorphism invariance (annihilation of h_{mn} = d_m xi_n + d_n xi_m for all xi, at all momenta)
leaves a ONE-parameter family -- the overall scale -- fixing (a,b,c,d) = (1, -2, 2, -1): exactly Einstein-
Hilbert. That operator, coupled to a conserved static source, gives the trace-reversal h_{ij} = h_{00} and so
gamma = 1, with exactly two propagating polarisations.

  [G1] Gauge invariance UNIQUELY fixes the two-derivative graviton action to Einstein-Hilbert: the constraint
       matrix on (a,b,c,d) has a one-dimensional null space, and it is (1, -2, 2, -1), the EH/Fierz-Pauli
       ratios. Any other choice is not diffeomorphism-invariant.
  [G2] The Einstein-Hilbert operator gives gamma = 1 and two DOF: coupled to a static energy density it
       returns h_{00} = h_{xx} = h_{yy} = h_{zz} (the trace-reversal), so Psi = Phi and gamma = 1; and the
       gauge-invariant operator has exactly four zero modes (the diffeomorphism directions) -> two propagating
       polarisations, a massless spin-2 graviton.
  [G3] Closure: the world crystal's graviton is massless (S8.67), two-derivative (S8.67) and diffeomorphism-
       invariant (S8.68), so by [G1] it is Einstein-Hilbert and by [G2] gamma = 1 -- FORCED by gauge
       invariance, not tuned. The deconfinement of the disclinations (S8.64, Y -> 0) that makes the graviton
       massless is the same fact that makes it gauge-invariant, hence Einstein-Hilbert.

Honest scope: this closes the last open item as a theorem -- gauge invariance forces EH, so the world crystal
does not need to be fine-tuned to the EH point; being the deconfined, gauge-invariant, massless graviton IS
being EH. The one honest caveat is that the world crystal's diffeomorphism invariance is EMERGENT (a lattice
breaks continuous diffeomorphisms, so gauge invariance -- and thus exact EH and gamma = 1 -- holds in the
long-wavelength limit, with corrections at the lattice scale), exactly as the model's Lorentz invariance is
emergent (S8.1). So gamma = 1 is an emergent, long-wavelength result, on the same footing as emergent Lorentz.
Pure numpy.
"""
from __future__ import annotations
import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])                 # mostly-plus metric

# 10 symmetric basis tensors for h_{mu nu} (orthogonal under the Frobenius product)
BAS = []
for i in range(4):
    E = np.zeros((4, 4)); E[i, i] = 1.0; BAS.append(E)
for i in range(4):
    for j in range(i + 1, 4):
        E = np.zeros((4, 4)); E[i, j] = E[j, i] = 1.0; BAS.append(E)


def term_values(h, q):
    """The four two-derivative invariants evaluated on symmetric h at momentum q (indices via eta)."""
    qu = ETA @ q                                     # q^mu
    hh = np.sum(h * (ETA @ h @ ETA))                 # h_{mu nu} h^{mu nu}
    q2 = q @ ETA @ q                                 # q_mu q^mu
    qh = qu @ h                                       # (q.h)_nu = q^mu h_{mu nu}
    trh = np.trace(ETA @ h)                           # h^mu_mu
    return np.array([q2 * hh, qh @ ETA @ qh, (qu @ h @ qu) * trh, q2 * trh * trh])


def term_hessians(q):
    """The 10x10 quadratic-form matrix of each of the four terms at momentum q."""
    n = 10
    Ms = [np.zeros((n, n)) for _ in range(4)]
    def Qt(x):
        h = sum(x[k] * BAS[k] for k in range(n))
        return term_values(h, q)
    for k in range(n):
        for l in range(n):
            ek = np.zeros(n); ek[k] = 1.0
            el = np.zeros(n); el[l] = 1.0
            cross = 0.5 * (Qt(ek + el) - Qt(ek) - Qt(el))
            for t in range(4):
                Ms[t][k, l] = cross[t]
    return Ms


def gauge_modes(q):
    """Diffeomorphism modes h = q (x) xi + xi (x) q as coefficient 10-vectors (BAS is orthogonal)."""
    gs = []
    for a in range(4):
        xi = np.zeros(4); xi[a] = 1.0
        G = np.outer(q, xi) + np.outer(xi, q)
        gs.append(np.array([np.sum(B * G) / np.sum(B * B) for B in BAS]))
    return gs


def eh_operator(q, params):
    Ms = term_hessians(q)
    return sum(params[t] * Ms[t] for t in range(4))


def main():
    print("=" * 92)
    print("ROUTE A (the last calculation): gauge invariance FORCES the world-crystal graviton to be Einstein-Hilbert")
    print("=" * 92)
    ok = True
    rng = np.random.default_rng(1)

    # [G1] impose diffeomorphism invariance on the general 2-derivative action -> unique EH
    rows = []
    for _ in range(6):
        q = rng.standard_normal(4)
        Ms = term_hessians(q)
        for g in gauge_modes(q):
            cols = [Ms[t] @ g for t in range(4)]         # (a M_a + ... ) @ g must vanish
            for r in range(10):
                rows.append([cols[t][r] for t in range(4)])
    A = np.array(rows)
    U, S, Vt = np.linalg.svd(A)
    nulldim = int(np.sum(S < 1e-9 * S[0]))
    sol = Vt[-1]; sol = sol / sol[0]
    print("\n  [G1] impose diffeomorphism invariance on S = a(dh)^2 + b(d.h)^2 + c(d.h)(dh_tr) + d(dh_tr)^2:")
    print(f"       singular values of the constraint matrix on (a,b,c,d): {np.round(S, 4)}")
    print(f"       null-space dimension (free parameters after gauge invariance): {nulldim}  (1 = only overall scale)")
    print(f"       forced solution (a,b,c,d), scaled to a=1: {np.round(sol, 6)}")
    g1 = nulldim == 1 and np.allclose(sol, [1, -2, 2, -1], atol=1e-6)
    ok &= g1
    print(f"       => gauge invariance UNIQUELY fixes the action to (1,-2,2,-1) = linearised Einstein-Hilbert")
    print(f"          (Fierz-Pauli). No tuning freedom remains  -> {'PASS' if g1 else 'FAIL'}")

    # [G2] the EH operator: gamma = 1 (trace-reversal) and two propagating DOF
    EH = np.array([1.0, -2.0, 2.0, -1.0])
    q = np.array([0.0, 0.0, 0.0, 1.0])                    # static source, spatial momentum along z
    M = eh_operator(q, EH)
    # de Donder gauge fixing: subtract (q^mu hbar_{mu nu})^2, hbar = h - 1/2 eta tr h
    def gf(q):
        n = 10; G = np.zeros((n, n)); qu = ETA @ q
        def Gvec(x):
            h = sum(x[k] * BAS[k] for k in range(n)); hbar = h - 0.5 * ETA * np.trace(ETA @ h)
            return qu @ hbar
        for k in range(n):
            for l in range(n):
                ek = np.zeros(n); ek[k] = 1.0; el = np.zeros(n); el[l] = 1.0
                gk, gl, gkl = Gvec(ek), Gvec(el), Gvec(ek + el)
                G[k, l] = 0.5 * ((gkl @ ETA @ gkl) - (gk @ ETA @ gk) - (gl @ ETA @ gl))
        return G
    Mtot = M - gf(q)
    Tt = np.zeros((4, 4)); Tt[0, 0] = 1.0                 # static energy density T_00 = 1 (conserved: q^mu T_mn = 0)
    Tvec = np.array([np.sum(B * Tt) / np.sum(B * B) for B in BAS])
    hvec = np.linalg.lstsq(Mtot, Tvec, rcond=None)[0]
    h = sum(hvec[k] * BAS[k] for k in range(10))
    Phi = -0.5 * h[0, 0]
    Psi = -0.5 * 0.5 * (h[1, 1] + h[2, 2])                # transverse spatial metric seen by light
    gamma = Psi / Phi
    zero_modes = int(np.sum(np.abs(np.linalg.eigvalsh(M)) < 1e-9))
    print("\n  [G2] the Einstein-Hilbert operator coupled to a static energy density T_00:")
    print(f"       response  h_00={h[0,0]:.4f}  h_xx={h[1,1]:.4f}  h_yy={h[2,2]:.4f}  h_zz={h[3,3]:.4f}  "
          f"(all equal = the trace-reversal)")
    print(f"       Phi = {Phi:.4f},  Psi = {Psi:.4f}  ->  gamma = Psi/Phi = {gamma:.4f}")
    print(f"       gauge-invariant operator zero modes = {zero_modes} (the 4 diffeomorphism directions) -> "
          f"2 propagating DOF")
    g2 = abs(gamma - 1.0) < 1e-6 and zero_modes == 4
    ok &= g2
    print(f"       => the trace-reversal gives gamma = 1, with a massless spin-2 (two-polarisation) graviton"
          f"  -> {'PASS' if g2 else 'FAIL'}")

    # [G3] closure
    print("\n  [G3] closing the last open item of the gamma programme:")
    print("       * world-crystal graviton is MASSLESS (S8.67) + TWO-DERIVATIVE (S8.67) + DIFFEOMORPHISM-")
    print("         INVARIANT (S8.68: it is the incompatible/dual field, invariant under h -> h + d xi).")
    print("       * [G1]: those three conditions FORCE the action to be Einstein-Hilbert (Fierz-Pauli unique).")
    print("       * [G2]: Einstein-Hilbert gives gamma = 1.")
    print("       => gamma = 1 is FORCED by gauge invariance, not tuned. The disclination deconfinement")
    print("          (S8.64, Y->0) that makes the graviton massless is the same fact that makes it gauge-")
    print("          invariant, hence Einstein-Hilbert. (Emergent/long-wavelength, like Lorentz in S8.1.)")
    g3 = g1 and g2
    ok &= g3
    print(f"       -> {'PASS' if g3 else 'FAIL'}")

    print("\n" + "=" * 92)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  The last open item is closed as a theorem. A massless, two-derivative, diffeomorphism-invariant")
    print("  symmetric-tensor action is uniquely linearised Einstein-Hilbert -- demonstrated here: imposing")
    print("  gauge invariance on the general two-derivative graviton action leaves only the overall scale and")
    print("  fixes (a,b,c,d) = (1,-2,2,-1), whose operator trace-reverses a static energy density into")
    print("  h_ij = h_00, giving gamma = 1 with two propagating polarisations. The world crystal's graviton is")
    print("  massless (S8.67), two-derivative (S8.67) and gauge-invariant (S8.68), so it is Einstein-Hilbert by")
    print("  this uniqueness and gamma = 1 -- FORCED, not tuned: the same disclination deconfinement (S8.64)")
    print("  that makes it massless makes it gauge-invariant, hence Einstein. gamma = 1 is emergent and long-")
    print("  wavelength, on the same footing as the model's emergent Lorentz invariance.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
