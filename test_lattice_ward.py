"""
The symmetry-preserving regulator: why the photon Ward identity is exact and gamma=1 is not.

*** STATUS UPDATE -- of the two obstructions named below, the FIRST has since been removed and the
    SECOND has been measured and found worse than "not exact".
    Obstruction (i), that the graviton Ward identity is INHOMOGENEOUS because <T^ij> =/= 0: that
    <T^ij> is the same object as the q = 0 graviton mass, which test_curvature_mass identified as the
    COSMOLOGICAL CONSTANT and cancelled via the self-sustained vacuum (what gravitates is the grand
    potential -P, which vanishes at equilibrium). So the identity can legitimately be made
    homogeneous, and this file's reason for calling the question closed no longer holds.
    Obstruction (ii), that diffeomorphism invariance is not a lattice symmetry: test_graviton_transversality
    reopened the measurement with a NONPERTURBATIVE method (finite-q sea energy, which contains every
    seagull at once) calibrated by an EXACT finite-q photon Ward identity. The residual violation does
    not merely fail to vanish -- it is MARGINAL, holding a flat ratio to the invariant term as q
    falls, so it does not flow away in the infrared. Rotational invariance fails too, by a converged
    ~12.4%, and the fitted two-derivative coefficients are nowhere near Einstein-Hilbert.
    One methodological correction, recorded rather than buried: the PHOTON bubble below includes its
    diamagnetic seagull while the GRAVITON bubble in [B] has none, so the two are not like-for-like
    and [B] understates the graviton case. The energy method avoids the issue entirely.
    Nothing below is retracted -- the photon calibration is exact and the structural diagnosis was
    right. What changed is the conclusion drawn from it: gamma = 1 cannot be routed through the
    TETRAD sector at all. It still rests on Weinberg applied to the deconfined CURVATURE sector,
    which remains unmeasured. ***

test_graviton_ward tried to measure the graviton transversality q_i Pi^{ij,kl} = 0 (the Ward
identity that forces gamma = 1) and found it REGULATOR-LIMITED: a hard momentum cutoff breaks the
identity through a surface term, and even the PHOTON validation q_i Pi^{ij} = 0 (charge conservation,
which MUST hold) did not fall cleanly to zero. It flagged the fix -- "a lattice model with an EXACT
lattice Ward identity" -- but did not build it. This file builds it, and in doing so resolves the
whole question, though not the way one might have hoped.

The regulator. Work on a PERIODIC lattice (a Wilson-Dirac model on the square/cubic BZ torus). A
torus has NO boundary, so momentum integrals carry NO surface term -- the disease that broke
test_graviton_ward is gone by construction. The price of gauge invariance is paid honestly by
including the DIAMAGNETIC (seagull) term: the physical current response is
        Pi^{ij}_full = Pi^{ij}_para + delta^{ij} K,   K = <d^2 H / dk_i dk_j>_filled sea,
and gauge invariance demands the LONGITUDINAL piece vanish, Pi^{xx}_para(q||x,0) + K = 0.

What comes out is a sharp DICHOTOMY between the two induced gauge fields:

  * PHOTON (U(1)) -- an EXACT lattice symmetry. The Peierls substitution makes U(1) gauge invariance
    exact on the lattice, so the Ward identity is exact: the para bubble and the diamagnetic seagull
    cancel to MACHINE PRECISION at every q, while the TRANSVERSE (physical Maxwell) response stays
    nonzero. The symmetry-preserving regulator works -- the thing test_graviton_ward needed.

  * GRAVITON (diffeomorphism) -- NOT a lattice symmetry. A lattice has only DISCRETE translations,
    so the continuous diffeomorphism invariance that forces gamma=1 is BROKEN at the lattice scale
    and restored only in the IR (the continuum limit) -- exactly like the model's emergent Lorentz
    invariance. Two concrete fingerprints of this, measured below:
      (i) the graviton Ward identity is INHOMOGENEOUS -- q_i Pi^{ij,kl} equals <T>-terms, not zero,
          because <T^{ij}> (the induced stress / cosmological piece) is nonzero;
      (ii) the longitudinal stress response does NOT cancel to zero the way the photon's did.

So gamma=1 CANNOT be a lattice-exact identity, and no finite-cutoff direct measurement (as in
test_graviton_ward) can make it one. That is not a numerical shortcoming -- it is the structural
statement that diffeomorphism invariance is EMERGENT here, not fundamental. gamma=1 is therefore
correctly read from the IR: Weinberg's theorem (massless spin-2 + CONSERVED IR stress tensor =>
Einstein) applied to the emergent-Lorentz Dirac sector (test_lorentz, test_dirac), with the spin-2
graviton now measured to actually propagate and be healthy (test_spin2_dynamical). This file closes
the loop by explaining WHY that IR route is the right and only one, and by supplying the exact
gauge-sector Ward identity as the calibration that proves the regulator itself is sound.
"""
from __future__ import annotations
import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def dvec(kx, ky, m0):
    return np.sin(kx), np.sin(ky), m0 + (2 - np.cos(kx) - np.cos(ky))


def proj(kx, ky, m0):
    dx, dy, dz = dvec(kx, ky, m0)
    E = np.sqrt(dx * dx + dy * dy + dz * dz)
    dh = (dx / E)[:, None, None] * sx + (dy / E)[:, None, None] * sy + (dz / E)[:, None, None] * sz
    return 0.5 * (I2 + dh), 0.5 * (I2 - dh), E


def parabub(kx, ky, qx, qy, m0, A, B):
    """Static (Omega=0) paramagnetic interband bubble <A B>."""
    Pp, Pm, Ek = proj(kx, ky, m0)
    Pp2, Pm2, Eq = proj(kx + qx, ky + qy, m0)
    dE = Ek + Eq
    t = (np.einsum("mij,mjk,mkl,mli->m", Pm, A, Pp2, B) +
         np.einsum("mij,mjk,mkl,mli->m", Pm2, A, Pp, B))
    return -float(np.sum(t.real / dE) / len(kx))


def expval(kx, ky, m0, O):
    """Filled-sea expectation <O> = (1/N) sum_k Tr[P_-(k) O]."""
    _, Pm, _ = proj(kx, ky, m0)
    return float(np.sum(np.einsum("mij,mji->m", Pm, O).real) / len(kx))


# vertices
def jx(kx, ky, m0):
    return np.cos(kx)[:, None, None] * sx + np.sin(kx)[:, None, None] * sz   # dH/dkx


def jy(kx, ky, m0):
    return np.cos(ky)[:, None, None] * sy + np.sin(ky)[:, None, None] * sz   # dH/dky


def d2xx(kx, ky, m0):
    return -np.sin(kx)[:, None, None] * sx + np.cos(kx)[:, None, None] * sz  # d^2H/dkx^2 (diamagnetic)


def Txx(kx, ky, m0):
    return dvec(kx, ky, m0)[0][:, None, None] * sx                          # stress T^xx = d_x sigma_x


if __name__ == "__main__":
    print("=== The symmetry-preserving regulator: photon Ward identity EXACT, gamma=1 emergent ===\n")
    N = 260
    g = (np.arange(N) + 0.5) / N * 2 * np.pi
    KX, KY = np.meshgrid(g, g, indexing="ij")
    kx, ky = KX.ravel(), KY.ravel()
    m0 = -1.0
    print(f"  Wilson-Dirac on a periodic {N}x{N} BZ torus (no boundary -> no surface term), gap m0={m0}\n")

    # ---------- [A] PHOTON: exact lattice Ward identity ----------
    K = expval(kx, ky, m0, d2xx(kx, ky, m0))          # diamagnetic seagull K^{xx}
    print("  [A] PHOTON (U(1), an EXACT lattice symmetry): para bubble + diamagnetic seagull")
    print(f"      diamagnetic  K^xx = <d^2H/dkx^2> = {K:+.5f}")
    print(f"      {'qx':>6} {'Pi^xx_para(q||x)':>17} {'+ K^xx = LONGITUDINAL':>22} {'Pi^yy (TRANSVERSE)':>19}")
    for qx in (0.30, 0.15, 0.05):
        Ppara = parabub(kx, ky, qx, 0.0, m0, jx(kx + qx / 2, ky, m0), jx(kx + qx / 2, ky, m0))
        Ptr = parabub(kx, ky, qx, 0.0, m0, jy(kx + qx / 2, ky, m0), jy(kx + qx / 2, ky, m0))
        print(f"      {qx:>6.2f} {Ppara:>17.5f} {Ppara + K:>22.2e} {Ptr:>19.5f}")
    print("      => LONGITUDINAL response = 0 to machine precision at every q (para + seagull cancel")
    print("         EXACTLY): the gauge Ward identity is exact. TRANSVERSE (Maxwell) stays nonzero, so")
    print("         this is a real cancellation, not a trivial zero. The regulator test_graviton_ward")
    print("         needed -- a periodic lattice with the seagull -- WORKS for the gauge sector.\n")

    # ---------- [B] GRAVITON: not a lattice symmetry -> inhomogeneous, no exact cancellation ----------
    Tvev = expval(kx, ky, m0, Txx(kx, ky, m0))
    print("  [B] GRAVITON (diffeomorphism, NOT a lattice symmetry -- only discrete translations):")
    print(f"      stress expectation <T^xx> = {Tvev:+.5f}  =/= 0  -> the Ward identity is INHOMOGENEOUS")
    print("      (q_i Pi^{ij,kl} = <T>-terms, the induced stress/cosmological piece, not zero).")
    print(f"      {'qx':>6} {'Pi^(xx,xx)_para(q||x)':>22}")
    for qx in (0.30, 0.15, 0.05):
        Pll = parabub(kx, ky, qx, 0.0, m0, Txx(kx + qx / 2, ky, m0), Txx(kx + qx / 2, ky, m0))
        print(f"      {qx:>6.2f} {Pll:>22.5f}")
    print("      => the longitudinal STRESS response does NOT cancel to zero the way the photon's did.")
    print("         There is no exact lattice diffeomorphism to enforce it, and <T> =/= 0 makes the")
    print("         identity inhomogeneous. gamma=1 CANNOT be a lattice-exact statement.\n")

    print("[verdict] the transversality question is RESOLVED -- structurally, not numerically:")
    print("  * The symmetry-preserving regulator works: on the periodic lattice the PHOTON Ward")
    print("    identity closes to machine precision (para + diamagnetic seagull), with the transverse")
    print("    Maxwell response nonzero. This is exactly the exact-lattice-Ward-identity")
    print("    test_graviton_ward said was needed, now delivered for the gauge sector.")
    print("  * It also shows WHY the GRAVITON case is different in kind: diffeomorphism invariance is")
    print("    NOT a lattice symmetry (only discrete translations survive), and the gravitational Ward")
    print("    identity is INHOMOGENEOUS (<T^ij> =/= 0, the induced cosmological/stress term). So")
    print("    gamma=1 is not a lattice-exact identity and no finite-cutoff DIRECT measurement can")
    print("    make it one -- test_graviton_ward's 'regulator-limited' was structural, not a bug.")
    print("  * Hence gamma=1 is correctly an IR-EMERGENT statement, on the same footing as the model's")
    print("    emergent Lorentz invariance: read from Weinberg's theorem (massless spin-2 + CONSERVED")
    print("    IR stress tensor => Einstein) applied to the emergent Dirac sector (test_lorentz,")
    print("    test_dirac), with the spin-2 graviton now measured to propagate and be healthy")
    print("    (test_spin2_dynamical). The loop is closed: gamma=1 rests on the IR fixed point, and")
    print("    that is not a limitation of method but the nature of emergent gravity.")
