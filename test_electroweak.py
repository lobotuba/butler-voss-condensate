"""
The electroweak Higgs mechanism: the medium's condensate breaks SU(2)xU(1) to a single massless photon.

test_yang_mills showed the fermion loop induces genuine non-Abelian Yang-Mills (SU(2), SU(3)); the
screening arc (test_screening_gauged) showed a gauged U(1) condensate gives its gauge boson a Meissner
mass (the Abelian Higgs mechanism, measured lambda_London ~ 1/e). The Standard Model's electroweak sector
is the non-Abelian marriage of the two: a Higgs DOUBLET condensate breaks SU(2)_L x U(1)_Y down to a
single U(1)_EM, giving the three weak bosons (W+, W-, Z) a mass while leaving the photon exactly massless.
This file shows that breaking works in the model's own terms -- the Higgs is the medium's condensate
amplitude, the gauge fields are the induced ones -- and reproduces the defining electroweak structure: a
mass spectrum with EXACTLY one massless boson, and the W/Z mass ratio fixed by the weak mixing angle.

The mechanism. The condensate is a doublet Phi in the fundamental 2 of SU(2)_L with hypercharge Y = 1/2.
Its equilibrium amplitude is a vacuum expectation value, which the model can put at <Phi> = (0, v/sqrt2).
The induced gauge kinetic term for the condensate is the covariant derivative |D_mu Phi|^2 with
        D_mu Phi = (d_mu - i g (tau^a/2) W^a_mu - i g' Y B_mu) Phi,
so a constant gauge field costs |(g tau^a/2 W^a + g' Y B) <Phi>|^2 -- a MASS for whichever gauge
directions move the vacuum. The four gauge fields (W1,W2,W3,B) mix into a 4x4 mass matrix; its spectrum
is the electroweak prediction.

  [A] the mass spectrum: diagonalise the induced mass matrix -> {0, m_W, m_W, m_Z}. Exactly one massless
      boson (the photon) and three massive (W+, W-, Z), for any couplings -- the signature of the
      SU(2)xU(1) -> U(1) breaking pattern, not two massless or none.
  [B] the weak mixing angle: m_W/m_Z = cos(theta_W) = g/sqrt(g^2+g'^2), verified across g'/g.
  [C] WHY exactly one photon: the unbroken generator is the electric charge Q = T3 + Y, and the vacuum is
      Q-neutral (Q<Phi> = 0), so its gauge field stays massless; the other three combinations move the
      vacuum and gain mass. The single surviving photon is forced by the doublet's rep and its Y = 1/2.
"""
from __future__ import annotations
import numpy as np

TX = np.array([[0, 1], [1, 0]], dtype=complex)
TY = np.array([[0, -1j], [1j, 0]], dtype=complex)
TZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def mass_matrix(g, gp, v=1.0, Y=0.5):
    """4x4 gauge-boson mass matrix in the basis (W1, W2, W3, B) from |sum_A V^A Q^A <Phi>|^2, with the
    charge matrices Q^a = g tau^a/2 (a=1,2,3) and Q^B = g' Y I acting on the doublet vacuum (0, v/sqrt2)."""
    vac = np.array([0.0, v / np.sqrt(2)], dtype=complex)
    Q = [g * TX / 2, g * TY / 2, g * TZ / 2, gp * Y * I2]         # generators coupling to (W1,W2,W3,B)
    M = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            M[a, b] = (vac.conj() @ Q[a] @ Q[b] @ vac).real       # <Phi| Q^a Q^b |Phi>
    # |D Phi|^2 = sum_AB V^A V^B <Phi|Q^A Q^B|Phi>; the physical mass^2 (coeff of (1/2) V^A M^2 V^B) is 2x:
    return (M + M.T)                                              # symmetric physical mass^2 matrix


if __name__ == "__main__":
    print("=== The electroweak Higgs mechanism: SU(2)xU(1) -> U(1)_EM in the model ===\n")
    print("  The medium's condensate is a doublet Phi (fundamental of SU(2)_L, hypercharge Y=1/2). Its")
    print("  VEV breaks the induced SU(2)xU(1); the covariant derivative |D Phi|^2 gives the mass matrix.\n")

    g, gp, v = 0.65, 0.35, 1.0                                    # SU(2) and U(1)_Y couplings, VEV
    # ---------- [A] the mass spectrum ----------
    M2 = mass_matrix(g, gp, v)
    m = np.sqrt(np.clip(np.linalg.eigvalsh(M2), 0, None))
    print("  [A] GAUGE-BOSON MASS SPECTRUM (diagonalising the 4x4 induced mass matrix):")
    print(f"      couplings g = {g}, g' = {gp}, VEV v = {v}")
    print(f"      masses (sorted) = [{', '.join(f'{x:.4f}' for x in m)}]")
    mW_pred, mZ_pred = g * v / 2, np.sqrt(g * g + gp * gp) * v / 2
    print(f"      => spectrum = {{0 (photon), m_W, m_W, m_Z}}:  one MASSLESS + three MASSIVE.")
    print(f"         massless boson (photon): {m[0]:.2e}  (machine zero)")
    print(f"         W+, W- (degenerate):     {m[1]:.4f}, {m[2]:.4f}   vs  g v/2 = {mW_pred:.4f}")
    print(f"         Z:                       {m[3]:.4f}   vs  sqrt(g^2+g'^2) v/2 = {mZ_pred:.4f}")
    print("         Exactly one massless boson -- the SU(2)xU(1) -> U(1) signature, not zero or two.\n")

    # ---------- [B] the weak mixing angle ----------
    print("  [B] THE WEAK MIXING ANGLE: m_W/m_Z = cos(theta_W) = g/sqrt(g^2+g'^2):")
    print(f"      {'g_prime/g':>12} {'m_W/m_Z (measured)':>20} {'cos(theta_W)':>14} {'sin^2(theta_W)':>16}")
    for gpr in (0.0, 0.3, 0.5364, 0.8):                          # 0.5364 ~ physical sin^2 ~ 0.223
        M2r = mass_matrix(g, gpr * g, v)
        mr = np.sqrt(np.clip(np.linalg.eigvalsh(M2r), 0, None))
        ratio = mr[1] / mr[3] if mr[3] > 0 else float("nan")
        cw = 1.0 / np.sqrt(1 + gpr ** 2)
        print(f"      {gpr:>12.4f} {ratio:>20.5f} {cw:>14.5f} {1-cw*cw:>16.5f}")
    print("      => the measured W/Z mass ratio is exactly cos(theta_W): the weak mixing angle is not a")
    print("         separate input but the geometry of the W3-B mixing that the condensate induces.\n")

    # ---------- [C] why exactly one photon ----------
    print("  [C] WHY EXACTLY ONE PHOTON: the unbroken generator is Q = T3 + Y (electric charge).")
    vac = np.array([0.0, v / np.sqrt(2)], dtype=complex)
    Q_charge = TZ / 2 + 0.5 * I2                                  # Q = T3 + Y on the doublet
    print(f"      electric charge of the vacuum:  Q<Phi> = {np.linalg.norm(Q_charge @ vac):.2e}  (zero)")
    print("      the vacuum is Q-neutral, so the photon (the gauge field of Q) leaves it invariant and")
    print("      stays massless; the three orthogonal generators move the vacuum and acquire mass.")
    for name, T in (("T1", TX / 2), ("T2", TY / 2), ("T3-Y (Z-like)", TZ / 2 - 0.5 * I2)):
        print(f"        {name:>14} moves the vacuum by |T<Phi>| = {np.linalg.norm(T @ vac):.3f}  -> massive")
    print("      => one massless photon is forced by the doublet's representation and its Y = 1/2; a")
    print("         different hypercharge would leave no massless U(1), or a different unbroken one.\n")

    print("[verdict] the electroweak Higgs mechanism works in the model, with the correct structure:")
    print("  * The medium's condensate, as a doublet with hypercharge 1/2, breaks the induced SU(2)xU(1)")
    print("    to a single U(1): the gauge spectrum is {0, m_W, m_W, m_Z} -- EXACTLY one massless photon")
    print("    and three massive weak bosons, for any couplings [A]. This is the same condensate that")
    print("    gaps the amplitude mode for gravity and self-tunes the cosmological constant, now giving")
    print("    the weak bosons their mass -- and the same Meissner mechanism measured for the Abelian")
    print("    gauge field in the screening arc (test_screening_gauged), here in its non-Abelian form.")
    print("  * The W/Z mass ratio is the cosine of the weak mixing angle [B], and the surviving massless")
    print("    photon is the electric-charge direction Q = T3 + Y, kept massless because the vacuum is")
    print("    Q-neutral [C]. The breaking PATTERN and the mass relations are reproduced, not fitted.")
    print("  * HONEST scope -- as in test_yang_mills, this is the MECHANISM, not a derivation of the")
    print("    Standard Model. The gauge group SU(2)xU(1), the doublet representation, and the hypercharge")
    print("    Y = 1/2 are inputs; what is shown is that the condensate breaks them the way the Standard")
    print("    Model does, with one massless photon and the W/Z/theta_W relations. WHY the condensate sits")
    print("    in a Y = 1/2 doublet -- the assignment that leaves electromagnetism unbroken and is fixed")
    print("    with the fermion charges by anomaly cancellation -- is the next question, not settled here.")
