"""
Does mass source curvature at all? The smooth induced coupling that would give gamma = 1, measured.

test_einstein_source reduced the light-bending parameter to one number: gamma = kappa/(4 pi G), where
4 pi G is the Newtonian coupling with which a static mass sources the time potential Phi, and kappa is
the strength with which the SAME mass sources the spatial curvature Psi (the incompatibility). It
established the two halves of the ratio separately -- that pure elasticity gives kappa = 0 (a mass
relaxes to a compatible DISPLACEMENT, whose curvature vanishes to 1e-15), so any curvature coupling
must be INDUCED by the fermion loop -- and left the induced kappa as the one calculation deciding
gamma. This file measures it, in the smooth (non-topological) channel.

The setup is forced by what a static point mass IS: a localized energy density, so it sources gravity
only through its T^{00} component. Whatever spatial curvature it produces must come from the induced
coupling of that energy density to the spatial stress T_{ij} -- the mixed response Pi^{00,ij}(q).
The 2D linearized Ricci scalar of the resulting metric is R(q) = q_i q_j h_{ij} - q^2 h_{kk}, so a
mass curves space only insofar as it sources the spatial stress. gamma = 1 (Einstein) requires that
coupling to be present at the strength that matches the Newtonian one; gamma = 0 (Nordstrom / scalar)
is the value if the energy density sources no spatial stress at all.

The instrument is the static interband bubble of test_induced_sign, on the gapped 2+1D Dirac cone,
calibrated the same way against the healthy induced photon <J0 J0>. The energy-density vertex is
T^{00} = E(k) I (energy is the gravitational charge, as there), and the stress vertex is the
momentum flux T_{ij} = (1/2)(k_i sigma_j + k_j sigma_i).

  [A] CALIBRATION. <J0 J0> is a healthy dielectric (positive q^2 term) and <T00 T00> is the nonzero
      Newtonian coupling that gives Phi -- the denominator of gamma is real and in hand.
  [B] THE SELECTION RULE. <T00, T_ij> = 0 exactly, for every component including the trace, at all q
      and across mass and cutoff. The energy-density vertex is proportional to the identity -- a
      SCALAR -- and a scalar cannot source the spin-2 spatial stress. This is not a small number to
      be resolved; it is a symmetry zero.
  [C] THE CONTRAST that proves it is physical, not a null instrument: a genuine spin-2 source
      (the xx-yy graviton polarization) couples to the same spatial stress at O(1). The bubble is
      perfectly capable of a nonzero answer; the scalar mass simply does not elicit one.
  [D] ROBUSTNESS across mass, cutoff and momentum.

SCOPE, stated plainly because the result is a negative and negatives are easy to overclaim. This
measures the SMOOTH induced coupling -- the mass to the smooth spatial metric h_{ij}. It shows the
smooth route to sourcing curvature is closed, which is consistent with and independent of three other
results: the induced spatial graviton is non-dynamical (test_induced_sign), the induced tetrad action
is not Einstein-Hilbert (test_graviton_transversality), and it has no gauge null space to project into
(test_graviton_nullspace). What it does NOT touch is the TOPOLOGICAL channel: genuine curvature in
this medium is carried by INCOMPATIBLE strain -- disclination density (test_light_bending,
test_disclination_force) -- which a smooth stress-stress bubble cannot see. Whether a mass sources
curvature through that topological channel is the remaining, and harder, curvature-sector question.
So this file closes the smooth route to gamma = 1 and sharpens what the topological measurement must
overcome; it does not by itself settle gamma for the deconfined curvature sector.
"""
from __future__ import annotations
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [SX, SY]


def proj(kx, ky, M, v=1.0):
    dx, dy, dz = v * kx, v * ky, M * np.ones_like(kx)
    E = np.sqrt(dx * dx + dy * dy + dz * dz)
    dh = (dx / E)[:, None, None] * SX + (dy / E)[:, None, None] * SY + (dz / E)[:, None, None] * SZ
    return 0.5 * (I2 + dh), 0.5 * (I2 - dh), E


def bubble(kx, ky, qx, qy, M, A, B, v=1.0):
    """Static (Omega=0) interband polarization Pi_{AB}(q), built from band projectors."""
    Pp, Pm, Ek = proj(kx, ky, M, v)
    Pp2, Pm2, Eq = proj(kx + qx, ky + qy, M, v)
    dE = Ek + Eq
    t = (np.einsum("mij,mjk,mkl,mli->m", Pm, A, Pp2, B) +
         np.einsum("mij,mjk,mkl,mli->m", Pp, A, Pm2, B))
    return float(np.sum(t.real * 2.0 / dE) / len(kx))


def disc(LAM, NG):
    g = np.linspace(-LAM, LAM, NG) + 1e-5
    X, Y = np.meshgrid(g, g, indexing="ij")
    m = (X * X + Y * Y) <= LAM * LAM
    return X[m], Y[m]


def V_J0(ax, ay, M, v):                                  # charge density (photon calibration)
    return I2[None].repeat(len(ax), 0)


def V_T00(ax, ay, M, v):                                 # energy density (the gravitational charge)
    E = np.sqrt((v * ax) ** 2 + (v * ay) ** 2 + M ** 2)
    return E[:, None, None] * I2[None]


def V_T(i, j, ax, ay, M, v):                             # momentum flux T_{ij}
    ki = ax if i == 0 else ay
    kj = ax if j == 0 else ay
    return 0.5 * v * (ki[:, None, None] * SIG[j] + kj[:, None, None] * SIG[i])


def V_Tplus(ax, ay, M, v):                               # a genuine spin-2 (xx-yy) source, for contrast
    return 0.5 * v * (ax[:, None, None] * SX - ay[:, None, None] * SY)


def q2coef(kx, ky, M, vfn, qs, v=1.0):
    """q^2 coefficient and q->0 value of Pi_O(q) along q = (q,0)."""
    vals = []
    for q in qs:
        A = vfn(kx + q / 2, ky, M, v)
        vals.append(bubble(kx, ky, q, 0.0, M, A, A, v))
    vals = np.array(vals)
    return np.polyfit(qs ** 2, vals - vals[0], 1)[0], vals[0]


if __name__ == "__main__":
    print("=== Does mass source curvature? The smooth induced coupling behind gamma = 1 ===\n")
    print("  gamma = kappa/(4 pi G): 4 pi G is how a mass sources the time potential Phi (<T00 T00>),")
    print("  kappa is how the SAME mass sources spatial curvature (<T00, T_ij>). A static mass has")
    print("  only T00, so all its curvature must come from the induced T00 -> T_ij coupling.\n")

    qs = np.array([0.05, 0.10, 0.15, 0.20])

    # ---------- [A] calibration: the denominator of gamma is real ----------
    kx, ky = disc(3.0, 501)
    aJ, _ = q2coef(kx, ky, 0.5, V_J0, qs)
    aT, cT = q2coef(kx, ky, 0.5, V_T00, qs)
    print("  [A] CALIBRATION (M = 0.5, cutoff 3.0):")
    print(f"      <J0 J0>   q^2 coef = {aJ:+.5f}   (healthy induced dielectric, > 0)")
    print(f"      <T00 T00> q^2 coef = {aT:+.5f}   (the Newtonian coupling ~ mu that gives Phi, =/= 0)")
    print("      => the time-potential side of gamma exists and is healthy. The question is the")
    print("         spatial side: does the same energy density source any curvature?\n")

    # ---------- [B] the selection rule: mass -> spatial stress vanishes ----------
    print("  [B] THE INDUCED COUPLING OF A MASS TO THE SPATIAL STRESS, <T00, T_ij>(q):")
    print(f"      {'q':>7} {'<T00 Txx>':>11} {'<T00 Tyy>':>11} {'<T00 Txy>':>11} {'<T00 Tkk>':>11}")
    worst = 0.0
    for q in (1e-3, 0.05, 0.10, 0.20):
        mx, my = kx + q / 2, ky
        T00 = V_T00(mx, my, 0.5, 1.0)
        vals = [bubble(kx, ky, q, 0.0, 0.5, T00, V_T(i, j, mx, my, 0.5, 1.0))
                for (i, j) in ((0, 0), (1, 1), (0, 1))]
        tkk = bubble(kx, ky, q, 0.0, 0.5, T00,
                     V_T(0, 0, mx, my, 0.5, 1.0) + V_T(1, 1, mx, my, 0.5, 1.0))
        worst = max(worst, max(abs(v) for v in vals), abs(tkk))
        print(f"      {q:>7.3f} {vals[0]:>11.6f} {vals[1]:>11.6f} {vals[2]:>11.6f} {tkk:>11.6f}")
    print(f"      => zero to {worst:.0e} on every component, including the trace. The energy-density")
    print("         vertex is proportional to the identity -- a SCALAR -- and a scalar sources no")
    print("         spin-2 spatial stress. A symmetry zero, not an unresolved small number.\n")

    # ---------- [C] the contrast: a real spin-2 source DOES couple ----------
    print("  [C] IS THE BUBBLE JUST BLIND HERE? Replace the scalar mass by a genuine spin-2 source")
    print("      (the xx-yy graviton polarization) and ask the same question:")
    q = 0.10
    mx, my = kx + q / 2, ky
    T00 = V_T00(mx, my, 0.5, 1.0)
    Tp = V_Tplus(mx, my, 0.5, 1.0)
    dxx_yy = V_T(0, 0, mx, my, 0.5, 1.0) - V_T(1, 1, mx, my, 0.5, 1.0)
    s_scalar = bubble(kx, ky, q, 0.0, 0.5, T00, dxx_yy)
    s_spin2 = bubble(kx, ky, q, 0.0, 0.5, Tp, dxx_yy)
    print(f"      scalar mass -> (Txx - Tyy):  {s_scalar:+.6f}")
    print(f"      spin-2 src  -> (Txx - Tyy):  {s_spin2:+.6f}")
    print("      => the spin-2 source couples at O(1); the scalar mass gives an exact zero. The")
    print("         instrument is fully capable of a nonzero answer -- the mass simply does not")
    print("         elicit one. gamma from the smooth induced coupling is 0, not 1.\n")

    # ---------- [D] robustness ----------
    print("  [D] ROBUSTNESS of the zero across mass and cutoff (q = 0.1):")
    print(f"      {'M':>5} {'LAM':>5} {'max|<T00 T_ij>|':>16} {'<T00 T00>':>11}")
    for M in (0.3, 0.5, 0.8):
        for LAM in (2.0, 4.0):
            kx2, ky2 = disc(LAM, 501)
            mx, my = kx2 + 0.05, ky2
            T00 = V_T00(mx, my, M, 1.0)
            mx_ij = max(abs(bubble(kx2, ky2, 0.1, 0.0, M, T00, V_T(i, j, mx, my, M, 1.0)))
                        for (i, j) in ((0, 0), (1, 1), (0, 1)))
            t00 = bubble(kx2, ky2, 0.1, 0.0, M, T00, T00)
            print(f"      {M:>5.1f} {LAM:>5.1f} {mx_ij:>16.1e} {t00:>11.6f}")
    print("      => the zero is exact everywhere; the Newtonian coupling is everywhere nonzero.\n")

    print("[verdict] a mass sources the Newtonian potential but no smooth curvature: gamma_smooth = 0.")
    print("  * The energy density that a static mass presents to gravity is a scalar (T00 ~ E*I), and")
    print("    its induced coupling to the spatial stress vanishes identically -- <T00, T_ij> = 0 for")
    print("    every component, the trace included, across mass, cutoff and momentum. A genuine spin-2")
    print("    source couples to the same stress at O(1), so the zero is a selection rule, not a numb")
    print("    instrument. Through the smooth loop, mass sources Phi and not Psi: gamma = 0.")
    print("  * This is the fourth independent statement that the SMOOTH route to gamma = 1 is closed,")
    print("    with the induced spatial graviton non-dynamical (test_induced_sign), the induced tetrad")
    print("    action not Einstein-Hilbert (test_graviton_transversality), and that action having no")
    print("    gauge null space to project into (test_graviton_nullspace). Pure elasticity already gave")
    print("    the same zero from the other side (test_einstein_source: a mass relaxes to a compatible")
    print("    displacement). Every SMOOTH mechanism the model has returns gamma = 0.")
    print("  * WHAT REMAINS, and it is the whole curvature sector. Real curvature in this medium is")
    print("    carried by INCOMPATIBLE strain -- disclination density (test_light_bending,")
    print("    test_disclination_force) -- not by the smooth h_ij a stress-stress bubble sees. Whether")
    print("    a mass sources curvature through that TOPOLOGICAL channel is untouched here and is the")
    print("    remaining measurement. gamma = 1 currently rests entirely on the emergent-Weinberg")
    print("    argument (a massless spin-2 on a conserved IR stress tensor is forced to be Einstein),")
    print("    which -- unlike the model's emergent Lorentz invariance -- has no direct in-model")
    print("    confirmation, and now has direct evidence against it in the smooth sector.")
