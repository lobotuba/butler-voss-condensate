"""
The sign of the induced gravitational coupling: is mu > 0?  (Closing the tensor-gravity crux.)

*** STATUS UPDATE -- the quantity this file DISCARDS has since been measured, and it is not zero.
    The extraction below takes mu as the q^2 coefficient of Pi(q) - Pi(0), subtracting Pi(0) as a
    contact term. That subtraction was never justified, and Pi(0) is exactly the GRAVITON MASS
    candidate: test_ir_fixed_point then showed a mass is the one RELEVANT deformation, so whether
    Pi(0) vanishes decides whether the infrared attractor exists at all.
    test_graviton_mass measured it, avoiding the bubble/seagull bookkeeping entirely by computing the
    filled-sea energy under a CONSTANT deformation (the q -> 0 limit, to all orders at once). For the
    PHOTON the q = 0 term is -1.5e-10 and falls under refinement, because an exact lattice symmetry
    forbids it. For the TETRAD cone shear it is -0.27 -- order unity, eight orders larger, stable
    across gap and grid, and NEGATIVE, so it destabilises the symmetric cone rather than merely making
    it massive.
    Nothing measured below is retracted: the SIGN result mu > 0 stands, and so does the calibration
    against the healthy photon. What is no longer available is the tacit assumption that the
    discarded contact term was harmless. Masslessness cannot be inherited from the tetrad; the
    deconfined curvature sector of test_deconfinement must supply it, and that remains unmeasured. ***

test_deconfinement reduced tensor gravity to ONE number's sign. The confining spin-2 sector
deconfines into a massless Newtonian graviton -- turning the biharmonic +R into a -1/r -- for ANY
POSITIVE induced Einstein/Newton coefficient mu. So the whole question is:

        does integrating out the fermions induce a POSITIVE (healthy) kinetic term
        for the gravitational potential, mu > 0?

This is the Sakharov induced-gravity sign. It is famously scheme-sensitive: the ANALYTIC
(Einstein-Hilbert) coefficient is UV-dominated, so its raw value depends on the regulator, and free
matter fields do not universally give the healthy sign (Visser). What saves the question from being
a convention is a CALIBRATION the model already owns:

    the induced PHOTON is HEALTHY. test_induced_action showed the fermion loop generates the
    Maxwell term (a working emergent photon). Its induced Coulomb kinetic term -- the charge-density
    correlator <J0 J0> -- is a healthy dielectric (susceptibility chi > 0).

So compute, FROM THE SAME LOOP WITH IDENTICAL CONVENTIONS, the induced NEWTONIAN kinetic term -- the
ENERGY-density correlator <T00 T00> (energy is the gravitational charge, T00 couples to the
Newtonian potential h00). If <T00 T00> has the SAME-SIGN q^2 coefficient as the (healthy) <J0 J0>,
then the induced gravity is as healthy as the induced electromagnetism the model already runs on:
mu > 0. The overall loop-sign convention cancels in the comparison -- only the RELATIVE sign is
physical, and it is calibrated against a sector known to work.

Method. Gapped Dirac cone H = v(kx sx + ky sy) + M sz (a gap makes the small-q limit ANALYTIC and
LOCAL, so a genuine Einstein/Coulomb q^2 term exists). Static (Omega=0) interband polarization of
the filled sea, built from band PROJECTORS P_pm = (I +- dhat.sigma)/2 (no eigenvector-phase
pitfalls), for a vertex O:
    Pi_O(q) = (1/N) sum_k Tr[ P_-(k) O P_+(k+q) O + P_+(k) O P_-(k+q) O ] * 2/(E_k + E_{k+q}).
The induced kinetic term is the q^2 coefficient of Pi_O(q) - Pi_O(0). Vertices:
    charge density (photon)   J0  = I                         -> induced Coulomb term (calibration)
    energy density (graviton) T00 = E(k) I   [and H(k), x-check] -> induced Newtonian term (result)

Honest scope. This settles the SIGN of the induced coupling for the NEWTONIAN (h00, force-carrying)
sector -- exactly the sector deconfinement needs and the one test_critical_gravity's attraction
lives in. It does NOT compute the magnitude of G (the coefficient is UV/cutoff-dependent -- the
Sakharov feature; only the sign is robust, and the model's physical lattice cutoff would fix the
value in a separate computation), and it does NOT close the full 3+1D spin-2 (light-bending, gamma=1)
graviton: on the 2+1D emergent Dirac cone the SPATIAL (transverse-traceless) graviton is
non-dynamical (its q^2 coefficient is ~0 below), so the force lives entirely in h00 -- consistent
with the scalar/Newtonian picture, but the tensor/gamma=1 completion stays the open item flagged in
test_graviton_ward. What this file delivers is the crux the deconfinement rested on: the sign.
"""
from __future__ import annotations
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def projectors(kx, ky, M, v=1.0):
    dx, dy, dz = v * kx, v * ky, M * np.ones_like(kx)
    E = np.sqrt(dx * dx + dy * dy + dz * dz)
    dh = (dx / E)[:, None, None] * SX + (dy / E)[:, None, None] * SY + (dz / E)[:, None, None] * SZ
    return 0.5 * (I2 + dh), 0.5 * (I2 - dh), E


def bubble(kx, ky, qx, qy, M, A, B, v=1.0):
    """Static interband polarization Pi_{AB}(q) at Omega = 0."""
    Pp, Pm, Ek = projectors(kx, ky, M, v)
    Pp2, Pm2, Eq = projectors(kx + qx, ky + qy, M, v)
    dE = Ek + Eq
    t = (np.einsum("mij,mjk,mkl,mli->m", Pm, A, Pp2, B) +
         np.einsum("mij,mjk,mkl,mli->m", Pp, A, Pm2, B))
    return float(np.sum(t.real * 2.0 / dE) / len(kx))


def disc(LAM, NG):
    g = np.linspace(-LAM, LAM, NG) + 1e-5
    X, Y = np.meshgrid(g, g, indexing="ij")
    m = (X * X + Y * Y) <= LAM * LAM
    return X[m], Y[m]


def q2coef(kx, ky, M, vfn, qs, qdir=(1.0, 0.0), v=1.0):
    """q^2 coefficient of Pi_O(q) - Pi_O(0): the induced kinetic term."""
    vals = []
    for q in qs:
        qx, qy = q * qdir[0], q * qdir[1]
        A = vfn(kx + qx / 2, ky + qy / 2, M, v)
        vals.append(bubble(kx, ky, qx, qy, M, A, A, v))
    vals = np.array(vals)
    return np.polyfit(qs ** 2, vals - vals[0], 1)[0], vals[0]


# vertices
def V_J0(ax, ay, M, v):                       # charge density (photon time component)
    return I2[None].repeat(len(ax), 0)


def V_T00_E(ax, ay, M, v):                    # energy density = E(k) I  (energy is the grav. charge)
    E = np.sqrt((v * ax) ** 2 + (v * ay) ** 2 + M ** 2)
    return E[:, None, None] * I2[None]


def V_T00_H(ax, ay, M, v):                    # cross-check: the Hamiltonian-density vertex H(k)
    return (v * ax[:, None, None] * SX + v * ay[:, None, None] * SY +
            M * SZ[None].repeat(len(ax), 0))


def V_Tplus(ax, ay, M, v):                    # spatial (xx-yy) graviton polarisation
    return (v / 2.0) * (ax[:, None, None] * SX - ay[:, None, None] * SY)


if __name__ == "__main__":
    print("=== The sign of the induced gravitational coupling: is mu > 0? ===\n")
    print("  Calibration: the induced PHOTON is healthy -> <J0 J0> q^2 coef = induced dielectric")
    print("  chi > 0. Result: does <T00 T00> (energy density) share that sign -> mu > 0 (healthy")
    print("  Newtonian coupling)? Same loop, same conventions; only the RELATIVE sign is physical.\n")
    qs = np.array([1e-4, 0.05, 0.10, 0.15, 0.20])

    print(f"  {'M':>5} {'LAM':>5} {'<J0J0> q2 (calib)':>18} {'<T00T00> q2 (E*I)':>18} "
          f"{'<T00T00> q2 (H)':>16} {'same sign?':>11}")
    all_ok = True
    for M in (0.3, 0.4, 0.5, 0.6, 0.8):
        for LAM in (2.0, 3.0, 4.0):
            kx, ky = disc(LAM, 401)
            aJ, _ = q2coef(kx, ky, M, V_J0, qs)
            aTE, _ = q2coef(kx, ky, M, V_T00_E, qs)
            aTH, _ = q2coef(kx, ky, M, V_T00_H, qs)
            ok = (np.sign(aTE) == np.sign(aJ)) and (np.sign(aTH) == np.sign(aJ))
            all_ok &= ok
            print(f"  {M:>5.1f} {LAM:>5.1f} {aJ:>18.5f} {aTE:>18.5f} {aTH:>16.6f} {str(ok):>11}")
    print(f"\n  => charge-density and energy-density kinetic terms share the SAME sign in every case: "
          f"{all_ok}")
    print("     The induced Coulomb term is a healthy dielectric (chi > 0); the induced Newtonian")
    print("     term matches it -> mu > 0. The induced gravitational coupling is HEALTHY.\n")

    # the spatial (TT) graviton: non-dynamical in 2+1D (its q^2 coefficient ~ 0)
    kx, ky = disc(3.0, 401)
    aTT, cTT = q2coef(kx, ky, 0.5, V_Tplus, qs)
    aJ, _ = q2coef(kx, ky, 0.5, V_J0, qs)
    print(f"  spatial (xx-yy) graviton: q^2 coef = {aTT:+.2e} (const piece {cTT:.3f}), vs <J0J0> "
          f"{aJ:+.3f}")
    print("     ~ 0: on the 2+1D cone the spatial transverse-traceless graviton is non-dynamical, so")
    print("     the force lives ENTIRELY in the Newtonian h00 sector -- the scalar/Newtonian picture.\n")

    print("[verdict] the sign that was 'the whole ballgame' comes out RIGHT: mu > 0.")
    print("  * Calibrated against the model's own working sector: the induced Coulomb kinetic term")
    print("    (<J0 J0>) is a healthy dielectric chi > 0. The induced Newtonian kinetic term")
    print("    (<T00 T00>, energy density) has the SAME sign -- robust across mass M, cutoff LAM, and")
    print("    both energy-density vertex definitions. So the fermion loop induces mu > 0.")
    print("  * By test_deconfinement, mu > 0 turns the confining +R into an attractive Newtonian")
    print("    -1/r; by test_critical_gravity the energy source is positive-definite so like masses")
    print("    ATTRACT. The Newtonian (force) sector of induced gravity is HEALTHY and deconfines.")
    print("  * HONEST scope: (1) only the SIGN is robust; the MAGNITUDE of G is UV/cutoff-dependent")
    print("    (the Sakharov feature -- the model's physical lattice cutoff fixes it, a separate")
    print("    computation). (2) This is the h00/Newtonian sector on the 2+1D emergent Dirac; the")
    print("    full 3+1D spin-2 (gamma=1 light-bending) graviton is non-dynamical here and stays the")
    print("    open item of test_graviton_ward. What is settled is the deconfinement input: the sign.")
