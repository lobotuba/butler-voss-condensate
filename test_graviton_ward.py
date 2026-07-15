"""
Is the induced graviton transverse? The Ward identity that forces gamma = 1.

test_einstein_source reduced the spin-2 wall to ONE question: does mass source the curvature
sector with the Einstein strength (gamma = 1)? Weinberg's answer: a massless spin-2 field coupled
to a CONSERVED stress tensor is forced to be Einstein (gamma = 1). Conservation shows up as
TRANSVERSALITY of the induced graviton self-energy:
        q_i Pi^{ij,kl}(q) = 0        (the Ward identity of linearised diffeomorphism invariance)
exactly as gauge invariance of the induced PHOTON is q_i Pi^{ij}(q) = 0 (charge conservation).

The honest way to compute this without lattice stress-tensor vertex ambiguities: work in the IR
theory the lattice FLOWS to. The project has established repeatedly (test_lorentz, test_dirac)
that the medium's low-energy sector is a Lorentz-invariant DIRAC fermion. A continuum Dirac
fermion has an exactly conserved, symmetric stress tensor, so IF its induced graviton is
transverse, then the medium -- flowing to this theory in the IR -- inherits gamma = 1. The
vertices there are unambiguous:
        current   J^i  = sigma_i                       (photon,  spin-1)
        stress    T^ij = 1/2 (k_i sigma_j + k_j sigma_i)   (graviton, spin-2)

Method: the static (omega = 0) interband polarisation of the filled Dirac sea (a cutoff Lambda
regularises the UV), built entirely from band PROJECTORS P_pm = (I +/- dhat.sigma)/2 so there are
no eigenvector-phase pitfalls:
    Pi^{AB}(q) = (1/M) sum_k [ Tr(P_-(k) A P_+(k+q) B) / (E_-(k) - E_+(k+q))
                             + Tr(P_+(k) A P_-(k+q) B) / (E_+(k) - E_-(k+q)) ]
with the vertices A, B evaluated at the symmetric midpoint k + q/2.

Checks:
  * CURRENT (validation): q_i Pi^{ij}_JJ must be ~ 0 -- emergent gauge invariance. If the method
    gets this, it is trusted.
  * STRESS (the result): q_i Pi^{ij,kl}_TT ~ 0 would mean the induced graviton is transverse =
    diffeomorphism-invariant = coupled to a conserved source = gamma = 1.

HONEST OUTCOME (recorded rather than hidden). The numerical check did NOT land cleanly. The
CURRENT validation -- which MUST be exactly transverse by charge conservation -- stays at the
few-to-tens-of-percent level and does not fall cleanly to zero; hard cutoff, soft cutoff, larger
Lambda, and Pauli-Villars subtraction all leave it erratic. This is the classic disease: a
momentum-cutoff fermion bubble breaks the Ward identity through a cutoff surface term, and a
symmetry-preserving regulator (proper dimensional regularisation, or a lattice model with an
EXACT lattice Ward identity) is needed to see transversality cleanly -- beyond this quick setup.
Because the validation is not clean, the stress ratio here is NOT a trustworthy measurement of
graviton transversality, and is not reported as one. What survives is the ANALYTIC argument,
resting on results the project HAS established: the IR theory is a Lorentz-invariant Dirac fermion
(test_lorentz, test_dirac); its stress tensor is conserved by Noether (translation + Lorentz
symmetry); a conserved stress tensor induces a transverse spin-2 (Weinberg) => gamma = 1. The
numbers below are an honest record of a regulator-limited attempt, not the proof.
"""
from __future__ import annotations
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [SX, SY]


def dvec(kx, ky, m):
    return kx, ky, m * np.ones_like(kx)


def projectors(kx, ky, m):
    dx, dy, dz = dvec(kx, ky, m)
    nrm = np.sqrt(dx * dx + dy * dy + dz * dz)
    ax, ay, az = dx / nrm, dy / nrm, dz / nrm
    dhat_sigma = (ax[:, None, None] * SX + ay[:, None, None] * SY + az[:, None, None] * SZ)
    Pp = 0.5 * (I2 + dhat_sigma)
    Pm = 0.5 * (I2 - dhat_sigma)
    return Pp, Pm, nrm


def tr4(A, B, C, D):
    return np.einsum("mij,mjk,mkl,mli->m", A, B, C, D)


def polarization(kx, ky, qx, qy, m, Amats, Bmats, Lam):
    """Return Pi^{AB} for lists of vertex matrices A (len nA) and B (len nB): shape (nA,nB)."""
    Pp_k, Pm_k, Ek = projectors(kx, ky, m)
    Pp_q, Pm_q, Eq = projectors(kx + qx, ky + qy, m)
    d1 = (-Ek) - (Eq)          # E_-(k) - E_+(k+q)  (negative)
    d2 = (Ek) - (-Eq)          # E_+(k) - E_-(k+q)  (positive)
    out = np.zeros((len(Amats), len(Bmats)), dtype=complex)
    for ia, A in enumerate(Amats):
        for ib, B in enumerate(Bmats):
            t1 = tr4(Pm_k, A, Pp_q, B) / d1
            t2 = tr4(Pp_k, A, Pm_q, B) / d2
            out[ia, ib] = (t1 + t2).sum()
    return out.real / len(kx)


def make_grid(Lam, Nk):
    ks = np.linspace(-Lam, Lam, Nk)
    KX, KY = np.meshgrid(ks, ks, indexing="ij")
    kx, ky = KX.ravel(), KY.ravel()
    inside = (kx ** 2 + ky ** 2) <= Lam ** 2      # rotationally symmetric cutoff (no lattice anisotropy)
    return kx[inside], ky[inside]


if __name__ == "__main__":
    print("=== Is the induced graviton transverse? The Ward identity forcing gamma = 1 ===\n")
    Lam, Nk, m = 4.0, 320, 0.03
    kx, ky = make_grid(Lam, Nk)
    print(f"  IR Dirac fermion, filled sea, rotational cutoff Lambda = {Lam}, {len(kx)} k-points, m = {m}\n")

    for (qx, qy) in [(0.30, 0.0), (0.15, 0.0), (0.075, 0.0)]:
        qn = np.hypot(qx, qy)
        kb_x, kb_y = kx + qx / 2, ky + qy / 2      # symmetric midpoint for the vertices

        # ---- current (photon): vertices sigma_x, sigma_y ----
        J = [SX[None].repeat(len(kx), 0), SY[None].repeat(len(kx), 0)]
        PiJ = polarization(kx, ky, qx, qy, m, J, J, Lam)      # 2x2
        longJ = np.array([qx * PiJ[0, j] + qy * PiJ[1, j] for j in range(2)])
        transJ = np.linalg.norm(longJ) / (np.linalg.norm(PiJ) + 1e-30)

        # ---- stress (graviton): vertices T^{ij} = 1/2(k_i s_j + k_j s_i), midpoint k ----
        def Tmat(i, j):
            ki = kb_x if i == 0 else kb_y
            kj = kb_x if j == 0 else kb_y
            return 0.5 * (ki[:, None, None] * SIG[j] + kj[:, None, None] * SIG[i])
        comps = [(0, 0), (1, 1), (0, 1)]           # xx, yy, xy
        T = [Tmat(i, j) for (i, j) in comps]
        PiT = polarization(kx, ky, qx, qy, m, T, T, Lam)      # 3x3 in (xx,yy,xy)

        # longitudinal: q_i Pi^{ij,kl}. assemble the 2x2x2x2 tensor from the 3 independent comps
        idx = {(0, 0): 0, (1, 1): 1, (0, 1): 2, (1, 0): 2}
        def PiT4(i, j, k, l):
            return PiT[idx[(i, j)], idx[(k, l)]]
        long_norm, full_norm = 0.0, 0.0
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    val = sum(q * PiT4(i, j, k, l) for i, q in ((0, qx), (1, qy)))
                    long_norm += val ** 2
                    full_norm += PiT4(0, j, k, l) ** 2 + PiT4(1, j, k, l) ** 2
        transT = np.sqrt(long_norm) / (np.sqrt(full_norm) + 1e-30)

        print(f"  |q| = {qn:.3f}:   current  |q.Pi|/|Pi| = {transJ:.2e}   "
              f"stress  |q.Pi|/|Pi| = {transT:.2e}")

    print("\n  NOTE: the CURRENT ratio is the validation, and it does NOT fall cleanly to zero")
    print("  (charge conservation demands it must). So this momentum-cutoff bubble breaks the")
    print("  Ward identity -- a regulator artifact -- and the STRESS ratio is therefore NOT a")
    print("  trustworthy measurement of graviton transversality. Reported honestly, not spun.\n")

    print("[verdict] the direct numerical Ward identity is REGULATOR-LIMITED (inconclusive), but")
    print("  the conclusion follows analytically from what the project already established:")
    print("  * The numerics fail cleanly because a sharp momentum cutoff violates the Ward identity")
    print("    via a surface term; even the photon (current) validation stays at percent-to-tens-")
    print("    of-percent, and Pauli-Villars / soft cutoffs do not fix it. A symmetry-preserving")
    print("    scheme (dimensional regularisation, or an exact lattice Ward identity) is required --")
    print("    a separate, careful computation, not this quick bubble.")
    print("  * The ANALYTIC argument stands on established results: the medium's IR sector is a")
    print("    Lorentz-invariant Dirac fermion (test_lorentz, test_dirac); its stress tensor is")
    print("    conserved by Noether (translation + Lorentz symmetry); a massless spin-2 coupled to")
    print("    a conserved stress tensor is forced by Weinberg's theorem to be Einstein => gamma=1,")
    print("    light bending x2. So the single missing coupling of test_einstein_source is supplied,")
    print("    in the IR, by the emergent-Lorentz flow -- but that is an ARGUMENT, not (yet) a clean")
    print("    in-model measurement. The honest status of the wall: crossed in principle by the IR")
    print("    fixed point; a direct numerical confirmation of graviton transversality remains open,")
    print("    as do the UV/lattice obstructions (Weinberg-Witten, the induced cosmological term).")
