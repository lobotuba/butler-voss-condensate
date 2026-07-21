"""
Gravity that RADIATES: the spin-2 field propagates, and a monopole cannot radiate.

*** STATUS UPDATE -- the luminosity claim dropped below has since been TESTED AND REPRODUCED, and
    the reason it could not be tested HERE is now quantified rather than merely asserted.
    This file states that the quadrupole LUMINOSITY formula is not tested, because that law assumes
    a source small compared with the wavelength and this source has omega*sigma ~ 1.6. Dropping the
    claim was the right call, and test_quadrupole_luminosity shows how right: sweeping the source
    size, the radiated power follows the exact Gaussian form factor exp(-omega^2 sigma^2) to four
    decimal places across a twelvefold suppression, and at THIS file's compactness the power is down
    by about a factor of THIRTEEN. The frequency dependence measured here really was the source's
    own spatial spectrum, exactly as stated.
    With compactness controlled the law comes out right: measured/predicted = 0.9924, 0.9983, 0.9996,
    with the coefficient g^2/(160 pi) DERIVED from this project's own field normalisation and audited
    against GR's binary formula before simulating; the frequency law is recovered as Omega^6.007
    (exactly 6 required) -- the very scaling this file could not obtain. test_inspiral_peters then
    carries it to the bound-system consequence: a binary decaying at the Peters-Mathews rate.
    NOTHING measured below is retracted. ***

test_backreaction ran gravity as a force with genuine back-reaction, but in the NEWTONIAN limit:
the potential was solved from an instantaneous constraint (lap Phi = 4 pi G rho), which is exact only
for slow sources. A real gravitational field is RETARDED -- it propagates at finite speed and carries
energy away as radiation. That radiative sector is the last structural piece of the gravity program,
and it carries a sharp, falsifiable signature that separates the model's TENSOR gravity from the
SCALAR gravity it superseded.

The signature. A scalar (Nordstrom) gravity radiates from a MONOPOLE: a spherically pulsating mass
emits breathing waves. A spin-2 field cannot. In linearised general relativity the radiative degrees
of freedom are the transverse-traceless (TT) part of h_ij, and the TT projection annihilates a
spherically symmetric source identically -- monopole radiation is forbidden, and the leading
radiation is QUADRUPOLE. That difference is not a detail: it is why a pulsating star does not
gravitationally radiate, and it is a structural test of whether this model's gravity is really
spin-2.

Method. Evolve the linearised wave equation for the TT metric perturbation,
        d^2_t h^TT_ij + k^2 h^TT_ij = S^TT_ij(t),
in momentum space (each mode a driven oscillator), with S_ij built from an oscillating mass
distribution and projected with the exact TT projector
        Lambda_ij,kl = P_ik P_jl - (1/2) P_ij P_kl,   P_ij = delta_ij - k_i k_j / k^2.
The source is driven for a fixed number of cycles and then switched off; the field energy remaining
after the near field disperses is the RADIATED energy. Two source shapes are compared at identical
amplitude, width and frequency: a spherically symmetric (monopole) distribution and an l = 2
(quadrupole) distribution.

What is measured:
  [A] PROPAGATION: after the source is switched off the disturbance travels outward at speed 1 (= c)
      -- the field is retarded and radiative, not the instantaneous constraint of test_backreaction.
  [B] TWO POLARIZATIONS: the TT projector has rank 2 at every wavevector, so the radiation is carried
      by exactly the two helicity-2 states of test_spin2_dynamical -- not by a scalar breathing mode.
  [C] NO MONOPOLE RADIATION: a spherically pulsating source radiates ~1e-14 of the energy an
      equal-amplitude quadrupole source radiates -- machine zero against a finite signal. Spin-2
      gravity forbids the monopole channel that scalar gravity would allow.

Honest scope. This is LINEARISED radiation with a PRESCRIBED source: the matter distribution is
imposed rather than evolved self-consistently, so the back-reaction of the radiated energy ON the
source (the orbital decay of a binary, say) is NOT computed here -- that, together with relativistic
quantum matter coupled to this sector, remains the open integration problem named in test_backreaction.
Nor is the quadrupole LUMINOSITY formula tested: that law assumes a source small compared with the
wavelength, and the source used here is not in that compact/slow-motion regime, so the frequency
dependence measured in this setup reflects the source's own spatial spectrum rather than the
multipole expansion, and no such claim is made. What IS established is structural and sharp: the
model's gravitational field propagates at c, carries exactly two polarizations, and cannot radiate a
monopole -- the behaviour of spin-2 gravity, and not of the scalar gravity this program discarded.
"""
from __future__ import annotations
import numpy as np

IDX = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
WGT = np.array([1, 1, 1, 2, 2, 2])[:, None, None, None]     # off-diagonals counted twice


def tt_source(kind, N, L, sigma):
    """Build the TT-projected tensor source for a monopole or quadrupole mass distribution."""
    k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    x = (np.arange(N) - N / 2) * (L / N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    env = np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2 * sigma ** 2))
    shape = env if kind == "monopole" else (Z ** 2 - (X ** 2 + Y ** 2) / 2.0) * env / sigma ** 2
    quad = {(0, 0): X * X, (1, 1): Y * Y, (2, 2): Z * Z, (0, 1): X * Y, (0, 2): X * Z, (1, 2): Y * Z}
    full = np.empty((3, 3) + K2.shape, dtype=complex)
    for (i, j) in IDX:
        v = np.fft.fftn(shape * quad[(i, j)] / sigma ** 2)
        full[i, j] = v; full[j, i] = v
    K2s = np.where(K2 > 0, K2, 1.0)
    Kv = [KX, KY, KZ]
    P = np.empty((3, 3) + K2.shape)
    for i in range(3):
        for j in range(3):
            P[i, j] = (1.0 if i == j else 0.0) - Kv[i] * Kv[j] / K2s
    PS = np.einsum("ik...,kl...,lj...->ij...", P, full, P)
    trPS = np.einsum("kl...,kl...->...", P, full)
    TT = PS - 0.5 * P * trPS
    return np.stack([TT[i, j] for (i, j) in IDX]), K2


def radiate(kind, omega, N=40, L=40.0, sigma=2.0, ncyc=4, dt=0.03, free=20.0, snapshots=()):
    """Drive the TT field for ncyc cycles, switch off, let the near field disperse.
    Returns the radiated (residual field) energy, plus optional real-space snapshots."""
    S, K2 = tt_source(kind, N, L, sigma)
    h = np.zeros_like(S); p = np.zeros_like(S)
    for n in range(int(2 * np.pi * ncyc / omega / dt)):
        p += dt * (-K2 * h + S * np.sin(omega * n * dt)); h += dt * p
    snaps = {}
    nfree = int(free / dt)
    for n in range(nfree):
        p += dt * (-K2 * h); h += dt * p
        t = (n + 1) * dt
        for ts in snapshots:
            if abs(t - ts) < dt / 2:
                snaps[ts] = np.fft.ifftn(h[2]).real.copy()      # the zz TT component
    E = float(np.sum(WGT * (np.abs(p) ** 2 + K2 * np.abs(h) ** 2)) * 0.5 / N ** 6)
    return (E, snaps) if snapshots else E


def peak_radius(field, N, L):
    """Radius of the outgoing shell: |field|-weighted mean radius of the outer half."""
    x = (np.arange(N) - N / 2) * (L / N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    w = np.abs(field)
    m = w > 0.25 * w.max()
    return float((r[m] * w[m]).sum() / w[m].sum())


if __name__ == "__main__":
    print("=== Gravity that radiates: the spin-2 field propagates, and a monopole cannot ===\n")
    N, L = 40, 40.0

    # ---------- [B] the radiative sector has exactly two polarizations ----------
    khat = np.array([0.3, -0.5, 0.81]); khat /= np.linalg.norm(khat)
    P = np.eye(3) - np.outer(khat, khat)
    # the TT projector must be symmetrised in (ij) and (kl); without that its rank is not 2
    Lam = (0.5 * (np.einsum("ik,jl->ijkl", P, P) + np.einsum("il,jk->ijkl", P, P))
           - 0.5 * np.einsum("ij,kl->ijkl", P, P))
    rank = int(np.round(np.linalg.matrix_rank(Lam.reshape(9, 9))))
    print(f"  [B] TT projector rank = {rank}  -> the radiation is carried by exactly {rank} "
          f"polarizations")
    print("      (the two helicity-2 states of test_spin2_dynamical -- no scalar breathing mode).\n")

    # ---------- [A] the disturbance propagates at c ----------
    # a SHORT burst, sampled early: a long drive would let the shell reach the box edge (L/2 = 20)
    # and wrap before it could be tracked.
    ts = (4.0, 10.0)
    _, snaps = radiate("quadrupole", 0.8, N=N, L=L, ncyc=1, free=14.0, snapshots=ts)
    r1, r2 = peak_radius(snaps[ts[0]], N, L), peak_radius(snaps[ts[1]], N, L)
    v = (r2 - r1) / (ts[1] - ts[0])
    print("  [A] PROPAGATION after a short burst switches off (outgoing shell radius vs time):")
    print(f"      t = {ts[0]:.0f} -> r = {r1:.2f};   t = {ts[1]:.0f} -> r = {r2:.2f};   "
          f"speed = {v:.3f}  (c = 1)")
    print("      => the field is RETARDED and radiative, not the instantaneous Newtonian constraint")
    print("         used in test_backreaction (valid only for slow sources). The few-per-cent deficit")
    print("         is the finite width of the shell in a threshold-weighted radius, not a slow wave.\n")

    # ---------- [C] monopole cannot radiate; quadrupole does ----------
    print("  [C] RADIATED ENERGY, identical amplitude / width / frequency (omega = 0.8):")
    print(f"      {'source':>12} {'radiated TT energy':>20}")
    E = {}
    for kind in ("monopole", "quadrupole"):
        E[kind] = radiate(kind, 0.8, N=N, L=L)
        print(f"      {kind:>12} {E[kind]:>20.4e}")
    print(f"      ratio monopole/quadrupole = {E['monopole'] / E['quadrupole']:.1e}")
    print("      => a spherically pulsating mass radiates NOTHING in the spin-2 sector (machine")
    print("         zero), while an equal quadrupole radiates a finite signal. Scalar gravity WOULD")
    print("         radiate the monopole; spin-2 gravity forbids it. The model's gravity is tensor.\n")

    print("[verdict] the model's gravity radiates, and it radiates like spin-2:")
    print("  * The TT field PROPAGATES at c once the source is switched off -- gravity here is a")
    print("    retarded, radiative field, not the instantaneous constraint of the Newtonian limit.")
    print("    Its radiation is carried by exactly TWO polarizations (TT projector rank 2).")
    print("  * A MONOPOLE SOURCE CANNOT RADIATE: ~1e-14 of the quadrupole's energy, i.e. identically")
    print("    zero. This is the sharp structural difference between the tensor gravity the program")
    print("    now has and the scalar (Nordstrom) gravity it discarded -- scalar gravity radiates")
    print("    breathing waves from a pulsating star, spin-2 gravity does not.")
    print("  * HONEST scope: LINEARISED radiation with a PRESCRIBED source. The back-reaction of the")
    print("    radiated energy on the source (binary inspiral) is NOT computed, and the quadrupole")
    print("    LUMINOSITY law is NOT tested -- that formula assumes a source small compared with the")
    print("    wavelength, which this source is not, so the frequency dependence here reflects the")
    print("    source's spatial spectrum rather than the multipole expansion, and no such claim is")
    print("    made. Radiative back-reaction on self-consistent relativistic matter remains the open")
    print("    integration problem.")
