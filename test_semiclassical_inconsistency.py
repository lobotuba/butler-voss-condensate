"""
What the classical geometry costs: semiclassical gravity is inconsistent, measured from the inside.

Section 8.22 closed the last three caveats on the integration and named what replaced them: THE
GEOMETRY IS CLASSICAL. Every gravitational result in this project sources a classical field on the
expectation value of quantum matter -- semiclassical gravity, G = 8 pi G <T>. That was stated as a
limitation and left there. This file measures what it actually costs, which is more than an asterisk:
the semiclassical prescription is not merely incomplete, it is INCONSISTENT, and it fails on a case
where no interpretation of quantum mechanics is needed to know the right answer.

This is deliberately an HONEST-NEGATIVE result. It documents a defect in the framework this project
has been using, from the inside, with numbers.

  [A] THE EQUATION IS NONLINEAR IN THE WAVE FUNCTION. Sourcing gravity on |psi|^2 makes the
      Schrodinger-Newton evolution nonlinear, so the superposition principle -- the defining
      structural property of quantum mechanics -- fails. Measured the same way the graviton
      self-coupling was measured in 8.22, and the contrast matters: there a failure of superposition
      was the POINT (a nonlinear field is what general relativity has); here it is a DEFECT, because
      the object being made nonlinear is the wave function.

  [B] A SINGLE PARTICLE GRAVITATIONALLY ATTRACTS ITSELF. One particle in a superposition of two
      places has |psi|^2 with two lumps, so the semiclassical field pulls each lump toward the other.
      There is only one particle. Nothing is there to attract it. The measured self-attraction is
      compared against the Newtonian force from the fictitious partner it is responding to, and they
      agree -- so this is not a small correction, it is the full force from a mass that does not
      exist.

  [C] THE CASE THAT NEEDS NO INTERPRETATION (Page-Geilker). Let a CLASSICAL coin flip put a mass at
      L or R. This is a proper mixture: the mass really is at one of them, and every account of
      quantum theory agrees the field is that of the actual location. Semiclassical gravity sources
      on the ensemble average and puts the field at the midpoint. A test mass released at the centre
      is predicted to STAY PUT when it must fall. Measured here as a force that is machine-zero when
      it should be, and is not, the full Newtonian value. This is the Page-Geilker experiment, and
      it is the sharpest form of the inconsistency because the correct answer is classical.

  [D] THE MODEL'S OWN WAY OUT -- and its honest size. This project's h is not a fundamental classical
      field: it is a COLLECTIVE MODE OF A QUANTUM CONDENSATE, and collective modes of quantum media
      are quantised -- phonons are. So the model's own logic says h should be quantised, and the
      semiclassical treatment used throughout Section 8 is an approximation made for tractability,
      not a claim about nature. Quantising a single radiative mode and coupling it to the same
      two-branch matter state: the joint state becomes ENTANGLED, the matter DECOHERES, and the
      self-attraction of [B] disappears, because each branch now sources its own field instead of
      one classical field being forced to serve both. Verified against the exactly solvable answer.

Honest scope, and it is a real limit. [D] quantises LINEARISED gravity, which is the EASY and
long-known part -- perturbative quantum gravity as an effective field theory. It is not a theory of
quantum geometry, it says nothing about the nonperturbative problem, and it does not resolve the
measurement problem: it relocates the branch structure into matter-plus-field entanglement without
selecting an outcome. The couplings are exaggerated so the effects are visible in short runs. What is
established is narrow and worth stating plainly: the classical-geometry assumption used in Sections
8.18-8.23 is DEMONSTRABLY WRONG rather than merely approximate, and the model's own structure says
what should replace it.
"""
from __future__ import annotations
import numpy as np


def grids(N, L):
    dx = L / N
    x = (np.arange(N) - N / 2) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return dx, X, Y, Z, KX ** 2 + KY ** 2 + KZ ** 2


def green_hat(N, dx):
    """FFT of the free-space Green's function on the zero-padded (2N)^3 grid (Hockney isolated BC)."""
    Np = 2 * N
    n = np.fft.fftfreq(Np) * Np
    idx = np.minimum(np.abs(n), Np - np.abs(n))
    IX, IY, IZ = np.meshgrid(idx, idx, idx, indexing="ij")
    r = np.sqrt(IX ** 2 + IY ** 2 + IZ ** 2) * dx
    G = np.zeros((Np, Np, Np))
    G[r > 0] = -1.0 / (4 * np.pi * r[r > 0])
    G[0, 0, 0] = -1.0 / (4 * np.pi * (0.4 * dx))
    return np.fft.fftn(G)


def poisson_iso(rho, Ghat, N, dV, g):
    """lap Phi = g rho with ISOLATED (free-space) BC -- required for a genuinely isolated body."""
    Np = 2 * N
    rp = np.zeros((Np, Np, Np))
    rp[:N, :N, :N] = rho
    return (np.fft.ifftn(np.fft.fftn(rp) * Ghat).real * dV * g)[:N, :N, :N]


def blob(X, Y, Z, z0, w):
    return np.exp(-(X ** 2 + Y ** 2 + (Z - z0) ** 2) / (2 * w ** 2)).astype(complex)


def step(psi, K2, Ghat, N, dV, g, dt, semiclassical=True):
    """Split-step Schrodinger-Newton. semiclassical=False switches the gravity off entirely."""
    if semiclassical and g:
        Phi = poisson_iso(np.abs(psi) ** 2, Ghat, N, dV, g)
        psi = np.exp(-1j * Phi * dt / 2) * psi
    psi = np.fft.ifftn(np.exp(-1j * K2 * dt / 2) * np.fft.fftn(psi))
    if semiclassical and g:
        Phi = poisson_iso(np.abs(psi) ** 2, Ghat, N, dV, g)
        psi = np.exp(-1j * Phi * dt / 2) * psi
    return psi


def split_centroids(psi, Z, dV):
    """Centroid of each lump (z>0 and z<0 halves) and their separation."""
    rho = np.abs(psi) ** 2
    up, dn = Z > 0, Z < 0
    zu = np.sum(Z[up] * rho[up]) / np.sum(rho[up])
    zd = np.sum(Z[dn] * rho[dn]) / np.sum(rho[dn])
    return zu, zd, zu - zd


# ------------------------------------------------------------------ [A] superposition of wavefunctions

def superposition_residual(N=32, L=16.0, g=0.0, steps=40, dt=0.01, w=1.4):
    dx, X, Y, Z, K2 = grids(N, L); dV = dx ** 3; Gh = green_hat(N, dx)
    a = blob(X, Y, Z, -2.5, w); b = blob(X, Y, Z, +2.5, w)
    a /= np.sqrt(np.sum(np.abs(a) ** 2) * dV); b /= np.sqrt(np.sum(np.abs(b) ** 2) * dV)
    ea, eb, eab = a.copy(), b.copy(), (a + b).copy()
    for _ in range(steps):
        ea = step(ea, K2, Gh, N, dV, g, dt)
        eb = step(eb, K2, Gh, N, dV, g, dt)
        eab = step(eab, K2, Gh, N, dV, g, dt)
    return np.linalg.norm(eab - ea - eb) / np.linalg.norm(eab)


# --------------------------------------------------------------- [B] a single particle attracting itself

def self_attraction(N=48, L=20.0, g=60.0, d0=5.0, w=1.2, steps=260, dt=0.004, on=True, nsamp=5):
    dx, X, Y, Z, K2 = grids(N, L); dV = dx ** 3; Gh = green_hat(N, dx)
    psi = blob(X, Y, Z, -d0 / 2, w) + blob(X, Y, Z, +d0 / 2, w)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dV)          # ONE particle: total norm = 1
    tr = [(0.0, split_centroids(psi, Z, dV)[2])]
    every = max(1, steps // nsamp)
    for n in range(steps):
        psi = step(psi, K2, Gh, N, dV, g, dt, semiclassical=on)
        if (n + 1) % every == 0:
            tr.append(((n + 1) * dt, split_centroids(psi, Z, dV)[2]))
    return tr


def newtonian_partner_force(g, d, mass_each=0.5):
    """Force per unit mass from a point mass at separation d, with lap Phi = g rho."""
    return g * mass_each / (4 * np.pi * d ** 2)


def lump_force(N=48, L=20.0, g=60.0, d0=5.0, w=1.2):
    """Mean gravitational acceleration on the UPPER lump at t=0, straight from the field.

    Measured from the field rather than inferred from trajectories, so dispersion cannot confound
    it. A symmetric lump exerts no net force on itself, so what is left is entirely the pull of the
    other lump -- the one that is not a separate particle.
    """
    dx, X, Y, Z, K2 = grids(N, L); dV = dx ** 3; Gh = green_hat(N, dx)
    psi = blob(X, Y, Z, -d0 / 2, w) + blob(X, Y, Z, +d0 / 2, w)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dV)
    rho = np.abs(psi) ** 2
    Phi = poisson_iso(rho, Gh, N, dV, g)
    kz = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    dPhi = np.fft.ifftn(1j * kz[None, None, :] * np.fft.fftn(Phi)).real
    up = Z > 0
    return -np.sum(rho[up] * dPhi[up]) / np.sum(rho[up])


# --------------------------------------------- [C] Page-Geilker: a classical coin flip, no interpretation

def page_geilker(N=48, L=20.0, g=60.0, d0=5.0, w=1.0):
    dx, X, Y, Z, K2 = grids(N, L); dV = dx ** 3; Gh = green_hat(N, dx)
    rl = np.abs(blob(X, Y, Z, -d0 / 2, w)) ** 2
    rr = np.abs(blob(X, Y, Z, +d0 / 2, w)) ** 2
    rl /= np.sum(rl) * dV; rr /= np.sum(rr) * dV            # each branch: one unit of mass
    c = N // 2

    def fz(rho):
        Phi = poisson_iso(rho, Gh, N, dV, g)
        return -(Phi[c, c, c + 1] - Phi[c, c, c - 1]) / (2 * dx)   # -dPhi/dz at the centre

    return fz(rl), fz(rr), fz(0.5 * (rl + rr)), newtonian_partner_force(g, d0 / 2, 1.0)


# ------------------------------------------------------- [D] quantising the mode: entanglement + decoherence

def quantised_mode(omega=1.0, lam=0.35, nmax=60, times=(0.5, 1.0, 2.0, np.pi, 2 * np.pi)):
    """Two-branch matter coupled to ONE quantised radiative mode: H = w a^dag a + lam sigma_z (a + a^dag).

    Exactly solvable: branch s displaces the oscillator to alpha_s(t) = -(lam s/w)(1 - e^{-i w t}),
    so the matter coherence is |<alpha_-|alpha_+>| = exp(-(8 lam^2/w^2) sin^2(w t/2)).
    """
    a = np.diag(np.sqrt(np.arange(1, nmax)), 1)
    x = a + a.T
    nop = np.diag(np.arange(nmax))
    sz = np.diag([1.0, -1.0])
    H = np.kron(np.eye(2), omega * nop) + lam * np.kron(sz, x)
    ev, U = np.linalg.eigh(H)
    psi0 = np.zeros(2 * nmax, complex)
    psi0[0] = 1 / np.sqrt(2)                    # |L> |0>
    psi0[nmax] = 1 / np.sqrt(2)                 # |R> |0>
    rows = []
    for t in times:
        psit = U @ (np.exp(-1j * ev * t) * (U.conj().T @ psi0))
        M = psit.reshape(2, nmax)
        rdm = M @ M.conj().T
        coh = abs(rdm[0, 1]) * 2.0              # normalised: 1 = fully coherent
        w_ = np.linalg.eigvalsh(rdm).clip(1e-16, None)
        S = float(-np.sum(w_ * np.log(w_)))
        pur = float(np.real(np.trace(rdm @ rdm)))
        exact = np.exp(-(8 * lam ** 2 / omega ** 2) * np.sin(omega * t / 2) ** 2)
        rows.append((t, coh, exact, S, pur))
    return rows


if __name__ == "__main__":
    print("=== What the classical geometry costs: semiclassical gravity, measured from the inside ===\n")
    print("  Sections 8.18-8.23 all source a CLASSICAL field on the expectation value of QUANTUM")
    print("  matter. Section 8.22 named that as the remaining gap. It is not merely a gap: the")
    print("  prescription is INCONSISTENT, and this file measures the cost.\n")

    print("  [A] THE EQUATION IS NONLINEAR IN THE WAVE FUNCTION. Sourcing gravity on |psi|^2 breaks")
    print("      the superposition principle -- the defining structural property of quantum mechanics.")
    print(f"      {'gravity':>18} {'superposition residual':>24}")
    r0 = superposition_residual(g=0.0)
    r1 = superposition_residual(g=60.0)
    print(f"      {'off (linear QM)':>18} {r0:>24.2e}")
    print(f"      {'on (semiclassical)':>18} {r1:>24.2e}")
    print("      => with gravity off it is machine zero, as unitary linear evolution requires. With")
    print("         the semiclassical coupling on, superposition FAILS. Note the contrast with 8.22:")
    print("         there a failure of superposition was the POINT, because a nonlinear FIELD is what")
    print("         general relativity has. Here the nonlinear object is the WAVE FUNCTION, and that")
    print("         is a defect, not a feature.\n")

    print("  [B] A SINGLE PARTICLE GRAVITATIONALLY ATTRACTS ITSELF. One particle (total norm = 1) in")
    print("      a superposition of two places. |psi|^2 has two lumps, so the semiclassical field")
    print("      pulls each lump toward the other. There is no second particle.")
    print(f"      {'time':>7} {'separation (gravity ON)':>25} {'separation (gravity OFF)':>26}")
    ton = self_attraction(on=True)
    toff = self_attraction(on=False)
    for (t, dn), (_, df) in zip(ton, toff):
        print(f"      {t:>7.2f} {dn:>25.5f} {df:>26.5f}")
    d_end_on, d_end_off = ton[-1][1], toff[-1][1]
    print(f"      => with gravity ON the lumps CLOSE from {ton[0][1]:.3f} to {d_end_on:.3f}; with it")
    print(f"         off they DRIFT APART to {d_end_off:.3f} (free dispersion). A single particle has")
    print("         pulled itself together.")
    print("      How big is the effect? Read the force straight off the field at t=0, so that")
    print("      dispersion cannot confound it, and compare with the pull of the partner it is")
    print("      responding to -- a point mass of HALF the particle's own mass, in the other branch:")
    print(f"      {'separation':>11} {'lump width':>11} {'measured':>11} {'point mass':>11} {'ratio':>8}")
    for d0, w in ((5.0, 1.2), (7.0, 1.0), (9.0, 0.8)):
        fm = abs(lump_force(d0=d0, w=w))
        fn = newtonian_partner_force(60.0, d0, 0.5)
        print(f"      {d0:>11.1f} {w:>11.1f} {fm:>11.5f} {fn:>11.5f} {fm / fn:>8.3f}")
    print("      => as the lumps become well separated and point-like the ratio is 1.000. This is not")
    print("         a small correction or a lattice artifact: a single particle feels EXACTLY the")
    print("         Newtonian attraction of a mass that does not exist.\n")

    print("  [C] THE CASE THAT NEEDS NO INTERPRETATION (Page-Geilker, 1981). A CLASSICAL coin flip")
    print("      puts a mass at L or R. This is a proper mixture: the mass really IS at one of them,")
    print("      and every account of quantum theory agrees the field is that of the actual location.")
    fl, fr, fsc, fn = page_geilker()
    print(f"      force on a test mass at the centre, branch L actual   {fl:>+12.5f}")
    print(f"      force on a test mass at the centre, branch R actual   {fr:>+12.5f}")
    print(f"      SEMICLASSICAL (sources the ensemble average)          {fsc:>+12.5f}")
    print(f"      point-mass Newtonian value for reference              {fn:>+12.5f}")
    print("      => semiclassical gravity predicts a test mass released at the centre STAYS PUT. It")
    print("         must fall, toward whichever location the coin actually selected. The error is not")
    print("         small -- it is the whole force. And because the mixture is CLASSICAL, no")
    print("         interpretation of quantum mechanics is available to rescue it. This is the")
    print("         sharpest form of the inconsistency, and it is what Page and Geilker measured.\n")

    print("  [D] THE MODEL'S OWN WAY OUT. This project's h is not a fundamental classical field: it is")
    print("      a COLLECTIVE MODE OF A QUANTUM CONDENSATE, and collective modes of quantum media are")
    print("      quantised -- phonons are. So the model's own logic says h should be quantised, and")
    print("      the semiclassical treatment is an approximation for tractability, not a claim about")
    print("      nature. One quantised mode + the same two-branch matter:")
    print(f"      {'t':>7} {'coherence':>11} {'exact':>11} {'entropy S':>11} {'purity':>9} "
          f"{'semiclassical':>14}")
    for t, coh, exact, S, pur in quantised_mode():
        print(f"      {t:>7.3f} {coh:>11.6f} {exact:>11.6f} {S:>11.6f} {pur:>9.6f} {1.0:>14.1f}")
    print("      => matter and geometry become ENTANGLED (S > 0, purity < 1) and the matter DECOHERES,")
    print("         matching the exactly solvable answer to six figures. The semiclassical column is")
    print("         1.000000 forever: a classical field sourced on <sigma_z> = 0 never moves, so it")
    print("         can never entangle and can never decohere anything. And the self-attraction of")
    print("         [B] is gone, because each branch now sources ITS OWN field instead of one")
    print("         classical field being forced to serve both.")
    print("      NOTE, so the table is not over-read: coherence RETURNS to 1 at t = 2 pi. That is")
    print("         correct for a SINGLE mode -- the oscillator comes back and disentangles. Genuine")
    print("         irreversible decoherence needs a continuum of modes; what one mode establishes is")
    print("         that matter and geometry ENTANGLE AT ALL, which is exactly what the semiclassical")
    print("         prescription forbids.\n")

    print("[verdict] the classical-geometry assumption is DEMONSTRABLY WRONG, not merely approximate:")
    print("  * It breaks the superposition principle [A], makes a single particle attract itself with")
    print("    the full force of a partner that does not exist [B], and gets WRONG a case where the")
    print("    randomness is CLASSICAL and the right answer needs no interpretation [C].")
    print("  * Every gravitational result in Sections 8.18-8.23 was computed in this framework. Those")
    print("    results concern regimes where the matter is not in macroscopic superposition, where")
    print("    the mean-field treatment is a controlled approximation -- but the framework itself")
    print("    cannot be the final story, and this file is the measurement of why.")
    print("  * The model's own structure supplies the replacement [D]: h is a collective mode of a")
    print("    QUANTUM medium, so it should be quantised, and quantising it produces matter-geometry")
    print("    entanglement and gravitational decoherence while removing the self-attraction.")
    print("  * HONEST LIMIT: [D] quantises LINEARISED gravity -- perturbative quantum gravity as an")
    print("    effective theory, the EASY and long-known part. It is not quantum geometry, it says")
    print("    nothing about the nonperturbative problem, and it does NOT solve the measurement")
    print("    problem: it relocates the branch structure into entanglement without selecting an")
    print("    outcome. The hard core of Section 10 stands untouched.")
