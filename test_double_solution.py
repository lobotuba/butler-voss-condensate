"""
Frontier 3, step 3 (the last depth): de Broglie's DOUBLE SOLUTION -- deriving the
guidance equation v = grad(S)/m from the medium, instead of inserting it.

Steps 1-2 established: the emergent field's slow envelope obeys Schrodinger with hbar a
material property (the "pilot wave", carrying its phase S emerges for free), and |psi|^2
is the stochastic equilibrium of a diffusion (Born rule). BUT the link between the two --
the guidance rule that a PARTICLE rides the pilot wave's phase, v = grad(S)/m -- was still
put in by hand (test_born_rule set the drift from |psi|).

de Broglie's double-solution program (1927, 1950s): there is ONE field. Far from the
particle it is the smooth pilot wave psi_bg = a e^{i k x} (phase S = hbar k x, so
grad(S)/m = hbar k/m). AT the particle it has a localized SOLITON core. The claim: the
soliton is dragged by the background's phase gradient and drifts at exactly v = grad(S)/m,
so the guidance equation is DERIVED from the field's own (nonlinear) dynamics.

This is decisive and genuinely uncertain. In LINEAR theory superposition means the
background just slides through a resting particle -- no guidance (control below, g=0).
Guidance can only come from the NONLINEAR self-coupling of the medium coupling the soliton
to the background phase. Whether the resulting drift is *exactly* grad(S)/m is the open
question de Broglie never closed.

Medium model (units hbar = m = 1): the emergent massive field in its NR limit is the
nonlinear Schrodinger fluid
    i psi_t = -1/2 psi_xx - g |psi|^2 psi
(dispersion omega = k^2/2, so group velocity dω/dk = k = grad(S) -- the de Broglie value).
Its bright soliton eta*sech(eta x) is a legitimate localized "particle" of the same field.

Test: seed a RESTING soliton (no phase of its own) plus a smooth pilot wave of wavenumber
k; evolve; measure the soliton's drift velocity vs the prediction v = k. Sweep k for the
slope (de Broglie => slope 1). Control g=0 (linear) => no guidance.
"""
from __future__ import annotations
import numpy as np


def centroid(x, psi, L, remove=None):
    """Locate the particle robustly. If a separate pilot wave is present, remove its
    carrier (a single on-grid Fourier mode) so |residue|^2 is the clean particle density
    with NO interference cross term. Then take a periodic (circular) centroid so box
    wrap-around is handled exactly."""
    if remove is not None:
        ph = np.fft.fft(psi)
        ph[remove] = 0.0                              # delete the pilot-wave mode e^{i k_bg x}
        w = np.abs(np.fft.ifft(ph)) ** 2              # particle-only density
    else:
        w = np.abs(psi) ** 2
    theta = 2 * np.pi * x / L
    z = (w * np.exp(1j * theta)).sum()
    ang = np.angle(z) % (2 * np.pi)
    return float(ang * L / (2 * np.pi))


def evolve(n_wave, a_bg=0.0, k_particle=0.0, g=1.0, eta=1.0,
           L=200.0, N=2048, dt=0.005, T=30.0, x0=100.0):
    """Split-step (Strang) NLS. n_wave = integer harmonic; k = 2 pi n_wave / L (on-grid).
    a_bg=0 => SELF-GUIDED: soliton carries its own phase k_particle (no separate wave).
    a_bg>0 => PILOT-GUIDED: resting soliton + a separate pilot wave of wavenumber k."""
    x = np.linspace(0.0, L, N, endpoint=False)
    dx = x[1] - x[0]
    kf = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    k = 2 * np.pi * n_wave / L
    particle = (eta / np.cosh(eta * (x - x0))).astype(complex) * np.exp(1j * k_particle * x)
    background = a_bg * np.exp(1j * k * x)                        # the separate pilot wave, grad(S) = k
    psi = particle + background
    remove = n_wave if a_bg > 0 else None
    Lin = np.exp(-0.5j * kf ** 2 * dt)                           # full linear (kinetic) step
    ts, cs = [], []
    for n in range(int(T / dt)):
        psi = psi * np.exp(1j * g * np.abs(psi) ** 2 * dt / 2)   # half nonlinear
        psi = np.fft.ifft(Lin * np.fft.fft(psi))                # full linear
        psi = psi * np.exp(1j * g * np.abs(psi) ** 2 * dt / 2)   # half nonlinear
        if n % 20 == 0:
            ts.append(n * dt)
            cs.append(centroid(x, psi, L, remove))
    # unwrap the periodic centroid so a steady drift is a straight line
    cs = np.unwrap(np.array(cs) * 2 * np.pi / L) * L / (2 * np.pi)
    return np.array(ts), cs


def drift(ts, cs, t_settle=6.0):
    """Least-squares slope of centroid(t) after a short settling time."""
    m = ts > t_settle
    A = np.vstack([ts[m], np.ones(m.sum())]).T
    v, _ = np.linalg.lstsq(A, cs[m], rcond=None)[0]
    return float(v)


if __name__ == "__main__":
    print("=== Frontier 3, step 3: de Broglie double solution -- deriving guidance v = grad(S)/m ===\n")
    print("  ONE nonlinear field: smooth pilot wave (phase gradient k) + a RESTING soliton particle.")
    print("  Q: is the particle dragged to drift at exactly v = grad(S)/m = k, from the field alone?\n")

    L = 200.0
    ns = [3, 6, 9]
    ks = [2 * np.pi * n / L for n in ns]

    def slope_of(runner):
        vs = np.array([drift(*runner(n)) for n in ns])
        s = float(np.linalg.lstsq(np.array(ks)[:, None], vs[:, None], rcond=None)[0][0, 0])
        return vs, s

    print("  [A] SELF-GUIDED: the soliton carries its OWN phase e^{ikx} (its own de Broglie wave)")
    print(f"      {'k = grad(S)/m':>14} {'v_measured':>12} {'v/k':>8}")
    vA, sA = slope_of(lambda n: evolve(n, a_bg=0.0, k_particle=2 * np.pi * n / L))
    for k, v in zip(ks, vA):
        print(f"      {k:>14.3f} {v:>12.4f} {v / k:>8.3f}")
    print(f"      => slope = {sA:.3f}   (de Broglie lambda = h/p for the particle's OWN wave)\n")

    print("  [B] PILOT-GUIDED: a RESTING soliton + a SEPARATE pilot wave of phase gradient k")
    print(f"      {'k = grad(S)/m':>14} {'v_measured':>12} {'v/k':>8}")
    vB, sB = slope_of(lambda n: evolve(n, a_bg=0.15, k_particle=0.0))
    for k, v in zip(ks, vB):
        print(f"      {k:>14.3f} {v:>12.4f} {v / k:>8.3f}")
    print(f"      => slope = {sB:.3f}   (guidance by a SEPARATE wave => would be 1.000)\n")

    print("  [C] CONTROL g=0 (linear): resting bump + pilot wave, superposition only")
    print(f"      {'k':>14} {'v_measured':>12}")
    for n, k in zip(ns, ks):
        print(f"      {k:>14.3f} {drift(*evolve(n, a_bg=0.15, k_particle=0.0, g=0.0)):>12.4f}")
    print("      => no self-coupling, no motion: the resting particle stays put.\n")

    print("[verdict] a sharp split -- half of de Broglie's guidance is derived, half is not:")
    print(f"  * DERIVED (A, slope {sA:.2f} = 1): a localized excitation moves at the phase gradient of")
    print("    its OWN carrier, exactly. v = grad(S)/m = hbar k/m is built into the medium (the NLS")
    print("    envelope is Galilean-covariant), so de Broglie's lambda = h/p is a THEOREM here, not a")
    print("    postulate -- the physical content of 'the particle rides its wave' emerges for free.")
    print(f"  * NOT derived (B, slope {sB:.2f} ~ 0): a resting particle is NOT dragged to v = grad(S)/m")
    print("    by a SEPARATE pilot wave. At the particle core the total phase is dominated by the")
    print("    particle's own (flat) phase, so there is nothing to guide it; the nonlinear coupling")
    print("    only scatters it weakly and erratically. de Broglie's full 'double solution' -- a")
    print("    soliton phase-locked to and steered by a distinct pilot wave -- is not realized by")
    print("    naive nonlinear superposition (the very step he never rigorously closed).")
    print("\n  => Frontier 3 summit. Emergent from the condensate: the wave equation and hbar (step 1),")
    print("     the Born rule as a stochastic attractor (step 2), and de Broglie v = grad(S)/m for a")
    print("     particle's own wave (step 3A). NOT emergent: guidance by a separate pilot wave and the")
    print("     selection of a single definite outcome -- the hard core of the measurement problem")
    print("     stays a postulate. An honest boundary: most of QM is medium mechanics; the residue is")
    print("     exactly the piece that is unsolved for everyone.")
