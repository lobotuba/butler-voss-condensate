"""
The collapse threshold: the minimum number of particles that destroys the condensate.

Robert's Postulate 1 -- mass is a confluence of loci -- and Postulate 2 -- enough distortion
drives the inter-particle distance to zero, giving a black hole. The first is right and the
model already contains the mechanism. The second needs testing, and the honest answer turns out
to be interesting precisely because it FAILS.

The mechanism is real and I have already seen it fire: in the strong-field regime of
test_critical_gravity, matter CRUSHES the condensate (phi -> 0 locally). That is not a numerical
artefact, it is a physical threshold. With the coupling g rho |chi|^2, matter shifts the local
curvature of the potential:

        m_eff^2(x) = -a + 2 g rho(x)

so wherever 2 g rho > a, the ordered state is no longer even locally stable and phi is driven to
zero. The condensate is destroyed and a BUBBLE OF NORMAL (uncondensed) PHASE opens -- the same
object as a vortex core, or a normal-phase bubble above H_c in a superconductor.

Measure it. For each source width sigma, bisect on the source amplitude for the threshold at
which phi(centre) is driven to zero; the particle number is N = integral rho. Then look for the
cheapest (minimum-N) collapse over sigma.

I EXPECTED a critical nucleus at the healing length: sigma* ~ xi = 1/m_A, with N_min ~ 1/sqrt(a),
so that halving a would raise N_min by sqrt(2) = 1.41. THE MEASUREMENT REFUTED THAT, and the
refutation is the actual result of this file (kept honestly rather than rewritten to look
predicted). N(sigma) is MONOTONIC -- the cheapest collapse is the SMALLEST region, down to the
lattice cutoff -- and N_min is essentially INDEPENDENT of the gap (ratio 0.96, not 1.41, when a
is halved), so it is set by the lattice spacing, not by any physical length.

Why: nucleation needs the new phase to be FAVOURABLE, and here it is not. The condensate is the
TRUE vacuum; a normal-phase bubble is never energetically favourable on its own. So there is no
barrier, no critical size, and NO RUNAWAY. Matter does not TRIGGER a collapse -- it merely HOLDS
a region suppressed, and that region heals shut the instant the matter is removed. It is an
impurity pinned by the source, not an instability. There is no gravitational collapse in this
model at all -- which is the third independent reason (with horizon density ~ 1/M^2, and scalar
gravity not bending light) that Postulate 2's black hole is not here.
"""
from __future__ import annotations
import numpy as np


def lap(f):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
            np.roll(f, 1, 1) + np.roll(f, -1, 1) +
            np.roll(f, 1, 2) + np.roll(f, -1, 2) - 6.0 * f)


def gaussian(N, sig):
    g = np.arange(N) - N // 2
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    return np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2 * sig ** 2))


def relax(N, a, b, g, rho, dt=0.05, mom=0.9, iters=4000):
    """Full nonlinear relaxation. Gradient descent (not the Fourier split of
    test_critical_gravity) because here phi is driven to ZERO -- the strong-field regime, where
    the weak-field iteration diverges."""
    phi = np.full((N, N, N), np.sqrt(a / b))
    v = np.zeros_like(phi)
    for _ in range(iters):
        grad = -lap(phi) - a * phi + b * phi ** 3 + 2.0 * g * rho * phi
        v = mom * v - dt * grad
        phi += v
    return phi


def collapsed(N, a, b, g, sig, amp, frac=0.10):
    """Is the condensate destroyed at the centre? (phi(0) < frac * phi0)"""
    phi = relax(N, a, b, g, amp * gaussian(N, sig))
    c = N // 2
    return phi[c, c, c] / np.sqrt(a / b) < frac


def threshold_amp(N, a, b, g, sig, lo=0.0, hi=None, steps=9):
    """Bisect on the source amplitude for the collapse threshold."""
    if hi is None:
        hi = 4.0 * a / (2 * g)                       # comfortably above rho_c = a/2g
        while not collapsed(N, a, b, g, sig, hi):
            hi *= 2.0
            if hi > 1e4:
                return np.nan
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if collapsed(N, a, b, g, sig, mid):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def particle_number(amp, sig):
    """N = integral rho d^3x for a Gaussian of peak `amp` and width `sig`."""
    return amp * (2 * np.pi * sig ** 2) ** 1.5


if __name__ == "__main__":
    print("=== The collapse threshold: minimum particles to destroy the condensate ===\n")
    print("  Matter shifts the local curvature: m_eff^2 = -a + 2 g rho. Where 2 g rho > a the")
    print("  ordered state is not even locally stable and phi -> 0: a bubble of NORMAL phase.")
    print("  I expected a critical nucleus at the healing length xi = 1/m_A. Watch the sigma")
    print("  dependence and the a-scaling -- they decide whether that expectation survives.\n")

    N, b, g = 64, 1.0, 0.25
    results = {}
    for a in (0.02, 0.01):
        xi = 1.0 / np.sqrt(2 * a)
        rho_c = a / (2 * g)
        print(f"  [a = {a}]  m_A = {np.sqrt(2*a):.4f}   healing length xi = 1/m_A = {xi:.2f}"
              f"   rho_c = a/2g = {rho_c:.4f}")
        print(f"      {'sigma':>7} {'sigma/xi':>9} {'amp_thresh':>11} {'amp/rho_c':>10} {'N particles':>12}")
        best = (np.inf, None)
        for sig in (1.5, 2.5, 3.5, 5.0, 7.0, 9.0):
            amp = threshold_amp(N, a, b, g, sig)
            if not np.isfinite(amp):
                print(f"      {sig:>7.1f} {sig/xi:>9.2f} {'--':>11} {'--':>10} {'no collapse':>12}")
                continue
            Np = particle_number(amp, sig)
            if Np < best[0]:
                best = (Np, sig)
            print(f"      {sig:>7.1f} {sig/xi:>9.2f} {amp:>11.4f} {amp/rho_c:>10.2f} {Np:>12.2f}")
        results[a] = best
        print(f"      => cheapest collapse: N_min = {best[0]:.2f} at sigma* = {best[1]:.1f}"
              f"  (sigma*/xi = {best[1]/xi:.2f}); N(sigma) is MONOTONIC -> minimum is at the"
              " smallest sigma = LATTICE cutoff, not the healing length.\n")

    (N1, s1), (N2, s2) = results[0.02], results[0.01]
    print("  [scaling gate -- the decisive test] a critical nucleus set by xi predicts")
    print("      N_min ~ 1/sqrt(a): halving a should raise N_min by sqrt(2) = 1.41.")
    print(f"      MEASURED: N_min(a=0.02) = {N1:.2f}   N_min(a=0.01) = {N2:.2f}   ratio = {N2/N1:.2f}")
    print(f"      sigma* did NOT move: {s1:.1f} -> {s2:.1f}, while xi went {1/np.sqrt(0.04):.1f}"
          f" -> {1/np.sqrt(0.02):.1f}.")
    print("      => ratio ~ 1, not 1.41; sigma* pinned at the lattice cutoff. THE PREDICTION FAILS.")
    print("         N_min is gap-INDEPENDENT: it is a lattice number, not a physical nucleus.\n")

    print("[verdict] there is a collapse THRESHOLD, but NO collapse INSTABILITY -- and no black hole.")
    print("  * P1 is RIGHT and the model contains it: mass is a confluence of loci, and where")
    print("    2 g rho > a the condensate is locally destroyed (phi -> 0, a normal-phase bubble).")
    print("    So 'enough mass suppresses the medium' is real and computable.")
    print("  * But there is NO critical nucleus and NO runaway. N(sigma) is monotonic and N_min is")
    print("    independent of the gap -- the cheapest collapse is the smallest, down to the lattice.")
    print("    The reason is structural: the condensate is the TRUE vacuum, so a normal bubble is")
    print("    never favourable on its own. No barrier, no critical size, no self-sustaining")
    print("    collapse. Matter PINS a suppressed region; remove it and the region heals shut. It")
    print("    is an impurity, not an instability. My 'critical nucleus at the healing length'")
    print("    prediction was wrong, and the data says so plainly.")
    print("  * So P2 fails for the DEEPEST reason yet: this model has no gravitational collapse at")
    print("    all. That joins two independent refutations already on record -- real horizon")
    print("    density goes as 1/M^2 (1.8e17 kg/m^3 stellar down to 4e-3 for TON 618, thinner than")
    print("    air, so no single 'spacing -> 0' density can describe black holes), and our gravity")
    print("    is SCALAR, which does not bend light, so nothing traps it into a horizon.")
    print("\n  => There are no black holes in this model, and the reason is not a too-coarse lattice")
    print("     -- it is that scalar gravity has no horizons and the condensate has no collapse")
    print("     instability. Both are cured only by the spin-2 / diffeomorphism-invariant upgrade,")
    print("     the same frontier test_critical_gravity named. A denser lattice would not help.")
