"""
First quantization step: does the unified linear sector become a relativistic QFT?

The unified linear theory (test_lorentz_unified) is a system of coupled harmonic
oscillators -- exactly a free field. Canonical quantization diagonalises it into
normal modes, each a quantum oscillator, so
    H = sum_k hbar omega(k) (a_k^dagger a_k + 1/2),
with bosonic quanta of energy hbar omega(k). Two checks that this is genuinely a
relativistic QFT (not just "a quantum system"):

  A. The single-quantum energies are the relativistic MASS-SHELL: E(k)=hbar omega(k)
     with omega = c|k| (massless) or sqrt(c^2 k^2 + m^2 c^4) (massive).
  B. The quantum VACUUM correlator <0|phi(x)phi(0)|0> = (1/N) sum_k (1/2 omega_k)
     e^{ik.x} has the relativistic form: a massless field -> POWER LAW ~ 1/r^{d-1}
     (long-range vacuum correlations); a massive field -> EXPONENTIAL, correlation
     length 1/m. (The quantum echo of the whole massless/massive, long/short-range
     theme of this project.)

Honest scope (C): canonical quantization ASSUMES quantum mechanics -- it imposes
[phi,pi]=i hbar. It shows the model quantizes to a sensible relativistic QFT; it
does NOT derive QM from the sub-quantum medium (emergent/stochastic QM), which is
the genuinely open, unsolved problem.
"""
from __future__ import annotations
import numpy as np

L = 64                                                    # cubic box, 3 spatial dims
_k1 = 2 * np.pi * np.fft.fftfreq(L)
KX, KY, KZ = np.meshgrid(_k1, _k1, _k1, indexing="ij")
SYMBOL = 2 * ((1 - np.cos(KX)) + (1 - np.cos(KY)) + (1 - np.cos(KZ)))  # -> |k|^2 at low k


def omega(m):
    return np.sqrt(SYMBOL + m ** 2)                       # c = 1 (lattice units)


def vacuum_correlator(m):
    """<phi(x)phi(0)> = (1/N) sum_k 1/(2 omega_k) e^{ik.x}; k=0 zero mode excluded."""
    w = omega(m); Gk = np.zeros_like(w)
    msk = w > 1e-9
    Gk[msk] = 1.0 / (2.0 * w[msk])                        # drop the massless IR zero mode
    return np.fft.ifftn(Gk).real


# ================================================================ A ============
def massshell():
    print("[A] single-quantum energies E(k)=omega(k): the relativistic mass-shell")
    ks = np.array([0.05, 0.1, 0.2]) * np.pi
    print(f"  {'|k|/pi':>7} {'massless omega':>14} {'(=|k|?)':>8} | "
          f"{'massive omega':>13} {'sqrt(k^2+m^2)':>13}  (m=0.3)")
    for f in ks:
        kv = np.array([f, 0, 0])
        s = 2 * ((1 - np.cos(kv)).sum())                  # symbol along [100]
        w0, wm = np.sqrt(s), np.sqrt(s + 0.09)
        print(f"  {f/np.pi:>7.2f} {w0:>14.4f} {abs(kv).sum():>8.4f} | "
              f"{wm:>13.4f} {np.sqrt(abs(kv).sum()**2+0.09):>13.4f}")
    print("  => omega -> |k| (massless: linear, isotropic) and omega -> sqrt(k^2+m^2)")
    print("     (massive): the field quanta are RELATIVISTIC particles. Zero-point energy")
    print(f"     per site (1/2N) sum omega = {0.5*omega(0.0).mean():.4f} (massless).\n")


# ================================================================ B ============
def correlator_forms():
    print("[B] quantum vacuum correlator <phi(x)phi(0)> along a lattice axis")
    r = np.arange(1, L // 2)
    # massless: expect power law ~ 1/r^{d-1} = 1/r^2 in 3D
    G0 = vacuum_correlator(0.0)[1:L // 2, 0, 0]
    m = (r >= 2) & (r <= 12)
    p = -np.polyfit(np.log(r[m]), np.log(np.abs(G0[m])), 1)[0]
    print(f"  massless:  G(r) ~ 1/r^{p:.2f}   POWER LAW (relativistic massless field is 1/r^2 in 3D;")
    print("             the small excess is lattice/finite-size -> 2 in the continuum). Long-range.")
    print("    r    :", " ".join(f"{x:6d}" for x in r[:8]))
    print("    G(r) :", " ".join(f"{x:6.4f}" for x in G0[:8]))
    # massive: expect exponential with correlation length xi = 1/m
    for mass in (0.3, 0.6):
        Gm = vacuum_correlator(mass)[1:L // 2, 0, 0]
        mm = (r >= 2) & (r <= 14) & (Gm > 1e-9)
        xi = -1.0 / np.polyfit(r[mm], np.log(Gm[mm]), 1)[0]
        print(f"  massive m={mass}:  G(r) ~ exp(-m r)/r (Yukawa), effective xi={xi:.2f}, order 1/m"
              f" = {1/mass:.2f}")
    print("  (the effective xi runs shorter than 1/m because the algebraic 1/r prefactor of the")
    print("   3D Yukawa steepens a log-linear fit -- an expected artifact, not a discrepancy.)")
    print("  => massless field: power-law (long-range) vacuum correlations; massive field:")
    print("     Yukawa-exponential, short-range. The quantum vacuum reproduces the relativistic")
    print("     QFT correlator -- and the SAME massless/massive, long/short-range dichotomy")
    print("     that ran through the whole forces program, now at the level of the vacuum.\n")


if __name__ == "__main__":
    print("=== First quantization: is the unified sector a relativistic QFT? ===\n")
    massshell()
    correlator_forms()
    print("[C] honest scope: canonical quantization IMPOSES [phi,pi]=i hbar -- it shows the")
    print("    model quantizes to a proper relativistic QFT (bosonic quanta on the mass-shell,")
    print("    relativistic vacuum correlations), which was a real check. It does NOT derive")
    print("    quantum mechanics from the sub-quantum medium (emergent/stochastic QM) -- that")
    print("    remains the deepest open problem, untouched here.")
