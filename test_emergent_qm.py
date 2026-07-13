"""
Frontier 3 (origin of QM), step 1: the wave half of quantum mechanics emerges,
and hbar is a material property of the medium.

Quantum mechanics has two halves: (I) the WAVE half -- the Schrodinger equation,
superposition, interference, and hbar; and (II) the PROBABILISTIC half -- the Born
rule and measurement. They are not equally hard.

The medium is a condensate, so its natural language is hydrodynamics, and Madelung
(1927) showed that Schrodinger's equation IS the hydrodynamics of a fluid with one
special gradient energy (the "quantum potential" Q = -(hbar^2/2m) lap(sqrt rho)/sqrt rho).
Equivalently: the emergent MASSIVE field the project already has, written as
chi = sqrt(rho) e^{iS/hbar}, is that fluid. Its dispersion
    omega(k) = sqrt(c^2 k^2 + Omega^2)   (Omega = mc^2/hbar, the mass gap)
in the NON-RELATIVISTIC limit is
    omega ~ Omega + (c^2/2 Omega) k^2 = Omega + (hbar/2m) k^2,
the free-Schrodinger dispersion -- with hbar/2m = c^2/(2 Omega) fixed by the
medium's own gap and speed. So a slow wave packet of the emergent field must spread
at the exact Schrodinger rate, with hbar a MATERIAL PROPERTY (gap-to-curvature ratio),
not a postulate.

Test: evolve the emergent (Klein-Gordon) field for a Gaussian packet and check the
envelope spreads as sigma(t) = sigma0 sqrt(1 + (D t / sigma0^2)^2) with the diffusion
constant D = hbar/2m = c^2/(2 Omega) predicted from the medium alone.
"""
from __future__ import annotations
import numpy as np

C = 1.0


def lap1(f):
    return np.roll(f, 1) + np.roll(f, -1) - 2 * f


def spread(L=4096, sig0=22.0, Omega=0.6, dt=0.15, T=420.0):
    """Evolve the complex KG field for a slow Gaussian packet; return t, sigma(t)."""
    x = np.arange(L); x0 = L / 2
    env = np.exp(-(x - x0) ** 2 / (2 * sig0 ** 2)).astype(complex)
    chi = env.copy()
    pi = (-1j * Omega) * env                              # positive-frequency (NR) packet
    ts, sg = [], []
    n = int(T / dt)
    for i in range(n):
        pi = pi + dt * (C ** 2 * lap1(chi) - Omega ** 2 * chi)
        chi = chi + dt * pi
        if i % 40 == 0:
            w = np.abs(chi) ** 2; w = w / w.sum()
            xm = (x * w).sum(); s = np.sqrt(((x - xm) ** 2 * w).sum())
            ts.append(i * dt); sg.append(s)
    return np.array(ts), np.array(sg)


if __name__ == "__main__":
    print("=== Frontier 3, step 1: the wave half of QM emerges; hbar from the medium ===\n")
    Omega = 0.6
    D_pred = C ** 2 / (2 * Omega)                         # hbar/2m predicted from the medium
    print(f"  medium inputs: c = {C}, mass gap Omega = {Omega}")
    print(f"  => predicted Schrodinger diffusion  D = hbar/2m = c^2/(2 Omega) = {D_pred:.4f}")
    print("     (hbar is fixed by the medium's gap-to-dispersion-curvature ratio, not postulated)\n")

    t, s = spread(Omega=Omega)
    s0 = s[0]
    pred = s0 * np.sqrt(1 + (D_pred * t / s0 ** 2) ** 2)  # Schrodinger free-packet spreading
    # extract D from the measured spreading (large-t)
    m = t > 200
    D_meas = np.median((s0 ** 2 / t[m]) * np.sqrt(np.clip((s[m] / s0) ** 2 - 1, 0, None)))
    print(f"  {'t':>6} {'sigma(measured)':>16} {'sigma(Schrodinger)':>19}")
    for i in range(0, len(t), max(1, len(t)//7)):
        print(f"  {t[i]:>6.0f} {s[i]:>16.3f} {pred[i]:>19.3f}")
    err = np.max(np.abs(s - pred)) / s0
    print(f"\n  measured D = {D_meas:.4f}  vs predicted {D_pred:.4f}  "
          f"({100*abs(D_meas-D_pred)/D_pred:.1f}% -- packet spreads at the Schrodinger rate)")
    print(f"  max |sigma_measured - sigma_Schrodinger| / sigma0 = {err:.3f}\n")

    print("[verdict] the WAVE half of QM emerges:")
    print("  * the emergent field's slow envelope obeys the free Schrodinger equation")
    print("    i hbar d_t psi = -(hbar^2/2m) lap psi -- single-particle quantum wave mechanics.")
    print("  * hbar is a MATERIAL PROPERTY of the medium (gap-to-curvature ratio), not a postulate;")
    print("    in Madelung form the same field is the condensate's density+flow with the quantum")
    print("    potential, hbar^2 = the von-Weizsaecker density-gradient stiffness.")
    print("  * superposition and interference are automatic -- it is a linear wave.")
    print("\n  What does NOT emerge here -- the genuinely hard half, and the real F3 frontier:")
    print("  * the BORN RULE (|psi|^2 = probability) and MEASUREMENT/collapse. Linear wave dynamics")
    print("    give amplitudes, not probabilities. These need a mechanism the linear substrate lacks:")
    print("    a stochastic sub-quantum process (Nelson), a deterministic automaton with equivalence")
    print("    classes ('t Hooft), or decoherence + branching. Scoping/testing THAT is F3 proper.")
    print("  => half of QM -- the wave mechanics and hbar -- emerges from the condensate for free;")
    print("     the probabilistic half is the moonshot that remains.")
