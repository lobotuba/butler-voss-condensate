"""
The off-equilibrium cosmological constant: what the expansion lag predicts, and why it is not the answer.

The CC thread so far: the equilibrium vacuum has P = 0, so rho_Lambda = -P = 0 exactly, self-tuned and
robust to Weinberg's no-go. That is the whole 10^122 fine-tuning, dissolved -- but it makes Lambda EXACTLY
zero, and the observed value is small but NONzero (rho_Lambda ~ 10^-122 M_Planck^4). The only way to get a
nonzero Lambda from this mechanism is for the vacuum to sit slightly OFF equilibrium. This file builds the
one physical reason it would: an EXPANDING universe drags the medium off equilibrium faster than it can
relax back. It measures what that predicts, on the medium's own equation of state, and reports honestly
where it lands.

The dynamics. The universe expanding at Hubble rate H dilutes the medium (density drops as it is stretched);
the medium relaxes back toward its equilibrium density at a finite rate 1/tau (node creation/annihilation
is not instantaneous). For the fractional density departure delta = (n - n_eq)/n_eq,
        d(delta)/dt = -3 H - delta / tau        =>   steady state   delta_ss = -3 H tau,
so expansion holds the medium a fraction ~H tau BELOW equilibrium. Rarefied means under tension: P < 0, so
the gravitating rho_Lambda = -P > 0 -- a positive, dark-energy-like term. Its size is set by the medium's
bulk modulus B (measured here) and the departure: rho_Lambda = -P(delta_ss) ~ B * 3 H tau, i.e.

        rho_Lambda  is proportional to  H          -- a DYNAMICAL dark energy, not a constant.

  [A] the medium's equation of state: measure the bulk modulus B = dP/d(delta) at equilibrium (real LJ
      condensate), so the off-equilibrium pressure is the model's, not a guess.
  [B] the relaxation steady state: integrate the dilution-vs-relaxation ODE, confirm delta_ss = -3 H tau,
      and show rho_Lambda scales LINEARLY with H (dynamical dark energy).
  [C] the magnitude, honestly: with Planck-scale medium parameters (B ~ M_Planck^4, tau ~ t_Planck),
      rho_Lambda ~ M_Planck^3 H ~ 10^-61 M_Planck^4 at H = H_0 -- the GEOMETRIC MEAN of the problem, still
      ~61 orders too big. The equilibrium mechanism removed all 122 orders; the expansion lag puts ~61 back.
  [D] the equation of state: rho_Lambda ~ H gives an evolving w (about -1/2 in the matter era, -> -1 in de
      Sitter) -- a real, falsifiable signature, but the matter-era value is in tension with the observed
      w ~ -1. So the expansion lag is the wrong size AND the wrong w: the observed dark energy is probably
      not the medium's expansion lag, and the honest boundary of the CC solution is the equilibrium result.
"""
from __future__ import annotations
import numpy as np

from test_vacuum_gravitates import per_atom_energy, pressure, equilibrium

A0, _, _ = equilibrium()                                    # equilibrium spacing of the LJ condensate


def P_of_delta(delta_n, sigma=1.0, eps=1.0):
    """Pressure at a fractional DENSITY departure delta_n = (n - n_eq)/n_eq. In 2D n ~ 1/a^2, so the
    spacing at density n_eq(1 + delta_n) is a = a0 / sqrt(1 + delta_n)."""
    a = A0 / np.sqrt(1.0 + delta_n)
    return pressure(a, sigma, eps)


def bulk_modulus():
    """B = dP/d(delta_n) at equilibrium -- the medium's stiffness against a uniform density change."""
    d = 1e-4
    return (P_of_delta(d) - P_of_delta(-d)) / (2 * d)


def steady_state(H, tau, dt=None, nsteps=20000):
    """Integrate d(delta)/dt = -3 H - delta/tau to steady state; returns delta_ss (should be -3 H tau)."""
    dt = dt or 0.01 * tau
    delta = 0.0
    for _ in range(nsteps):
        delta += dt * (-3 * H - delta / tau)
    return delta


if __name__ == "__main__":
    print("=== The off-equilibrium cosmological constant: what the expansion lag predicts ===\n")
    print("  Equilibrium gives Lambda = 0 exactly. A nonzero Lambda needs the vacuum off equilibrium.")
    print("  Cosmic expansion does that: it dilutes the medium faster than it relaxes back, holding it a")
    print("  fraction ~H*tau below equilibrium -> a tension -> a positive, dark-energy-like rho_Lambda.\n")

    # ---------- [A] the medium's equation of state ----------
    B = bulk_modulus()
    print("  [A] THE MEDIUM'S EQUATION OF STATE (real LJ condensate):")
    print(f"      equilibrium spacing a0 = {A0:.4f},  bulk modulus B = dP/d(delta_n) = {B:.4f}")
    print(f"      near equilibrium P(delta) ~ B*delta: P(+1%) = {P_of_delta(0.01):+.5f}, "
          f"P(-1%) = {P_of_delta(-0.01):+.5f}  (compress -> pressure, rarefy -> tension)\n")

    # ---------- [B] relaxation steady state and the H-scaling ----------
    tau = 1.0
    print("  [B] RELAXATION STEADY STATE (dilution vs relaxation), tau = 1 in medium units:")
    print(f"      {'H':>10} {'delta_ss (ODE)':>16} {'-3 H tau (analytic)':>20} {'rho_Lambda = -P':>16}")
    Hs = [1e-3, 1e-4, 1e-5, 1e-6]
    rls = []
    for H in Hs:
        dss = steady_state(H, tau)
        rl = -P_of_delta(dss)
        rls.append(rl)
        print(f"      {H:>10.0e} {dss:>16.3e} {-3*H*tau:>20.3e} {rl:>16.3e}")
    slope = np.polyfit(np.log(Hs), np.log(rls), 1)[0]
    print(f"      => delta_ss = -3 H tau confirmed, and rho_Lambda scales as H^{slope:.2f} -- LINEAR in H.")
    print("         The off-equilibrium dark energy is DYNAMICAL (rho_Lambda ~ H), not a constant.\n")

    # ---------- [C] the magnitude, at Planck scale ----------
    H0_over_MP = 1.2e-61                                    # H_0 / M_Planck (observed)
    rhoL_over_MP4 = 3 * H0_over_MP                          # rho_Lambda/M_P^4 ~ 3 B tau H with B~M_P^4, tau~1/M_P
    observed = 1e-122
    print("  [C] THE MAGNITUDE, at Planck-scale medium parameters (B ~ M_Planck^4, tau ~ t_Planck):")
    print(f"      rho_Lambda ~ 3 B tau H ~ 3 M_Planck^3 H  =>  rho_Lambda/M_Planck^4 ~ 3 H0/M_Planck ~ {rhoL_over_MP4:.0e}")
    print(f"      observed rho_Lambda/M_Planck^4 ~ {observed:.0e}")
    print(f"      predicted / observed ~ {rhoL_over_MP4/observed:.0e}  -- still ~61 orders too big.")
    print("      This is the GEOMETRIC MEAN of the problem: the equilibrium mechanism removed all 122")
    print("      orders (Lambda = 0 exactly); the expansion lag puts sqrt of them -- ~61 -- back. To land")
    print("      at 10^-122 the vacuum must sit ~10^-122 off equilibrium, i.e. ~61 orders CLOSER than the")
    print("      naive H*tau lag: the smallness is relocated to 'why so near equilibrium', not derived.\n")

    # ---------- [D] the equation of state, and the honest verdict ----------
    print("  [D] THE EQUATION OF STATE: rho_Lambda ~ H gives w = -1 - Hdot/(3H^2). In matter domination")
    print("      (Hdot = -3/2 H^2) that is w = -1/2; as dark energy takes over (H -> const) w -> -1.")
    print("      So the lag predicts an EVOLVING w near -1/2 in the matter era -- a real signature, but in")
    print("      tension with the observed w ~ -1. Wrong size AND wrong w.\n")

    print("[verdict] the off-equilibrium residual is a real prediction, but not the observed dark energy:")
    print("  * FORM (right): an expanding universe holds the medium a fraction ~H*tau off equilibrium, so")
    print("    the residual gravitating energy is DYNAMICAL, rho_Lambda proportional to H -- measured on the")
    print("    medium's own equation of state [A][B], not assumed.")
    print("  * MAGNITUDE (wrong): at Planck-scale stiffness and relaxation, rho_Lambda ~ M_Planck^3 H ~")
    print("    10^-61 M_Planck^4 today -- the geometric mean of the CC problem, ~61 orders above the")
    print("    observed value [C]. The equilibrium mechanism dissolves all 122 orders; the lag restores ~61.")
    print("  * EQUATION OF STATE (wrong): w ~ -1/2 in the matter era, versus the observed w ~ -1 [D].")
    print("  * HONEST BOUNDARY OF THE CC SOLUTION: what the model settles is the EQUILIBRIUM cosmological")
    print("    constant -- exactly zero, self-tuned, Weinberg-robust (test_vacuum_gravitates, test_cc_weinberg).")
    print("    The observed small NONzero Lambda is NOT derived: the simplest off-equilibrium source (the")
    print("    cosmic expansion lag) is both too large and the wrong w, so the dark energy is probably a")
    print("    different, smaller departure -- why the universe sits so extraordinarily near the medium's")
    print("    equilibrium is the question this thread relocates the coincidence to, and does not close.")
