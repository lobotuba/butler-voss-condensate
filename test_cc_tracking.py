"""
Deriving the tracking efficiency -- and finding the (H tau)^2 lead does NOT survive it.

test_cc_attractor floated a lead: if the medium creates nodes to track cosmic expansion, the linear
dilution -3H cancels and the residual departure is the second-order acceleration lag, delta ~ (H tau)^2,
giving rho_Lambda ~ M_P^2 H^2 = rho_crit -- the observed dark energy. That was explicitly CONDITIONAL on
"efficient tracking", left underived. This file derives the tracking efficiency from the medium's actual
node-creation physics, and the honest result is that the lead does not hold: the efficiency the condensate
provides gives the LINEAR residual (10^-61), not the second-order one (10^-122).

The efficiency sets the power. Write the node-creation rate as a part that offsets a fraction epsilon of
the dilution. The density departure then obeys
        d(delta)/dt = -3 H (1 - epsilon) - delta/tau     =>   delta_ss = -3 (1 - epsilon) H tau,
so rho_Lambda ~ B * (1 - epsilon) * H. The observed rho_Lambda ~ H^2 needs (1 - epsilon) ~ H tau -- a
tracking that is near-perfect but short by EXACTLY a fractional H tau. epsilon = 0 (no tracking) gives the
linear 10^-61; epsilon = 1 (perfect growth) gives zero; only the marginal epsilon = 1 - H tau gives 10^-122.

What the medium actually does. Node creation in a condensate is driven by the CHEMICAL POTENTIAL -- how
far the local density sits from equilibrium -- i.e. by the tension P = B*delta, a STATE (deficit) variable.
It is not driven by the expansion RATE H: the medium senses its own local density, not the global Hubble
rate. A deficit-driven creation gives epsilon = 0 (the medium only creates nodes once it is already
rarefied), hence the linear residual delta ~ -3 H tau and rho_Lambda ~ 10^-61. The alternative -- a
condensate that grows perfectly at its equilibrium density -- is epsilon = 1 and gives rho_Lambda = 0. The
observed value sits between them, at the marginal epsilon = 1 - H tau, which neither natural mechanism
selects. So the tracking efficiency, once derived, does not deliver the observed dark energy.

  [A] the efficiency-to-power map: delta_ss = -3(1-epsilon)H tau, and rho_Lambda at H_0 for epsilon = 0,
      1 - H tau, and 1 -- only the fine-tuned marginal value lands on 10^-122.
  [B] the medium's node-creation drive is DEFICIT-driven (the tension P = B*delta, a state variable),
      not RATE-driven (H): measured on the real condensate. So epsilon = 0.
  [C] therefore delta ~ -3 H tau (linear), rho_Lambda ~ 10^-61 -- the lead does not survive.
"""
from __future__ import annotations
import numpy as np

from test_cc_offequilibrium import bulk_modulus, P_of_delta

B = bulk_modulus()
TAU = 1.0
H0 = 1.2e-61                                                # observed H_0 / M_Planck
OBSERVED = 1e-122                                           # rho_Lambda / M_Planck^4


def delta_ss(one_minus_eps, H):
    """Steady-state density departure vs the tracking SHORTFALL (1-eps): delta = -3(1-eps)H tau.
    Parametrised by (1-eps) directly so the marginal case 1-eps = H tau ~ 1e-61 does not underflow
    against 1 in float64."""
    return -3 * one_minus_eps * H * TAU


if __name__ == "__main__":
    print("=== Deriving the tracking efficiency: does the (H tau)^2 lead survive? ===\n")
    print("  The residual dark energy is delta_ss = -3(1-epsilon)H tau, so rho_Lambda ~ B(1-epsilon)H.")
    print("  Only a tracking short by exactly a fractional H tau (epsilon = 1 - H tau) gives rho ~ H^2.\n")

    # ---------- [A] efficiency -> predicted dark energy ----------
    print("  [A] EFFICIENCY -> DARK ENERGY at H_0 (which epsilon reproduces the observed 10^-122?):")
    print(f"      {'shortfall 1-epsilon':>22} {'delta_ss':>14} {'rho_Lambda/M_P^4':>18} {'vs observed':>14}")
    for sf, tag in [(1.0, "1     (no tracking)"), (H0 * TAU, "H*tau (marginal)"), (0.0, "0     (perfect)")]:
        d = delta_ss(sf, H0)
        rho = B * abs(d)
        ratio = f"{rho/OBSERVED:.0e}x" if rho > 0 else "0 (< obs)"
        print(f"      {tag:>22} {d:>14.1e} {rho:>18.1e} {ratio:>14}")
    print("      => only the FINE-TUNED marginal shortfall 1-epsilon = H*tau lands near 10^-122. No")
    print("         tracking overshoots by ~61 orders; perfect growth gives exactly zero. Not generic.\n")

    # ---------- [B] the medium's node-creation drive: deficit- or rate-driven? ----------
    print("  [B] WHAT DRIVES NODE CREATION in the condensate -- the deficit (state) or the rate H?")
    print("      A condensate adds nodes in response to its CHEMICAL POTENTIAL = the tension P = B*delta,")
    print("      a function of the local density departure. Measured on the real medium:")
    for d in (-0.02, -0.01, 0.01, 0.02):
        P = P_of_delta(d)
        print(f"        delta = {d:>+5.2f}  ->  tension P = {P:>+8.4f}   (drive to create nodes ~ -P = {-P:>+8.4f})")
    print(f"      P/delta ~ B = {B:.1f} > 0: the creation drive is DEFICIT-driven (proportional to delta),")
    print("      a STATE variable. The medium senses its own local density, NOT the global Hubble rate H.")
    print("      A deficit-driven creation only fires once the medium is ALREADY rarefied => epsilon = 0.\n")

    # ---------- [C] the verdict ----------
    rho_passive = B * abs(delta_ss(1.0, H0))               # 1-eps = 1: no tracking (deficit-driven)
    print("  [C] THEREFORE the derived efficiency is epsilon = 0 (deficit-driven), not 1 - H*tau:")
    print(f"      delta ~ -3 H tau (linear),  rho_Lambda/M_P^4 ~ {rho_passive:.0e}  -- ~61 orders too big.")
    print("      The (H tau)^2 lead required epsilon = 1 - H*tau (marginal rate-tracking), which the")
    print("      condensate's deficit-driven node creation does NOT provide. The lead does not survive.\n")

    print("[verdict] the tracking efficiency, derived, does NOT deliver the observed dark energy:")
    print("  * The residual dark energy is rho_Lambda ~ B(1-epsilon)H, so the OBSERVED rho ~ H^2 requires")
    print("    a tracking efficiency short by exactly a fractional H*tau (epsilon = 1 - H*tau) -- a")
    print("    fine-tuned, non-generic value [A].")
    print("  * The medium's node creation is DEFICIT-driven: it responds to the local chemical potential")
    print("    (tension P = B*delta, measured), not to the global expansion rate H [B]. That gives")
    print("    epsilon = 0 and the LINEAR residual delta ~ -3 H tau, hence rho_Lambda ~ 10^-61 -- 61 orders")
    print("    too big [C]. The alternative, a condensate growing perfectly at equilibrium, is epsilon = 1")
    print("    and gives rho_Lambda = 0. Neither natural mechanism lands on 10^-122.")
    print("  * So test_cc_attractor's [C] tracking lead was too optimistic: derived rather than assumed,")
    print("    the efficiency gives 10^-61 (deficit-driven) or 0 (perfect growth), bracketing but not")
    print("    hitting the observed 10^-122. The coincidence is not explained by node tracking.")
    print("  * WHAT STANDS (unchanged): the EQUILIBRIUM cosmological constant is exactly zero, self-tuned")
    print("    and Weinberg-robust (test_vacuum_gravitates, test_cc_weinberg), and initial conditions are")
    print("    erased by the attractor in a Planck time (test_cc_attractor [A][B]). The observed NONzero")
    print("    Lambda remains underived -- the honest boundary of the CC solution, now confirmed against")
    print("    the one mechanism that looked like it might cross it.")
