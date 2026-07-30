"""
Why does the universe sit so near the medium's equilibrium? An attractor, an aging universe -- and the
coincidence rho_Lambda ~ rho_crit, if the medium tracks the expansion.

test_cc_offequilibrium left the honest gap: the equilibrium cosmological constant is zero, but the
observed nonzero value needs the vacuum a fraction ~10^-122 off equilibrium, while the naive expansion
lag gives ~10^-61 (61 orders too big). This file pushes on WHY the universe is so extraordinarily near
equilibrium -- the initial-conditions and relaxation-dynamics question -- and finds three things.

  [A] EQUILIBRIUM IS A DYNAMICAL ATTRACTOR that ERASES initial conditions. The relaxation
      d(delta)/dt = -3H - delta/tau has homogeneous solution delta ~ e^(-t/tau): any initial departure
      -- any initial vacuum energy, however large -- decays in a few relaxation times tau ~ t_Planck.
      So there is NO initial cosmological constant to fine-tune: whatever the medium started with is
      forgotten in a Planck time, and today's Lambda is the driven ATTRACTOR value, not a relic. The
      "why is the initial Lambda tuned to 10^-122" version of the problem simply does not arise here.

  [B] Lambda DECAYS as the universe AGES. The attractor value tracks H: rho_Lambda ~ H, so as the
      universe expands and H falls, the cosmological constant falls with it. Lambda is small NOW because
      the universe is OLD (H << M_Planck). Its smallness is the smallness of H/M_Planck -- the age of
      the universe -- not a tuning. This is a decaying, dynamical dark energy, large in the early
      universe and tiny today.

  [C] THE COINCIDENCE rho_Lambda ~ rho_crit, if node creation tracks the expansion. The passive attractor
      gives rho_Lambda/M_P^4 ~ H/M_P ~ 10^-61 (ONE power of H/M_P) -- still 61 orders too big. But a
      medium that CREATES nodes to fill expanding space cancels the leading dilution (-3H); what it
      cannot track is the CHANGE in the expansion rate, so the residual departure is set by Hdot over a
      relaxation time -- SECOND order, delta ~ (H tau)^2. Then rho_Lambda ~ M_P^4 (H/M_P)^2 = M_P^2 H^2 =
      rho_crit: exactly the observed "coincidence" that dark energy sits at the critical density. The
      extra power of H/M_P -- the whole remaining 61 orders -- is the acceleration lag of a tracking
      medium. This is conditional (efficient tracking is assumed, the coefficient is order-unity), but
      it is the first thing in the arc that lands on 10^-122 rather than 10^-61.
"""
from __future__ import annotations
import numpy as np

from test_vacuum_gravitates import equilibrium
from test_cc_offequilibrium import bulk_modulus

B = bulk_modulus()                                          # medium stiffness (from the real LJ condensate)
MP = 1.0                                                    # Planck units: M_Planck = 1, tau = 1
TAU = 1.0
H0 = 1.2e-61                                                # observed H_0 / M_Planck


def relax(delta0, H, tau=TAU, dt=0.02, T=30.0):
    """Integrate d(delta)/dt = -3H - delta/tau from delta0; return the trajectory sampled at a few t/tau."""
    n = int(T / dt)
    d = delta0
    out = {0.0: d}
    marks = [1.0, 3.0, 10.0, 30.0]
    for i in range(1, n + 1):
        d += dt * (-3 * H - d / tau)
        t = i * dt
        for m in marks:
            if abs(t - m) < dt / 2:
                out[m] = d
    return out


if __name__ == "__main__":
    print("=== Why so near equilibrium? Attractor, cosmic age, and rho_Lambda ~ rho_crit ===\n")

    # ---------- [A] the attractor erases initial conditions ----------
    H = 1e-3
    dss = -3 * H * TAU
    print("  [A] EQUILIBRIUM IS AN ATTRACTOR that erases initial conditions (H = 1e-3, tau = 1):")
    print(f"      {'t / tau =':>16} {'0':>10} {'1':>10} {'3':>10} {'10':>10} {'30':>10}")
    for d0 in (-0.5, 0.0, 0.3, 0.8):
        tr = relax(d0, H)
        row = " ".join(f"{tr[m]:>10.3e}" for m in (0.0, 1.0, 3.0, 10.0, 30.0))
        print(f"      delta0 = {d0:>+5.2f}   {row}")
    print(f"      => every initial departure (any initial vacuum energy) decays to the SAME attractor")
    print(f"         delta_ss = -3H*tau = {dss:.2e} within a few tau ~ a few Planck times. There is no")
    print("         initial cosmological constant to fine-tune: the medium forgets it in a Planck time.\n")

    # ---------- [B] Lambda decays as the universe ages ----------
    print("  [B] Lambda DECAYS as the universe AGES (attractor rho_Lambda ~ B*3H*tau ~ H):")
    print(f"      {'H / M_Planck':>16} {'rho_Lambda/M_P^4 (passive ~ H)':>32}")
    for Hf in (1e-1, 1e-3, 1e-6, 1e-9):
        print(f"      {Hf:>16.0e} {B * 3 * Hf * TAU:>32.2e}")
    print("      => the cosmological constant falls with H: it is small NOW because the universe is OLD")
    print("         (H << M_Planck). The smallness of Lambda is the smallness of H/M_Planck, not a tuning.\n")

    # ---------- [C] the coincidence: rho_Lambda ~ rho_crit if the medium tracks the expansion ----------
    print("  [C] THE COINCIDENCE rho_Lambda ~ rho_crit, if node creation tracks the expansion:")
    print("      passive medium (no tracking): residual delta ~ H*tau (linear), rho_Lambda ~ H")
    print("      tracking medium: node creation cancels -3H; residual is the acceleration lag,")
    print("      delta ~ (H*tau)^2 (second order), so rho_Lambda ~ H^2.")
    rho_crit = 3 * MP ** 2 * H0 ** 2                        # ~ critical density today (8piG/3 folded into O(1))
    print(f"      {'model':>18} {'delta_ss':>14} {'rho_Lambda/M_P^4 at H_0':>26} {'vs observed 1e-122':>20}")
    d_passive = 3 * H0 * TAU
    d_track = 3 * (H0 * TAU) ** 2
    print(f"      {'passive (~H)':>18} {d_passive:>14.1e} {B * d_passive:>26.1e} {B*d_passive/1e-122:>18.0e}x")
    print(f"      {'tracking (~H^2)':>18} {d_track:>14.1e} {B * d_track:>26.1e} {B*d_track/1e-122:>18.0e}x")
    print("      => the passive attractor overshoots by ~61 orders (one power of H/M_P); the tracking")
    print("         medium, whose residual is the SECOND-order acceleration lag, gives rho_Lambda ~ H^2 =")
    print("         M_P^2 H^2 = rho_crit -- the observed coincidence that dark energy sits at the critical")
    print("         density, to order unity. The missing 61 orders ARE the extra power of H/M_P.\n")

    print("[verdict] 'why so near equilibrium' has a three-part answer, two solid and one conditional:")
    print("  * ATTRACTOR (solid): the relaxation drives ANY initial state to delta ~ H*tau in a few")
    print("    Planck times, so there is no initial cosmological constant to tune -- the initial-conditions")
    print("    version of the CC problem does not arise. Whatever we observe is the attractor, not a relic.")
    print("  * COSMIC AGE (solid): the attractor value tracks H, so Lambda decays as the universe expands")
    print("    and is small today because the universe is old (H << M_Planck). The smallness IS H/M_Planck.")
    print("  * COINCIDENCE (conditional, and the new lead): a medium that creates nodes to track the")
    print("    expansion cancels the linear dilution; its residual is the acceleration lag, delta ~")
    print("    (H*tau)^2, giving rho_Lambda ~ M_P^2 H^2 = rho_crit -- exactly the observed dark-energy")
    print("    density, with the whole 61-order gap supplied by the second power of H/M_Planck. This")
    print("    assumes efficient tracking and fixes only the scaling (order-unity coefficient), but it is")
    print("    the first mechanism in the arc that reaches 10^-122, and it predicts DYNAMICAL dark energy")
    print("    tracking the critical density -- the natural next thing to pin down (the tracking")
    print("    efficiency, and the coefficient) and to confront with the evolving-w data.")
