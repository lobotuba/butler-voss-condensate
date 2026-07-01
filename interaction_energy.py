"""
interaction_energy.py -- quantify the gravitational CONTACT STRENGTH via the
adiabatic interaction energy U(d) (the clean static force method).

The dynamic probes (drift of dispersing lumps) all failed -- swamped by
dispersion. The right tool is static: freeze two gentle masses at separation d,
relax the medium, take the total energy; subtract the FROZEN-medium baseline
(which is only the direct field overlap) to isolate the medium-mediated
(gravitational) interaction U(d).  Then F(d) = -dU/dd.

Result (N=400, beta=60): a short-range attractive WELL --
  contact binding |U| ~= 0.92 at d=4, peak force ~= 0.5 at d~=5, vanishing by
  d~=6 (screening range ~5-6 = the mass size). NOT a long-range power law.
Saves figures/fig4_interaction_energy.png.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bvc_core import relax_medium
import integration_phase3_variational as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)


def two_lumps(X, d, amp=0.8, w=2.5):
    x = X
    return amp * (np.exp(-(((x[:, 0] - d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)) +
                  np.exp(-(((x[:, 0] + d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)))


def E_config(Xref, d, relax, beta=60.0, steps=5000, dt=0.002):
    f = V.VariationalCoupled(Xref, beta=beta, m2=1.0, g_min=0.02, damping=1.0, dt=dt)
    f.u = two_lumps(Xref, d); f.pi = np.zeros(f.N)
    if relax:                                    # relax medium (field frozen) to minimum
        for _ in range(steps):
            _, FX, _ = f.forces(); f.Wn += 0.5*FX*dt; f.X += f.Wn*dt
            _, FX, _ = f.forces(); f.Wn += 0.5*FX*dt; f.Wn *= 0.97
    return f.energy()


def main():
    Xref = relax_medium(N=400, seed=3, steps=6000)
    ds = [4, 5, 6, 7, 8, 10]
    resp = {}
    print("d   E_frozen    E_relaxed   response")
    for d in ds:
        Ef = E_config(Xref, d, False); Er = E_config(Xref, d, True)
        resp[d] = Er - Ef
        print(f"{d:2d}  {Ef:10.3f} {Er:10.3f}   {resp[d]:+.4f}")
    dfar = ds[-1]
    U = {d: resp[d] - resp[dfar] for d in ds}    # medium-mediated interaction energy
    print("\ngravitational interaction U(d):")
    for d in ds:
        print(f"  d={d:2d}: U={U[d]:+.4f}")
    print("\nforce F(d) = -dU/dd:")
    for i in range(1, len(ds)-1):
        F = -(U[ds[i+1]] - U[ds[i-1]]) / (ds[i+1] - ds[i-1])
        print(f"  d={ds[i]:2d}: F={F:+.4f}")
    print(f"\ncontact binding |U| at d={ds[0]}: {abs(U[ds[0]]):.4f};  "
          f"U(d>=6) ~ 0 => screening range ~5-6")

    # plot the potential well
    dd = np.array(ds); uu = np.array([U[d] for d in ds])
    plt.figure(figsize=(7, 5))
    plt.axhline(0, color="gray", lw=0.8)
    plt.plot(dd, uu, "o-", color="darkviolet", lw=2)
    plt.xlabel("separation d"); plt.ylabel("interaction energy U(d)")
    plt.title("Gravity is short-range: attractive well, screened by d≈6")
    plt.annotate(f"contact binding ≈ {abs(uu[0]):.2f}", xy=(dd[0], uu[0]),
                 xytext=(dd[0]+1.2, uu[0]+0.15),
                 arrowprops=dict(arrowstyle="->", color="black"))
    plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"{OUT}/fig4_interaction_energy.png", dpi=110); plt.close()
    print(f"saved {OUT}/fig4_interaction_energy.png")


if __name__ == "__main__":
    main()
