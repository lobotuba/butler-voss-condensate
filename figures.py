"""
figures.py -- render the project's headline results to PNGs (figures/).

Regenerates the communicable highlights:
  fig1  self-assembly of the medium (H10): disordered cloud -> hexagonal lattice
  fig2  a topological vortex (H6): |psi| core + phase winding = charge
  fig3  GRAVITY (3d): two masses drift together with coupling, apart without it
"""
from __future__ import annotations
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bvc_core import lj_forces_energy, relax_medium, coordination
import integration_phase3_variational as V
from prototype_complex import ComplexFabric

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)


# ---- fig1: self-assembly (H10) --------------------------------------------
def relax_snapshots(N=300, steps=6000, dt=0.005, cool=0.999, seed=0):
    rng = np.random.default_rng(seed); Rd = math.sqrt(N * 1.7 / math.pi); X = []
    while len(X) < N:
        p = (rng.random(2) * 2 - 1) * Rd
        if math.hypot(*p) <= Rd and all(math.hypot(*(p - q)) > 0.95 for q in X):
            X.append(p)
    X = np.array(X); X0 = X.copy(); Vv = np.zeros_like(X); F, _ = lj_forces_energy(X)
    for k in range(steps):
        X += Vv * dt + 0.5 * F * dt ** 2
        Fn, _ = lj_forces_energy(X); Vv += 0.5 * (F + Fn) * dt; Vv *= cool; F = Fn
    return X0 - X0.mean(0), X - X.mean(0)


def fig1():
    X0, Xf = relax_snapshots()
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    ax[0].scatter(X0[:, 0], X0[:, 1], s=16, c="steelblue")
    ax[0].set_title("initial: disordered cloud")
    sc = ax[1].scatter(Xf[:, 0], Xf[:, 1], s=16, c=coordination(Xf), cmap="viridis")
    ax[1].set_title("self-assembled: hexagonal, spacing set by the medium (H10)")
    for a in ax:
        a.set_aspect("equal"); a.axis("off")
    fig.colorbar(sc, ax=ax[1], fraction=0.046, label="neighbours")
    plt.tight_layout(); plt.savefig(f"{OUT}/fig1_self_assembly.png", dpi=110); plt.close()


# ---- fig2: topological vortex = charge (H6) -------------------------------
def fig2():
    cf = ComplexFabric(rows=72, cols=72, potential="mexicanhat")
    cf.set_vortices([(0.5, 0.5, 1)])
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    s0 = ax[0].scatter(cf.pos[:, 0], cf.pos[:, 1], c=np.abs(cf.psi), s=7, cmap="inferno")
    ax[0].set_title("|psi|  (amplitude dips to 0 at the vortex core)")
    s1 = ax[1].scatter(cf.pos[:, 0], cf.pos[:, 1], c=np.angle(cf.psi), s=7, cmap="hsv")
    ax[1].set_title("phase arg(psi): the 2pi winding = quantized charge (H6)")
    for a in ax:
        a.set_aspect("equal"); a.axis("off")
    fig.colorbar(s0, ax=ax[0], fraction=0.046)
    fig.colorbar(s1, ax=ax[1], fraction=0.046, label="phase")
    plt.tight_layout(); plt.savefig(f"{OUT}/fig2_vortex_charge.png", dpi=110); plt.close()


# ---- fig3: GRAVITY -- two masses drift together (3d) ----------------------
def seed_two(f, d, amp=1.5, w=2.0):
    x = f.X
    f.u = amp * (np.exp(-(((x[:, 0] - d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)) +
                 np.exp(-(((x[:, 0] + d/2) ** 2) + x[:, 1] ** 2) / (2 * w ** 2)))
    f.pi = np.zeros(f.N)


def separation(f):
    e = 0.5 * f.pi ** 2 + 0.5 * f.m2 * f.u ** 2
    L = f.X[:, 0] < 0; R = ~L
    return float(np.linalg.norm((f.X[R] * e[R, None]).sum(0) / e[R].sum() -
                                (f.X[L] * e[L, None]).sum(0) / e[L].sum()))


def fig3():
    cloud = relax_medium(N=300, seed=3)
    plt.figure(figsize=(7.5, 5))
    for beta, col, lab in [(0, "crimson", "β=0  (no coupling)"),
                           (60, "steelblue", "β=60 (gravity)")]:
        f = V.VariationalCoupled(cloud, beta=beta, m2=1.0, g_min=0.02, damping=1.0, dt=0.001)
        seed_two(f, 9.0)
        ts, ss = [0.0], [separation(f)]
        for k in range(8000):
            f.step()
            if (k + 1) % 400 == 0:
                ts.append(f.time); ss.append(separation(f))
        plt.plot(ts, ss, color=col, lw=2.2, label=lab)
    plt.xlabel("time"); plt.ylabel("separation of the two masses")
    plt.title("Emergent gravity: two masses drift together (Phase 3d)")
    plt.legend(); plt.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(f"{OUT}/fig3_gravity_separation.png", dpi=110); plt.close()


if __name__ == "__main__":
    print("rendering figures ->", OUT)
    fig1(); print("  fig1 self-assembly done")
    fig2(); print("  fig2 vortex/charge done")
    fig3(); print("  fig3 gravity done")
    print("done.")
