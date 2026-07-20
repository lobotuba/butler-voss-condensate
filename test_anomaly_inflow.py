"""
Chirality without inconsistency: the anomaly is quantized, and the bulk supplies it (Callan-Harvey).

test_domain_wall showed the model can carry a SINGLE chiral fermion on a domain wall, evading
Nielsen-Ninomiya. test_yang_mills showed the fermion loop induces genuine non-Abelian Yang-Mills.
The Standard Model needs both AT ONCE, and that is where the real obstruction lives: the SM's SU(2)_L
couples CHIRALLY, and a chiral gauge theory is INCONSISTENT unless its anomalies cancel. A chiral
fermion's gauge current is not conserved in a background field -- the anomaly -- so a lone chiral
fermion coupled to a gauge field is not a theory at all.

The domain-wall construction survives this, but only by a specific mechanism: CALLAN-HARVEY ANOMALY
INFLOW. Each wall is INDIVIDUALLY anomalous; the charge it appears to lose is supplied by the BULK,
whose Chern-Simons response pumps charge from one wall to the other. The full system is consistent
because the lattice as a whole is vector-like. The content is quantitative and QUANTIZED:

    (bulk Chern number C)  =  (number of chiral modes per wall)  =  (the anomaly coefficient)
                           =  (charge pumped per flux quantum, Laughlin).

This file measures all three and checks they agree, with a trivial-phase control.

Model (the same Wilson-Dirac / QWZ system as test_domain_wall):
    H(k) = sin(kx) sx + sin(ky) sy + (M - cos kx - cos ky) sz,
topological (|C| = 1) for 0 < |M| < 2, trivial (C = 0) for |M| > 2. On a strip (kx a good quantum
number, y open) each edge is a domain wall to the trivial vacuum and binds one chiral branch.

What is measured:
  [A] the bulk anomaly coefficient, as the Chern number C = (1/4pi) int dhat . (d_kx dhat x d_ky dhat):
      -1 in the topological phase, 0 in the trivial phase.
  [B] the wall content, as the NET SIGNED zero-energy crossings of the edge-localized branches across
      the Brillouin zone -- the spectral flow, i.e. the charge pumped onto that wall per flux quantum:
      +1 on one wall, -1 on the other; 0,0 in the trivial control.
  [C] consistency: the two wall anomalies SUM TO ZERO (the lattice is vector-like, as
      Nielsen-Ninomiya requires), while each wall separately is anomalous with |coefficient| = |C|.
      That is anomaly inflow: neither wall is a consistent theory alone; the pair plus the bulk is.

HONEST SCOPE. This shows the model can host a CHIRAL fermion consistently -- its anomaly is real,
quantized, and exactly supplied by bulk inflow -- which is the mechanism Standard-Model chirality
requires. It does NOT derive the Standard Model's own anomaly cancellation. The SM is a standalone
FOUR-dimensional chiral gauge theory: its anomalies cancel among its OWN fermion content (the
quark/lepton hypercharge conspiracy, sum of Y = 0 and sum of Y^3 = 0 per generation), with no bulk to
lean on. Here the bulk is doing the work, so the wall theory is anomaly-free only together WITH the
bulk. Getting a standalone anomaly-free chiral spectrum -- the actual SM fermion content and
hypercharges -- is not attempted and is not fixed by the medium. What is established is the weaker
but necessary statement: chirality plus gauge fields is CONSISTENTLY realizable here, by inflow, with
a quantized coefficient.
"""
from __future__ import annotations
import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def chern(M, N=200):
    """Bulk Chern number = the anomaly coefficient, C = (1/4pi) int dhat.(dx dhat x dy dhat)."""
    k = np.linspace(0, 2 * np.pi, N, endpoint=False)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    dx, dy, dz = np.sin(KX), np.sin(KY), (M - np.cos(KX) - np.cos(KY))
    n = np.sqrt(dx * dx + dy * dy + dz * dz)
    d = np.stack([dx / n, dy / n, dz / n], -1)
    dX = np.gradient(d, k, axis=0)
    dY = np.gradient(d, k, axis=1)
    integ = np.einsum("ijk,ijk->ij", d, np.cross(dX, dY))
    return float(integ.sum() * (k[1] - k[0]) ** 2 / (4 * np.pi))


def strip_H(kx, Ny, M):
    """Wilson-Dirac strip at momentum kx (y open, Ny sites) -- as in test_domain_wall."""
    on = np.sin(kx) * sx + (M - np.cos(kx)) * sz
    hop = -0.5 * sz - 0.5j * sy
    H = np.zeros((2 * Ny, 2 * Ny), dtype=complex)
    for j in range(Ny):
        H[2 * j:2 * j + 2, 2 * j:2 * j + 2] = on
        if j + 1 < Ny:
            H[2 * j:2 * j + 2, 2 * j + 2:2 * j + 4] = hop
            H[2 * j + 2:2 * j + 4, 2 * j:2 * j + 2] = hop.conj().T
    return H


def edge_branches(M, Ny=80, NK=400, gap=0.9):
    """Energies of the in-gap edge-localized branches vs kx, for the bottom and top walls.
    The zone is scanned symmetrically about kx = 0 so the crossing is interior (the branch merges
    into the bulk near the zone boundary, so a [0, 2pi) scan would split it across the wrap-around)."""
    ks = np.linspace(-np.pi, np.pi, NK)
    Eb = np.full(NK, np.nan)
    Et = np.full(NK, np.nan)
    ys = np.arange(Ny)
    for i, kx in enumerate(ks):
        w, v = np.linalg.eigh(strip_H(kx, Ny, M))
        for idx in np.where(np.abs(w) < gap)[0]:
            dens = (np.abs(v[:, idx].reshape(Ny, 2)) ** 2).sum(1)
            if dens.max() < 0.05:                      # bulk state, not wall-bound
                continue
            ybar = (ys * dens).sum()
            if ybar < Ny / 3:
                if np.isnan(Eb[i]) or abs(w[idx]) < abs(Eb[i]):
                    Eb[i] = w[idx]
            elif ybar > 2 * Ny / 3:
                if np.isnan(Et[i]) or abs(w[idx]) < abs(Et[i]):
                    Et[i] = w[idx]
    return ks, Eb, Et


def spectral_flow(ks, E):
    """Net SIGNED zero-energy crossings of a branch = charge pumped onto that wall per flux quantum."""
    net = 0
    for i in range(len(ks) - 1):
        a, b = E[i], E[i + 1]
        if np.isnan(a) or np.isnan(b):
            continue
        if a < 0 <= b:
            net += 1
        elif a > 0 >= b:
            net -= 1
    return net


if __name__ == "__main__":
    print("=== Chirality without inconsistency: a quantized anomaly, supplied by bulk inflow ===\n")
    print("  A chiral gauge theory is inconsistent unless its anomalies cancel. On a domain wall the")
    print("  cancellation is Callan-Harvey INFLOW: each wall is anomalous, the bulk supplies the")
    print("  charge. The check: (bulk Chern C) = (chiral modes per wall) = (charge pumped/flux).\n")

    print(f"  {'phase':>18} {'M':>5} {'bulk Chern C':>13} {'wall A flow':>12} {'wall B flow':>12} "
          f"{'sum':>5} {'|C| = |flow|?':>14}")
    for name, M in (("TOPOLOGICAL", 1.0), ("trivial (control)", 3.0)):
        C = chern(M)
        ks, Eb, Et = edge_branches(M)
        nb, nt = spectral_flow(ks, Eb), spectral_flow(ks, Et)
        ok = abs(round(C)) == abs(nb) == abs(nt)
        print(f"  {name:>18} {M:>5.1f} {C:>13.4f} {nb:>+12d} {nt:>+12d} {nb+nt:>+5d} {str(ok):>14}")

    print("\n  [A] the bulk Chern number is the ANOMALY COEFFICIENT: -1.000 in the topological phase,")
    print("      0.000 in the trivial phase -- quantized, and it is what the bulk can pump.")
    print("  [B] each wall carries exactly ONE chiral branch: its net signed zero-crossing across the")
    print("      zone is +1 on one wall and -1 on the other. That spectral flow IS the charge pumped")
    print("      onto that wall per flux quantum -- i.e. that wall's gauge charge is NOT conserved.")
    print("      Each wall, by itself, is ANOMALOUS. (Trivial control: no crossings, no anomaly.)")
    print("  [C] the two flows SUM TO ZERO. The lattice as a whole is vector-like and anomaly-free,")
    print("      exactly as Nielsen-Ninomiya requires -- so the charge a wall loses is not destroyed,")
    print("      it is pumped through the BULK to the other wall. Anomaly inflow, and |C| equals the")
    print("      per-wall coefficient, so the three numbers agree.\n")

    print("[verdict] chirality is CONSISTENTLY realizable in this medium, by anomaly inflow:")
    print("  * The anomaly is not a pathology here but a quantized bookkeeping identity: the bulk")
    print("    Chern number, the number of chiral modes per wall, and the charge pumped per flux")
    print("    quantum are ONE integer. Each wall is individually anomalous; the wall pair plus the")
    print("    bulk is exactly anomaly-free. This is the mechanism Standard-Model chirality needs,")
    print("    and it works in the model's own domain-wall construction (test_domain_wall) alongside")
    print("    the induced non-Abelian gauge fields (test_yang_mills).")
    print("  * HONEST ceiling: this is NOT the Standard Model's anomaly cancellation. The SM is a")
    print("    standalone FOUR-dimensional chiral gauge theory whose anomalies cancel among its OWN")
    print("    fermion content -- the quark/lepton hypercharge conspiracy (sum Y = 0, sum Y^3 = 0 per")
    print("    generation) -- with no bulk to lean on. Here the BULK does the cancelling, so the wall")
    print("    theory is anomaly-free only together with it. Producing a standalone anomaly-free")
    print("    chiral spectrum -- the actual SM fermion content, representations and hypercharges --")
    print("    is not attempted and is not fixed by the medium. Chirality is shown CONSISTENT; the")
    print("    Standard Model's particular chiral content remains an input.")
