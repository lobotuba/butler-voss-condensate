"""
Escaping the chirality wall: a single chiral fermion on a domain wall.

test_dirac.py hit Nielsen-Ninomiya: lattice fermions come as opposite-chirality
pairs, so a single chiral (Standard-Model) fermion seems forbidden. The standard
escape (Kaplan domain-wall fermions; Callan-Harvey) is to regularise the fermion
as a 2-band Wilson-Dirac (Chern) insulator and put it on a DOMAIN WALL: a single
chiral mode binds to the wall, while its opposite-chirality partner lives on the
OTHER wall, spatially separated. Locally, on one wall, there is one chiral fermion.

Wilson-Dirac (QWZ) Bloch Hamiltonian:
  H(k) = sin(k_x) sx + sin(k_y) sy + (M - cos k_x - cos k_y) sz,
topological (Chern +/-1) for 0<|M|<2, trivial for |M|>2. On a strip (k_x good, y
open) the topological phase has ONE chiral branch crossing the gap on EACH edge,
with OPPOSITE velocities (chiralities). An edge is a domain wall to the trivial
vacuum outside, so each edge hosts a single chiral fermion.
"""
from __future__ import annotations
import numpy as np

sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)


def strip_H(kx, Ny, M):
    """2Ny x 2Ny Wilson-Dirac strip Hamiltonian at momentum kx (y open, 1..Ny)."""
    on = np.sin(kx) * sx + (M - np.cos(kx)) * sz          # on-site
    hop = -0.5 * sz - 0.5j * sy                           # y-hop j -> j+1
    H = np.zeros((2 * Ny, 2 * Ny), complex)
    for j in range(Ny):
        H[2*j:2*j+2, 2*j:2*j+2] = on
        if j + 1 < Ny:
            H[2*j:2*j+2, 2*j+2:2*j+4] = hop
            H[2*j+2:2*j+4, 2*j:2*j+2] = hop.conj().T
    return H


def near_zero(Ny, M, kx, n=2):
    """The n states nearest E=0: list of (E, <y>, peak-density, dens)."""
    w, v = np.linalg.eigh(strip_H(kx, Ny, M))
    ys = np.arange(Ny); out = []
    for i in np.argsort(np.abs(w))[:n]:
        dens = (np.abs(v[:, i].reshape(Ny, 2)) ** 2).sum(1)
        out.append((w[i], (ys * dens).sum(), dens.max(), dens))
    return out


def report(name, M, Ny=60):
    print(f"\n[{name}]  M = {M}  (Wilson-Dirac strip, Ny={Ny})")
    em = near_zero(Ny, M, 0.3); emL = near_zero(Ny, M, 0.3 + 1e-3)
    print(f"  {'E(kx=0.3)':>10} {'velocity':>9} {'<y>':>6} {'peak dens':>10} {'classification':>16}")
    for (E, y, pk, dens), (E2, *_ ) in zip(em, emL):
        localized = pk > 0.1 and (y < 4 or y > Ny - 5)   # concentrated near a boundary
        if localized:
            vel = (E2 - E) / 1e-3
            cls = f"CHIRAL edge ({'bottom' if y < Ny/2 else 'top'})"
            print(f"  {E:>+10.4f} {vel:>+9.3f} {y:>6.1f} {pk:>10.3f} {cls:>16}")
        else:
            print(f"  {E:>+10.4f} {'--':>9} {y:>6.1f} {pk:>10.3f} {'bulk (no edge)':>16}")


if __name__ == "__main__":
    print("=== A single chiral fermion on a domain wall (Callan-Harvey / Kaplan) ===")
    report("TOPOLOGICAL", M=1.0)
    report("trivial", M=3.0)
    print("\nVerdict:")
    print("  TOPOLOGICAL (0<|M|<2): two in-gap branches CROSS E=0, one localized on each edge")
    print("  with OPPOSITE velocity (chirality). Each edge = a domain wall to the trivial")
    print("  vacuum, and hosts a SINGLE chiral fermion; the Nielsen-Ninomiya partner is the")
    print("  opposite-chirality mode on the OTHER edge -- spatially separated, not on this wall.")
    print("  trivial (|M|>2): gapped, NO in-gap crossing -- no chiral fermion.")
    print("  => the chirality wall is EVADED the standard way: a lattice that is vector-like")
    print("     overall can still carry a single chiral fermion on a domain wall. This is the")
    print("     mechanism a fundamental version of the model would use for Standard-Model chirality.")
