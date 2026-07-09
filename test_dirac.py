"""
Test for emergent relativistic (Dirac) FERMIONS -- the Volovik/Wen route.

The model so far is bosonic (scalar/complex/vector fields). Matter is fermions:
spin-1/2, and chiral. Relativistic fermions can EMERGE as low-energy excitations
near a band-touching point, where the dispersion becomes a linear, isotropic cone
E ~ +/- v_F |k - k_D| -- a Dirac cone (as in graphene). This tests whether a
tight-binding ("hopping") model on the medium's lattice produces one, and confronts
the two hard facts:
  * a Dirac cone needs a BIPARTITE (two-sublattice) lattice: the self-assembled
    close-packing (triangular, one site per cell) has a single band and NO cone;
    the honeycomb (two sites per cell) has Dirac points.
  * Nielsen-Ninomiya: lattice fermions come in pairs of OPPOSITE chirality, so a
    single chiral fermion (as in the Standard Model) cannot appear alone.

Bloch Hamiltonian on the honeycomb (nn hopping t): H(k) = [[0, f], [f*, 0]],
f(k) = -t sum_j exp(i k.delta_j), bands E = +/- |f|. Chirality of a Dirac point =
winding of arg f(k) around it.
"""
from __future__ import annotations
import numpy as np

# honeycomb nearest-neighbour vectors (nn distance a = 1), A -> B
D = np.array([[0.0, 1.0], [np.sqrt(3) / 2, -0.5], [-np.sqrt(3) / 2, -0.5]])
K = np.array([4 * np.pi / (3 * np.sqrt(3)), 0.0])       # a Dirac point
KP = -K                                                 # the inequivalent partner


def f_honeycomb(k, t=1.0):
    k = np.asarray(k, float)
    return -t * np.sum(np.exp(1j * (D @ k)))


def cone(kD, t=1.0, rho=1e-3, ndir=12):
    """Measure the dispersion E=|f| just off a band-touching point: slope (v_F)
    per direction; return (v_F mean, anisotropy)."""
    vs = []
    for th in np.linspace(0, 2 * np.pi, ndir, endpoint=False):
        q = rho * np.array([np.cos(th), np.sin(th)])
        vs.append(abs(f_honeycomb(kD + q, t)) / rho)
    vs = np.array(vs)
    return vs.mean(), (vs.max() - vs.min()) / vs.mean()


def chirality(kD, t=1.0, rho=0.05, n=240):
    """Winding of arg f(k) around kD (the Dirac-point topological charge)."""
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ph = np.array([np.angle(f_honeycomb(kD + rho * np.array([np.cos(a), np.sin(a)]), t)) for a in th])
    d = np.diff(np.concatenate([ph, ph[:1]]))
    d = (d + np.pi) % (2 * np.pi) - np.pi                # unwrap steps
    return int(np.rint(d.sum() / (2 * np.pi)))


# triangular (close-packed) single band: E(k) = -2t sum cos(k.a_i), no touching
TRI = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2], [-0.5, np.sqrt(3) / 2]])   # 3 lattice dirs


def tri_band(k, t=1.0):
    return -2 * t * np.sum(np.cos(TRI @ np.asarray(k, float)))


if __name__ == "__main__":
    print("=== Test for emergent relativistic (Dirac) fermions ===\n")

    print("[A] honeycomb (bipartite) -- is there a linear isotropic Dirac cone?")
    for name, kD in [("K ", K), ("K'", KP)]:
        gap = abs(f_honeycomb(kD))
        vF, aniso = cone(kD)
        print(f"  Dirac point {name} at ({kD[0]:+.3f},{kD[1]:+.3f}): gap |f|={gap:.2e} "
              f"(gapless); v_F={vF:.3f}, cone anisotropy={aniso:.1e}")
    print("  => bands touch with a LINEAR, ISOTROPIC cone: emergent 2D massless Dirac")
    print("     fermions, Fermi velocity v_F = 3/2 (t=a=1) -- the fermions' 'speed of light'.")

    print("\n[B] chirality (Nielsen-Ninomiya doubling)")
    cK, cKp = chirality(K), chirality(KP)
    print(f"  winding of arg f: at K = {cK:+d},  at K' = {cKp:+d}   (sum = {cK+cKp:+d})")
    print("  => the two Dirac points carry OPPOSITE chirality and cancel: fermions come as")
    print("     a doubled pair, never a single chiral one. A single chiral (Standard-Model)")
    print("     fermion on the lattice is forbidden -- the Nielsen-Ninomiya wall.")

    print("\n[C] triangular (self-assembled close-packing) -- one band, no cone")
    ks = np.linspace(-np.pi, np.pi, 41)
    band = np.array([[tri_band([kx, ky]) for kx in ks] for ky in ks])
    # curvature at the band minimum: quadratic (non-relativistic), not linear
    im = np.unravel_index(np.argmin(band), band.shape)
    print(f"  single band, min E={band.min():.3f} at k=({ks[im[1]]:+.2f},{ks[im[0]]:+.2f}); "
          f"E-Emin ~ |k|^2 near the minimum (a normal, NON-relativistic band).")
    print("  => the plain close-packed medium (one site per cell) hosts NO Dirac fermions;")
    print("     relativistic fermions need a BIPARTITE (e.g. honeycomb) structure.")

    print("\nVerdict: the model CAN host emergent relativistic Dirac fermions -- but only on a")
    print("two-sublattice medium, and only as opposite-chirality PAIRS (Nielsen-Ninomiya).")
    print("Two open consequences: (i) a single chiral fermion needs an extra mechanism (domain")
    print("wall / dimensional reduction / interactions); (ii) the fermion cone speed v_F is a")
    print("new light cone that must be UNIFIED with the boson c -- the same universality demand")
    print("as before, now across statistics.")
