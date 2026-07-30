"""
Where the gauge group and the generations come from: discrete topology of the emergent band structure.

The anomaly result (test_anomaly_hypercharge) fixed the hypercharges but left the deepest inputs open:
why the gauge group SU(3)xSU(2)xU(1), why the fermion representations, why THREE generations. Deriving
the specific Standard-Model values from a specific lattice is a genuine open problem and is not attempted
here. What this file establishes is the structural fact that reframes the question: in an emergent-fermion
medium these are not continuous parameters to be fine-tuned but DISCRETE data -- symmetry and topology of
the band structure -- so they are quantized (small integers, specific compact groups) by their nature.

Two pieces, both computed on the medium's own kind of band structure:

  [A] THE FLAVOUR SYMMETRY IS A COUNT OF DIRAC POINTS. On a bipartite lattice the low-energy fermions
      live at isolated Dirac points, and the number of them is fixed by the lattice (fermion doubling),
      not tunable. The honeycomb has exactly two (the valleys K, K'), so its emergent low-energy theory
      carries a two-fold flavour multiplet and an internal symmetry rotating it -- the seed of a
      non-Abelian flavour/gauge symmetry, which the fermion loop (test_yang_mills) then gauges. The
      GROUP is the symmetry of the degenerate multiplet; the multiplet size is a lattice output.

  [B] THE GENERATION NUMBER IS A TOPOLOGICAL INDEX. The number of chiral fermion families localised on a
      domain wall equals the jump in the bulk Chern number across it (bulk-boundary correspondence, the
      mechanism of test_domain_wall). The Chern number is an INTEGER, computed here from the Berry
      curvature of a Wilson-Dirac band, that jumps only when the gap closes -- so the number of chiral
      generations is quantized and robust, not a continuous parameter. Why exactly three is then a
      discrete topological property of the physical medium's band structure; that it is a small integer
      at all is because it is a Chern index.
"""
from __future__ import annotations
import numpy as np

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


# ------------------------------------------------------------- [A] Dirac points of the honeycomb ------
def honeycomb_f(kx, ky):
    """Nearest-neighbour structure factor f(k) = sum_delta exp(i k.delta); Dirac points are its zeros."""
    d1 = np.array([1.0, 0.0])
    d2 = np.array([-0.5, np.sqrt(3) / 2])
    d3 = np.array([-0.5, -np.sqrt(3) / 2])
    return (np.exp(1j * (kx * d1[0] + ky * d1[1])) + np.exp(1j * (kx * d2[0] + ky * d2[1]))
            + np.exp(1j * (kx * d3[0] + ky * d3[1])))


def count_dirac_points(N=600):
    """Count INEQUIVALENT zeros of |f(k)| in one Brillouin zone. Sample exactly one reciprocal-lattice
    primitive cell (U, V in [0,1) along b1, b2) so each inequivalent Dirac point appears once."""
    a1 = np.array([1.5, np.sqrt(3) / 2])                             # honeycomb Bravais vectors
    a2 = np.array([1.5, -np.sqrt(3) / 2])                            # (a1 = d1-d3, a2 = d1-d2)
    B = 2 * np.pi * np.linalg.inv(np.array([a1, a2])).T              # reciprocal vectors (rows b1, b2)
    b1, b2 = B[0], B[1]
    u = np.linspace(0, 1, N, endpoint=False)
    U, V = np.meshgrid(u, u, indexing="ij")
    KX = U * b1[0] + V * b2[0]
    KY = U * b1[1] + V * b2[1]
    lo = np.abs(honeycomb_f(KX, KY)) < 0.05
    # cluster the sub-threshold grid points (in crystal coords, with periodic wrap) into distinct valleys
    coords = [(U[i, j], V[i, j]) for i, j in zip(*np.where(lo))]
    clusters = []
    for (cu, cv) in coords:
        if not any(min(abs(cu - qu), 1 - abs(cu - qu)) < 0.15 and min(abs(cv - qv), 1 - abs(cv - qv)) < 0.15
                   for (qu, qv) in clusters):
            clusters.append((cu, cv))
    return len(clusters)


# ------------------------------------------------------------- [B] Chern number of a Wilson-Dirac band -
def dvec(kx, ky, m0):
    return np.stack([np.sin(kx), np.sin(ky), m0 + 2 - np.cos(kx) - np.cos(ky)], axis=-1)


def chern_number(m0, N=200):
    """C = (1/4 pi) integral d_hat . (d_kx d_hat x d_ky d_hat) over the BZ, for H = d.sigma."""
    g = np.linspace(0, 2 * np.pi, N, endpoint=False)
    KX, KY = np.meshgrid(g, g, indexing="ij")
    d = dvec(KX, KY, m0)
    dh = d / np.linalg.norm(d, axis=-1, keepdims=True)
    dx = np.gradient(dh, g, axis=0)
    dy = np.gradient(dh, g, axis=1)
    integrand = np.sum(dh * np.cross(dx, dy), axis=-1)
    return float(np.sum(integrand) * (g[1] - g[0]) ** 2 / (4 * np.pi))


if __name__ == "__main__":
    print("=== Gauge group and generations from the discrete topology of the band structure ===\n")

    # ---------- [A] flavour symmetry = number of Dirac points ----------
    nD = count_dirac_points()
    print("  [A] THE FLAVOUR MULTIPLET IS A COUNT OF DIRAC POINTS (fixed by the lattice, not tunable):")
    print(f"      honeycomb (two sublattices): number of Dirac points in the BZ = {nD}  (the valleys K, K')")
    print("      => the low-energy theory carries a two-fold flavour multiplet; the internal symmetry")
    print("         rotating the degenerate valleys is the seed of a non-Abelian flavour/gauge symmetry,")
    print("         which the fermion loop (test_yang_mills) gauges. A richer lattice (more sublattices")
    print("         or orbitals) yields more Dirac points and a larger multiplet: the GROUP is the")
    print("         symmetry of the multiplet, and the multiplet SIZE is a discrete lattice output.\n")

    # ---------- [B] generation number = Chern index ----------
    print("  [B] THE GENERATION NUMBER IS A TOPOLOGICAL INDEX (a Chern integer of the band):")
    print("      Wilson-Dirac band d = (sin kx, sin ky, m0 + 2 - cos kx - cos ky); Chern number C(m0):")
    print(f"      {'m0':>8} {'Chern C':>12} {'phase':>28}")
    for m0 in (1.0, -1.0, -3.0, -5.0):
        C = chern_number(m0)
        phase = ("trivial (C=0)" if abs(C) < 0.3 else f"topological (C={round(C):+d})")
        print(f"      {m0:>8.1f} {C:>12.3f} {phase:>28}")
    print("      => C is an INTEGER that jumps only when the band gap closes. By bulk-boundary")
    print("         correspondence (test_domain_wall) the number of chiral fermion families on a domain")
    print("         wall equals the jump in C across it -- so a wall between the C=+1 and C=-1 regions")
    print("         carries |(+1)-(-1)| = 2 chiral families, and a wall against a trivial region carries")
    print("         one. Higher-Chern bands (longer-range hopping) give more. The generation count is")
    print("         therefore a quantized topological integer, not a continuous parameter.\n")

    print("[verdict] the gauge group and the generation number are discrete band-structure data, not")
    print("  continuous fine-tunings -- which is why they are specific compact groups and small integers:")
    print("  * The emergent flavour multiplet is the set of Dirac points of the lattice (two for the")
    print("    honeycomb), and the emergent gauge group is the symmetry that rotates that multiplet,")
    print("    gauged by the fermion loop. Both are outputs of the lattice's band structure [A].")
    print("  * The number of chiral generations is a Chern index -- an integer, robust, changing only at")
    print("    a gap closing -- so it is quantized by topology [B]. 'Why a small integer' is answered:")
    print("    because it is a topological invariant; 'why exactly three' is a discrete property of the")
    print("    physical medium's band structure.")
    print("  * HONEST scope: this does NOT derive SU(3)xSU(2)xU(1), the specific representations, or the")
    print("    number three. Those require the physical medium's actual lattice, which is not fixed here.")
    print("    What is established is the CATEGORY of the answer: the Standard Model's gauge group and")
    print("    generation count are symmetry and topology of the emergent fermion band structure --")
    print("    discrete, quantized data -- so they are the right KIND of object (compact groups, integer")
    print("    families) rather than arbitrary continuous inputs. With the hypercharges fixed by anomaly")
    print("    cancellation (test_anomaly_hypercharge), what remains input is which lattice the medium")
    print("    realizes -- and that is now a question about one discrete structure, not a list of tunings.")
