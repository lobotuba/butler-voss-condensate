"""
Route 1: the elasticity-fracton duality -- the medium's hidden long-range gravity sector.

test_graviton.py found gravity-by-density screens (Bitter-Crum) because it couples
to energy density = a center of DILATATION, and the medium relaxes it. But that is
not the medium's only defect sector. In 2D linear elasticity every defect sources
the same BIHARMONIC (Airy) equation
    nabla^4 chi = source,
and the defects differ only by how many derivatives hit the delta -- i.e. their
multipole order -- which fixes their range:

    disclination (curvature)   source = s delta        -> chi ~ r^2 ln r   LONG-RANGE
    dislocation  (torsion)     source = b . grad delta -> ~ ln r           LONG-RANGE
    dilatation   (compression) source = nabla^2 delta  -> contact          SCREENED (Bitter-Crum)

So ENERGY density couples as a dilatation -- the MOST screened multipole -- which is
exactly why gravity-by-density screened. CURVATURE couples as a disclination -- long
-range. In the elasticity-fracton duality (Pretko-Radzihovsky) these defects are the
charges of a rank-2 SYMMETRIC-TENSOR gauge theory (the disclination is a 'fracton'),
which is the structure of linearized gravity. And in 2D gravity a point MASS is a
conical deficit is a DISCLINATION. So the medium's gravity sector is already present
and long-range -- we were coupling to the wrong multipole.

Computed from the biharmonic Green's function G4 = FFT^-1(1/k^4); the defect
interactions are its multipole derivatives: 1/k^4 (disclination), 1/k^2 (dislocation),
1 (dilatation).
"""
from __future__ import annotations
import numpy as np

L = 256
_k = 2 * np.pi * np.fft.fftfreq(L)
KX, KY = np.meshgrid(_k, _k, indexing="ij")
K2 = 2 * ((1 - np.cos(KX)) + (1 - np.cos(KY)))            # lattice Laplacian symbol ~ |k|^2


def green(power):
    """FFT^-1 of 1/K2^power (k=0 excluded): the defect interaction kernel."""
    G = np.zeros_like(K2); m = K2 > 1e-9
    G[m] = 1.0 / K2[m] ** power
    return np.fft.ifft2(G).real


def main():
    print("=== Route 1: elasticity-fracton duality -- the medium's long-range gravity sector ===\n")
    Gdisc = green(2)     # 1/k^4  disclination (curvature)
    Gdisl = green(1)     # 1/k^2  dislocation  (torsion)
    Gdilat = green(0)    # 1      dilatation   (compression = energy density)
    r = np.arange(1, L // 4)
    gc, gl, gd = Gdisc[r, 0], Gdisl[r, 0], Gdilat[r, 0]

    print(f"  {'r':>4} {'disclination':>13} {'dislocation':>12} {'dilatation':>11}")
    for i in (1, 2, 4, 8, 16, 32, 48):
        print(f"  {r[i]:>4} {gc[i]:>13.4f} {gl[i]:>12.4f} {gd[i]:>11.2e}")

    m = (r >= 4) & (r <= 60)
    sl = np.polyfit(np.log(r[m]), gl[m], 1)[0]       # dislocation: G ~ slope * ln r
    dil = np.abs(gd[r > 0]).max()                    # dilatation: contact magnitude
    # disclination place in the hierarchy: Laplacian of its kernel = the dislocation kernel
    lapc = np.fft.ifft2(K2 * np.fft.fft2(Gdisc)).real
    rel = np.abs((lapc - Gdisl)[r, 0]).max() / (np.abs(Gdisl[r, 0]).max() + 1e-12)

    print(f"\n  DILATATION (energy density, source nabla^2 delta): |G(r>0)| <= {dil:.1e} ~ 0 --")
    print("    a CONTACT term only: SCREENED (Bitter-Crum). This IS the gravity-by-density coupling.")
    print(f"  DISLOCATION (torsion, source grad delta): G ~ {sl:+.3f} ln r + const (2D-Coulomb log;")
    print(f"    1/2pi = {1/(2*np.pi):.3f}): LONG-RANGE.")
    print("  DISCLINATION (curvature, source delta): the biharmonic 1/k^4 kernel -- the LEAST-")
    print("    screened multipole (free-space interaction ~ r^2 ln r, even longer-range). Check:")
    print(f"    Laplacian(disclination kernel) = dislocation kernel to rel-err {rel:.0e}, confirming")
    print("    it sits one integration BEYOND the (already long-range) dislocation. (Its magnitude")
    print("    is periodic-box-saturated here, so the clean quantitative contrast is the two above.)")
    print("\n  => the multipole order of the source sets the range: energy density couples as a")
    print("     dilatation (nabla^2 delta, the most-screened multipole) -- exactly why gravity-by-")
    print("     density screened. CURVATURE couples as a disclination (delta) -- LONG-RANGE and")
    print("     unscreened. Same medium, a different (tensor-gauge) sector.")
    print("\n  Interpretation (Pretko-Radzihovsky elasticity-fracton duality):")
    print("   * these defects are the charges of a rank-2 SYMMETRIC-TENSOR gauge theory -- the")
    print("     structure of LINEARIZED GRAVITY; the disclination is an immobile 'fracton'.")
    print("   * in 2D gravity a point MASS = a conical deficit = a DISCLINATION. So the medium")
    print("     already contains gravitational 'masses' that curve space around them, with a")
    print("     genuinely LONG-RANGE interaction -- the sector gravity-by-density missed.")
    print("\n  => Route 1 overcomes the screening: couple gravity to CURVATURE (disclinations,")
    print("     the tensor-gauge/fracton sector) rather than to energy density (dilatation).")
    print("     Open next: 3D (where the tensor gauge field is a dynamical spin-2 graviton), and")
    print("     wiring matter energy to disclination density so 'mass curves the medium' long-range.")


if __name__ == "__main__":
    main()
