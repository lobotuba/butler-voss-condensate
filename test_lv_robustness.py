"""
Sharpening the one prediction: how robust is the Lorentz-violation coefficient?

The model's single live falsifiable prediction (Frontier 1, test_lv_prediction; confronted with data in
test_lv_confrontation) is a quadratic, subluminal, species-universal Lorentz violation with an effective
scale E_QG,2 = E_Planck / sqrt(zeta_boost), zeta_boost ~ 0.245 read off the emergent fcc dispersion. But
0.245 came from ONE microscopic choice -- an fcc lattice, nearest-neighbour graph Laplacian. If E_QG,2
is to be a prediction rather than a fit, the honest question is which parts of it are firm and which are
soft: does the coefficient survive changing the lattice, the neighbour range, the direction? This file
stress-tests exactly that, separating the falsifiable structural content from the model-dependent detail.

For a graph-Laplacian dispersion omega^2(k) = SUM_delta (1 - cos k.delta), expanding in k gives an
isotropic quadratic speed and a quartic Lorentz-violating correction, so along a direction u,
    v(k)/c = 1 - zeta_boost(u) (k/k_max)^2 - ... ,   k_max = pi / a (a = nn spacing ~ Planck length),
with zeta_boost the direction-averaged coefficient and its spread the crystallographic anisotropy.

  [G1] BASELINE: the fcc value reproduces test_lv_prediction (zeta_boost ~ 0.245 -> E_QG,2 ~ 2.5e19 GeV).
  [G2] LATTICE ROBUSTNESS: fcc, bcc, sc all give zeta_boost of the same order and the same (subluminal)
       sign, so E_QG,2 stays within a small factor of a few x Planck -- the SCALE is robust; the exact
       0.245 is lattice-specific.
  [G3] NEIGHBOUR-RANGE ROBUSTNESS: adding the second-neighbour shell shifts zeta but keeps its order and
       sign -- the prediction does not hinge on the coupling stopping at nearest neighbours.
  [G4] THE FIRM FALSIFIERS: zeta_boost > 0 (SUBLUMINAL) for every lattice, and the anisotropy is a
       subleading fraction of it that averages down over directions -- so 'quadratic + subluminal +
       one-cone, E_QG ~ few x Planck' is the robust, falsifiable content; the coefficient and the
       crystallographic pattern are the soft, model-dependent detail.
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_fcc, R0

E_PLANCK = 1.22e19          # GeV


# ------------------------------------------------------------- lattice builders ---
def cubic(kind, m=6, a=1.0):
    """Return centered points of a simple/bcc/fcc cubic lattice in a box of half-width m."""
    r = range(-m, m + 1)
    pts = [(i, j, k) for i in r for j in r for k in r]
    P = np.array(pts, float)
    if kind == "sc":
        sel = np.ones(len(P), bool)
    elif kind == "bcc":
        sel = ((P.astype(int).sum(1)) % 2 == 0)            # even-sum simple-cubic sublattice = bcc
    elif kind == "fcc":
        sel = ((P[:, 0] + P[:, 1]) % 2 == 0) & ((P[:, 1] + P[:, 2]) % 2 == 0)
    return P[sel] * a


def shells(P, n_shells=1, tol=1e-6):
    """Neighbour vectors of the central point, out to n_shells distinct distances."""
    c = P[np.argmin((P ** 2).sum(1))]
    d = P - c
    r = np.sqrt((d ** 2).sum(1))
    d = d[r > tol]; r = r[r > tol]
    uniq = np.sort(np.unique(np.round(r, 6)))
    cut = uniq[n_shells - 1] + tol
    return d[r <= cut], uniq[0]                            # (neighbour vectors, nn distance a)


# ----------------------------------------------------- dispersion / coefficient ---
def zeta(delta, a, dirs, fracs=(0.05, 0.1, 0.15, 0.2)):
    """zeta_boost (direction-averaged) and zeta_aniso (spread) from v(k)/c = 1 - zeta (k/kmax)^2."""
    kmax = np.pi / a
    fr = np.array(fracs)

    def vc(frac, u):
        u = np.array(u, float); u = u / np.linalg.norm(u)
        k = frac * kmax * u
        k0 = 1e-4 * kmax * u
        c0 = np.sqrt((1 - np.cos(delta @ k0)).sum()) / np.linalg.norm(k0)
        return (np.sqrt((1 - np.cos(delta @ k)).sum()) / np.linalg.norm(k)) / c0

    iso = np.array([np.mean([1 - vc(f, u) for u in dirs]) for f in fr])
    zb = np.polyfit(fr ** 2, iso, 1)[0]
    an = np.array([np.ptp([vc(f, u) for u in dirs]) for f in fr])
    za = np.polyfit(fr ** 2, an, 1)[0]
    return zb, za


CUBIC_DIRS = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]


def main():
    print("=== Sharpening the one prediction: robustness of the LV coefficient ===\n")
    ok = True

    # ---- [G1] baseline: reproduce test_lv_prediction on fcc ----
    Xf = perfect_fcc(radius=9.0)
    df, af = shells(Xf, 1)
    zb_fcc, za_fcc = zeta(df, af, CUBIC_DIRS)
    Eqg_fcc = E_PLANCK / np.sqrt(zb_fcc)
    g1 = abs(zb_fcc - 0.245) < 0.03
    ok &= g1
    print("  [G1] baseline (fcc, nn) -- reproduce test_lv_prediction:")
    print(f"       zeta_boost = {zb_fcc:.3f} (expect ~0.245),  E_QG,2 = E_Pl/sqrt(zeta) = {Eqg_fcc:.2e} GeV"
          f"  -> {'PASS' if g1 else 'FAIL'}\n")

    # ---- [G2] robustness across lattice type ----
    print("  [G2] robustness across microscopic lattice (nn graph Laplacian):")
    print(f"       {'lattice':>8} {'zeta_boost':>11} {'sign':>6} {'E_QG,2 (GeV)':>14} {'/E_Planck':>10}")
    band = []
    for kind in ("sc", "bcc", "fcc"):
        d, a = shells(cubic(kind), 1)
        zb, _ = zeta(d, a, CUBIC_DIRS)
        band.append(zb)
        Eqg = E_PLANCK / np.sqrt(zb)
        print(f"       {kind:>8} {zb:>11.3f} {'sub' if zb>0 else 'SUPER':>6} {Eqg:>14.2e} {Eqg/E_PLANCK:>10.2f}")
    band = np.array(band)
    g2 = np.all(band > 0) and (band.max() / band.min() < 4.0)
    ok &= g2
    print(f"       zeta stays same sign (subluminal) and within a factor {band.max()/band.min():.1f}; "
          f"E_QG,2 in [{E_PLANCK/np.sqrt(band.max()):.1e}, {E_PLANCK/np.sqrt(band.min()):.1e}] GeV"
          f"  -> {'PASS' if g2 else 'FAIL'}\n")

    # ---- [G3] robustness to neighbour range ----
    print("  [G3] robustness to neighbour range (fcc, nn vs nn+2nn):")
    rows = []
    for ns in (1, 2):
        d, a = shells(Xf, ns)
        zb, _ = zeta(d, a, CUBIC_DIRS)
        rows.append(zb)
        print(f"       shells={ns}: zeta_boost = {zb:.3f}  (sign {'sub' if zb>0 else 'SUPER'})")
    g3 = (rows[0] > 0) and (rows[1] > 0) and (0.25 < rows[1] / rows[0] < 4.0)
    ok &= g3
    print(f"       same sign and order across neighbour range -> {'PASS' if g3 else 'FAIL'}\n")

    # ---- [G4] the firm falsifiers vs the soft detail ----
    print("  [G4] firm structural content vs model-dependent detail:")
    all_sub = np.all(band > 0)
    aniso_frac = abs(za_fcc) / zb_fcc
    g4 = all_sub and (aniso_frac < 0.6)
    ok &= g4
    print(f"       SUBLUMINAL (zeta>0) for every lattice: {bool(all_sub)}  (a firm, lattice-independent falsifier)")
    print(f"       QUADRATIC (n=2) by construction of the k^4 term: firm")
    print(f"       anisotropy zeta_aniso/zeta_boost = {aniso_frac:.2f} (fcc): subleading, and averages down")
    print(f"         over directions in a poly-domained medium -> the SOFT, model-dependent part"
          f"  -> {'PASS' if g4 else 'FAIL'}\n")

    print("=" * 80)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  The prediction separates cleanly into firm and soft. FIRM (falsifiable, lattice-independent):")
    print("  the Lorentz violation is QUADRATIC (n=2), SUBLUMINAL (zeta>0 on every lattice tested), and --")
    print("  with one-cone universality (test_emergent_tetrad) -- species-universal, at an effective scale")
    print(f"  E_QG,2 of a few x the Planck energy ({E_PLANCK/np.sqrt(band.max()):.1e}-{E_PLANCK/np.sqrt(band.min()):.1e} GeV across lattices). SOFT (model-")
    print("  dependent): the exact coefficient 0.245 and the crystallographic anisotropy pattern, which")
    print("  move with the microscopic lattice and average down in a poly-domained medium. So what nature")
    print("  would actually test -- n=2 vs n=1, subluminal vs superluminal, one cone vs species-dependent,")
    print("  and a near-Planckian scale -- is exactly the robust part; the fit-like number is not load-bearing.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
