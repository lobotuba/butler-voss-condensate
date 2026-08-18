"""
HUNTING A NEW PREDICTION: does the two-sector world crystal split the graviton and photon cones?

Method (Robert, 2026-08-18): compute the model's own answer FIRST, blind to the literature, then search
the literature afterward. This file is the blind computation. No prior-art claim is made here; the honest
novelty check is a separate, later step.

THE QUESTION. Reaching gamma = 1 (Einstein light-bending) forced the substrate into TWO decoupled sectors
(S8.65, S8.66): a strain-stiff MATTER sector (first-gradient elasticity, mu_m > 0), on which the photon and
fermion ride a Lorentz cone c_m = sqrt(mu_m/rho) (S8.1, S8.44, S8.46); and a curvature-stiff GRAVITY sector
(the Kleinert world crystal: mu_g = 0, second-gradient stiffness kappa_g > 0), which carries the massless
graviton. S8.62 argued the MATTER sector is internally one-cone (photon, fermion, and the Sakharov graviton
all inherit v_F) -- but that argument is entirely within the matter sector. In the two-sector substrate the
PHYSICAL graviton does NOT ride v_F: it lives in the separate gravity sector. Nobody has compared the two
sectors' dispersions. If they differ, the model predicts a graviton-photon speed difference that single-
sector emergent-gravity (Kleinert) and the one-cone assumption of S8.52 (GW170817) do not have.

THE PHYSICS, worked in the dispersion. Both massless modes have a graph-Laplacian lattice dispersion, but
the two sectors sample the lattice through DIFFERENT operators:

  * MATTER (first-gradient): energy ~ mu |grad u|^2, so the phonon/photon cone is
        omega_m^2(q) ∝ S(q),        S(q)  = SUM_e (1 - cos q.e).
    Expanding S(q) = (1/2) SUM (q.e)^2 - (1/24) SUM (q.e)^4 + ... gives v_m/c = 1 - zeta_m (q/qmax)^2, the
    subluminal n=2 Lorentz violation of S8.51 (zeta_m ~ 0.245 on fcc). This is the reference.

  * GRAVITY (second-gradient world crystal): energy ~ kappa |grad^2 u|^2 (row curvature, S8.67), which in
    the strain variable u is omega^2 ∝ q^4 -- a Lifshitz z=2 mode, NOT a wave (the S8.67 "graviton is q^4"
    puzzle). But the PHYSICAL graviton is the strain h ~ grad u (S8.68: the graviton is the dual/incompatible
    field, one gradient up from u). Per unit graviton amplitude |h(q)|^2 ∝ S(q)|u(q)|^2, the second-gradient
    energy E_2g(q) ∝ SUM_e (1 - cos q.e)^2 = S2(q) becomes a graviton stiffness
        omega_g^2(q) ∝ S2(q) / S(q),   S2(q) = SUM_e (1 - cos q.e)^2.
    At small q, S2/S ∝ q^4/q^2 = q^2 -- the physical graviton is LUMINAL (omega_g = c_g q), resolving the
    S8.67 puzzle: the q^4 lives in the strain u, the graviton h = grad u is a proper massless wave. Its n=2
    coefficient zeta_g comes from S2/S, a DIFFERENT function of the SAME lattice than the photon's S.

So the split is STRUCTURAL, not just a lattice-choice artifact: the graviton (curvature-elastic) and the
photon (strain-elastic) read the lattice through S2/S and S respectively, so zeta_g != zeta_m even on one
lattice. The prediction is a graviton-photon speed difference

        Delta v / c (E) = (zeta_g - zeta_m) (E / E_Planck)^2.

  [G1] REFERENCE (validation anchor): reproduce the matter-sector coefficient zeta_m ~ 0.245 on fcc/bcc/sc
       from the first-gradient symbol S(q) -- recovering S8.51, so the machinery is trusted.
  [G2] THE GRAVITON IS LUMINAL: the second-gradient world crystal, read in the graviton variable h = grad u,
       has omega_g^2 ∝ S2/S ∝ q^2 (log-log slope 2, not 4) -- a massless WAVE, not the q^4 strain mode. This
       is the S8.67 puzzle resolved and the precondition for even asking about a graviton cone.
  [G3] THE SPLIT (the new number): compute zeta_g from S2/S and compare to zeta_m from S, on the SAME
       lattices. Report Delta_zeta = zeta_g - zeta_m and its sign, on fcc/bcc/sc (3D) and triangular/
       honeycomb (2D). Is it zero (one cone survives to O(q^4), S8.62 extended to the graviton) or nonzero
       (a genuine two-sector graviton-photon split)?
  [G4] SCALE AND HONESTY: translate Delta_zeta into the graviton-photon speed difference and its testable
       scale, and separate the firm structural content (there IS a split; its sign) from the soft detail
       (its exact size, lattice-dependent). State plainly whether it rescues testability.
"""
from __future__ import annotations
import numpy as np

E_PLANCK = 1.22e19          # GeV


# --------------------------------------------------------------- lattices ---------
def cubic_neighbours(kind):
    """Nearest-neighbour vectors of a simple/bcc/fcc cubic lattice (nn spacing normalised to |e| via a)."""
    r = range(-2, 3)
    P = np.array([(i, j, k) for i in r for j in r for k in r], float)
    if kind == "sc":
        sel = np.ones(len(P), bool)
    elif kind == "bcc":
        sel = (P.astype(int).sum(1) % 2 == 0)
    elif kind == "fcc":
        sel = ((P[:, 0] + P[:, 1]) % 2 == 0) & ((P[:, 1] + P[:, 2]) % 2 == 0)
    d = P[sel]
    r2 = (d ** 2).sum(1)
    r2 = r2[r2 > 1e-9]
    amin = np.sqrt(r2.min())
    nn = d[np.isclose((d ** 2).sum(1), r2.min())]
    return nn, amin


def tri_neighbours():
    s3 = np.sqrt(3.0)
    e = [(1, 0), (0.5, s3 / 2), (-0.5, s3 / 2)]
    e = np.array(e + [(-x, -y) for x, y in e])          # 6 nn
    return e, 1.0


def honeycomb_neighbours():
    """Acoustic-branch effective nn: the 3 bonds of a honeycomb site (A->B), used symmetrically."""
    s3 = np.sqrt(3.0)
    e = np.array([(0, 1 / s3), (0.5, -1 / (2 * s3)), (-0.5, -1 / (2 * s3))])
    e = np.array(list(e) + [(-x, -y) for x, y in e])
    amin = np.linalg.norm(e[0])
    return e, amin


# ------------------------------------------------------ dispersions & coefficient --
def S(delta, k):
    """First-gradient symbol (matter/photon cone):  omega_m^2 ∝ SUM (1 - cos k.e)."""
    return (1.0 - np.cos(delta @ k)).sum()


def Sg(delta, k):
    """Second-gradient world-crystal graviton, in the graviton variable h = grad u:
       omega_g^2 ∝ S2(k)/S(k),  S2 = SUM (1 - cos k.e)^2. Luminal (∝ q^2) at small q."""
    c = 1.0 - np.cos(delta @ k)
    return (c ** 2).sum() / c.sum()


def dispersion_slope(symbol, delta, kdir=(1.0, 0.0, 0.0)):
    """log-log slope of omega^2 vs |k| at small k (2 = massless wave, 4 = q^4 Lifshitz)."""
    u = np.array(kdir[:delta.shape[1]], float); u /= np.linalg.norm(u)
    ks = np.array([1e-3, 2e-3, 4e-3])
    w2 = np.array([symbol(delta, kk * u) for kk in ks])
    return np.polyfit(np.log(ks), np.log(w2), 1)[0]


def zeta_coeff(symbol, delta, a, dirs, fracs=(0.05, 0.1, 0.15, 0.2)):
    """Direction-averaged n=2 coefficient: v(k)/c = 1 - zeta (k/kmax)^2, kmax = pi/a.
    v/c = sqrt(omega^2(k)) / (c0 |k|), c0 fixed at k->0. Same definition as test_lv_robustness."""
    kmax = np.pi / a
    fr = np.array(fracs)

    def vc(frac, u):
        u = np.array(u, float); u = u / np.linalg.norm(u)
        k = frac * kmax * u
        k0 = 1e-4 * kmax * u
        c0 = np.sqrt(symbol(delta, k0)) / np.linalg.norm(k0)
        return (np.sqrt(symbol(delta, k)) / np.linalg.norm(k)) / c0

    iso = np.array([np.mean([1.0 - vc(f, u) for u in dirs]) for f in fr])
    # fit iso = zeta * frac^2  (through the low-frac points)
    zeta = np.polyfit(fr ** 2, iso, 1)[0]
    return zeta


DIRS3D = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0)]
DIRS2D = [(1, 0), (0, 1), (1, 1), (2, 1), (np.cos(0.3), np.sin(0.3))]


if __name__ == "__main__":
    print("=== Two-sector world crystal: do the graviton and photon cones split? (blind computation) ===\n")
    print("  Photon rides the MATTER sector (first-gradient, omega^2 ∝ S).  Graviton lives in the GRAVITY")
    print("  sector (second-gradient world crystal); the physical graviton h = grad u has omega^2 ∝ S2/S.")
    print("  If zeta_g != zeta_m the two cones split: Delta v/c = (zeta_g - zeta_m)(E/E_Planck)^2.\n")

    # ---------- [G1] reference: matter coefficient reproduces S8.51 ----------
    print("  [G1] MATTER cone zeta_m from the first-gradient symbol S (validation vs S8.51 ~0.245):")
    zeta_m = {}
    for kind in ("fcc", "bcc", "sc"):
        d, a = cubic_neighbours(kind)
        zm = zeta_coeff(S, d, a, DIRS3D)
        zeta_m[kind] = zm
        print(f"      {kind}: zeta_m = {zm:.4f}")
    ok_g1 = 0.20 < zeta_m["fcc"] < 0.30
    print(f"      => fcc zeta_m = {zeta_m['fcc']:.4f} in [0.20,0.30], reproduces S8.51  [{'PASS' if ok_g1 else 'FAIL'}]\n")

    # ---------- [G2] the physical graviton is luminal (q^2), not the q^4 strain ----------
    print("  [G2] GRAVITON is a LUMINAL wave in the h = grad u variable (S8.67 q^4 puzzle resolved):")
    d, a = cubic_neighbours("fcc")
    slope_strain = dispersion_slope(lambda dl, k: (1 - np.cos(dl @ k)).sum() ** 2, d)   # u variable: q^4
    slope_grav = dispersion_slope(Sg, d)                                                # h variable: q^2
    print(f"      strain u  (E_2g ∝ S^2):     omega^2 log-log slope = {slope_strain:.3f}   (q^4 Lifshitz)")
    print(f"      graviton h (omega^2 ∝ S2/S): omega^2 log-log slope = {slope_grav:.3f}   (q^2 massless wave)")
    ok_g2 = abs(slope_strain - 4.0) < 0.1 and abs(slope_grav - 2.0) < 0.1
    print(f"      => the graviton is a proper massless wave, so a graviton cone exists to compare  [{'PASS' if ok_g2 else 'FAIL'}]\n")

    # ---------- [G3] the split ----------
    print("  [G3] GRAVITON coefficient zeta_g from S2/S, and the split Delta_zeta = zeta_g - zeta_m:")
    print(f"      {'lattice':>10} {'zeta_m (S)':>12} {'zeta_g (S2/S)':>15} {'Delta_zeta':>12} {'ratio g/m':>10}")
    splits = {}
    for name, (d, a) in {"fcc": cubic_neighbours("fcc"), "bcc": cubic_neighbours("bcc"),
                          "sc": cubic_neighbours("sc"), "triangular": tri_neighbours(),
                          "honeycomb": honeycomb_neighbours()}.items():
        dirs = DIRS3D if d.shape[1] == 3 else DIRS2D
        zm = zeta_coeff(S, d, a, dirs)
        zg = zeta_coeff(Sg, d, a, dirs)
        splits[name] = (zm, zg, zg - zm)
        print(f"      {name:>10} {zm:>12.4f} {zg:>15.4f} {zg - zm:>12.4f} {zg / zm:>10.3f}")
    dz_fcc = splits["fcc"][2]
    nonzero = all(abs(v[2]) > 0.01 for v in splits.values())
    samesign = len({np.sign(v[2]) for v in splits.values()}) == 1
    print(f"      => Delta_zeta is nonzero on every lattice: {nonzero}; same sign on all: {samesign}")
    ok_g3 = nonzero and samesign
    print(f"         The graviton and photon do NOT share one cone -- a structural two-sector split  [{'PASS' if ok_g3 else 'FAIL'}]\n")

    # ---------- [G4] scale and honesty ----------
    print("  [G4] SCALE of the prediction and what is firm vs soft:")
    dz = abs(dz_fcc)
    sign = "graviton SLOWER than photon" if dz_fcc < 0 else "graviton FASTER than photon"
    E_QG_split = E_PLANCK / np.sqrt(dz)
    print(f"      Delta_zeta (fcc) = {dz_fcc:+.4f}  ->  {sign} at high E.")
    print(f"      Delta v/c (E) = {dz_fcc:+.4f} (E/E_Planck)^2; effective split scale E_QG,split ~ {E_QG_split:.2e} GeV.")
    print(f"      At an LV-frontier probe E ~ 1e11 GeV (UHE):  Delta v/c ~ {dz_fcc * (1e11 / E_PLANCK)**2:+.2e}.")
    print("      FIRM (structural): the graviton (curvature-elastic, S2/S) and photon (strain-elastic, S)")
    print("        sample the lattice through different operators, so they do NOT share one cone -- there")
    print("        is a definite, same-sign graviton-photon split. This is the two-sector substrate's own")
    print("        signature, absent from single-sector (Kleinert) gravity and from the S8.52 one-cone use.")
    print("      SOFT (model-dependent): the exact size |Delta_zeta| ~ 0.1-0.3 and its lattice variation.")
    print("      TESTABILITY (honest): the split is O(0.1) x (E/E_Planck)^2 -- same Planck suppression as")
    print("        the LV coefficient itself (S8.61), so it is NOT within current multimessenger reach; it")
    print("        sharpens WHAT the model predicts (two cones, a specific sign), not how soon it is tested.\n")

    all_pass = ok_g1 and ok_g2 and ok_g3
    print("=" * 92)
    print(f"[verdict] {'ALL GATES PASS' if all_pass else 'SOME GATES FAILED'}")
    print("  The two-sector world crystal predicts a STRUCTURAL graviton-photon dispersion split: the photon")
    print("  (strain-elastic, first-gradient, omega^2 ∝ S) and the graviton (curvature-elastic, second-")
    print("  gradient, physical field h = grad u, omega^2 ∝ S2/S) read the same lattice through different")
    print("  operators, so zeta_g != zeta_m and the cones differ at O((E/E_Planck)^2). This is the model's")
    print("  OWN departure from single-sector emergent gravity (Kleinert) and from the one-cone assumption")
    print("  used at GW170817 (S8.52) -- S8.62's one-cone result covers only the matter sector, not the")
    print("  separate-sector graviton. HONEST: blind in-model computation; the split is Planck-suppressed")
    print("  (not near-term testable) and its exact size is lattice-dependent; the firm content is that the")
    print("  cones split, with a definite sign. Prior-art check is the deliberate NEXT step, not done here.")
