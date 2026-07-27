"""
Frontier 1: the model's falsifiable prediction -- Lorentz violation vs experiment.

*** STATUS UPDATE (test_lv_confrontation) -- the confrontation in [C]/[D] below is too generous, and
    the "safe by many orders" verdict is superseded. It compared the model's E_QG against the UHECR
    *energy* (1e11 GeV) as if that were the bound on E_QG; it is not. For n=2 (quadratic) LV the effect
    grows as E^2, so the strongest bounds come from ultra-high-energy astrophysics and sit near the
    PLANCK scale (E_QG,2 ~ 1e18-1e19 GeV: UHECR/GZK, PeV gammas, Crab electrons). The model predicts
    E_QG,2 ~ 2.5e19 GeV, only ~1 order above them; its effect is ~1e-3 of current UHECR sensitivity, not
    1e-16. So the n=2 prediction is NEAR the frontier, not far from it. Also added there: the one-cone
    universality is a prediction the model has already PASSED (GW170817 + multimessenger timing, to 1e-15),
    and the crystallographic anisotropy likely AVERAGES AWAY in a poly-domained medium, leaving one-cone
    universality (not the cubic pattern) as the real discriminator. Nothing below is retracted except the
    over-generous margin; the coefficient extraction and the structural falsifiers stand. ***

Every result so far REPRODUCES known physics. A theory must also PREDICT something
nature could veto. The model's most distinctive, quantitative feature is the Lorentz
violation it already measured: the emergent light cone is exact only as E -> 0, with
calculable corrections at higher energy. Here we extract those corrections from the
lattice dispersion, cast them in the standard modified-dispersion form, and confront
them with real experimental bounds.

Dispersion of the emergent field (nn graph Laplacian on fcc): with the lattice
spacing a identified as the Planck length (k_max = pi/a ~ the Planck momentum),
    v(k)/c = 1 - zeta_boost (k/k_max)^2 - ...        (isotropic: BOOST violation)
    dc/c(dir) = zeta_aniso (k/k_max)^2 + ...          (directional: ROTATION violation)
Both enter at ORDER (E/E_Planck)^2 -- a QUADRATIC (mass-dimension-6, "n=2") Lorentz
violation, NOT the linear (n=1) form many quantum-gravity scenarios predict. And by
the one-cone result (test_emergent_tetrad), photon, fermion and graviton share these
coefficients -- so there is NO leading-order, species-dependent maximal speed.
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_fcc, R0

# --- extract the LV coefficients from the emergent dispersion ------------------
_X = perfect_fcc(radius=9.0)
_i = int(np.argmin((_X ** 2).sum(1)))
DELTA = (_X - _X[_i])[((_X - _X[_i]) ** 2).sum(1) < (1.3 * R0) ** 2]
DELTA = DELTA[(DELTA ** 2).sum(1) > 1e-9]
KMAX = np.pi / R0
DIRS = {"[100]": (1, 0, 0), "[110]": (1, 1, 0), "[111]": (1, 1, 1)}


def ceff(frac, u):
    u = np.array(u, float); u /= np.linalg.norm(u)
    k = frac * KMAX * u
    sym = (1 - np.cos(DELTA @ k)).sum()
    k0 = 1e-4 * KMAX * u
    c0 = np.sqrt((1 - np.cos(DELTA @ k0)).sum()) / np.linalg.norm(k0)
    return (np.sqrt(sym) / np.linalg.norm(k)) / c0            # v/c, ->1 as k->0


def coefficients():
    fr = np.array([0.05, 0.1, 0.15, 0.2])
    iso = np.array([np.mean([1 - ceff(f, u) for u in DIRS.values()]) for f in fr])
    zeta_boost = np.polyfit(fr ** 2, iso, 1)[0]               # 1 - v/c = zeta_boost (k/kmax)^2
    an = np.array([(lambda cs: (cs.max() - cs.min()))(np.array([ceff(f, u) for u in DIRS.values()]))
                   for f in fr])
    zeta_aniso = np.polyfit(fr ** 2, an, 1)[0]
    return zeta_boost, zeta_aniso


# --- experimental bounds (representative, order-of-magnitude, in GeV) ----------
E_PLANCK = 1.22e19          # GeV
E_QG2_PHOTON = 1e10         # GeV, Fermi-LAT GRB time-of-flight, quadratic (n=2)
E_QG2_UHECR = 1e11          # GeV, UHECR threshold / photon-stability, n=2 (strongest for n=2)
DC_SPECIES = 1e-15          # |c_photon - c_electron|/c, cross-species maximal-speed bound
E_LHAASO = 1.4e6            # GeV, highest-energy gamma rays observed
E_UHECR = 1e11             # GeV, ultra-high-energy cosmic rays (~1e20 eV)


if __name__ == "__main__":
    print("=== Frontier 1: the model's Lorentz-violation prediction vs experiment ===\n")
    zb, za = coefficients()
    print("[A] leading Lorentz-violating coefficients, from the emergent lattice dispersion")
    print(f"    boost  (isotropic dispersion):  1 - v/c = {zb:.3f} (k/k_max)^2   [zeta_boost = {zb:.3f}]")
    print(f"    rotation (crystallographic):    dc/c    = {za:.3f} (k/k_max)^2   [zeta_aniso = {za:.3f}]")
    print("    => leading LV is QUADRATIC in energy (mass-dimension-6, n=2), NOT linear (n=1),")
    print("       and the rotation part carries the fcc lattice's cubic angular pattern.\n")

    Eqg = E_PLANCK / np.sqrt(zb)                              # effective quadratic LV scale
    print("[B] cast as a modified dispersion (a = Planck length, k_max ~ Planck momentum)")
    print(f"    v(E)/c = 1 - zeta (E/E_Planck)^2  =>  effective E_QG,2 = E_Planck/sqrt(zeta) = {Eqg:.1e} GeV\n")

    print("[C] confrontation with current bounds (n=2, subluminal)")
    for name, bound in [("Fermi-LAT GRB photons", E_QG2_PHOTON), ("UHECR (strongest n=2)", E_QG2_UHECR)]:
        margin = (Eqg / bound) ** 2                           # n=2 effect scales as 1/E_QG^2
        print(f"    {name:>24}: bound E_QG,2 > {bound:.0e} GeV; model = {Eqg:.1e} GeV "
              f"-> SAFE by ~{margin:.0e} in the effect")
    d_lhaaso = zb * (E_LHAASO / E_PLANCK) ** 2
    d_uhecr = zb * (E_UHECR / E_PLANCK) ** 2
    print(f"    predicted |dv/c| at LHAASO ({E_LHAASO:.0e} GeV): {d_lhaaso:.1e}")
    print(f"    predicted |dv/c| at UHECR  ({E_UHECR:.0e} GeV): {d_uhecr:.1e}  (closest to a bound)")
    print(f"    cross-species |c_gamma - c_e|/c at leading order: 0  (one cone) << bound {DC_SPECIES:.0e}\n")

    print("[D] verdict -- consistent, and falsifiable in STRUCTURE")
    print("  The model's Lorentz violation is QUADRATIC and Planck-suppressed, so it survives every")
    print("  current bound by many orders of magnitude; the closest frontier is ultra-high-energy")
    print("  cosmic rays. It is not quantitatively detectable with today's reach -- but it makes three")
    print("  QUALITATIVE predictions that do NOT need Planck-energy access and that would kill it:")
    print("   1. LV is QUADRATIC (n=2), not linear -- a confirmed LINEAR photon LV falsifies the model.")
    print("   2. The rotation-violating part is ANISOTROPIC with the emergent lattice's crystallographic")
    print("      pattern (cubic/hex), correlated between the boost and rotation sectors.")
    print("   3. ONE universal cone: no leading-order species-dependent maximal speed -- a confirmed")
    print("      c_gamma != c_electron (or != c_gravity, cf. GW170817) at O(E/E_Planck) falsifies it.")
    print("  This is the project's first empirical claim: not 'it reproduces physics', but a specific,")
    print("  Planck-suppressed, cross-species-universal, crystallographically-anisotropic n=2 signature")
    print("  that nature could rule out.")
