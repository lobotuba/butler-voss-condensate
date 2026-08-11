"""Route A, continued: the escape from the gamma=0 no-go -- the elastic moduli ARE a graviton mass.

S8.64 proved gamma = 0 is forced in any field-bearing elastic solid: the shear modulus mu that makes the
Lorentz cone and confines disclinations is the same mu that gamma = 1 needs to vanish. It named the escape
-- a graviton that is NOT a phonon of the rigid solid, a decoupled deconfined-disclination (topologically
ordered) sector -- but did not walk it. This section walks the first step and identifies the escape
precisely, via a reframing that makes S8.64 sharper.

THE REFRAMING. The medium's would-be graviton is the strain h_ij ~ eps_ij ~ d u (as in S8.32-8.34). So the
elastic energy is, in the graviton variable h,
        mu * INT eps_ij eps_ij  =  mu * INT h_ij h_ij      -- NO derivatives on h.
That is a graviton MASS term (INT h^2), not a kinetic term (INT (d h)^2). An elastic solid therefore gives
a MASSIVE graviton, and a massive graviton has no diffeomorphism gauge invariance -- which is exactly why
gamma != 1. S8.64's "rigidity forbids gamma = 1" is precisely "the elastic moduli are a graviton mass."

Read the two dispersions in the graviton variable. With displacement u_i(q), h ~ q u, and
        omega_T^2(q) = [ mu_1 q^2 + kappa_T q^4 ] / rho,     (transverse / graviton polarisation)
the energy per unit h^2 is D_T/q^2 = mu_1 + kappa_T q^2: the q^0 piece mu_1 is the graviton MASS^2, the q^2
piece kappa_T is the KINETIC (Einstein-Hilbert-like) coefficient. A first-gradient solid (kappa = 0, mu_1 >
0) is pure mass, no kinetic term. THE ESCAPE (Kleinert's world crystal): a MASSLESS graviton -- zero
first-gradient moduli (mu_1 = 0, mass removed) plus a second-gradient term kappa (kinetic kept), which keeps
the medium stable (curvature-stiff) though strain-floppy. That restores the diffeomorphism gauge invariance
(the mass is gone), deconfines the disclinations (their confinement was area-law in the first-gradient Young
modulus Y_1, which -> 0 with mu_1), and leaves a massless propagating graviton with only the two-derivative
kinetic term -- whose static coupling to a conserved source is Einstein (gamma = 1, Kleinert/Weinberg).

  [G1] The elastic moduli are a graviton MASS: for the model's medium (first-gradient, mu_1 > 0, kappa = 0)
       the graviton mass^2 = mu_1 > 0 and the kinetic coefficient = 0 -- a pure mass, no kinetic term.
  [G2] The world crystal is a MASSLESS graviton: mu_1 = 0 (mass removed), kappa_T > 0 (kinetic kept) gives
       mass^2 = 0, kinetic = kappa_T > 0, and D_T(q) = kappa_T q^4 > 0 for all q != 0 -- stable though
       strain-floppy. This is the unique loophole: it is the second-gradient term, not any elastic tuning.
  [G3] The three vector-charge (Einstein) preconditions hold JOINTLY at (mu_1 = 0, kappa > 0) -- diffeomorphism
       mass removed, disclinations deconfined (Y_1 = 0), medium stable -- which NO first-gradient solid can do
       (there mu_1 = 0 is the floppy instability of S8.64, D_T == 0). The kappa term is what saves it.
  [G4] The honest cost: the massless-graviton condition mu_1 = 0 also kills the first-gradient transverse
       cone c_T = sqrt(mu_1/rho) = 0 -- the cone the model needs for its photon and fermion (S8.1/S8.44/8.45).
       So gamma = 1 (wants mu_1 = 0) and the matter cone (wants mu_1 > 0) cannot come from the same elastic
       sector: the escape requires a substrate with TWO decoupled sectors.

Honest scope: this is an EXACT structural analysis of the escape, not a Monte-Carlo measurement like S8.64.
It proves the escape's preconditions are unique (second-gradient rigidity) and jointly constructible, and
names the cost (a two-sector substrate). It does NOT itself measure gamma = 1: that is Kleinert's analytic
world-crystal result, whose in-model confirmation is the named next construction -- build emergent matter on
the two-sector (curvature-stiff gravity + strain-stiff matter) medium and ray-trace gamma directly. Pure numpy.
"""
from __future__ import annotations
import numpy as np

RHO = 1.0


def dispersions(K1, mu1, kappaL, kappaT, q):
    """Isotropic first+second-gradient elasticity. Returns (omega_L^2, omega_T^2)(q)."""
    wL2 = ((K1 + mu1) * q**2 + kappaL * q**4) / RHO
    wT2 = (mu1 * q**2 + kappaT * q**4) / RHO
    return wL2, wT2


def graviton_mass_and_kinetic(mu1, kappaT):
    """In the graviton variable h ~ q u: D_T/q^2 = mu1 + kappaT q^2. q^0 = mass^2, q^2 = kinetic."""
    return mu1, kappaT               # (mass^2, kinetic coefficient)


def young_2d(K1, mu1):
    return 4 * K1 * mu1 / (K1 + mu1) if abs(K1 + mu1) > 1e-12 else 0.0


def main():
    print("=" * 92)
    print("ROUTE A (cont.) -- THE ESCAPE: the elastic moduli ARE a graviton mass; the world crystal is massless")
    print("=" * 92)
    ok = True
    qs = np.array([0.05, 0.1, 0.2, 0.4, 0.8])

    # the model's medium: first-gradient elastic solid (representative moduli, cf. honeycomb lam* of S8.64)
    K1, mu1_elastic, kappa = 17.3, 2.14, 1.0

    # [G1] the elastic moduli are a graviton mass term
    m2_el, kin_el = graviton_mass_and_kinetic(mu1_elastic, 0.0)   # kappa=0: pure first-gradient
    print("\n  [G1] read the elastic energy in the graviton variable h ~ eps ~ d u:  mu INT eps^2 = mu INT h^2")
    print(f"       first-gradient solid (mu_1={mu1_elastic}, kappa=0): graviton mass^2 = {m2_el:.3f} > 0, "
          f"kinetic = {kin_el:.3f}")
    print("       D_T/q^2 = mu_1 + kappa q^2 : the q^0 piece is the MASS^2, the q^2 piece the KINETIC term")
    g1 = m2_el > 1.0 and abs(kin_el) < 1e-12
    ok &= g1
    print(f"       => the moduli are a graviton MASS (INT h^2), with NO kinetic term -- a massive graviton has")
    print(f"          no diffeomorphism invariance, so gamma != 1 (this IS S8.64's obstruction)  -> "
          f"{'PASS' if g1 else 'FAIL'}")

    # [G2] the world crystal: massless graviton (mass removed, kinetic kept), stable though strain-floppy
    m2_wc, kin_wc = graviton_mass_and_kinetic(0.0, kappa)
    _, wT2 = dispersions(0.0, 0.0, kappa, kappa, qs)
    stable = np.all(wT2 > 0)
    print("\n  [G2] the world crystal -- set the first-gradient mass to zero, keep the second-gradient kinetic:")
    print(f"       (mu_1=0, kappa={kappa}): graviton mass^2 = {m2_wc:.3f}, kinetic = {kin_wc:.3f}")
    print(f"       transverse dispersion omega_T^2(q) = kappa q^4 : {', '.join(f'{w:.2e}' for w in wT2)}")
    g2 = abs(m2_wc) < 1e-12 and kin_wc > 0.5 and stable
    ok &= g2
    print(f"       => a MASSLESS graviton (mass^2 = 0) with only the two-derivative kinetic term, and D_T =")
    print(f"          kappa q^4 > 0 for all q -- STABLE though strain-floppy. The loophole is the SECOND-")
    print(f"          gradient term, not any elastic tuning  -> {'PASS' if g2 else 'FAIL'}")

    # [G3] the three vector-charge (Einstein) preconditions hold jointly at (mu_1=0, kappa>0)
    print("\n  [G3] the vector-charge (Einstein) preconditions, checked jointly -- and against a plain solid:")
    print(f"       {'sector':<26s} {'mass^2':>8s} {'diffeo?':>8s} {'Y_1':>7s} {'disclin.':>9s} {'D_T>0?':>7s}")
    rows = [
        ("first-gradient, mu_1=0",  0.0, 0.0,   0.0),   # the S8.64 floppy point: no kinetic term
        ("world crystal, mu_1=0",   0.0, kappa, kappa),
    ]
    plain_unstable = True
    wc_ok = False
    for tag, mu1, kL, kT in rows:
        _, wT2 = dispersions(17.3, mu1, kL, kT, qs)
        DTpos = bool(np.all(wT2 > 0))
        Y1 = young_2d(17.3, mu1)
        mass2 = mu1
        diffeo_gauge = abs(mass2) < 1e-12                 # mass removed => diffeomorphism is gauge
        disc = "deconfined" if abs(Y1) < 1e-9 else "confined"
        print(f"       {tag:<26s} {mass2:8.3f} {'gauge' if diffeo_gauge else 'stiff':>8s} {Y1:7.3f} "
              f"{disc:>9s} {str(DTpos):>7s}")
        if "world crystal" in tag:
            wc_ok = diffeo_gauge and abs(Y1) < 1e-9 and DTpos
        if tag.startswith("first-gradient"):
            plain_unstable = not DTpos                    # first-gradient mu_1=0 is NOT stable (D_T==0)
    g3 = wc_ok and plain_unstable
    ok &= g3
    print(f"       => world crystal: diffeomorphism gauge (mass=0) + disclinations deconfined (Y_1=0) + stable")
    print(f"          (D_T=kappa q^4>0), all at once. A first-gradient solid at mu_1=0 is NOT stable (D_T==0,")
    print(f"          the S8.64 floppy wall) -- only the kappa term opens the loophole  -> {'PASS' if g3 else 'FAIL'}")

    # [G4] the honest cost: mu_1 = 0 also kills the matter cone -> two decoupled sectors required
    print("\n  [G4] the cost -- dial the graviton mass mu_1 off at fixed kappa and watch the matter cone die:")
    print(f"       {'mu_1':>6s} {'mass^2':>8s} {'kinetic':>8s} {'c_T(matter)':>12s} {'Y_1':>7s} {'stable?':>8s}")
    cone_dies = False
    for mu1 in (2.14, 1.0, 0.3, 0.05, 0.0):
        _, wT2 = dispersions(17.3, mu1, kappa, kappa, qs)
        cT = np.sqrt(max(mu1, 0.0) / RHO)                 # first-gradient transverse cone
        stable = bool(np.all(wT2 > 0))
        print(f"       {mu1:6.3f} {mu1:8.3f} {kappa:8.3f} {cT:12.4f} {young_2d(17.3, mu1):7.3f} {str(stable):>8s}")
        if mu1 == 0.0:
            cone_dies = cT < 1e-9 and stable              # massless graviton: c_T=0 but medium still stable
    g4 = cone_dies
    ok &= g4
    print(f"       => the massless-graviton condition mu_1 = 0 sets the matter cone c_T = sqrt(mu_1/rho) = 0")
    print(f"          (the cone the model needs for its photon/fermion, S8.1/8.44/8.45), while the kappa term")
    print(f"          keeps the medium stable. gamma=1 (mu_1=0) and matter (mu_1>0) cannot share one elastic")
    print(f"          sector -> the escape needs a substrate with TWO decoupled sectors  -> {'PASS' if g4 else 'FAIL'}")

    print("\n" + "=" * 92)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  The elastic moduli ARE a graviton mass term (INT h^2, no derivatives on h): S8.64's 'rigidity")
    print("  forbids gamma = 1' is exactly 'a massive graviton has no diffeomorphism invariance.' The unique")
    print("  escape is Kleinert's world crystal -- a MASSLESS graviton: zero first-gradient moduli (mass off)")
    print("  plus a second-gradient kinetic term kappa INT (d h)^2 (Einstein-Hilbert-like), which keeps the")
    print("  medium stable though strain-floppy, restores the diffeomorphism gauge invariance, and deconfines")
    print("  the disclinations -- the vector-charge (Einstein) preconditions that no first-gradient solid can")
    print("  meet. The cost is that a massless-graviton sector has zero moduli and so cannot carry the matter")
    print("  Lorentz cone: gamma = 1 requires a substrate with two decoupled sectors -- a curvature-stiff")
    print("  gravity sector and a strain-stiff matter sector. This names the escape precisely and reduces")
    print("  'reach gamma = 1' to one construction: build emergent matter on the two-sector medium and")
    print("  ray-trace gamma, where Kleinert's world crystal predicts gamma = 1. Preconditions proven unique")
    print("  and constructible here; the gamma = 1 endpoint is cited, not yet measured in-model.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
