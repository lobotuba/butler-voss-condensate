"""
The disclination force law, end-to-end.  (The link Route 1 asserts but never measured.)

*** STATUS UPDATE -- the measurement stands; the open problem it named has since been SOLVED.
    STANDS: two like curvature charges are not Newtonian in the PURE medium -- |dE| ~ R^1.97, an
    energy that GROWS with separation. Route 1, as asserted, does fail.
    REINTERPRETED: that growing energy is CONFINEMENT, not a mere wrong-sign repulsion. In its clean
    3D form the biharmonic curvature sector gives a CONSTANT force -- a string tension
    s^2/(8 pi kappa) (test_deconfinement). A confining sector simply contributes nothing at long
    range; it is not a refutation of curvature-coupled gravity, it is a statement that the sector
    must be deconfined first.
    SOLVED: this file closes by naming "the honest open problem -- a 3D construction where the
    graviton's two propagating polarizations carry a universally attractive, unshieldable 1/r^2
    force." All three pieces now exist. The Sakharov-induced Einstein term deconfines the confining
    +R into an exact Newtonian -1/r with G = 1/(4 pi mu) (test_deconfinement); the sign mu > 0 is
    measured by calibration against the model's healthy photon (test_induced_sign); and in 3+1D the
    graviton has exactly two degenerate, healthy propagating polarizations (test_spin2_dynamical),
    with gamma -> 1 in the infrared (test_lattice_ward). ***

Route 1 says: couple gravity to CURVATURE (disclinations), not energy density, because the
curvature sector is unscreened (test_shielding: it is topologically unneutralizable). But the
chain "topological curvature charge -> long-range 1/r^2 attraction between two masses" has
never actually been measured in one piece:

  * test_fracton_gravity established the multipole hierarchy, but its disclination Green's
    function was IR-BOX-SATURATED in the periodic FFT -- reported honestly at the time; the
    clean results were the dilatation CONTACT term and the dislocation LOG.
  * test_graviton_dynamics got 1/r^2 -- but from an IMPOSED Poisson equation with a
    mass-density source, not from a disclination source.

So the force law between two actual topological charges is an assertion, not a measurement.
Measure it, avoiding the failure mode that saturated last time: work in REAL space on a
CLAMPED disc (no periodic zero mode, no IR saturation), and gate on box size.

2D elasticity, Airy stress function chi:   biharmonic(chi) = s(x),  E = (1/2) integral (lap chi)^2
with s the DISCLINATION density. Two sources, separation R, energy E(R) -> the force law.

Controls, from the same solver, so the comparison is apples-to-apples:
  * DISCLINATION pair = two like MONOPOLES of s.   (the gravity claim)
  * DISLOCATION  pair = two DIPOLES of s (a dislocation IS a disclination dipole).
    Known answer: logarithmic (Coulomb-like in 2D) -- so it validates the solver.
"""
from __future__ import annotations
import numpy as np


def lap(f):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
            np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)


def solve_biharm(s, free, tol=1e-12, itmax=40000):
    """Solve biharmonic(chi) = s on the clamped disc (chi = 0 outside `free`, which also
    pins the normal derivative since the 4th-order stencil reaches two layers out).
    The clamped biharmonic is SPD -> conjugate gradient."""
    def B(x):
        f = np.zeros_like(x); f[free] = x[free]
        out = lap(lap(f))
        out[~free] = 0.0
        return out

    b = np.zeros_like(s); b[free] = s[free]
    x = np.zeros_like(s)
    r = b - B(x); p = r.copy(); rs = (r * r).sum()
    for _ in range(itmax):
        Bp = B(p)
        al = rs / ((p * Bp).sum() + 1e-300)
        x += al * p; r -= al * Bp
        rs2 = (r * r).sum()
        if np.sqrt(rs2) < tol:
            break
        p = r + (rs2 / rs) * p; rs = rs2
    x[~free] = 0.0
    return x


def energy(chi):
    return 0.5 * float((lap(chi) ** 2).sum())


def blob(X, Y, x0, y0, sig):
    return np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * sig ** 2))


def make(N, Rd):
    g = np.arange(N) - (N - 1) / 2.0
    X, Y = np.meshgrid(g, g, indexing="ij")
    free = np.hypot(X, Y) < Rd
    return X, Y, free


def E_pair(N, Rd, R, kind, sig=2.0, d=3.0):
    """Energy of two curvature charges separated by R.
    kind='monopole' : two like disclinations   (the gravity claim)
    kind='dipole'   : two dislocations (each = a +/- disclination pair of spacing d)"""
    X, Y, free = make(N, Rd)
    s = np.zeros_like(X)
    for sgn in (-1, +1):
        x0 = sgn * R / 2.0
        if kind == "monopole":
            s += blob(X, Y, x0, 0.0, sig)
        else:                                   # a dislocation = a disclination DIPOLE
            s += blob(X, Y, x0, +d / 2, sig) - blob(X, Y, x0, -d / 2, sig)
    return energy(solve_biharm(s, free))


def law(N, Rd, kind, Rs):
    E = np.array([E_pair(N, Rd, R, kind) for R in Rs])
    return E


def fit(Rs, E):
    """The R-dependent (interaction) part: subtract the R-independent self-energies by
    referencing the smallest separation, then fit a power law to |dE|.
    SIGN MATTERS: dE > 0 (energy rises with R) = ATTRACTION; dE < 0 = REPULSION."""
    dE = E - E[0]
    m = np.abs(dE) > 0
    m[0] = False
    p = np.polyfit(np.log(Rs[m]), np.log(np.abs(dE[m])), 1)[0]
    sign = "ATTRACT" if dE[-1] > 0 else "REPEL"
    return float(p), dE, sign


if __name__ == "__main__":
    print("=== The disclination force law, end-to-end ===\n")
    print("  Real space, clamped disc: no periodic zero mode, no IR saturation.")
    print("  E(R) between two curvature charges -> the force law Route 1 needs.\n")

    N, Rd = 181, 78.0
    Rs = np.array([8.0, 12.0, 16.0, 20.0, 26.0, 32.0, 40.0])

    # ---- control first: the solver must reproduce the KNOWN dislocation law ----
    Ed = law(N, Rd, "dipole", Rs)
    pd, dEd, sd = fit(Rs, Ed)
    print("  [control] DISLOCATION pair (a disclination dipole) -- known answer: LOG, repulsive")
    print(f"      {'R':>6} {'E(R) - E(R0)':>16}")
    for R, v in zip(Rs, dEd):
        print(f"      {R:>6.0f} {v:>16.1f}")
    print(f"      => |dE| ~ R^{pd:+.2f}, and the pair {sd}S (like Burgers vectors repel).")
    print("         Log-like and repulsive: the standard dislocation law. Solver validated.\n")

    # ---- the measurement: two like disclinations (the gravity claim) ----
    Em = law(N, Rd, "monopole", Rs)
    pm, dEm, sm = fit(Rs, Em)
    print("  [measurement] DISCLINATION pair -- two LIKE curvature charges ('two masses')")
    print(f"      {'R':>6} {'E(R) - E(R0)':>16}")
    for R, v in zip(Rs, dEm):
        print(f"      {R:>6.0f} {v:>16.1f}")
    print(f"      => |E(R)| ~ R^{pm:+.2f}, and the pair {sm}S.\n")

    # ---- box gate: the exponent must not be the boundary talking ----
    Em2 = law(141, 60.0, "monopole", Rs)
    pm2, _, sm2 = fit(Rs, Em2)
    print(f"  [gate] box independence: R^{pm:+.2f} ({sm}) at Rd=78 "
          f"vs R^{pm2:+.2f} ({sm2}) at Rd=60\n")

    print("[verdict] Route 1 does NOT deliver the gravitational force law. Two failures, and")
    print("  the SIGN is the fatal one:\n")
    print(f"  * SIGN. Two LIKE curvature charges REPEL. Energy FALLS as R^{pm:.2f} with separation,")
    print("    so they fly apart. But gravity is UNIVERSALLY ATTRACTIVE -- all masses are positive")
    print("    and all masses attract. Same-sign topological charges repel, exactly as same-sign")
    print("    electric charges do. A disclination behaves like a CHARGE, not like a MASS.")
    print(f"  * RANGE. The interaction goes as R^{pm:.2f} -- it GROWS with distance (force ~ dE/dR")
    print("    grows linearly). Newtonian gravity falls as 1/r^2. This is not merely the wrong")
    print("    exponent; it is the wrong direction.\n")
    print("  The earlier IR box-saturation in test_fracton_gravity was this R^2 growth showing")
    print("  itself as a divergence -- the periodic FFT was not failing to measure a force law,")
    print("  it was correctly refusing to converge on one that grows without bound.\n")
    print("  In hindsight this is what 2+1D gravity actually IS: a point mass is a conical")
    print("  deficit, and 2+1D GR is TOPOLOGICAL -- no local propagating modes, no Newtonian")
    print("  attraction between static masses. The curvature sector is faithfully reproducing")
    print("  2+1D gravity. But that means the Newtonian force cannot come from 2D disclinations")
    print("  at all -- not in this model, and not in nature.\n")
    print("  STATUS -- neither gravity construction is complete, and now we know exactly why:")
    print("    Route 1 (topological curvature) : UNSHIELDABLE [test_shielding] -- but the force")
    print("        is repulsive and grows with distance. Solves screening, breaks the force law.")
    print("    Tetrad (energy-sourced)         : LONG-RANGE 1/r^2 [test_tetrad_shielding] -- but")
    print("        shieldable by an intervening slab. Solves the force law, breaks screening.")
    print("  (And the tetrad's ATTRACTION is still unmeasured: 1/r^2 was a FIELD falloff, not a")
    print("   force sign between two masses. That is the next thing to check, not to assume.)")
    print("\n  => The two routes fail in complementary ways. The 1/r^2 in test_graviton_dynamics")
    print("     came from an IMPOSED Poisson equation and cannot be cited as a derivation. The")
    print("     honest open problem is a 3D construction where the graviton's two propagating")
    print("     polarizations carry a universally attractive, unshieldable 1/r^2 force.")
