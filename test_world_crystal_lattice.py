"""Route A, the build -- rung 1: the world crystal as an EXPLICIT lattice (not a continuum dispersion).

S8.65 identified the escape from the gamma=0 no-go (the elastic moduli are a graviton mass; the massless
graviton needs a strain-floppy, curvature-stiff medium), and S8.66 measured gamma=1 on it at the linear-
response level. Both wrote the world-crystal dispersion down by hand (omega_T^2 = kappa q^4). This section
BUILDS that medium as an actual lattice and measures its graviton on the lattice, so the substrate is a real
thing, not an assumed form.

THE LATTICE. A triangular lattice of nodes with displacement u(n), whose energy penalises the CURVATURE of
atomic rows -- for every collinear triple (n-e, n, n+e) along each of the 3 lattice directions e,
        E_grav = (kappa/2) SUM_n SUM_e | u(n+e) - 2 u(n) + u(n-e) |^2,
and NO nearest-neighbour stretching springs. A uniform (affine) displacement u(n) = A x_n has vanishing
second difference, so it costs ZERO energy: every homogeneous strain is a zero mode -- the moduli vanish, so
by S8.65 the graviton is MASSLESS. But any non-affine (curved) deformation costs kappa > 0, so the medium is
still a stable solid -- strain-floppy, curvature-stiff. The control is the ordinary first-gradient medium: the
same triangular lattice with central (stretching) springs, E_el = (k/2) SUM |[u(n+e)-u(n)].e|^2, whose
homogeneous strains cost the moduli (a massive graviton, gamma=0 -- S8.64).

  [G1] Massless graviton: the second-gradient lattice has omega^2(q) ~ q^4 as q->0 (a massless, two-derivative
       graviton dispersion), versus omega^2 ~ q^2 for the first-gradient control. Measured by the log-log slope
       of the dispersion at small q: 4 vs 2.
  [G2] The moduli vanish (the graviton mass is zero): omega^2(q)/q^2 -> 0 for the second-gradient lattice (no
       q^2 stiffness -- zero shear and bulk modulus), while it -> a finite modulus for the first-gradient
       control. This is S8.65's "moduli = graviton mass" measured directly on the lattice: zero moduli = massless.
  [G3] Yet stable: the second-gradient dynamical matrix is positive semi-definite across the whole Brillouin
       zone (no negative modes), with zero modes only the affine ones -- a strain-floppy but curvature-stiff
       stable solid. Corroborated in real space: a uniform shear applied to a finite patch costs ~0 on the
       second-gradient lattice and a finite energy on the first-gradient one.

Honest scope: this builds the GRAVITY sector of the two-sector substrate (S8.66) as an explicit lattice, and
verifies its defining properties -- massless (q^4), zero moduli, stable, strain-floppy/curvature-stiff. The
isotropic |d^2 u|^2 energy gives both graviton polarisations the same q^4 dispersion; matching the exact
Einstein-Hilbert tensor structure on the lattice (so that gamma=1 follows on the lattice as it did in S8.66's
linear response) is a refinement for the next rung, as is adding the strain-stiff matter sector and its back-
reaction. What is shown here: the curvature-stiff, massless-graviton substrate is a concrete lattice. Pure numpy.
"""
from __future__ import annotations
import numpy as np

SQ3 = np.sqrt(3.0)
# triangular lattice: 3 bond directions (each with +/- neighbours -> 6 nearest neighbours)
DIRS = [np.array([1.0, 0.0]),
        np.array([0.5, SQ3 / 2]),
        np.array([-0.5, SQ3 / 2])]
RHO = 1.0


def D_secondgradient(q, kappa=1.0):
    """Bloch dynamical matrix of the row-curvature (second-gradient) lattice. Isotropic in polarisation."""
    s = 0.0
    for e in DIRS:
        f = 2.0 * (1.0 - np.cos(q @ e))          # |second difference| factor along e
        s += f * f
    return kappa * s * np.eye(2)                 # |d^2 u|^2 couples both components equally


def D_firstgradient(q, k=1.0):
    """Bloch dynamical matrix of the central-force (first-gradient) triangular lattice -- the massive control."""
    D = np.zeros((2, 2))
    for e in DIRS:
        f = 2.0 * (1.0 - np.cos(q @ e))          # |first difference| factor along e
        D += k * f * np.outer(e, e)
    return D


def min_omega2(Dfun, q):
    return np.linalg.eigvalsh(Dfun(q)).min() / RHO


def slope_loglog(Dfun, direction, qs):
    """Log-log slope of min omega^2 vs |q| at small q (the dispersion exponent)."""
    d = direction / np.linalg.norm(direction)
    w = np.array([min_omega2(Dfun, qq * d) for qq in qs])
    return np.polyfit(np.log(qs), np.log(w), 1)[0]


def modulus_q2(Dfun, direction, q0=1e-3):
    """omega^2 / q^2 at small q -- the q^2 stiffness (a modulus; the graviton mass^2 of S8.65)."""
    d = direction / np.linalg.norm(direction)
    return min_omega2(Dfun, q0 * d) / q0**2


def bz_min_eig(Dfun, ngrid=60):
    """Minimum eigenvalue of the dynamical matrix over the Brillouin zone (stability)."""
    b1 = 2 * np.pi * np.array([1.0, -1 / SQ3])
    b2 = 2 * np.pi * np.array([0.0, 2 / SQ3])
    lo = np.inf
    for i in range(ngrid):
        for j in range(ngrid):
            q = (i / ngrid) * b1 + (j / ngrid) * b2
            if np.linalg.norm(q) < 1e-6:
                continue
            lo = min(lo, np.linalg.eigvalsh(Dfun(q)).min())
    return lo


def uniform_shear_energy(second_gradient, npatch=6, amp=0.02):
    """Real-space: energy of a uniform shear u=(gamma*y, 0) on a finite triangular patch (open boundaries)."""
    pts = {}
    for i in range(npatch):
        for j in range(npatch):
            pts[(i, j)] = i * DIRS[0] + j * DIRS[1]
    A = np.array([[0.0, amp], [0.0, 0.0]])                    # simple shear
    u = {k: A @ x for k, x in pts.items()}
    E = 0.0
    for (i, j), x in pts.items():
        for e_idx, (di, dj) in enumerate([(1, 0), (0, 1), (-1, 1)]):   # the 3 lattice directions in index space
            p, m = (i + di, j + dj), (i - di, j - dj)
            if p in pts and m in pts:
                if second_gradient:
                    d2 = u[p] - 2 * u[(i, j)] + u[m]           # second difference
                    E += 0.5 * (d2 @ d2)
                else:
                    e = DIRS[e_idx]
                    d1 = (u[p] - u[(i, j)]) @ e                # central-force first difference
                    E += 0.5 * d1 * d1
    return E


def main():
    print("=" * 92)
    print("ROUTE A (the build) -- RUNG 1: the world crystal as an EXPLICIT lattice (massless graviton)")
    print("=" * 92)
    ok = True
    qs = np.array([0.02, 0.03, 0.05, 0.08, 0.12])
    dirs = {"Gamma-K (1,0)": DIRS[0], "Gamma-M (0,1)": np.array([0.0, 1.0])}

    # [G1] massless graviton: omega^2 ~ q^4 (slope 4) vs first-gradient q^2 (slope 2)
    print("\n  [G1] dispersion exponent (log-log slope of omega^2 vs |q| at small q):")
    print(f"       {'direction':<16s} {'2nd-gradient (world crystal)':>30s} {'1st-gradient (control)':>26s}")
    g1 = True
    for name, d in dirs.items():
        s2 = slope_loglog(D_secondgradient, d, qs)
        s1 = slope_loglog(D_firstgradient, d, qs)
        print(f"       {name:<16s} {s2:>30.3f} {s1:>26.3f}")
        g1 &= abs(s2 - 4.0) < 0.15 and abs(s1 - 2.0) < 0.15
    ok &= g1
    print(f"       => world crystal: omega^2 ~ q^4 (massless two-derivative graviton); control: omega^2 ~ q^2")
    print(f"          -> {'PASS' if g1 else 'FAIL'}")

    # [G2] the moduli vanish (graviton mass = 0)
    print("\n  [G2] the q^2 stiffness omega^2/q^2 at small q -- the modulus, i.e. the graviton mass^2 (S8.65):")
    m2_wc = modulus_q2(D_secondgradient, DIRS[0])
    m2_el = modulus_q2(D_firstgradient, DIRS[0])
    print(f"       2nd-gradient (world crystal): omega^2/q^2 = {m2_wc:.3e}  -> ~0 : NO q^2 stiffness (moduli=0)")
    print(f"       1st-gradient (control):       omega^2/q^2 = {m2_el:.3e}  -> finite modulus (massive graviton)")
    g2 = m2_wc < 1e-4 and m2_el > 0.1
    ok &= g2
    print(f"       => zero moduli = massless graviton, directly on the lattice (S8.65 'moduli = graviton mass')")
    print(f"          -> {'PASS' if g2 else 'FAIL'}")

    # [G3] stable despite being strain-floppy
    print("\n  [G3] stability -- minimum dynamical-matrix eigenvalue over the Brillouin zone (q != 0):")
    lo_wc = bz_min_eig(D_secondgradient)
    lo_el = bz_min_eig(D_firstgradient)
    Esh_wc = uniform_shear_energy(second_gradient=True)
    Esh_el = uniform_shear_energy(second_gradient=False)
    print(f"       2nd-gradient (world crystal): min eig over BZ = {lo_wc:.4f}  (>= 0: STABLE, no negative modes)")
    print(f"       1st-gradient (control):       min eig over BZ = {lo_el:.4f}")
    print(f"       real-space uniform-shear energy on a finite patch: world crystal = {Esh_wc:.2e}, "
          f"control = {Esh_el:.2e}")
    g3 = lo_wc > -1e-9 and Esh_wc < 1e-9 and Esh_el > 1e-6
    ok &= g3
    print(f"       => the world-crystal lattice is stable (min eig >= 0) yet a uniform shear costs ~0 -- a")
    print(f"          strain-floppy but curvature-stiff stable solid, exactly S8.65's medium  -> {'PASS' if g3 else 'FAIL'}")

    print("\n" + "=" * 92)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  The world crystal is a concrete lattice, not an assumed dispersion. A triangular lattice whose")
    print("  energy penalises row curvature (second differences) with no stretching springs has: a massless")
    print("  graviton (omega^2 ~ q^4, slope 4 vs the control's 2); vanishing moduli (omega^2/q^2 -> 0), i.e. a")
    print("  zero graviton mass measured directly on the lattice (S8.65); and full stability across the")
    print("  Brillouin zone despite a uniform shear costing zero -- strain-floppy, curvature-stiff. This is the")
    print("  gravity sector of S8.66's two-sector substrate, now built. Next rungs: match the exact Einstein-")
    print("  Hilbert tensor structure on the lattice, add the strain-stiff matter sector, and ray-trace gamma")
    print("  end-to-end with self-consistent back-reaction.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
