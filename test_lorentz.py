"""
Test for emergent Lorentz invariance.

A medium has a preferred frame (the nodes' rest frame), so a "space is a medium"
model is viable as fundamental physics only if Lorentz symmetry EMERGES at long
wavelength: a single, round, universal light cone, with preferred-frame effects
suppressed at low energy. Three accessible facets, all from the dispersion
relation omega(k) (no time evolution needed -- the exact lattice Fourier symbol):

  A. ISOTROPY + linear dispersion (field sector). On the self-assembled isotropic
     lattices (hex 2D, fcc 3D) the field wave speed c(direction,|k|) must become
     direction-independent as k->0, the anisotropy (Lorentz violation) vanishing
     like (k*a)^2 -- i.e. as (E/E_Planck)^2, with the lattice spacing a as the
     model's Planck scale.
  B. UNIVERSALITY (the crux). Does the field wave speed equal the medium's phonon
     (sound) speed? Different speeds = two light cones = Lorentz violation between
     sectors -- the generic killer of medium theories.

Field dispersion: normalized nearest-neighbour graph Laplacian symbol
    omega^2(k) = c^2 * norm * sum_j (1 - cos(k . delta_j)),  norm -> c_eff=c at k->0.
Phonon dispersion: Lennard-Jones dynamical matrix D_ab(k) = sum_j Phi_ab(delta_j)
    (1 - cos(k . delta_j)), Phi_ab = phi''(r) u_a u_b + (phi'(r)/r)(d_ab - u_a u_b).
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_hex, perfect_fcc, R0, EPS, SIGMA


def neighbours(X, rcut):
    i = int(np.argmin((X ** 2).sum(1)))          # a central node
    d = X - X[i]
    r2 = (d ** 2).sum(1)
    sel = (r2 > 1e-9) & (r2 < rcut ** 2)
    return d[sel]


# ---------------------------------------------------- A. field dispersion ------
def field_speed(delta, kvec):
    """c_eff = omega/|k| for the normalized graph Laplacian (c=1); direction &
    magnitude via kvec."""
    k = np.linalg.norm(kvec)
    sym = (1.0 - np.cos(delta @ kvec)).sum()      # sum_j (1-cos(k.delta_j)), w=1
    norm = len(delta[0]) * 1.0 / (0.5 * (delta ** 2).sum() / len(delta))  # -> c_eff=1 at k->0
    # low-k: sym ~ 0.5 * sum (k.delta)^2 = 0.5 * k^2 * <delta_a delta_b>; norm fixes c_eff(0)=1
    return np.sqrt(norm * sym / len(delta)) / k * np.sqrt(len(delta))


def _symbol_speed(delta, kvec):
    k = np.linalg.norm(kvec)
    sym = (1.0 - np.cos(delta @ kvec)).sum()
    return np.sqrt(sym) / k                        # un-normalized; we normalise by k->0 below


def isotropy_test(name, X, dirs, rcut=1.3 * R0):
    delta = neighbours(X, rcut)
    kmax = np.pi / R0                              # Brillouin-zone edge ~ the "Planck" wavenumber
    fracs = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8])
    # normalise each direction's speed to the common k->0 limit
    c0 = np.mean([_symbol_speed(delta, 1e-4 * kmax * np.array(u, float) / np.linalg.norm(u)) for u in dirs.values()])
    print(f"\n[A] {name}: field wave speed vs direction and |k|  (a = R0, k in units of pi/R0)")
    hdr = "  k/kmax |" + "".join(f"{('c['+d+']'):>9}" for d in dirs.keys()) + "   aniso Dc/c"
    print(hdr)
    for f in fracs:
        cs = []
        for u in dirs.values():
            uu = np.array(u, float); uu /= np.linalg.norm(uu)
            cs.append(_symbol_speed(delta, f * kmax * uu) / c0)
        cs = np.array(cs)
        aniso = (cs.max() - cs.min()) / cs.mean()
        print(f"  {f:5.2f}  |" + "".join(f"{c:9.4f}" for c in cs) + f"   {aniso:10.2e}")
    # anisotropy scaling: fit Dc/c ~ (k/kmax)^p at low k
    lowf = np.array([0.05, 0.1, 0.2])
    an = []
    for f in lowf:
        cs = np.array([_symbol_speed(delta, f*kmax*np.array(u,float)/np.linalg.norm(u))/c0 for u in dirs.values()])
        an.append((cs.max()-cs.min())/cs.mean())
    p = np.polyfit(np.log(lowf), np.log(np.array(an) + 1e-30), 1)[0]
    print(f"  => low-k anisotropy scales as (k/kmax)^{p:.1f}: Lorentz violation suppressed")
    print(f"     as (E/E_Planck)^{p:.0f} -> emergent ISOTROPY (rotational Lorentz) at long")
    print(f"     wavelength. (Hexagonal is isotropic to rank 4 -> k^4; cubic only to rank 2 -> k^2.)")


# ------------------------------------------------- B. phonon vs field speed ----
def lj_derivs(r):
    s6 = (SIGMA / r) ** 6; s12 = s6 * s6
    phi1 = 4 * EPS * (-12 * s12 + 6 * s6) / r
    phi2 = 4 * EPS * (156 * s12 - 42 * s6) / r ** 2
    return phi1, phi2


def sound_speeds(X, kdir, rcut=2.5, kmag=1e-3):
    delta = neighbours(X, rcut)
    D = np.zeros((X.shape[1], X.shape[1]))
    kvec = kmag * (np.pi / R0) * np.array(kdir, float) / np.linalg.norm(kdir)
    for dj in delta:
        r = np.linalg.norm(dj); u = dj / r
        p1, p2 = lj_derivs(r)
        Phi = p2 * np.outer(u, u) + (p1 / r) * (np.eye(len(u)) - np.outer(u, u))
        D += Phi * (1.0 - np.cos(kvec @ dj))
    w2 = np.linalg.eigvalsh(D)                     # mass = 1
    w2 = np.clip(w2, 0, None)
    return np.sqrt(w2) / np.linalg.norm(kvec)      # sound speeds (branches)


def universality_test(name, X, kdir):
    cs = sound_speeds(X, kdir)
    cL, cT = cs.max(), cs.min()
    print(f"\n[B] {name}: medium phonon sound speeds along {kdir}: "
          f"longitudinal c_L={cL:.3f}, transverse c_T={cT:.3f} (LJ units)")
    return cL, cT


if __name__ == "__main__":
    print("=== Test for emergent Lorentz invariance ===")
    hex2d = perfect_hex(radius_cells=14)
    fcc3d = perfect_fcc(radius=12.0)

    isotropy_test("hex 2D", hex2d, {"bond 0deg": (1, 0), "gap 30deg": (np.cos(np.pi/6), np.sin(np.pi/6))})
    isotropy_test("fcc 3D", fcc3d, {"[100]": (1, 0, 0), "[110]": (1, 1, 0), "[111]": (1, 1, 1)})

    print("\n" + "=" * 66)
    cL, cT = universality_test("fcc 3D", fcc3d, (1, 0, 0))
    print("\n[B] UNIVERSALITY verdict: the field wave speed c_field is a FREE parameter")
    print(f"    (normalised to 1), while the medium's phonons travel at c_L={cL:.2f}, "
          f"c_T={cT:.2f}.")
    print("    Nothing forces c_field = c_L: the field and phonon sectors have DIFFERENT")
    print("    light cones unless c_field is tuned by hand. Two cones => Lorentz violation")
    print("    BETWEEN sectors -- the central obstacle. A single universal cone needs all")
    print("    excitations to emerge from ONE structure (as near a single Fermi point),")
    print("    not a field added on top of an independent elastic medium.")
