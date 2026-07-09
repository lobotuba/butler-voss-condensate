"""
Unified prototype: one operator for medium AND field -> a single universal cone.

test_lorentz.py found the obstacle: the field wave speed is a free parameter and
the medium's phonons travel at the Lennard-Jones sound speeds, so the two sectors
have different light cones (Lorentz violation between sectors). Worse, a stable
CENTRAL-force solid always has c_L > c_T (c_L^2 - c_T^2 = (K + mu/3)/rho > 0), so
the medium alone already carries TWO cones.

The fix is structural: govern all low-energy dynamics by ONE isotropic operator.
Keep central (LJ) forces to SET the isotropic geometry (self-assembly, H10), but
let the DYNAMICS of both the medium displacement and the matter field be a single
VECTOR-HOOKE graph Laplacian,
    E = 1/2 K sum_bonds |U_i - U_j|^2,   U in R^{D+1}  (D medium + 1 field),
which penalises the full relative displacement, not just the bond length. Its
dynamical matrix is diagonal and isotropic, D_ab(k) = delta_ab * K * S(k) with
S(k) = sum_delta (1 - cos k.delta), so every polarisation -- longitudinal,
transverse, and the internal 'field' component -- has the identical dispersion:
one universal cone, c_L = c_T = c_field, by construction. The residual lattice
anisotropy is the same emergent (k a)^2..4 suppression measured in test_lorentz.
"""
from __future__ import annotations
import numpy as np

from bvc_core import perfect_hex, perfect_fcc, R0
from test_lorentz import neighbours, lj_derivs


def central_speeds(X, kdir, rcut=2.5, kmag=1e-3):
    """Longitudinal & transverse sound speeds of the CENTRAL-force (LJ) medium."""
    delta = neighbours(X, rcut); dim = X.shape[1]
    kvec = kmag * (np.pi / R0) * np.array(kdir, float) / np.linalg.norm(kdir)
    D = np.zeros((dim, dim))
    for dj in delta:
        r = np.linalg.norm(dj); u = dj / r
        p1, p2 = lj_derivs(r)
        D += (p2 * np.outer(u, u) + (p1 / r) * (np.eye(dim) - np.outer(u, u))) * (1 - np.cos(kvec @ dj))
    w = np.sqrt(np.clip(np.linalg.eigvalsh(D), 0, None)) / np.linalg.norm(kvec)
    return w.max(), w.min()                       # c_L, c_T


def unified_speeds(X, kdir, rcut=1.3 * R0, kmag=1e-3):
    """All polarisation speeds of the VECTOR-HOOKE unified network (medium+field)."""
    delta = neighbours(X, rcut); dim = X.shape[1]
    kvec = kmag * (np.pi / R0) * np.array(kdir, float) / np.linalg.norm(kdir)
    S = (1 - np.cos(delta @ kvec)).sum()
    D = S * np.eye(dim)                            # diagonal & isotropic in polarisation
    w = np.sqrt(np.clip(np.linalg.eigvalsh(D), 0, None)) / np.linalg.norm(kvec)
    return w.max(), w.min()                        # c_L, c_T (equal by construction)


def report(name, X, dirs):
    print(f"\n{name}: light-cone speeds by polarisation (LJ units), along symmetry directions")
    print(f"  {'dir':>7} | {'CENTRAL c_L':>11} {'c_T':>7} {'c_L/c_T':>8} | "
          f"{'UNIFIED c_L':>11} {'c_T':>7} {'c_L/c_T':>9}")
    for d in dirs:
        cL, cT = central_speeds(X, d); uL, uT = unified_speeds(X, d)
        print(f"  {str(d):>7} | {cL:>11.3f} {cT:>7.3f} {cL/cT:>8.3f} | "
              f"{uL:>11.3f} {uT:>7.3f} {uL/uT:>9.6f}")


if __name__ == "__main__":
    print("=== Unified prototype: one operator -> one universal light cone ===")
    fcc = perfect_fcc(radius=12.0)
    report("fcc 3D", fcc, [(1, 0, 0), (1, 1, 0), (1, 1, 1)])

    # the crux, quantified: c_L/c_T is the split between medium cones; for the field
    # sector, the central medium leaves c_field a FREE parameter, the unified one ties it.
    cL, cT = central_speeds(fcc, (1, 0, 0)); uL, uT = unified_speeds(fcc, (1, 0, 0))
    print("\nUniversality verdict:")
    print(f"  CENTRAL-force medium:  c_L/c_T = {cL/cT:.3f}  (two cones already), and the matter")
    print(f"                         field's speed is an INDEPENDENT free parameter -> 3rd cone.")
    print(f"  UNIFIED (vector-Hooke): c_L/c_T = {uL/uT:.6f}  (one cone), and the matter field is")
    print(f"                         the (D+1)-th component of the SAME operator -> c_field = c_L")
    print(f"                         EXACTLY. One universal light cone, no tuning.")
    print("  => Universality is recovered by governing medium AND field with a single isotropic")
    print("     operator on the self-assembled (central-force) geometry: matter and medium are")
    print("     one structure, not a field bolted onto an elastic solid. The low-k anisotropy is")
    print("     the same emergent (E/E_Planck)^2..4 suppression measured in test_lorentz -- so the")
    print("     single cone is also ROUND. (Cost: bonds need shear stiffness; bare central forces")
    print("     give c_L>c_T and cannot. Geometry still comes from central-force self-assembly.)")
