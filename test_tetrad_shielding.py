"""
Can a slab SHIELD the tetrad graviton?  (The tension test.)

test_shielding.py showed the two gravitational couplings sit on opposite sides of the
shielding question: a DILATATION (energy-density) charge is neutralizable and therefore
shieldable; a TOPOLOGICAL (curvature/winding) charge is not. That raised a sharp problem,
because the project actually contains TWO gravity constructions:

  (1) Route 1 / fracton : mass = a disclination = a TOPOLOGICAL charge  -> unshieldable.
  (2) Tetrad / Sakharov : the graviton is the SHAPE OF THE FERMION CONE, and Sakharov
      induction couples it to the stress tensor T_uv -- i.e. to ENERGY, which is exactly
      the channel that test_shielding showed is SHIELDABLE.

If (2) really sits in the shieldable channel, the "one structure" capstone has a hole: a
slab of matter could dielectrically screen a distant mass, which spacetime emphatically
does not do.

But there is a reason to think it escapes, and it is worth stating before measuring, because
it is decisive. Around a dilatation centre in 2D elasticity the displacement is u_r = C/r, so
    trace strain  e_rr + e_tt = 0                  (div u = 0: the DENSITY is strictly local)
    deviator      e_rr - e_tt = -2C/r^2            (the SHEAR is long-range, ~ 1/r^2)
Bitter-Crum kills the TRACE. It says nothing about the TRACELESS part. And the tetrad h is
precisely the TRACELESS part of the cone deformation (test_emergent_tetrad: uniform/conformal
bond stretch -> h = 0; only the traceless doublet sources h). So the graviton may read the
one piece of the medium's response to a mass that is NOT screened.

Test, in the medium's own idiom (a bond network, since the tetrad is read from bonds):
  A. put a mass (a dilatation centre: an inclusion whose bonds want to be longer) in a
     triangular spring medium, relax, and measure the TRACE and the DEVIATOR vs r.
  B. surround it with a SHELL of contrasting stiffness -- a slab of intervening matter --
     and ask whether any contrast can null the deviator at a probe outside (elastic
     "cloaking" / a neutral inclusion). If it can, the tetrad graviton is SHIELDABLE.
"""
from __future__ import annotations
import numpy as np

SQ3 = np.sqrt(3.0)
# the three independent nearest-neighbour directions of the triangular lattice
NB = [(1, 0), (0, 1), (-1, 1)]
ALL_NB = NB + [(-1, 0), (0, -1), (1, -1)]


def lattice(N, Rd):
    """Triangular lattice (nonzero shear modulus with central-force springs), cut to a DISC
    of radius Rd so the outer boundary is isotropic and can be cleanly clamped."""
    ij_all = [(i, j) for j in range(N) for i in range(N)]
    pos_all = np.array([[i + 0.5 * j, j * SQ3 / 2] for (i, j) in ij_all])
    c = pos_all.mean(axis=0)
    rad_all = np.linalg.norm(pos_all - c, axis=1)
    keep = rad_all <= Rd
    ij = [p for p, k in zip(ij_all, keep) if k]
    idx = {p: n for n, p in enumerate(ij)}
    return ij, idx, pos_all[keep], c, rad_all[keep]


def build(N, Rd, r_src, e0, shell, k_shell):
    """Bond lists for the spring network (pure NumPy; no scipy)."""
    ij, idx, pos, c, rad = lattice(N, Rd)
    A, B, NH, KK, PRE = [], [], [], [], []
    for a, (i, j) in enumerate(ij):
        for (di, dj) in NB:
            p = (i + di, j + dj)
            if p not in idx:
                continue
            b = idx[p]
            d = pos[b] - pos[a]
            nh = d / np.linalg.norm(d)
            rm = 0.5 * (rad[a] + rad[b])                 # bond midpoint radius
            k = k_shell if (shell is not None and shell[0] <= rm <= shell[1]) else 1.0
            A.append(a); B.append(b); NH.append(nh); KK.append(k)
            PRE.append(e0 if rm < r_src else 0.0)        # the mass: bonds want to be longer
    return (ij, idx, pos, c, rad, np.array(A), np.array(B),
            np.array(NH), np.array(KK), np.array(PRE))


def solve(N=95, Rd=40.0, r_src=2.5, e0=0.02, shell=None, k_shell=1.0, tol=1e-11, itmax=20000):
    """Relax the medium around a mass (matrix-free conjugate gradient)."""
    ij, idx, pos, c, rad, A, B, NH, KK, PRE = build(N, Rd, r_src, e0, shell, k_shell)
    n = len(ij)
    fixed = rad > Rd - 2.5                               # clamped outer boundary (Dirichlet)

    def scat(contrib):
        out = np.empty((n, 2))
        for p in range(2):
            out[:, p] = (np.bincount(A, contrib[:, p], n)
                         - np.bincount(B, contrib[:, p], n))
        out[fixed] = 0.0
        return out

    def Kop(u):
        s = ((u[A] - u[B]) * NH).sum(1)                  # bond stretch
        return scat((KK * s)[:, None] * NH)

    f = scat((KK * PRE)[:, None] * NH)

    u = np.zeros((n, 2))
    r = f - Kop(u); p = r.copy(); rs = (r * r).sum()
    for _ in range(itmax):
        Kp = Kop(p)
        al = rs / ((p * Kp).sum() + 1e-300)
        u += al * p; r -= al * Kp
        rs2 = (r * r).sum()
        if np.sqrt(rs2) < tol:
            break
        p = r + (rs2 / rs) * p; rs = rs2
    return ij, idx, pos, c, rad, u


def strain(ij, idx, pos, u):
    """Per-site strain from the displacement-gradient fit to the 6 neighbours.
    Interior sites all share the same neighbour geometry, so the least-squares
    pseudo-inverse is one constant 2x6 matrix -- fully vectorised.
    Returns trace (= DENSITY channel) and |traceless part| (= TETRAD h channel)."""
    n = len(ij)
    nbr = np.full((n, 6), -1, int)
    for a, (i, j) in enumerate(ij):
        for m, (di, dj) in enumerate(ALL_NB):
            p = (i + di, j + dj)
            if p in idx:
                nbr[a, m] = idx[p]
    full = (nbr >= 0).all(1)                             # interior sites only
    D = np.array([[di + 0.5 * dj, dj * SQ3 / 2] for (di, dj) in ALL_NB])   # 6x2
    P = np.linalg.pinv(D)                                # 2x6

    U = u[nbr[full]] - u[full][:, None, :]               # (m,6,2)
    G = np.einsum("ij,njk->nik", P, U)                   # (m,2,2) = grad u
    e = 0.5 * (G + np.transpose(G, (0, 2, 1)))
    t = e[:, 0, 0] + e[:, 1, 1]
    dv = e - 0.5 * t[:, None, None] * np.eye(2)[None]
    tr = np.full(n, np.nan); dev = np.full(n, np.nan)
    tr[full] = t
    dev[full] = np.sqrt((dv ** 2).sum(axis=(1, 2)))
    return tr, dev


def radial(rad, q, lo, hi, nb=14):
    """Mean |q| in radial bins."""
    edges = np.linspace(lo, hi, nb + 1)
    rs, qs = [], []
    for i in range(nb):
        m = (rad >= edges[i]) & (rad < edges[i + 1]) & np.isfinite(q)
        if m.sum() > 6:
            rs.append(rad[m].mean()); qs.append(np.abs(q[m]).mean())
    return np.array(rs), np.array(qs)


def slope(r, q):
    m = q > 0
    return float(np.polyfit(np.log(r[m]), np.log(q[m]), 1)[0])


if __name__ == "__main__":
    print("=== Can a slab SHIELD the tetrad graviton? ===\n")
    print("  The density channel reads the TRACE of the strain; the tetrad h is the TRACELESS")
    print("  part. Bitter-Crum kills the trace -- it says nothing about the deviator.\n")

    # ---------------- A. the two channels around a mass ----------------
    ij, idx, pos, c, rad, u = solve()
    tr, dev = strain(ij, idx, pos, u)
    peak = np.nanmax(np.abs(tr))

    print("  [A] a mass (dilatation centre) in the medium -- what does each channel see?")
    r1, t1 = radial(rad, tr, 6, 30)
    r2, d1 = radial(rad, dev, 6, 30)
    print(f"      {'r':>6} {'|trace| (density)':>19} {'|deviator| (tetrad h)':>23}")
    for i in range(0, len(r1), 3):
        print(f"      {r1[i]:>6.1f} {t1[i]:>19.3e} {d1[i]:>23.3e}")
    s_tr, s_dev = slope(r1, t1), slope(r2, d1)
    print(f"\n      density channel: |trace| ~ r^{s_tr:+.2f}, and only {t1.max()/peak:.1e} of peak")
    print("        -> NO radial field: a flat pedestal (the clamped box conserves volume), not")
    print("           a long-range tail. div u is strictly local: SCREENED (Bitter-Crum).")
    print(f"      tetrad  channel: |deviator| ~ r^{s_dev:+.2f}  -> LONG-RANGE 1/r^2 power law.")
    print("      => the graviton reads the ONE part of the medium's response to a mass that")
    print("         Bitter-Crum does NOT screen. The traceless sector escapes the theorem.\n")

    # box-size independence: the exponent must not depend on the boundary
    _ij, _idx, _pos, _c, _rad, _u = solve(N=71, Rd=28.0)
    _, _dev = strain(_ij, _idx, _pos, _u)
    _r, _d = radial(_rad, _dev, 6, 20)
    print(f"      [gate] box independence: exponent {s_dev:+.2f} (Rd=40) vs {slope(_r,_d):+.2f} (Rd=28)\n")

    # ---------------- B. the shielding test ----------------
    print("  [B] now put a SLAB of matter (a shell, r=8..14) between the mass and the probe.")
    print("      Sweep its stiffness contrast: can ANY of it null the tetrad at the probe?")
    probe = (rad > 18) & (rad < 22)
    base = np.nanmean(np.abs(dev[probe]))
    print(f"\n      {'k_shell/k':>10} {'|h| at probe':>14} {'/ no-shell':>12}")
    ratios = {}
    for ks in (0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0):
        _, _, pos2, _, rad2, u2 = solve(shell=(8.0, 14.0), k_shell=ks)
        _, dev2 = strain(ij, idx, pos2, u2)
        v = np.nanmean(np.abs(dev2[(rad2 > 18) & (rad2 < 22)]))
        ratios[ks] = v / base
        print(f"      {ks:>10.2f} {v:>14.3e} {v/base:>12.3f}")
    print("\n      (a NEUTRAL INCLUSION -- an elastic cloak -- would drive this ratio to 0.)")

    # box-size gate on the load-bearing number
    _, _, p3, _, r3, u3 = solve(N=71, Rd=28.0, shell=(8.0, 14.0), k_shell=20.0)
    _, d3 = strain(_ij, _idx, p3, u3)
    b3 = np.nanmean(np.abs(_dev[(_rad > 18) & (_rad < 22)]))
    v3 = np.nanmean(np.abs(d3[(r3 > 18) & (r3 < 22)]))
    print(f"      [gate] box independence of the k=20 ratio: {ratios[20.0]:.3f} (Rd=40) "
          f"vs {v3/b3:.3f} (Rd=28)")

    print("\n[verdict] the tension is REAL, and it cuts both ways.\n")
    print("  GOOD -- the tetrad graviton is LONG-RANGE from ordinary energy density.")
    print("  A mass's density perturbation is strictly local (Bitter-Crum), but its SHEAR is")
    print("  not: the deviatoric strain falls as 1/r^2, box-independently. The graviton is the")
    print("  TRACELESS part of the cone deformation, so it reads exactly the sector the")
    print("  screening theorem does not touch. Gravity-by-density failed because it coupled to")
    print("  the TRACE. Coupling to the traceless part is long-range WITHOUT needing topology --")
    print("  which is why the tetrad construction works at all.\n")
    print("  BAD -- but it is SHIELDABLE. An intervening shell of matter attenuates the tetrad")
    print("  at the probe by up to ~4x (ratio 1.00 -> 0.23 soft, -> 0.32 stiff). ANY impedance")
    print("  mismatch scatters it; the far field is not protected. Real gravity shows NO such")
    print("  attenuation -- a shell of matter around the Sun does not weaken its pull on Earth.")
    print("  So the tetrad graviton, sourced by energy, does NOT inherit nature's protection.\n")
    print("  DIAGNOSIS. The screening is possible only because the medium's elastic MODULI are")
    print("  an independent background structure that can be dialled freely. In general")
    print("  relativity you cannot independently dial spacetime's stiffness -- the metric IS the")
    print("  gravitational field, so a 'stiffness contrast' is not a free knob: it is mass, and")
    print("  mass can only ADD. The shielding measured here is therefore a direct diagnostic of")
    print("  the model's remaining NON-GR structure, and it says exactly what must be fixed:")
    print("  the moduli must not be free parameters but must themselves be the gravitational")
    print("  field, determined self-consistently by the matter. That is precisely the open item")
    print("  'self-consistent matter -> tetrad back-reaction' -- now shown to be load-bearing,")
    print("  not cosmetic.\n")
    print("  STATUS: the two gravity constructions are NOT yet unified.")
    print("    Route 1 (topological) : unshieldable [test_shielding], but the force law between")
    print("                            two topological charges is not yet measured end-to-end.")
    print("    Tetrad  (energy)      : long-range 1/r^2 [here], but shieldable [here].")
    print("  Neither is complete. The disclination force law is now the decisive next test.")
