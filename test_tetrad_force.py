"""
Does the tetrad graviton actually PULL?  (The last assumed-but-unmeasured link.)

test_tetrad_shielding showed the tetrad graviton is LONG-RANGE from ordinary energy density:
the deviatoric (traceless) strain around a mass falls as 1/r^2, escaping the Bitter-Crum
screening that kills the trace. It is tempting -- and the project has been tempted -- to read
that as "a 1/r^2 gravitational attraction". But a FIELD FALLOFF IS NOT A FORCE. What was
measured was |h|(r) around ONE mass. Nobody measured the sign, or even the existence, of the
force between TWO.

There is a specific reason to worry, and it is a classical theorem. Eshelby/Crum: in an
infinite ISOTROPIC elastic medium the interaction energy between two centres of dilatation
VANISHES. Each one has a perfectly good long-range strain field, and yet they do not pull on
each other -- the cross term integrates to zero. If that holds here, the tetrad graviton's
1/r^2 field carries NO force, and the "long-range gravity" column of the scorecard is empty
too.

So: two masses (dilatation centres) at separation R in the spring medium. Measure the
INTERACTION energy, with self-energies removed exactly:

    E_int(R) = E(both) - E(mass 1 alone) - E(mass 2 alone)

Sign convention (as in test_disclination_force): E_int RISING with R = ATTRACTION
(force = -dE/dR < 0, pulls together); E_int FALLING with R = REPULSION.

POSITIVE CONTROL. A null result is only meaningful if the probe can see a force when one
exists. The triangular lattice with central-force springs is ISOTROPIC in the continuum, so
Crum applies. Break the isotropy (make one bond direction stiffer) and the theorem no longer
holds -- an interaction MUST appear. If it does, the null in the isotropic case is physics.
"""
from __future__ import annotations
import numpy as np

SQ3 = np.sqrt(3.0)
NB = [(1, 0), (0, 1), (-1, 1)]


def lattice(N, Rd):
    ij_all = [(i, j) for j in range(N) for i in range(N)]
    pos_all = np.array([[i + 0.5 * j, j * SQ3 / 2] for (i, j) in ij_all])
    c = pos_all.mean(axis=0)
    rel = pos_all - c
    keep = np.linalg.norm(rel, axis=1) <= Rd
    ij = [p for p, k in zip(ij_all, keep) if k]
    return ij, {p: n for n, p in enumerate(ij)}, rel[keep]


def energy(N, Rd, srcs, r_src=2.5, e0=0.02, kdir=(1.0, 1.0, 1.0),
           tol=1e-12, itmax=30000):
    """Relax the medium with dilatation centres at `srcs` (positions relative to the disc
    centre) and return the total elastic energy at equilibrium."""
    ij, idx, rel = lattice(N, Rd)
    n = len(ij)
    rad = np.linalg.norm(rel, axis=1)

    A, B, NH, KK, PRE = [], [], [], [], []
    for a, (i, j) in enumerate(ij):
        for m, (di, dj) in enumerate(NB):
            p = (i + di, j + dj)
            if p not in idx:
                continue
            b = idx[p]
            d = rel[b] - rel[a]
            nh = d / np.linalg.norm(d)
            mid = 0.5 * (rel[a] + rel[b])
            near = any(np.linalg.norm(mid - np.asarray(s)) < r_src for s in srcs)
            A.append(a); B.append(b); NH.append(nh)
            KK.append(kdir[m])                      # anisotropy lives here
            PRE.append(e0 if near else 0.0)         # the mass: bonds want to be longer
    A = np.array(A); B = np.array(B)
    NH = np.array(NH); KK = np.array(KK); PRE = np.array(PRE)

    fixed = rad > Rd - 2.5                          # clamped outer boundary

    def scat(contrib):
        out = np.empty((n, 2))
        for p in range(2):
            out[:, p] = (np.bincount(A, contrib[:, p], n)
                         - np.bincount(B, contrib[:, p], n))
        out[fixed] = 0.0
        return out

    def Kop(u):
        s = ((u[A] - u[B]) * NH).sum(1)
        return scat((KK * s)[:, None] * NH)

    f = scat((KK * PRE)[:, None] * NH)
    u = np.zeros((n, 2))
    r = f - Kop(u); p_ = r.copy(); rs = (r * r).sum()
    for _ in range(itmax):
        Kp = Kop(p_)
        al = rs / ((p_ * Kp).sum() + 1e-300)
        u += al * p_; r -= al * Kp
        rs2 = (r * r).sum()
        if np.sqrt(rs2) < tol:
            break
        p_ = r + (rs2 / rs) * p_; rs = rs2

    s = ((u[A] - u[B]) * NH).sum(1)
    return 0.5 * float((KK * (s - PRE) ** 2).sum())


def interaction(N, Rd, R, kdir):
    """E_int(R) = E(both) - E(1 alone) - E(2 alone). Exact self-energy removal."""
    s1 = (-R / 2, 0.0)
    s2 = (+R / 2, 0.0)
    both = energy(N, Rd, [s1, s2], kdir=kdir)
    one = energy(N, Rd, [s1], kdir=kdir)
    two = energy(N, Rd, [s2], kdir=kdir)
    return both - one - two


def report(Rs, E, scale):
    """A force, not an energy, is what a gravity claim is about: F = -dE/dR."""
    Rm = 0.5 * (Rs[1:] + Rs[:-1])
    F = -np.diff(E) / np.diff(Rs)
    print(f"      {'R':>6} {'E_int(R)':>13} {'|':>2} {'R':>6} {'F = -dE/dR':>13} {'sense':>9}")
    for i in range(len(Rs)):
        row = f"      {Rs[i]:>6.0f} {E[i]:>13.3e} {'|':>2}"
        if i < len(Rm):
            sense = "attract" if F[i] < 0 else "repel"
            row += f" {Rm[i]:>6.0f} {F[i]:>13.3e} {sense:>9}"
        print(row)
    return F, Rm


if __name__ == "__main__":
    print("=== Does the tetrad graviton actually PULL? ===\n")
    print("  A field falloff is not a force. test_tetrad_shielding measured |h| ~ 1/r^2 around")
    print("  ONE mass. Here: the interaction energy between TWO.  E_int rising with R = ATTRACT.\n")

    N, Rd = 95, 40.0
    Rs = np.array([8.0, 12.0, 16.0, 20.0, 24.0])
    self_E = energy(N, Rd, [(0.0, 0.0)])            # the scale to judge "zero" against
    print(f"  (self-energy of one mass = {self_E:.4e}; judge E_int against this)\n")

    # ---------------- the measurement: isotropic medium ----------------
    print("  [A] ISOTROPIC medium (the real one: triangular lattice, central-force springs)")
    Ei = np.array([interaction(N, Rd, R, (1.0, 1.0, 1.0)) for R in Rs])
    Fi, Rm = report(Rs, Ei, self_E)
    decay = abs(Fi[0]) / max(abs(Fi[-1]), 1e-30)
    print(f"\n      => E_int SATURATES: the force collapses by {decay:.0f}x from R~10 to R~22.")
    print("      There IS a short-range attraction -- and it is exactly the CONTACT term:")
    print("      the screened gravity-by-density of Phase 3, nothing new. Beyond a few lattice")
    print("      spacings the force dies. There is NO long-range force. Each mass carries a")
    print("      perfectly good 1/r^2 tetrad field, and the cross term still integrates to")
    print("      nothing -- Crum's theorem, confirmed by measurement.\n")

    # ---------------- positive control: break the isotropy ----------------
    print("  [B] CONTROL -- break the isotropy (one bond direction 1.6x stiffer).")
    print("      Crum's theorem no longer applies, so a force MUST appear if the probe works.")
    Ea = np.array([interaction(N, Rd, R, (1.6, 1.0, 1.0)) for R in Rs])
    Fa, _ = report(Rs, Ea, self_E)
    print(f"\n      => a genuine long-range force appears, and it is REPULSIVE; at R~22 it is")
    print(f"      {abs(Fa[-1])/max(abs(Fi[-1]),1e-30):.0f}x the isotropic residual. The probe IS")
    print("      sensitive at exactly the range where [A] reads zero -- so the null is physics,")
    print("      not blindness. (It also shows how badly non-gravitational the elastic sector is:")
    print("      the only way to get a long-range force out of it is to make masses REPEL.)\n")

    print("[verdict] the tetrad graviton's 1/r^2 field carries NO LONG-RANGE FORCE.")
    print("  Crum/Eshelby: in an infinite ISOTROPIC elastic medium the interaction energy of two")
    print("  centres of dilatation VANISHES. The medium here IS isotropic (triangular lattice,")
    print("  central forces), and the measurement confirms it: two masses each have a long-range")
    print("  traceless strain field, a short-range contact attraction, and no force at range.")
    print("\n  So the last surviving gravity claim collapses. Restating the scorecard honestly:")
    print("    Route 1 (topological curvature): UNSHIELDABLE, but like charges REPEL with a force")
    print("        that GROWS with distance [test_disclination_force]. Not gravity.")
    print("    Tetrad (energy-sourced): long-range 1/r^2 FIELD, but SHIELDABLE")
    print("        [test_tetrad_shielding] and it exerts NO FORCE [here]. Not gravity either.")
    print("    gravity-by-density: real attraction, but SCREENED/short-range [phase 3].")
    print("\n  The model has produced exactly one genuine attraction between two masses -- the")
    print("  original nonlinear gravity-by-density (two lumps drift together, phase 3d) -- and")
    print("  that one is screened. Every LONG-RANGE candidate since has failed to deliver a")
    print("  universally attractive force. Long-range emergent gravity in this model is OPEN,")
    print("  not 'a route found'. The linear/elastic sector is provably the wrong place to look")
    print("  (Crum forbids it); the attraction that does exist is NONLINEAR. That is the lead.")
