"""
test_lattice_selection -- which lattice does the medium select, and can it be
the bipartite HONEYCOMB that hosts Dirac fermions?
=============================================================================

Context (the frontier this sits on).  The Standard-Model arc (test_sm_structure,
test_anomaly_hypercharge, test_electroweak) compressed the remaining SM input to
ONE discrete choice: *which lattice the medium realises*.  The gauge group, the
fermion reps and the generation count are all band-structure data of that
lattice (Dirac-point count -> multiplet; Chern index -> generations).  The
sharpest fork inside that choice: the medium's self-assembly (bvc_core.relax_medium)
is driven by a CENTRAL Lennard-Jones pair force, whose ground state is the
close-packed TRIANGULAR lattice (coordination 6; fcc in 3D).  A close-packed
lattice has NO Dirac cone -> NO emergent fermions (test_dirac needs a BIPARTITE
honeycomb, 2 sites/cell, coordination 3).  So: can the medium select honeycomb,
and what does that cost in the interaction?

What is measured (all T=0 lattice sums, cohesive energy per site minimised over
the NN spacing; ground state = lowest bound interior minimum; pure numpy):

  G0  Gate: plain LJ selects triangular; honeycomb is the HIGHEST-energy
      candidate (coordination gate + ranking).

  A   CENTRAL two-scale potentials (LJ + a repulsive Gaussian shoulder) can
      destabilise close packing to the SQUARE lattice, but NEVER to honeycomb
      (0 of the whole scan).  Reason: a close-packed lattice dodges a single
      isotropic barrier by DILATION (expand until the 2nd-neighbour shell slips
      off the barrier, paying only a little well depth).  No number of isotropic
      length scales forces coordination 3.

  B   ANGULAR route.  A Stillinger-Weber 3-body penalty lam*sum(cos th + 1/2)^2
      over each site's NN bonds (minimised at 120 deg) is ZERO for honeycomb
      (all bonds 120 deg) and large for triangular (60 deg bonds).  As lam rises
      the ground state goes triangular -> square -> HONEYCOMB.  The measured
      switch to the fermion lattice needs angular stiffness lam* ~ 0.57 of the
      NN bond energy -- an O(1) amount.

  C   DYNAMICS.  The same LJ + angular force (analytic, gated against finite
      differences to 1e-10) is annealed.  A heated honeycomb seed is preserved
      under the angular force (coordination 3, angles 120 deg) but densifies to
      coordination 6 under the central force alone -- honeycomb has z=3 < 4, the
      2D Maxwell rigidity threshold, so it is mechanically floppy without angular
      bonds.  A disordered droplet close-packs (z=6) under central forces and
      forms honeycomb (z=3, angles 120 deg) under the angular force.
      NB the bond-orientational psi6 is ~1 for BOTH honeycomb and triangular
      (120-deg bonds give exp(6i.theta)=1), so coordination -- not psi6 -- is the
      order parameter that separates them.

Conclusion: the fermion-hosting lattice is not reachable from a central pair
potential; it requires order-unity DIRECTIONAL (3-body) bonding -- the same
non-central/shear stiffness the emergent-Lorentz work already required
(test_lorentz_unified: bare central forces give c_L>c_T always).  This reframes
the last SM input "which lattice" as one physical property of the medium: does
its bonding carry O(1) angular rigidity?
"""
import numpy as np

SQ3 = np.sqrt(3.0)

# ================================================================ lattices ===
def _cluster(kind, nmax=16):
    if kind == 'triangular':
        a1, a2, basis = np.array([1, 0.]), np.array([0.5, SQ3/2]), [np.zeros(2)]
    elif kind == 'square':
        a1, a2, basis = np.array([1, 0.]), np.array([0, 1.]), [np.zeros(2)]
    elif kind == 'honeycomb':
        a1, a2 = np.array([1.5, SQ3/2]), np.array([1.5, -SQ3/2])
        basis = [np.zeros(2), np.array([1., 0.])]
    elif kind == 'kagome':
        a1, a2 = np.array([2., 0.]), np.array([1., SQ3])
        basis = [np.zeros(2), np.array([1., 0.]), np.array([0.5, SQ3/2])]
    elif kind == 'chain':                         # z=2 control (180 deg bonds)
        a1, a2, basis = np.array([1, 0.]), np.array([0., 40.]), [np.zeros(2)]
    else:
        raise ValueError(kind)
    pts = [n1*a1 + n2*a2 + b
           for n1 in range(-nmax, nmax + 1) for n2 in range(-nmax, nmax + 1)
           for b in basis]
    return np.array(pts)

def _ref(P):
    return int(np.argmin(((P - P.mean(0))**2).sum(1)))

def neighbour_dists(kind, smax=6.0):
    P = _cluster(kind); i = _ref(P)
    d = np.sqrt(((P - P[i])**2).sum(1))
    return np.sort(d[(d > 1e-9) & (d <= smax)])

def nn_bonds(kind, tol=1.05):
    P = _cluster(kind); i = _ref(P)
    v = P - P[i]; d = np.sqrt((v**2).sum(1))
    return v[(d > 1e-9) & (d < tol)]              # first shell at unit NN

def angular_penalty(kind):
    """sum over unordered NN-bond pairs of (cos th + 1/2)^2 ; 0 at 120 deg."""
    b = nn_bonds(kind)
    if len(b) < 2:
        return 0.0
    u = b / np.linalg.norm(b, axis=1, keepdims=True)
    C = u @ u.T
    iu = np.triu_indices(len(u), k=1)
    return float(((C[iu] + 0.5)**2).sum())

LATT = ['triangular', 'square', 'honeycomb', 'kagome', 'chain']
Z = {'triangular': 6, 'square': 4, 'honeycomb': 3, 'kagome': 4, 'chain': 2}
DISTS = {k: neighbour_dists(k) for k in LATT}
PANG = {k: angular_penalty(k) for k in LATT}

# ============================================================== potentials ===
def Vlj(r):
    return 4.0*(r**-12 - r**-6)

def Vcs(r, A, rs, w):
    return Vlj(r) + A*np.exp(-((r - rs)**2)/(2*w*w))

def min_energy(dists, Vfunc, Rc=4.5, dgrid=None):
    """(E_min, d_min, is_real): a real crystal has a BOUND (E<0) INTERIOR
    minimum (not pinned at the grid edge = wants to dilate/unbind)."""
    if dgrid is None:
        dgrid = np.linspace(0.85, 2.60, 501)
    E = np.array([0.5*Vfunc((dists*d)[dists*d <= Rc]).sum() for d in dgrid])
    i = int(np.argmin(E))
    return E[i], dgrid[i], (0 < i < len(dgrid) - 1 and E[i] < 0.0)

def ground_state(Vfunc, cands=LATT):
    res = {k: min_energy(DISTS[k], Vfunc) for k in cands}
    order = sorted([k for k in cands if res[k][2]], key=lambda k: res[k][0])
    return order, res


# ================================ dynamics: 3-body self-assembly (part C) =====
# LJ pair + Stillinger-Weber 120-degree angular term, both smoothly cut by a C1
# cosine window so the analytic force can be gated against finite differences.
R1D, RCD = 1.30, 1.60                 # window: g=1 below R1D, 0 above RCD
C0 = -0.5                             # cos(120 deg)

def _gwin(r):
    g = np.ones_like(r); gp = np.zeros_like(r)
    m = (r > R1D) & (r < RCD)
    t = (r[m] - R1D)/(RCD - R1D)
    g[m] = 0.5*(1 + np.cos(np.pi*t))
    gp[m] = -0.5*np.pi/(RCD - R1D)*np.sin(np.pi*t)
    g[r >= RCD] = 0.0
    return g, gp

def _phi(r):
    return 4.0*(r**-12 - r**-6)
def _phip(r):
    return 4.0*(-12*r**-13 + 6*r**-7)

def energy_force(X, lam):
    """Total energy and force for LJ + lam*(cos th + 1/2)^2 angular bonds."""
    N = len(X); F = np.zeros_like(X); E = 0.0
    d = X[:, None, :] - X[None, :, :]
    r = np.sqrt((d**2).sum(-1)) + np.eye(N)*1e9
    g, gp = _gwin(r)
    iu = np.triu_indices(N, 1)
    rp = r[iu]
    E += (g[iu]*_phi(rp)).sum()
    dV = gp[iu]*_phi(rp) + g[iu]*_phip(rp)
    u = d[iu]/rp[:, None]
    fp = -dV[:, None]*u
    np.add.at(F, iu[0], fp)
    np.add.at(F, iu[1], -fp)
    for i in range(N):                            # 3-body, vertex at i
        nb = np.where(r[i] < RCD)[0]
        if len(nb) < 2:
            continue
        A = X[nb] - X[i]
        ra = np.sqrt((A**2).sum(1))
        ga, gpa = g[i, nb], gp[i, nb]
        s, t = np.triu_indices(len(nb), 1)
        As, At, ras, rat = A[s], A[t], ra[s], ra[t]
        ca = (As*At).sum(1)/(ras*rat)
        p = ca - C0
        E += lam*(ga[s]*ga[t]*p**2).sum()
        dca_s = At/(ras*rat)[:, None] - (ca/ras**2)[:, None]*As
        dca_t = As/(ras*rat)[:, None] - (ca/rat**2)[:, None]*At
        dh_s = ga[t][:, None]*(gpa[s][:, None]*(As/ras[:, None])*(p**2)[:, None]
                               + ga[s][:, None]*(2*p)[:, None]*dca_s)
        dh_t = ga[s][:, None]*(gpa[t][:, None]*(At/rat[:, None])*(p**2)[:, None]
                               + ga[t][:, None]*(2*p)[:, None]*dca_t)
        np.add.at(F, nb[s], -lam*dh_s)
        np.add.at(F, nb[t], -lam*dh_t)
        F[i] += lam*(dh_s + dh_t).sum(0)
    return E, F

def _fd_force(X, lam, eps=1e-6):
    Fn = np.zeros_like(X)
    for i in range(len(X)):
        for c in range(2):
            Xp = X.copy(); Xp[i, c] += eps
            Xm = X.copy(); Xm[i, c] -= eps
            Fn[i, c] = -(energy_force(Xp, lam)[0] - energy_force(Xm, lam)[0])/(2*eps)
    return Fn

def honeycomb_patch(rings=3, a=1.12):
    a1, a2 = np.array([1.5, SQ3/2])*a, np.array([1.5, -SQ3/2])*a
    b = np.array([a, 0.])
    pts = [n1*a1 + n2*a2 + bb
           for n1 in range(-rings, rings+1) for n2 in range(-rings, rings+1)
           for bb in (np.zeros(2), b)]
    X = np.array(pts); X = X - X.mean(0)
    return X[np.linalg.norm(X, axis=1) <= rings*1.5*a]

NNCUT = 1.45                          # NN shell cut for coordination (bond ~1.12)

def order_params(X):
    d = X[:, None, :] - X[None, :, :]
    r = np.sqrt((d**2).sum(-1)) + np.eye(len(X))*1e9
    coord = (r < NNCUT).sum(1)
    rad = np.linalg.norm(X - X.mean(0), axis=1)
    interior = rad < 0.62*rad.max()
    th = np.arctan2(d[:, :, 1], d[:, :, 0])
    psi6, ang120 = [], []
    for i in np.where(interior)[0]:
        nb = np.where(r[i] < NNCUT)[0]
        if len(nb) < 2:
            continue
        psi6.append(np.mean(np.exp(6j*th[i, nb])))
        u = (X[nb] - X[i]); u /= np.linalg.norm(u, axis=1, keepdims=True)
        C = np.clip(u @ u.T, -1, 1)
        iu = np.triu_indices(len(nb), 1)
        ang120.append(np.mean(np.abs(np.degrees(np.arccos(C[iu])) - 120) < 20))
    return dict(zmean=float(coord[interior].mean()),
                z3=float(np.mean(coord[interior] == 3)),
                z6=float(np.mean(coord[interior] == 6)),
                psi6=float(np.abs(np.mean(psi6))) if psi6 else 0.0,
                ang120=float(np.mean(ang120)) if ang120 else 0.0)

def anneal(X, lam, steps=3500, dt=0.0025, T0=0.15, fcap=30.0, seed=0):
    rng = np.random.default_rng(seed)
    X = X.copy()
    for k in range(steps):
        _, F = energy_force(X, lam)
        n = np.linalg.norm(F, axis=1, keepdims=True)
        F = np.where(n > fcap, F*fcap/(n + 1e-12), F)
        T = T0*(1 - k/steps)
        X += dt*F + np.sqrt(2*T*dt)*rng.normal(size=X.shape)
    return X


def main():
    print('=' * 70)
    print('LATTICE SELECTION: can the medium pick the fermion (honeycomb) lattice?')
    print('=' * 70)

    # --- coordination gate ---
    for k in LATT:
        z = int((np.abs(DISTS[k] - 1.0) < 1e-6).sum())
        assert z == Z[k], f'{k}: coordination {z} != {Z[k]}'
    print('[gate] coordination:', {k: Z[k] for k in LATT}, 'OK')

    # ---------------------------------------------------------------- G0 ----
    order, res = ground_state(Vlj, cands=['triangular', 'square',
                                          'honeycomb', 'kagome'])
    print('\n[G0] plain LJ (single scale) ground-state ranking:')
    for k in order:
        print(f'     {k:11s} E={res[k][0]:+8.4f}  d={res[k][1]:.3f}')
    assert order[0] == 'triangular', 'LJ must select triangular'
    assert order[-1] == 'honeycomb', 'honeycomb must be highest-energy under LJ'
    print('     => triangular (close-packed); honeycomb is HIGHEST. Gate OK.')

    # ----------------------------------------------------- PART A: central --
    print('\n[A] central two-scale scan  V = LJ + A*exp(-(r-rs)^2/2w^2)')
    As = np.linspace(0.3, 3.0, 22)
    rss = np.linspace(1.7, 2.7, 16)
    ws = [0.30, 0.40, 0.50]
    tally, npts = {}, 0
    for w in ws:
        for A in As:
            for rs in rss:
                order, _ = ground_state(lambda r: Vcs(r, A, rs, w),
                                        cands=['triangular', 'square',
                                               'honeycomb', 'kagome'])
                if not order:
                    continue
                npts += 1
                tally[order[0]] = tally.get(order[0], 0) + 1
    print('     ground-state tally:', tally, f'  ({npts} bound points)')
    print(f'     honeycomb wins: {tally.get("honeycomb", 0)} / {npts}')
    assert tally.get('honeycomb', 0) == 0, \
        'no central 2-scale potential should select honeycomb'
    assert tally.get('square', 0) > 0, 'a shoulder should reach square'
    print('     => central potentials reach SQUARE but NEVER honeycomb.')
    print('        (close packing dodges an isotropic barrier by dilation.)')

    # ----------------------------------------------------- PART B: angular --
    print('\n[B] angular 3-body term  E = E_pair(LJ) + lam * (cos th + 1/2)^2')
    Epair = {k: min_energy(DISTS[k], Vlj)[0] for k in LATT}
    print('     lattice      z   E_pair    angular P(120deg)')
    for k in LATT:
        print(f'       {k:10s} {Z[k]:>2d}  {Epair[k]:+7.3f}   {PANG[k]:7.4f}')
    assert PANG['honeycomb'] < 1e-9, 'honeycomb must have zero 120deg penalty'

    def total(k, lam):
        return Epair[k] + lam*PANG[k]

    def ground_ang(lam):
        return min(LATT, key=lambda k: total(k, lam))

    # the two crossovers along increasing lam
    lam_ts = (Epair['square'] - Epair['triangular']) / \
             (PANG['triangular'] - PANG['square'])          # tri -> squ
    lam_sh = (Epair['honeycomb'] - Epair['square']) / \
             (PANG['square'] - PANG['honeycomb'])           # squ -> hex
    print(f'\n     crossovers:  triangular->square  lam={lam_ts:.3f}'
          f'   square->honeycomb  lam={lam_sh:.3f}')
    print('     lam     ground state')
    for lam in [0.0, 0.1, 0.3, lam_sh*0.98, lam_sh*1.02, 0.8, 1.2]:
        print(f'       {lam:5.3f}   {ground_ang(lam)}')

    assert ground_ang(0.0) == 'triangular'
    assert ground_ang(lam_sh*1.02) == 'honeycomb'
    assert 0.3 < lam_sh < 1.0, 'honeycomb crossover should be O(1) bond energy'
    print(f'\n     => honeycomb (the FERMION lattice) becomes the ground state')
    print(f'        above angular stiffness lam* = {lam_sh:.3f} x NN bond energy.')

    # ------------------------------------ PART C: dynamical self-assembly ----
    print('\n[C] self-assembly under LJ + lam*(cos th + 1/2)^2 angular force')
    Xh = honeycomb_patch(3)
    F = energy_force(Xh + 1e-3, 1.0)[1]           # nonzero-config force gate
    Xg = Xh + np.random.default_rng(4).normal(size=Xh.shape)*0.03
    err = np.abs(energy_force(Xg, 0.8)[1] - _fd_force(Xg, 0.8)).max()
    fmag = np.abs(energy_force(Xg, 0.8)[1]).max()
    print(f'     force gate: analytic vs finite-diff rel err = {err/fmag:.1e}')
    assert err/fmag < 1e-6, 'force gate FAILED'

    # basin: heat a honeycomb seed, relax under central-only vs angular
    seed = order_params(Xh)
    print(f'     seed honeycomb (N={len(Xh)}): z={seed["zmean"]:.2f} '
          f'z3={seed["z3"]:.2f} ang120={seed["ang120"]:.2f}')
    Xheat = Xh + np.random.default_rng(2).normal(size=Xh.shape)*0.12
    b0 = order_params(anneal(Xheat, 0.0, seed=3))
    b1 = order_params(anneal(Xheat, 1.0, seed=3))
    print(f'     basin  lam=0.0 : z={b0["zmean"]:.2f} z3={b0["z3"]:.2f} '
          f"z6={b0['z6']:.2f} ang120={b0['ang120']:.2f}  (collapses)")
    print(f'     basin  lam=1.0 : z={b1["zmean"]:.2f} z3={b1["z3"]:.2f} '
          f"z6={b1['z6']:.2f} ang120={b1['ang120']:.2f}  (honeycomb held)")
    assert b0['z6'] > 0.4 and b0['z3'] < 0.1, 'central should densify to z=6'
    assert b1['z3'] > 0.9 and b1['ang120'] > 0.9, 'angular should hold honeycomb'

    # nucleation: disordered droplet -> anneal
    rng = np.random.default_rng(7)
    N = len(Xh); Rd = np.sqrt(N*1.9/np.pi); P = []
    while len(P) < N:
        p = (rng.random(2)*2 - 1)*Rd
        if np.hypot(*p) <= Rd and all(np.hypot(*(p - q)) > 1.0 for q in P):
            P.append(p)
    Xdis = np.array(P)
    n0 = order_params(anneal(Xdis, 0.0, steps=4500, T0=0.25, seed=5))
    n1 = order_params(anneal(Xdis, 1.0, steps=4500, T0=0.25, seed=5))
    print(f'     nucl   lam=0.0 : z={n0["zmean"]:.2f} z6={n0["z6"]:.2f}  '
          f'(triangular)')
    print(f'     nucl   lam=1.0 : z={n1["zmean"]:.2f} z3={n1["z3"]:.2f} '
          f"ang120={n1['ang120']:.2f}  (honeycomb domains)")
    assert n0['z6'] > 0.4, 'central droplet should close-pack'
    assert n1['z3'] > 0.4 and n1['ang120'] > 0.8, 'angular droplet -> honeycomb'
    print('     => honeycomb is only DYNAMICALLY selected with angular forces;')
    print('        under central forces it is floppy (z=3<4) and densifies to z=6.')

    print('\n' + '=' * 70)
    print('RESULT: a central pair potential cannot select the fermion lattice')
    print('(triangular/square only). Honeycomb needs O(1) DIRECTIONAL bonding')
    print(f'(lam* ~ {lam_sh:.2f} bond-energy) -- the same non-central stiffness')
    print('the emergent-Lorentz cone already required. The last SM input')
    print('"which lattice" = "does the medium have O(1) angular rigidity?"')
    print('=' * 70)


if __name__ == '__main__':
    main()
