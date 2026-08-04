"""
test_chiral_protection -- the Dirac cone survives phonon fluctuations to all
orders, by chiral symmetry (the dynamical closure of the two-cone seam)
=============================================================================

Section 8.46 showed a static phonon acts as a pseudo-gauge field and cannot gap
the Dirac cone.  This closes the dynamical version: a fluctuating phonon cannot
gap it either, and for a symmetry reason that holds to all orders.

A phonon is a bond-length modulation -- an OFF-DIAGONAL hopping term.  The
honeycomb Dirac Hamiltonian has the sublattice (chiral) symmetry S H S = -H with
S = diag(+1 on A, -1 on B).  Any bond term is chiral-ODD (S dH S = -dH, it
transforms like H), so it CANNOT generate the chiral-EVEN mass (S M S = +M) that
would open a gap -- at any order, static or dynamical.  The only gapping term is
an on-site sublattice mass, which a bond phonon does not produce.

Gates:
  G1  clean honeycomb: S H S = -H to machine precision (exact chiral symmetry).
  G2  a bond (phonon) perturbation is chiral-odd (S dH S + dH = 0); a staggered
      on-site mass is chiral-even (S M S - M = 0) -- so only the mass can gap it.
  G3  DOS: random BOND disorder (a frozen-phonon ensemble = the fluctuation
      content) keeps the spectrum gapless at E=0 up to O(1) amplitude, while a
      small staggered mass opens a hard gap.  The cone survives phonon
      fluctuations.

Honest scope -- the velocity refinement is NOT settled here.  Whether the fermion
cone v_F and the mechanical cone c_T actually MERGE (v_F = c_T, full Lorentz)
rather than merely coexisting is a two-velocity RG question.  That calculation is
regulator-limited: with two velocities there is no single Lorentz-covariant
cutoff, so a hard momentum/frequency cutoff breaks the symmetry it is meant to
measure -- the one-loop result fails its own Lorentz gate (gamma_v - gamma_c != 0
at v = c), exactly as test_graviton_ward's Ward identity failed under a hard
cutoff.  A symmetry-preserving scheme (dimensional regularisation, or a lattice/
BZ regulator as in test_lattice_ward) is needed to settle it.  What is robust
here is the gaplessness, because it rests on an exact symmetry, not a loop.
"""
import numpy as np
SQ3 = np.sqrt(3.0)


def honeycomb(L, bond_disorder=0.0, mass=0.0, seed=0):
    """Periodic L x L honeycomb tight-binding; returns H and the chiral op S."""
    r = np.random.default_rng(seed)
    a1, a2 = np.array([1.5, SQ3/2]), np.array([1.5, -SQ3/2])
    b = np.array([1.0, 0.0])
    pos, sub = [], []
    for i in range(L):
        for j in range(L):
            R = i*a1 + j*a2
            pos.append(R); sub.append(0)          # A
            pos.append(R + b); sub.append(1)      # B
    pos = np.array(pos); sub = np.array(sub); N = len(pos)
    Hs = np.column_stack([L*a1, L*a2]); Hinv = np.linalg.inv(Hs)
    D = pos[:, None, :] - pos[None, :, :]
    f = D @ Hinv.T; f -= np.round(f); Dm = f @ Hs.T
    dist = np.sqrt((Dm**2).sum(-1))
    H = np.zeros((N, N))
    ii, jj = np.where(np.triu((dist > 0.5) & (dist < 1.2)))
    for a, c in zip(ii, jj):
        t = 1.0 + bond_disorder*(r.random() - 0.5)
        H[a, c] = H[c, a] = -t
    H[np.arange(N), np.arange(N)] = mass*np.where(sub == 0, 1.0, -1.0)
    S = np.diag(np.where(sub == 0, 1.0, -1.0))
    return H, S


def min_gap(L, kind, W, nreal=8):
    ev = []
    for s in range(nreal):
        H, _ = honeycomb(L, bond_disorder=W if kind == 'bond' else 0.0,
                         mass=W if kind == 'mass' else 0.0, seed=s)
        ev.append(np.linalg.eigvalsh(H))
    ev = np.concatenate(ev)
    return float(np.min(np.abs(ev)))


def main():
    print('=' * 70)
    print('CHIRAL PROTECTION: the Dirac cone survives phonon fluctuations')
    print('=' * 70)
    L = 18

    # G1 -- exact chiral symmetry of the clean lattice
    H, S = honeycomb(L)
    err = np.abs(S @ H @ S + H).max()
    print(f'\n[G1] clean honeycomb: |S H S + H| = {err:.1e}  (chiral symmetry exact)')
    assert err < 1e-12, 'clean lattice must be chiral-symmetric'

    # G2 -- a bond perturbation is chiral-odd; a mass is chiral-even
    Hb, _ = honeycomb(L, bond_disorder=0.5, seed=1)
    dH = Hb - H                                   # pure bond (phonon) perturbation
    odd = np.abs(S @ dH @ S + dH).max()
    Hm, _ = honeycomb(L, mass=0.3)
    M = Hm - H                                    # pure staggered mass
    even = np.abs(S @ M @ S - M).max()
    print(f'[G2] bond (phonon) term chiral-ODD:  |S dH S + dH| = {odd:.1e}')
    print(f'     staggered mass  chiral-EVEN:    |S M S  -  M | = {even:.1e}')
    print('     => a bond term cannot make the (chiral-even) mass that gaps the cone.')
    assert odd < 1e-12 and even < 1e-12, 'chiral parities must be exact'

    # G3 -- DOS: bond disorder stays gapless; a mass gaps
    print(f'\n[G3] min|E| (half-gap), periodic L={L} ({2*L*L} sites):')
    g_clean = min_gap(L, 'bond', 0.0, nreal=1)
    print(f'     clean                    {g_clean:.4f}')
    for W in (0.3, 0.6, 1.0):
        g = min_gap(L, 'bond', W)
        print(f'     bond disorder W={W:<4}       {g:.4f}   (phonon-like, gapless)')
        assert g < 0.03, 'bond disorder must not gap the cone'
    for m in (0.1, 0.3):
        g = min_gap(L, 'mass', m)
        print(f'     staggered mass m={m:<4}      {g:.4f}   (gap = 2m)')
        assert abs(g - m) < 1e-6, 'a mass must open a gap of 2m'

    print('\n' + '=' * 70)
    print('RESULT: a phonon is chiral-odd, so it cannot generate the mass that')
    print('gaps the Dirac cone -- at any order, static or dynamical. Random bond')
    print('(phonon) disorder keeps the cone gapless to O(1); only a sublattice')
    print('mass (which a phonon cannot make) gaps it. The fermion cone survives')
    print('phonon fluctuations by exact chiral symmetry. (Whether v_F and c_T')
    print('MERGE is a separate, regulator-limited two-velocity RG question.)')
    print('=' * 70)


if __name__ == '__main__':
    main()
