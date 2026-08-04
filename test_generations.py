"""
The generation number as a MEASURED topological output, not an assertion.

test_sm_structure (S8.43) reframed the family number as a Chern index but only ever
computed |C| = 1 (a Wilson band) and only ASSERTED the bulk-boundary count. This file
makes the statement a measurement, on the medium's own kind of band structure.

The emergent low-energy fermion is two-band (two sublattices), H(k) = d(k).sigma. Its
generation content is the Chern number of the occupied band -- and that Chern number is
exactly the number of times the inter-sublattice coupling winds around the Brillouin
zone. We take a coupling that winds n times, (sin kx + i sin ky)^n, gapped by a Wilson
mass, and MEASURE the family number two independent ways that must agree:

  [G1] BULK. The Fukui-Hatsugai-Suzuki lattice flux (a gauge-invariant integer on any
       grid) gives the Chern number C = -n for winding n = 1, 2, 3. The generation
       number is dialled by the winding of one structural function, and it comes out an
       exact integer -- quantized by topology, not tuned.

  [G2] BOUNDARY (measured, not asserted). On a cylinder (periodic in x, open in y) we
       count the in-gap edge branches crossing E=0. The left edge carries exactly n
       chiral families -- the bulk-boundary correspondence of S8.43, now measured. The
       generations are the chiral modes bound to a boundary/wall.

  [G3] NIELSEN-NINOMIYA. The right edge carries the same n branches with the OPPOSITE
       chirality, so the closed lattice nets to zero. n net-chiral generations therefore
       cannot be a bulk 2D property; they require a defect (edge/domain wall), and the
       opposite chirality lives elsewhere as a compensating (mirror) sector.

  [G4] ROBUSTNESS -- and the unification of S8.43's two mechanisms. A winding-n point is a
       MULTI-WEYL point: the off-diagonal coupling w = dx + i dy has a degree-n zero,
       protected only by crystalline symmetry and therefore fine-tuned. A generic
       perturbation (here a uniform A-B coupling, w -> w + delta^n) splits it into exactly
       n unit-winding Dirac points, and the TOTAL winding is conserved = n. So the robust
       content of "n generations" is n ordinary (unit) Dirac points -- a count of fermion
       doublers -- not a fragile n-fold degeneracy. This identifies S8.43's two separate
       statements as one integer: the flavour count [A] (number of Dirac points) and the
       generation index [B] (the Chern number, = the summed chirality/winding of those
       points) are the same conserved topological charge.

Honest scope. This does NOT derive the number three. It shows the family number is a
winding integer that can take any value; nothing in the medium's currently-known
structure pins it to three. "Why three" becomes "the inter-sublattice coupling winds
three times" -- a geometric restatement of the same input, and the count is now a
measured band-structure invariant rather than a free continuous parameter. The model is
a two-band Chern caricature, not the Standard Model's actual fermion content.
"""
from __future__ import annotations
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def hk(kx, ky, n, m, b=1.0):
    """2-band Bloch H = d.sigma; the inter-sublattice coupling winds n times:
    (sin kx + i sin ky)^n, gapped by the Wilson mass dz = m + b(cos kx + cos ky)."""
    w = (np.sin(kx) + 1j * np.sin(ky)) ** n
    dz = m + b * (np.cos(kx) + np.cos(ky))
    return w.real * SX + w.imag * SY + dz * SZ


# ------------------------------------------------------------- [G1] bulk Chern number ------------------
def chern_fhs(n, m, b=1.0, N=48):
    """Fukui-Hatsugai-Suzuki lattice Chern of the lower band -- an exact integer on any grid."""
    g = np.linspace(0, 2 * np.pi, N, endpoint=False)
    U = np.empty((N, N, 2), complex)
    for i, kx in enumerate(g):
        for j, ky in enumerate(g):
            e, v = np.linalg.eigh(hk(kx, ky, n, m, b))
            U[i, j] = v[:, 0]
    F = 0.0
    for i in range(N):
        for j in range(N):
            u00, u10 = U[i, j], U[(i + 1) % N, j]
            u11, u01 = U[(i + 1) % N, (j + 1) % N], U[i, (j + 1) % N]
            F += np.angle(np.vdot(u00, u10) * np.vdot(u10, u11)
                          / np.vdot(u01, u11) / np.vdot(u00, u01))
    return F / (2 * np.pi)


# ------------------------------------------------------------- [G2],[G3] edge modes on a cylinder ------
def _yhops(kx, n, m, b, kys):
    """exact real-space y-hoppings t[dy] = (1/Mky) sum_ky H(kx,ky) e^{i ky dy} (range |dy|<=n)."""
    return {dy: sum(hk(kx, ky, n, m, b) * np.exp(1j * ky * dy) for ky in kys) / len(kys)
            for dy in range(-n, n + 1)}


def edge_branches(n, m, b=1.0, W=60, Nkx=361, edge_frac=0.12, ewin=0.06):
    """Cylinder (periodic x, open y with W sites). For each edge, count the kx-clusters
    where an edge-localised in-gap state sits at E~0 (= chiral branches) and the net
    chirality (sum of branch slope signs). Returns (nL, chiL, nR, chiR)."""
    kys = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    kxs = np.linspace(0, 2 * np.pi, Nkx)
    ne = max(1, int(edge_frac * W))
    hitsL, hitsR = [], []   # (kx, E) of edge-localised in-gap states near E=0
    for kx in kxs:
        t = _yhops(kx, n, m, b, kys)
        Hr = np.zeros((2 * W, 2 * W), complex)
        for y in range(W):
            for dy in range(-n, n + 1):
                yp = y + dy
                if 0 <= yp < W:
                    Hr[2 * y:2 * y + 2, 2 * yp:2 * yp + 2] += t[dy]
        Hr = 0.5 * (Hr + Hr.conj().T)
        e, v = np.linalg.eigh(Hr)
        wg = (np.abs(v) ** 2).reshape(W, 2, -1).sum(axis=1)
        wL, wR = wg[:ne].sum(axis=0), wg[-ne:].sum(axis=0)
        for E in e[(np.abs(e) < ewin) & (wL > 0.6) & (wR < 0.2)]:
            hitsL.append((kx, E))
        for E in e[(np.abs(e) < ewin) & (wR > 0.6) & (wL < 0.2)]:
            hitsR.append((kx, E))
    return (*_clusters(hitsL), *_clusters(hitsR))


def _clusters(hits, gap=0.3):
    """Cluster the (kx, E~0) crossings ON THE CIRCLE (kx=0 == kx=2pi). Return (count,
    net chirality) where chirality = sign(dE/dkx) of each branch at its crossing."""
    if not hits:
        return 0, 0
    hits = sorted(hits)
    kx = np.array([h[0] for h in hits]); E = np.array([h[1] for h in hits])
    ncnt = len(kx)
    circ = (np.roll(kx, -1) - kx) % (2 * np.pi)     # gap to the next point around the circle
    s = int(np.argmax(circ))                        # cut right after the widest gap -> no wrap
    order = np.roll(np.arange(ncnt), -(s + 1))
    ku = kx[order].astype(float); es = E[order]
    for i in range(1, ncnt):                         # unwrap kx along the ordered sequence
        while ku[i] < ku[i - 1]:
            ku[i] += 2 * np.pi
    splits = [0] + [i + 1 for i in np.where(np.diff(ku) > gap)[0]] + [ncnt]
    count, chi = len(splits) - 1, 0
    for a, b in zip(splits[:-1], splits[1:]):
        if b - a >= 2 and np.ptp(ku[a:b]) > 1e-9:
            chi += int(np.sign(np.polyfit(ku[a:b], es[a:b], 1)[0]))
    return count, chi


# ------------------------------------------------------------- [G4] multi-Weyl fragmentation ----------
def _s(kx, ky):
    return np.sin(kx) + 1j * np.sin(ky)


def w_winding(wf, cx, cy, R, N=2000):
    """(1/2pi) * closed integral of d(arg w) around a circle radius R centred at (cx,cy)."""
    th = np.linspace(0, 2 * np.pi, N, endpoint=False)
    wv = wf(cx + R * np.cos(th), cy + R * np.sin(th))
    ang = np.unwrap(np.angle(wv))
    return (ang[-1] - ang[0] + np.angle(wv[0] / wv[-1])) / (2 * np.pi)


def find_zeros(wf, box=0.4, grid=500, tol=2e-3, merge=0.1):
    """locate zeros of complex w in [-box,box]^2 by sub-tol grid points, de-duplicated."""
    ax = np.linspace(-box, box, grid)
    KX, KY = np.meshgrid(ax, ax, indexing='ij')
    lo = np.abs(wf(KX, KY)) < tol
    zeros = []
    for i, j in zip(*np.where(lo)):
        p = (KX[i, j], KY[i, j])
        if not any((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 < merge ** 2 for q in zeros):
            zeros.append(p)
    return zeros


def main():
    print('=' * 70)
    print('GENERATIONS: the family number is a measured band-structure invariant')
    print('=' * 70)
    M = -0.5   # in-window Wilson mass: only the origin winding-n vortex is inverted

    print('\n[G1] BULK Chern number dials to the winding n of the sublattice coupling:')
    print(f"     {'winding n':>10} {'FHS Chern C':>13}")
    for n in (1, 2, 3):
        C = chern_fhs(n, M)
        print(f'     {n:>10} {C:>13.3f}')
        assert abs(C - (-n)) < 1e-6, 'Chern number must equal -n (exact integer)'

    print('\n[G2] BOUNDARY measured -- in-gap chiral branches on a cylinder edge = n:')
    print('[G3] NIELSEN-NINOMIYA -- the two edges carry opposite chirality (net zero):')
    print(f"     {'n':>3} {'|C| bulk':>9} {'left #':>7} {'right #':>8} {'chi_L':>6} {'chi_R':>6}")
    for n in (1, 2, 3):
        nL, chiL, nR, chiR = edge_branches(n, M)
        print(f'     {n:>3} {n:>9} {nL:>7} {nR:>8} {chiL:>+6} {chiR:>+6}')
        assert nL == n and nR == n, 'each edge must carry exactly n chiral families'
        assert chiL == -chiR and abs(chiL) == n, 'edges must have opposite chirality (NN)'

    print('\n[G4] ROBUSTNESS -- a winding-n (multi-Weyl) point is fine-tuned; a generic')
    print('     perturbation splits it into n unit Dirac points, total winding conserved:')
    print(f"     {'n':>3} {'charge':>7} {'# points after split':>21} {'sum winding':>12}")
    for n in (1, 2, 3):
        charge = w_winding(lambda kx, ky: _s(kx, ky) ** n, 0, 0, 0.3)
        assert abs(charge - n) < 1e-6, 'multi-Weyl charge must equal the winding n'
        d = 0.15
        wf = (lambda kx, ky, n=n, d=d: _s(kx, ky) ** n + d ** n)   # generic unfolding
        zs = find_zeros(wf)
        tot = sum(w_winding(wf, cx, cy, 0.04) for (cx, cy) in zs)
        print(f'     {n:>3} {charge:>+7.2f} {len(zs):>21} {tot:>+12.2f}')
        assert len(zs) == n, 'must split into exactly n Dirac points'
        assert abs(tot - n) < 1e-6, 'each split point is unit-winding; total is conserved = n'
    print('     => the robust content of "n generations" is n unit Dirac points (fermion')
    print('        doubling); S8.43\'s Dirac-point count [A] and Chern index [B] are one integer.')

    print('\n' + '=' * 70)
    print('RESULT: the generation number is the Chern number of the emergent two-band')
    print('fermion -- the winding of the medium\'s inter-sublattice coupling around the')
    print('Brillouin zone. It is a measured integer (bulk lattice flux AND the counted')
    print('edge branches agree, C = -n for n = 1,2,3), so it is quantized by topology,')
    print('not tuned. By Nielsen-Ninomiya the closed lattice nets to zero, so the n')
    print('net-chiral generations live on a defect (edge/domain wall) with a compensating')
    print('opposite-chirality mirror sector. The winding-n point is a fine-tuned multi-Weyl')
    print('degeneracy that fragments into n unit Dirac points, so the robust family number is')
    print('a count of ordinary Dirac points -- unifying the flavour count and the Chern index')
    print('of S8.43 as one conserved integer. HONEST scope: this does not derive three -- the')
    print('winding is a free integer; "why three" is now "the coupling winds three times", a')
    print('measured geometric restatement of the same input.')
    print('=' * 70)


if __name__ == '__main__':
    main()
