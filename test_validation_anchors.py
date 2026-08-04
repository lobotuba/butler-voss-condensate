"""
Validation anchors: the topological machinery reproduces UN-TUNED textbook results.

The band-structure tooling used to measure the generation number (test_generations, S8.48)
raises a fair authenticity question: a Chern number computed for a coupling we CHOSE to wind
n times could be a construction reflecting itself back rather than a real measurement. The
guard against that is to run the SAME machinery on independent models it was NOT built for,
and check it reproduces their known analytic numbers.

  [A] HALDANE MODEL (Haldane, PRL 61, 2015 (1988)) -- the canonical Chern insulator, a
      DIFFERENT Hamiltonian from our winding construction. Its topological phase boundary is
      the analytic result |M| = 3*sqrt(3)*t2*sin(phi), Chern C = +-1 inside, 0 outside. The
      identical Fukui-Hatsugai-Suzuki link-variable flux used in test_generations reproduces
      that boundary (the number 3*sqrt(3) = 5.196...) and the +-1 lobes -- so the code
      computes real topology, not the winding we put in.

  [B] GRAPHENE nearest-neighbour tight binding -- textbook honeycomb numbers: total bandwidth
      6t (max|f| = 3), a gapless Dirac point at K, and Fermi velocity v_F = (3/2) t a.

Scope: this validates the CODE (internal correctness), not the physics claims' contact with
experiment. It shows the tooling is honest, nothing more.
"""
from __future__ import annotations
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)

# honeycomb geometry (a = 1): NN vectors delta; 2nd-neighbour vectors v = delta_i - delta_j
_d1 = np.array([0.5, np.sqrt(3) / 2]); _d2 = np.array([0.5, -np.sqrt(3) / 2]); _d3 = np.array([-1.0, 0.0])
DELTA = [_d1, _d2, _d3]
VNN = [_d2 - _d3, _d3 - _d1, _d1 - _d2]                          # |v| = sqrt(3), fixed chirality
BREC = 2 * np.pi * np.linalg.inv(np.array([_d1 - _d3, _d2 - _d3])).T   # reciprocal vectors (rows)


def haldane_H(k, M, t2, phi, t1=1.0):
    """Haldane Bloch Hamiltonian: real NN hopping t1, complex NNN hopping t2 e^{i phi},
    sublattice mass M."""
    hx = t1 * sum(np.cos(k @ d) for d in DELTA)
    hy = t1 * sum(np.sin(k @ d) for d in DELTA)
    hz = M - 2 * t2 * np.sin(phi) * sum(np.sin(k @ v) for v in VNN)
    h0 = 2 * t2 * np.cos(phi) * sum(np.cos(k @ v) for v in VNN)
    return h0 * I2 + hx * SX + hy * SY + hz * SZ


def chern_fhs_general(Hfunc, N=60):
    """The SAME Fukui-Hatsugai-Suzuki lattice flux as test_generations.chern_fhs, fed an
    arbitrary 2-band Bloch Hamiltonian; lower-band Chern, an exact integer on any grid."""
    u = np.linspace(0, 1, N, endpoint=False)
    U = np.empty((N, N, 2), complex)
    for i, a in enumerate(u):
        for j, b in enumerate(u):
            _, v = np.linalg.eigh(Hfunc(a * BREC[0] + b * BREC[1]))
            U[i, j] = v[:, 0]
    F = 0.0
    for i in range(N):
        for j in range(N):
            u00, u10 = U[i, j], U[(i + 1) % N, j]
            u11, u01 = U[(i + 1) % N, (j + 1) % N], U[i, (j + 1) % N]
            F += np.angle(np.vdot(u00, u10) * np.vdot(u10, u11)
                          / np.vdot(u01, u11) / np.vdot(u00, u01))
    return F / (2 * np.pi)


def graphene_f(k):
    return sum(np.exp(1j * (k @ d)) for d in DELTA)


def main():
    print('=' * 70)
    print('VALIDATION ANCHORS: the machinery reproduces un-tuned textbook results')
    print('=' * 70)

    print('\n[A] HALDANE MODEL -- topological boundary must sit at |M| = 3*sqrt(3)*t2*sin(phi):')
    print(f"    3*sqrt(3) = {3*np.sqrt(3):.4f}")
    print(f"    {'phi':>8} {'predicted |M|_c':>16} {'C below':>9} {'C above':>9}")
    for phi in (np.pi / 2, np.pi / 4, np.pi / 6):
        Mc = 3 * np.sqrt(3) * np.sin(phi)                     # analytic boundary, t2=1
        c_below = chern_fhs_general(lambda k, M=Mc - 0.3: haldane_H(k, M, 1.0, phi))
        c_above = chern_fhs_general(lambda k, M=Mc + 0.3: haldane_H(k, M, 1.0, phi))
        print(f"    {phi:>8.4f} {Mc:>16.3f} {c_below:>+9.2f} {c_above:>+9.2f}")
        assert abs(round(c_below)) == 1, 'inside the boundary the Haldane phase must be C=+-1'
        assert abs(round(c_above)) == 0, 'outside the boundary the Haldane phase must be trivial'

    print('\n[B] GRAPHENE nearest-neighbour tight binding (t = a = 1):')
    g = np.linspace(0, 1, 240, endpoint=False)
    absf = np.array([abs(graphene_f(a * BREC[0] + b * BREC[1])) for a in g for b in g])
    bw = absf.max()
    K = (2 * BREC[0] + BREC[1]) / 3
    fK = abs(graphene_f(K))
    q = 1e-4
    vF = np.mean([abs(graphene_f(K + q * np.array([np.cos(t), np.sin(t)]))) / q
                  for t in np.linspace(0, 2 * np.pi, 12, endpoint=False)])
    print(f"    bandwidth/2t = max|f| = {bw:.4f}   (textbook 3.000 -> total width 6t)")
    print(f"    |f(K)|               = {fK:.2e}   (textbook 0 -> gapless Dirac point)")
    print(f"    v_F = d|f|/dq at K   = {vF:.4f}   (textbook (3/2) t a = 1.5)")
    assert abs(bw - 3.0) < 1e-3, 'graphene bandwidth must be 6t (max|f| = 3)'
    assert fK < 1e-6, 'graphene must be gapless at K'
    assert abs(vF - 1.5) < 1e-3, 'graphene Fermi velocity must be (3/2) t a'

    print('\n' + '=' * 70)
    print('RESULT: the same Fukui-Hatsugai-Suzuki Chern code reproduces the Haldane phase')
    print('boundary |M| = 3*sqrt(3) t2 sin(phi) and the +-1 lobes, and the honeycomb code')
    print('reproduces graphene\'s bandwidth 6t and v_F = (3/2) t a -- both un-tuned. The')
    print('S8.48 topological integers are the code computing real topology, not a')
    print('construction reflecting itself back. (Validates the code, not contact with')
    print('experiment.)')
    print('=' * 70)


if __name__ == '__main__':
    main()
