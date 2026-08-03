"""
test_cone_unification -- the medium's non-central stiffness does triple duty:
it selects the fermion lattice, rigidifies it, AND sources its acoustic cone
=============================================================================

Context.  test_lattice_selection (Section 8.44) showed the bipartite honeycomb
-- the lattice that carries the emergent Dirac fermions -- is reachable only
with an O(1) angular (three-body) stiffness lam* ~ 0.57 of the bond energy; a
central pair force close-packs to triangular/fcc.  Independently,
test_lorentz_unified (Section 8.1) showed a purely central medium has c_L != c_T
(two acoustic cones, a Lorentz-violation between polarisations) and needs a
non-central stiffness to give one universal cone.  This asks whether those are
the SAME non-central stiffness by using the SAME force model (LJ central +
lam*(cos th + 1/2)^2 angular, exactly as in test_lattice_selection) and reading
the acoustic cone off its elastic constants as a function of the SAME lam.

Method: periodic elastic constants by homogeneous strain (minimum image, atoms
relaxed internally at fixed cell), pure numpy.  2D isotropic energy density
w = (lambda_2D/2)(tr eps)^2 + mu tr(eps^2), so a dilation gives the bulk modulus
K_2D and a shear gives mu; c_L = sqrt((K_2D+mu)/rho), c_T = sqrt(mu/rho).

Gates and results:
  G1  triangular, lam=0: c_L/c_T = sqrt(3) = 1.732 and nu = 1/3 (the central-force
      Cauchy solid -- matches test_cc_weinberg's 1.7330 and test_gamma_elastic's
      0.330), validating the machinery.
  G2  honeycomb, lam=0: mu = 0 (c_T = 0) -- floppy, coordination 3 < the 2D
      Maxwell threshold 4, so a central-force honeycomb has no transverse cone.
  G3  the honeycomb cone is elastically isotropic (shear the same in two
      orientations), so "cone" is well defined.
  sweep  honeycomb vs lam: the bulk modulus K is lam-INDEPENDENT (compression is
      central), while the shear modulus mu -- and hence the whole transverse cone
      c_T -- is sourced ENTIRELY by the angular term (mu ~ 0 at lam=0, rising
      ~linearly; c_T ~ sqrt(lam)).  At lam* ~ 0.57 (where honeycomb becomes the
      ground state) c_T is already O(1) (c_T/c_L ~ 0.33).  Full universality
      c_L=c_T is the mu >> K limit reached as lam grows past lam* (the vector-
      Hooke / equal-stretch-and-bend point of test_lorentz_unified).

Conclusion: one ingredient -- O(1) non-central (angular) stiffness -- has three
consequences absent from a purely central medium: it selects the honeycomb, it
makes it mechanically stable, and it is the sole source of its transverse cone.
Compression stays central.

Scope: the cone here is the MECHANICAL (medium+field) cone of test_lorentz_unified.
The fermion cone v_F = (3/2) t a is set by the electronic hopping t and is an
independent scale, so locking it to the boson cone (cross-statistics universality)
remains the separate composite-boson question (test_cone_universality /
test_cone_lock), which this angular stiffness does not by itself settle.
"""
import numpy as np

SQ3 = np.sqrt(3.0)
R1D, RCD, C0 = 1.30, 1.60, -0.5                   # angular window, cos(120 deg)


def _gwin(r):
    g = np.ones_like(r); gp = np.zeros_like(r)
    m = (r > R1D) & (r < RCD)
    t = (r[m] - R1D)/(RCD - R1D)
    g[m] = 0.5*(1 + np.cos(np.pi*t))
    gp[m] = -0.5*np.pi/(RCD - R1D)*np.sin(np.pi*t)
    g[r >= RCD] = 0.0
    return g, gp

def _phi(r):  return 4.0*(r**-12 - r**-6)
def _phip(r): return 4.0*(-12*r**-13 + 6*r**-7)

def cell(kind, n, a):
    """n x n supercell at NN spacing a; returns positions and cell matrix Hs
    (supercell lattice vectors as columns)."""
    if kind == 'triangular':
        a1, a2, basis = np.array([1, 0.]), np.array([0.5, SQ3/2]), [np.zeros(2)]
    elif kind == 'honeycomb':
        a1, a2 = np.array([1.5, SQ3/2]), np.array([1.5, -SQ3/2])
        basis = [np.zeros(2), np.array([1., 0.])]
    else:
        raise ValueError(kind)
    a1, a2 = a1*a, a2*a; basis = [b*a for b in basis]
    pts = [i*a1 + j*a2 + b for i in range(n) for j in range(n) for b in basis]
    return np.array(pts), np.column_stack([n*a1, n*a2])

def _minimg(D, Hs, Hinv):
    f = D @ Hinv.T
    f -= np.round(f)
    return f @ Hs.T

def energy_force(X, Hs, lam):
    """Periodic energy and analytic force (minimum image): LJ + lam angular."""
    N = len(X); Hinv = np.linalg.inv(Hs); F = np.zeros_like(X)
    D = _minimg(X[:, None, :] - X[None, :, :], Hs, Hinv)   # D[i,j] = r_i - r_j
    r = np.sqrt((D**2).sum(-1)); np.fill_diagonal(r, 1e9)
    g, gp = _gwin(r)
    iu = np.triu_indices(N, 1); rp = r[iu]
    E = (g[iu]*_phi(rp)).sum()
    dV = gp[iu]*_phi(rp) + g[iu]*_phip(rp)
    u = D[iu]/rp[:, None]; fp = -dV[:, None]*u
    np.add.at(F, iu[0], fp); np.add.at(F, iu[1], -fp)
    for i in range(N):
        nb = np.where(r[i] < RCD)[0]
        if len(nb) < 2: continue
        A = -D[i, nb]
        ra = np.sqrt((A**2).sum(1)); ga = g[i, nb]; gpa = gp[i, nb]
        s, t = np.triu_indices(len(nb), 1)
        As, At, ras, rat = A[s], A[t], ra[s], ra[t]
        ca = (As*At).sum(1)/(ras*rat); p = ca - C0
        E += lam*(ga[s]*ga[t]*p**2).sum()
        dca_s = At/(ras*rat)[:, None] - (ca/ras**2)[:, None]*As
        dca_t = As/(ras*rat)[:, None] - (ca/rat**2)[:, None]*At
        dh_s = ga[t][:, None]*(gpa[s][:, None]*(As/ras[:, None])*(p**2)[:, None]
                               + ga[s][:, None]*(2*p)[:, None]*dca_s)
        dh_t = ga[s][:, None]*(gpa[t][:, None]*(At/rat[:, None])*(p**2)[:, None]
                               + ga[t][:, None]*(2*p)[:, None]*dca_t)
        np.add.at(F, nb[s], -lam*dh_s); np.add.at(F, nb[t], -lam*dh_t)
        F[i] += lam*(dh_s + dh_t).sum(0)
    return E, F

def relax(X, Hs, lam, steps=600, dt=0.01):
    X = X.copy(); V = np.zeros_like(X)
    for _ in range(steps):
        V = 0.9*V + dt*energy_force(X, Hs, lam)[1]
        X += dt*V
    return X

def equilibrium(kind, lam, n=3):
    aa = np.linspace(1.02, 1.20, 19)
    E = [energy_force(*cell(kind, n, a), lam)[0] for a in aa]
    a = aa[int(np.argmin(E))]
    X, Hs = cell(kind, n, a)
    return relax(X, Hs, lam), Hs, a

def moduli(kind, lam, n=3, amp=0.01):
    X0, Hs0, a = equilibrium(kind, lam, n)
    A0 = abs(np.linalg.det(Hs0)); rho = len(X0)/A0
    E0 = energy_force(X0, Hs0, lam)[0]

    def w(strain):
        out = []
        for d in (amp, -amp):
            eps = strain*d
            Hs = (np.eye(2) + eps) @ Hs0
            X = relax(X0 @ (np.eye(2) + eps).T, Hs, lam)
            out.append((energy_force(X, Hs, lam)[0] - E0)/A0)
        return 0.5*sum(out)                        # symmetric second order

    K = w(np.eye(2))/(2*amp**2)                     # dilation eps = d*I
    mu = w(np.diag([1., -1.]))/(2*amp**2)           # pure shear
    mu2 = w(np.array([[0., 1.], [1., 0.]]))/(2*amp**2)   # shear at 45 deg
    cL = np.sqrt(max(K + mu, 0)/rho); cT = np.sqrt(max(mu, 0)/rho)
    return dict(a=a, K=K, mu=mu, mu2=mu2, cL=cL, cT=cT,
                ratio=(cL/cT if cT > 1e-6 else np.inf),
                nu=(K - mu)/(K + mu))


def main():
    print('=' * 70)
    print('CONE UNIFICATION: does angular stiffness source the fermion lattice')
    print('AND its acoustic cone from one O(1) knob?')
    print('=' * 70)

    print('\n[G1] triangular, lam=0  (expect c_L/c_T = sqrt3 = 1.732, nu = 1/3)')
    m = moduli('triangular', 0.0)
    print(f'     K={m["K"]:.3f}  mu={m["mu"]:.3f}  c_L/c_T={m["ratio"]:.3f}  '
          f'nu={m["nu"]:.3f}')
    assert abs(m['ratio'] - SQ3) < 0.03, 'triangular c_L/c_T must be sqrt3'
    assert abs(m['nu'] - 1/3) < 0.02, 'triangular nu must be 1/3'

    print('\n[G2] honeycomb, lam=0  (expect mu ~ 0 -> c_T = 0, floppy)')
    m0 = moduli('honeycomb', 0.0)
    print(f'     K={m0["K"]:.3f}  mu={m0["mu"]:.4f}  c_T={m0["cT"]:.4f}  '
          f'c_L={m0["cL"]:.3f}')
    assert abs(m0['mu']) < 0.15, 'central-force honeycomb must be floppy (mu~0)'

    print('\n[G3] honeycomb cone is elastically isotropic (shear two ways equal)')
    mi = moduli('honeycomb', 0.568)
    print(f'     mu(pure)={mi["mu"]:.3f}  mu(45deg)={mi["mu2"]:.3f}  '
          f'rel diff={abs(mi["mu"]-mi["mu2"])/mi["mu"]:.1e}')
    assert abs(mi['mu'] - mi['mu2'])/mi['mu'] < 0.05, 'cone must be isotropic'

    print('\n[sweep] honeycomb cone vs lam   (lam* selection = 0.568)')
    print('     lam      K       mu      c_L     c_T    c_T/c_L')
    Ks, rows = [], []
    for lam in [0.0, 0.1, 0.2, 0.3, 0.4, 0.568, 0.7, 1.0, 1.5]:
        m = moduli('honeycomb', lam)
        Ks.append(m['K']); rows.append((lam, m))
        rr = m['cT']/m['cL'] if m['cL'] > 1e-9 else 0.0
        print(f'     {lam:5.3f}  {m["K"]:6.3f}  {m["mu"]:6.3f}  '
              f'{m["cL"]:6.3f}  {m["cT"]:6.3f}   {rr:5.3f}')

    # K is central (lam-independent); mu (hence c_T) is entirely angular
    Kspread = (max(Ks) - min(Ks))/np.mean(Ks)
    mu_star = dict(rows)[0.568]['mu']
    print(f'\n     bulk modulus K spread over the whole sweep = {Kspread:.1e}'
          f'  (K is central, lam-independent)')
    print(f'     shear modulus at lam*=0.568: mu = {mu_star:.2f} '
          f'(0 at lam=0) -- the transverse cone is entirely angular-sourced')
    assert Kspread < 0.02, 'bulk modulus should be lam-independent'
    assert mu_star > 0.5, 'at lam* the honeycomb should have an O(1) shear cone'

    print('\n' + '=' * 70)
    print('RESULT: one O(1) non-central (angular) stiffness does triple duty --')
    print('it selects the honeycomb (lam* ~ 0.57), makes it mechanically stable,')
    print('and is the SOLE source of its transverse cone (c_T: 0 -> O(1); K stays')
    print('central and fixed). A purely central medium has none of the three.')
    print('Full universality c_L=c_T is the mu>>K limit past lam* (vector-Hooke).')
    print('=' * 70)


if __name__ == '__main__':
    main()
