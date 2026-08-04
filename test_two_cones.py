"""
test_two_cones -- the medium has two cone sectors, and only the fermionic one is
the emergent Lorentz cone
=============================================================================

test_cone_unification (Section 8.45) showed an O(1) angular stiffness gives the
honeycomb a mechanical (acoustic) cone c_T and selects the lattice.  But the
matter that lives ON the honeycomb are the Dirac fermions of Section 8.2, whose
cone v_F = (3/2) t a is set by the electronic hopping t -- a different scale from
the mechanical stiffness.  This asks the provenance question directly: are the
mechanical phonon cone and the fermion cone the same cone, and which one is the
Lorentz cone of the emergent world?

Measurements (pure numpy):

  A  the two cones are INDEPENDENT.  c_T is fixed by the mechanical model (LJ +
     angular, Section 8.45); v_F = (3/2) t a needs the hopping t, which is
     nowhere in the mechanical model.  So v_F/c_T is proportional to t and takes
     any value -- equality is a tuning (t = 2 c_T / 3a), not a symmetry.  This is
     the cross-statistics two-cone problem of test_cone_universality, now at the
     mechanical-phonon-vs-fermion level.

  B  the PHYSICAL emergent boson rides v_F, not c_T.  The lower edge of the
     interband particle-hole continuum, omega_min(q) = min_k t(|f(k+q)|+|f(k)|),
     goes to v_F |q| as q -> 0 (test_cone_lock's signature): a massless collective
     mode whose speed is the fermion v_F, independent of the mechanical c_T.  So
     the composite boson (the emergent photon/graviton of Sections 8.5, 8.12) is
     on the fermion cone; the mechanical phonon is a separate sector.

  C  the mismatch is BENIGN.  A phonon is a bond-length modulation, i.e. a
     strain-dependent hopping u_j -- an OFF-DIAGONAL Bloch term.  It preserves
     the sublattice (chiral) symmetry, so it cannot produce a mass; it only
     MOVES the Dirac point (a pseudo-gauge field), keeping the cone gapless until
     an O(1) distortion annihilates it in a Lifshitz transition (u ~ 1, far above
     any thermal phonon).  The contrast is an on-site sublattice mass, which gaps
     the cone at once -- exactly the term a bond phonon is forbidden from making.
     So the mechanical sound cannot impose its cone on the fermions or gap them;
     v_F is protected.

Conclusion: the medium carries two cones -- a mechanical/acoustic one (c_T, the
angular-sourced sound of the node lattice, Section 8.45) and a fermionic one
(v_F, the cone all emergent relativistic matter and its composite gauge/gravity
excitations ride, Sections 8.2, 8.5, 8.12).  They are independent knobs, but the
mismatch is harmless: a phonon couples to the fermions as a pseudo-gauge field,
which protects the Dirac cone, so the emergent Lorentz cone v_F is insulated from
the mechanical sound, and c_T is a decoupled sub-quantum spectator (Volovik).
The physical photon and graviton ride v_F because their kinetic terms are INDUCED
by the fermion loop (Sakharov: test_induced_action, test_induced_gravity), not
inherited from the bare node-lattice elasticity.  Honest scope: the protection is
shown at the level of a static (frozen) phonon acting as a gauge field / tetrad;
a full dynamical phonon self-energy on the fermion cone is a further computation,
but the cone's survival is topological (the Dirac points persist to an O(1)
Lifshitz distortion), not perturbative.
"""
import numpy as np
from test_cone_unification import moduli

SQ3 = np.sqrt(3.0)
ACC = 1.12                                        # relaxed honeycomb bond length

# honeycomb tight-binding: 3 nearest-neighbour bond vectors and a Dirac point
_D = [ACC*np.array([1, 0.]), ACC*np.array([-0.5, SQ3/2]), ACC*np.array([-0.5, -SQ3/2])]
_K = (2*np.pi/(3*ACC))*np.array([1, 1/SQ3])

def fk(k):
    return sum(np.exp(1j*np.dot(k, x)) for x in _D)

def fermi_velocity(t=1.0, q=1e-5):
    """Dirac velocity at K (mean over directions) and its anisotropy."""
    vs = [t*abs(fk(_K + q*np.array([np.cos(a), np.sin(a)])))/q
          for a in np.linspace(0, 2*np.pi, 24, endpoint=False)]
    return float(np.mean(vs)), float(np.ptp(vs)/np.mean(vs))

def hk_bond(k, u):
    """off-diagonal Bloch element with per-bond hopping t_j = 1 + u_j."""
    return sum((1 + u[j])*np.exp(1j*np.dot(k, _D[j])) for j in range(3))

def dirac_gap(u, ng=240):
    """2*min_k |h(k;u)| over the BZ, refined near the minimum."""
    G = 2*np.pi/ACC
    def scan(cx, cy, half, n):
        ax = np.linspace(-half, half, n)
        best = (1e9, 0, 0)
        for kx in cx + ax:
            for ky in cy + ax:
                v = abs(hk_bond(np.array([kx, ky]), u))
                if v < best[0]:
                    best = (v, kx, ky)
        return best
    b = scan(0, 0, G, ng)
    for _ in range(3):                            # refine
        b = scan(b[1], b[2], (2*G/ng), 40); G = 2*G/ng*20
    return 2*b[0]

def gap_staggered(m, u=(0, 0, 0), ng=120):
    """gap with an on-site sublattice mass m: H=[[m,h],[h*,-m]] -> 2*sqrt(m^2+min|h|^2)."""
    G = 2*np.pi/ACC; ax = np.linspace(-G, G, ng)
    mn = min(abs(hk_bond(np.array([kx, ky]), u))
             for kx in ax for ky in ax)
    return 2*np.sqrt(m**2 + mn**2)

def ph_continuum_edge(q, t=1.0, ng=500):
    """Lower edge of the interband particle-hole continuum near K:
    omega_min(q) = min_k t(|f(k+q)| + |f(k)|).  -> v_F|q| as q->0.
    Window scales with |q| so small q stays resolved."""
    half = max(6.0*np.linalg.norm(q), 0.05)
    ax = np.linspace(-half, half, ng)
    KX, KY = np.meshgrid(_K[0] + ax, _K[1] + ax)
    K1 = np.stack([KX.ravel(), KY.ravel()], 1)
    fa = np.abs([fk(k) for k in K1])
    fb = np.abs([fk(k + q) for k in K1])
    return float(t*(fa + fb).min())


def main():
    print('=' * 70)
    print('TWO CONES: is the mechanical phonon cone the fermion (Lorentz) cone?')
    print('=' * 70)

    # --- fermion cone gate ---
    vF, aniso = fermi_velocity(t=1.0)
    print(f'\n[gate] honeycomb Dirac cone: v_F = {vF:.4f} (theory (3/2)t a = '
          f'{1.5*ACC:.4f}), anisotropy {aniso:.1e}')
    assert abs(vF - 1.5*ACC) < 1e-3 and aniso < 1e-4, 'Dirac cone gate failed'

    # --- mechanical cone from Section 8.45 (same force model, lam*) ---
    m = moduli('honeycomb', 0.568)
    cT = m['cT']
    print(f'[mech] mechanical cone at lam*=0.568: c_T = {cT:.3f} '
          f'(K={m["K"]:.2f} central, mu={m["mu"]:.2f} angular)')

    # ---------------------- A: two independent knobs (t moves v_F; lam moves c_T)
    print('\n[A] the two cones answer to DIFFERENT knobs -> independent')
    print('     electronic t -> v_F (mechanics fixed):')
    for t in [0.5, 1.0, 1.5]:
        print(f'        t={t:.1f}:  v_F = {1.5*t*ACC:6.3f}   (c_T unchanged)')
    print('     mechanical lam -> c_T (hopping fixed):')
    cTs = {}
    for lam in [0.3, 0.568, 1.0]:
        cTs[lam] = moduli('honeycomb', lam)['cT']
        print(f'        lam={lam:.3f}:  c_T = {cTs[lam]:6.3f}   (v_F unchanged)')
    assert cTs[1.0] - cTs[0.3] > 0.5, 'c_T should move with lam while v_F does not'
    t_match = cT/(1.5*ACC)
    print(f'     => v_F = c_T only at the tuned hopping t = {t_match:.3f}; '
          f'no symmetry sets t against the mechanical stiffness.')

    # ------------------------------------------- B: physical boson on v_F -----
    print('\n[B] the collective particle-hole boson rides v_F (test_cone_lock)')
    print('     |q|      omega_min/|q|   (-> v_F as q->0)')
    edges = []
    for q in [0.08, 0.05, 0.03, 0.02, 0.01]:
        w = ph_continuum_edge(np.array([q, 0.0]))/q
        edges.append((q, w))
        print(f'     {q:.3f}     {w:8.4f}')
    w_small = edges[-1][1]
    print(f'     => omega_min/|q| -> {w_small:.3f} = v_F = {vF:.3f}: the emergent')
    print(f'        boson is on the FERMION cone, set by t, independent of c_T.')
    assert abs(w_small - vF)/vF < 0.03, 'collective mode should ride v_F'

    # -------------------------- C: is the two-cone mismatch harmful? ----------
    # A phonon is a bond-length modulation -> a strain-dependent hopping u_j.
    # It enters the fermions as an OFF-DIAGONAL (bond) term, so it preserves the
    # sublattice/chiral symmetry: it cannot make a mass. It only MOVES the Dirac
    # point (a pseudo-gauge field) until an O(1) distortion annihilates the cone.
    print('\n[C] the phonon acts as a pseudo-gauge field: it cannot gap the cone')
    print('     at realistic (small) amplitude every pattern stays gapless:')
    print('     phonon pattern      u=0.1   u=0.2   u=0.3   (Dirac gap)')
    pats = [('uniform',       np.array([1., 1, 1])),
            ('shear doublet',  np.array([2., -1, -1])),
            ('single bond',    np.array([1., 0, 0]))]
    for name, pat in pats:
        pn = pat/np.abs(pat).max()
        gaps = [dirac_gap(a*pn) for a in (0.1, 0.2, 0.3)]
        print(f'     {name:16s}  ' + '  '.join(f'{g:6.3f}' for g in gaps))
        assert max(gaps) < 0.05, f'{name} phonon must stay gapless at small u'
    # only an O(1) distortion gaps it, in a Lifshitz transition (pattern-dependent
    # threshold from the bond triangle inequality); uniform NEVER gaps (conformal)
    print('     only an O(1) distortion gaps it (Lifshitz cone annihilation):')
    for name, pat in pats:
        pn = pat/np.abs(pat).max()
        gaps = [dirac_gap(a*pn) for a in (0.5, 1.0)]
        print(f'     {name:16s}  u=0.5:{gaps[0]:5.2f}  u=1.0:{gaps[1]:5.2f}')
    assert dirac_gap(1.0*np.array([1., 1, 1])) < 0.05, 'uniform must never gap'
    assert dirac_gap(1.5*np.array([1., 0, 0])) > 0.5, 'O(1) single bond gaps it'
    # contrast: an on-site SUBLATTICE mass (which a bond phonon cannot make) gaps
    # it immediately -- that is the term the phonon is forbidden from producing
    gm = gap_staggered(0.3)
    print(f'     contrast: an on-site sublattice mass m=0.3 gaps it at once, '
          f'gap = {gm:.3f}')
    print(f'        -- a bond phonon produces no such term (chiral protection).')
    assert gm > 0.5, 'a sublattice mass should gap the cone'

    print('\n' + '=' * 70)
    print('RESULT: the medium has TWO cones -- mechanical c_T (angular-sourced')
    print('node-lattice sound, Section 8.45) and fermionic v_F (electronic) -- and')
    print('they are independent knobs. BUT the mismatch is BENIGN: a phonon is a')
    print('bond (off-diagonal) modulation, so it enters the fermions as a pseudo-')
    print('gauge field that only MOVES the Dirac point; it preserves chiral symmetry')
    print('and cannot make a mass, so it cannot gap the cone or replace v_F with c_T')
    print('(only an O(1) Lifshitz distortion annihilates the cone). So v_F, the')
    print('emergent Lorentz cone, is protected from the mechanical sound. The physical')
    print('photon/graviton ride v_F (fermion-induced, test_induced_action/gravity);')
    print('c_T is a decoupled sub-quantum spectator (Volovik).')
    print('=' * 70)


if __name__ == '__main__':
    main()
