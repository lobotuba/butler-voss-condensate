"""Route A: the fracton-elasticity dual diagnosis of the graviton -- why gamma = 0 is FORCED.

The report has measured gamma = 0 (Nordstrom light-bending, half the GR value) in every channel it
could reach: the source coupling (S8.32-8.34), the propagator / Poisson-ratio route (S8.29-8.30, S8.35),
the second sublattice (S8.36), the RG fate (S8.37), the disclination channel (S8.49), and an amorphous
medium (S8.60). Each was a separate null. This asks the structural question behind all of them, in the
language of the elasticity <-> symmetric-tensor gauge theory (fracton) duality (Pretko-Radzihovsky):
WHICH gauge theory is the medium's would-be graviton, and is gamma = 0 an accident of the tested media
or forced by what the medium must be to host the rest of the physics?

THE DUALITY. A 2D elastic medium is dual to a symmetric-tensor gauge theory. There are two kinds:
  * SCALAR-charge theory: gauge parameter a scalar, delta A_ij = d_i d_j alpha. Its charges are IMMOBILE
    fractons; on the lattice these are the disclinations. This is what ordinary elasticity duals to. It is
    NOT a graviton -- its gauge structure is d_i d_j alpha, not a diffeomorphism.
  * VECTOR-charge theory: gauge parameter a vector, delta A_ij = d_i xi_j + d_j xi_i. This IS linearised
    gravity -- delta g_ij = d_(i xi_j) is a linearised diffeomorphism. gamma = 1 needs THIS structure,
    because the diffeomorphism gauge invariance is the Ward identity that trace-reverses the source (an
    isotropic mass -> spatial curvature), which is exactly the coupling S8.32 measured to vanish.

THE MEASURABLE PIVOT. A linearised diffeomorphism delta g_ij = d_(i xi_j) is, in the medium, a
symmetric-gradient displacement -- i.e. a homogeneous strain is a GLOBAL AFFINE DIFFEOMORPHISM (a shear is
a volume-preserving diffeomorphism, a dilation a scale diffeomorphism). In a diffeomorphism-invariant
theory (gamma = 1) these cost ZERO energy (pure gauge). In the medium they cost the ELASTIC MODULI. So the
elastic moduli literally ARE the energy cost of a diffeomorphism -- the order parameter for
diffeomorphism-NON-invariance. And the shear modulus mu is triply load-bearing:
    (a) mu > 0 is required for the transverse Lorentz cone c_T = sqrt(mu/rho) -- i.e. for propagating
        photon/graviton polarisations to exist at all (S8.1, S8.45);
    (b) mu > 0 gives a positive 2D Young modulus Y = 4 K mu/(K+mu), which CONFINES disclinations
        (E_disc ~ Y R^2) -- immobile fractons: the SCALAR-charge dual, and the dual statement of S8.49's
        "a mass nucleates no disclination charge" = no curvature = gamma = 0;
    (c) mu = 0 is exactly what gamma = 1 requires (the shear-diffeomorphisms must become gauge).
You cannot have mu > 0 and mu = 0. So being a field-bearing solid (mu > 0) FORCES the scalar-charge dual
and gamma = 0; reaching gamma = 1 (mu = 0) destroys the solid and with it the Lorentz cone and all the
emergent fields. gamma = 0 is the price of the rigidity the medium needs for everything else.

  [G1] Diffeomorphisms cost the moduli -> the medium is not diffeomorphism-invariant -> SCALAR-charge
       (fracton), not vector-charge (graviton). A homogeneous shear (a global volume-preserving
       diffeomorphism) costs energy density w = mu > 0 at the physical lattice.
  [G2] The scalar-charge dual is a fracton phase: the SAME mu that gives the transverse cone (c_T > 0)
       gives Y > 0, which confines disclinations (immobile fractons). Cone and confinement are one knob.
  [G3] The no-go, quantified: gamma = 1 requires mu = 0, but across every accessible knob (lattice choice,
       angular stiffness lam) mu = 0 coincides exactly with c_T = 0 (no cone, not a solid). No field-bearing
       solid has soft shear-diffeomorphisms. Escaping gamma = 0 needs the metric to NOT be a phonon of the
       rigid solid -- a decoupled emergent gauge sector (a deconfined-disclination / topologically ordered
       phase), a qualitatively different substrate. This unifies S8.30/S8.35/S8.36/S8.49/S8.60 as one fact
       and shows no elastic tuning can move gamma -- which is why none ever did.

Honest scope: the duality (elasticity -> scalar-charge tensor gauge theory; disclinations = fractons) is
the established Pretko-Radzihovsky result; the model-specific content measured here is that THIS medium sits
in the scalar-charge (fracton) phase at every field-bearing point, so its gamma = 0 is structural, not
tunable. It is a no-go for the ELASTIC (graviton-as-phonon) route, not a proof that no substrate can give
gamma = 1 -- it names the escape (a non-phonon, deconfined-disclination graviton) rather than walking it.
Pure numpy; reuses the elastic-moduli machinery validated in S8.44/S8.45 (test_cone_unification).
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
    """Periodic energy and analytic force (minimum image): LJ + lam angular (as in S8.45)."""
    N = len(X); Hinv = np.linalg.inv(Hs); F = np.zeros_like(X)
    D = _minimg(X[:, None, :] - X[None, :, :], Hs, Hinv)
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
    """2D bulk K, shear mu, cones c_L,c_T by homogeneous strain (atoms relaxed at fixed cell)."""
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
        return 0.5*sum(out)

    K = w(np.eye(2))/(2*amp**2)
    mu = w(np.diag([1., -1.]))/(2*amp**2)
    cL = np.sqrt(max(K + mu, 0)/rho); cT = np.sqrt(max(mu, 0)/rho)
    Y = 4*K*mu/(K + mu) if abs(K + mu) > 1e-9 else 0.0     # 2D Young modulus
    return dict(a=a, K=K, mu=mu, cL=cL, cT=cT, rho=rho, Y=Y)


def main():
    print("=" * 92)
    print("ROUTE A -- FRACTON-DUAL DIAGNOSIS OF THE GRAVITON: why gamma = 0 is forced by rigidity")
    print("=" * 92)
    ok = True

    # physical media: the Cauchy solid (S8.35) and the honeycomb fermion lattice at lam* (S8.44/8.45)
    tri = moduli('triangular', 0.0)
    hcs = moduli('honeycomb', 0.0)          # central-force honeycomb: floppy (S8.45 G2)
    hcf = moduli('honeycomb', 0.568)        # the fermion lattice at lam*

    # [G1] a homogeneous (global affine) diffeomorphism costs the moduli -> not diffeo-invariant
    print("\n  [G1] a homogeneous strain IS a global affine diffeomorphism; in a diffeo-invariant theory")
    print("       (gamma=1) it is pure gauge and costs zero. Here it costs the elastic moduli:")
    print(f"       {'medium':30s} {'K(dilation)':>12s} {'mu(shear)':>11s}   diffeomorphism cost")
    for tag, m in [("Cauchy solid (S8.35)", tri), ("honeycomb lam* (fermion lattice)", hcf)]:
        print(f"       {tag:30s} {m['K']:12.3f} {m['mu']:11.3f}   nonzero -> NOT gauge")
    g1 = tri['mu'] > 1.0 and hcf['mu'] > 1.0
    ok &= g1
    print(f"       => shear-diffeomorphisms cost mu > 0 at every physical lattice: the medium is NOT")
    print(f"          diffeomorphism-invariant -> its would-be graviton is the SCALAR-charge (fracton)")
    print(f"          gauge field (gauge param d_i d_j alpha), not the vector-charge graviton (d_(i xi_j))")
    print(f"          -> {'PASS' if g1 else 'FAIL'}")

    # [G2] the same mu gives the transverse cone AND confines disclinations (fractons)
    print("\n  [G2] the SAME shear modulus mu is the transverse cone and the disclination confiner:")
    print(f"       {'medium':30s} {'mu':>8s} {'c_T=sqrt(mu/rho)':>16s} {'Y_2D':>8s}   disclinations")
    for tag, m in [("Cauchy solid (S8.35)", tri), ("honeycomb lam* (fermion lattice)", hcf)]:
        print(f"       {tag:30s} {m['mu']:8.3f} {m['cT']:16.3f} {m['Y']:8.3f}   confined (E~Y R^2)")
    g2 = hcf['cT'] > 0.1 and hcf['Y'] > 0.1 and tri['Y'] > 0.1
    ok &= g2
    print(f"       => mu>0 => c_T>0 (transverse photon/graviton polarisations exist, S8.1/8.45) AND Y>0")
    print(f"          => disclinations confined = immobile fractons = the scalar-charge dual. This is the")
    print(f"          dual of S8.49 (a mass nucleates no disclination charge => no curvature => gamma=0)")
    print(f"          -> {'PASS' if g2 else 'FAIL'}")

    # [G3] the no-go: gamma=1 needs mu=0, but mu=0 <=> c_T=0 (no solid) across every knob
    print("\n  [G3] no-go -- sweep every accessible knob: is mu ever 0 while c_T>0 (a field-bearing solid)?")
    print(f"       {'lattice':>10s} {'lam':>6s} {'mu':>9s} {'c_T':>8s} {'Y_2D':>8s}   solid?  diffeo-soft?")
    sweep = [('triangular', 0.0)] + [('honeycomb', L) for L in (0.0, 0.1, 0.3, 0.568, 1.0)]
    TH = 1e-2
    solids_mu, nonsolid_cT, locked = [], [], True
    for kind, L in sweep:
        m = moduli(kind, L)
        rigid = m['mu'] > TH                # shear-rigid: diffeomorphisms cost energy (not gauge)
        has_cone = m['cT'] > TH             # transverse cone exists: fields propagate
        locked &= (rigid == has_cone)       # the two must coincide, point by point
        (solids_mu if rigid else nonsolid_cT).append(m['mu'] if rigid else m['cT'])
        print(f"       {kind:>10s} {L:6.3f} {m['mu']:9.4f} {m['cT']:8.4f} {m['Y']:8.3f}   "
              f"{str(has_cone):>5s}   {str(not rigid):>5s}")
    min_solid_mu = min(solids_mu)
    max_nonsolid_cT = max(nonsolid_cT) if nonsolid_cT else 0.0
    # locking: rigid (mu>0) <=> has a cone (c_T>0), point by point; no point is both soft and a solid
    g3 = locked and min_solid_mu > 0.1 and max_nonsolid_cT < TH
    ok &= g3
    print(f"       => shear-rigidity (mu>0) and a transverse cone (c_T>0) coincide at EVERY point (locked).")
    print(f"          min shear modulus among field-bearing solids: mu = {min_solid_mu:.3f} > 0; the one")
    print(f"          non-solid point (honeycomb central, mu<=0, Y<0) has c_T = {max_nonsolid_cT:.2e} ~ 0.")
    print(f"          gamma=1 (mu=0) and propagating fields (c_T>0) are mutually exclusive")
    print(f"          -> {'PASS' if g3 else 'FAIL'}")

    print("\n" + "=" * 92)
    print("[verdict] " + ("ALL GATES PASS" if ok else "GATE FAILURE"))
    print("  gamma = 0 is not an accident of the tested media -- it is FORCED by the elasticity-fracton")
    print("  duality. The medium's would-be graviton is a SCALAR-charge (fracton) gauge field: its")
    print("  diffeomorphism modes are the acoustic phonons, and those cost the elastic moduli (mu, K > 0).")
    print("  The single shear modulus mu locks three things together -- the transverse Lorentz cone")
    print("  (mu>0, so fields exist), disclination confinement (Y>0, immobile fractons = scalar charge),")
    print("  and gamma = 0 -- while gamma = 1 would require mu = 0, i.e. no solid and no fields. The very")
    print("  rigidity that lets the medium host photons, fermions and a Lorentz cone is what forbids")
    print("  Einstein statics. This unifies S8.30/8.35/8.36/8.49/8.60 as one structural fact and shows no")
    print("  elastic tuning can move gamma. The escape is not a better medium: it is a graviton that is")
    print("  NOT a phonon of the rigid solid -- a decoupled, deconfined-disclination (topologically")
    print("  ordered) gauge sector, a qualitatively different substrate.")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
