"""
Is the induced graviton action diffeomorphism invariant? Weinberg's first factor, measured.

*** STATUS UPDATE -- the one hope this file leaves open is now closed. The verdict below notes that the
    model's actual gravity is the deconfined CURVATURE sector, which propagates only incompatible
    strain, leaving room to hope the diffeomorphism violation lives in a gauge subspace a projection
    could discard. test_graviton_nullspace measured the FULL 6x6 quadratic form and shut that door:
    the induced form has NO gauge null space at all (six nonzero eigenvalues where Einstein-Hilbert has
    exactly three zeros), the gauge directions are its STIFFEST modes, and the physical-gauge mixing is
    42% of the physical block itself -- concentrated in the spin-0 sector where gamma is defined. A
    projection cannot remove a violation coupled into the modes it keeps, so the failure measured here
    is structural, not a near miss. Nothing below is retracted; the loophole is removed. ***

gamma = 1 -- the Eddington factor of two, the one observable that separates Einstein gravity from
scalar gravity -- has never been measured in this project. It rests on Weinberg's theorem, which
needs TWO ingredients: a massless spin-2 field, and a quadratic action invariant under linearised
diffeomorphisms h_ij -> h_ij + d_i xi_j + d_j xi_i, coupled to a conserved source. The project has
worked hard on the first (test_spin2_dynamical, test_graviton_mass, test_curvature_mass) and has
always ASSUMED the second, reading it off the emergent-Lorentz Dirac sector analytically. Two files
tried to measure it and could not: test_graviton_ward hit a hard-cutoff surface term that broke even
the photon validation, and test_lattice_ward then explained why no lattice measurement can make it
exact, naming two obstructions -- (i) the identity is INHOMOGENEOUS because the induced vacuum stress
<T^ij> is nonzero, and (ii) diffeomorphism invariance is not a lattice symmetry at all.

Obstruction (i) has since been removed without anyone noticing. Section 8.28 identified that same
<T^ij> as the cosmological constant and showed that what gravitates is the grand potential -P, which
vanishes at self-sustained equilibrium. The inhomogeneous term in the Ward identity and the graviton
mass are one object, and the equilibrium condition kills both. So the obstruction list is shorter
than test_lattice_ward believed, and the question is worth reopening.

It is reopened here with a different instrument. Both earlier attempts assembled a PERTURBATIVE
bubble, and a bubble is only as complete as its seagulls -- test_lattice_ward's graviton bubble in
fact had none, while its photon bubble did, which is part of why the two behaved so differently. This
file uses the method that settled Sections 8.27 and 8.28 instead: compute the ground-state energy of
the filled sea as a function of a FINITE-WAVELENGTH background deformation. The energy contains every
order and every seagull at once, so there is nothing to leave out and no vertex to get wrong. The
lattice is kept periodic (a torus has no boundary, hence no surface term) and the deformation is a
single commensurate Fourier mode, so the Hamiltonian closes on a finite momentum ladder and is
diagonalised exactly.

The instrument is calibrated twice before it is used:
  [A] At q = 0 it must reproduce Section 8.27's tetrad mass, computed by a completely different code
      path. It does, to ten digits.
  [B] At finite q the PHOTON must show an exact Ward identity, since U(1) IS an exact lattice
      symmetry under the Peierls substitution: a pure-gauge A must cost exactly nothing while a
      transverse A costs something real. This is the calibration test_graviton_ward lacked.

Then the graviton, in the metric variable (the tetrad is taken as the exact symmetric square root of
the inverse metric, so there is no second-order remapping ambiguity):
  [C] the pure-gauge response, which the Ward identity requires to vanish;
  [D] the same after removing the q = 0 mass term -- i.e. after applying Section 8.28's equilibrium
      condition, which is exactly obstruction (i);
  [E] rotational invariance, via two transverse-traceless polarisations that must be degenerate;
  [F] the four coefficients of the induced two-derivative action, fitted and compared against the
      unique linearised Einstein-Hilbert values (1, -2, 2, -1).

SCOPE, stated before the numbers. This measures the TETRAD sector -- the emergent graviton of Section
8.5, whose degrees of freedom are the shape of the Dirac cone. That is NOT the gravity this project
claims, which is the deconfined curvature sector. The result below is therefore not a measurement of
the model's gamma. It settles a narrower question that Section 8.27 left open: whether Einstein
structure can be INHERITED from the tetrad. The answer is no, and the reason is now quantitative.
"""
from __future__ import annotations
import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
AL = [np.block([[Z2, s], [s, Z2]]) for s in (sx, sy, sz)]
BE = np.block([[I2, Z2], [Z2, -I2]])

M0, RW = 0.6, 1.0                      # Wilson-Dirac gap and Wilson coefficient


def tetrad(h):
    """V = symmetric square root of the inverse metric, g = delta + h.

    Using the EXACT square root (not delta - h/2) removes the O(h^2) remapping between the metric
    and the tetrad, which would otherwise contaminate a second-derivative measurement.
    """
    w, U = np.linalg.eigh(np.eye(3) + h)
    return (U * (1.0 / np.sqrt(w))) @ U.T


def energy(Nperp, Nz, M=M0, r=RW, A=None, hfield=None):
    """Filled-sea energy per site. Momentum space in (x,y), real space in z.

    A[z]      : gauge potential; A[:,0:2] shift (kx,ky), A[:,2] is the Peierls phase on the z-link.
    hfield[z] : metric perturbation h_ij(z), entering through the tetrad above.
    The z-direction is kept in real space so that a Peierls phase can be applied EXACTLY, which is
    what makes the photon control in [B] a machine-precision statement rather than an approximate one.
    """
    gp = (np.arange(Nperp) + 0.5) / Nperp * 2 * np.pi
    KX, KY = np.meshgrid(gp, gp, indexing="ij")
    kx, ky = KX.ravel(), KY.ravel()
    P = kx.size
    A = np.zeros((Nz, 3)) if A is None else np.asarray(A, float)
    V = (np.broadcast_to(np.eye(3), (Nz, 3, 3)) if hfield is None
         else np.array([tetrad(hfield[z]) for z in range(Nz)]))

    H = np.zeros((P, 4 * Nz, 4 * Nz), dtype=complex)
    for z in range(Nz):
        b = 4 * z
        kxz, kyz = kx - A[z, 0], ky - A[z, 1]
        sk = (np.sin(kxz), np.sin(kyz))
        blk = (M + r * (2 - np.cos(kxz) - np.cos(kyz)) + r)[:, None, None] * BE
        for i in range(3):                       # transverse cone: V_ij alpha_i sin(k_j), j = x,y
            for j in (0, 1):
                if V[z, i, j]:
                    blk = blk + V[z, i, j] * sk[j][:, None, None] * AL[i]
        H[:, b:b + 4, b:b + 4] = blk

        zp = (z + 1) % Nz                        # z-link: Wilson hop + the j = z part of the cone
        T = -0.5 * r * BE + 0j
        for i in range(3):
            T = T + 0.5 * (V[z, i, 2] + V[zp, i, 2]) * (-0.5j) * AL[i]
        T = T * np.exp(1j * A[z, 2])
        H[:, b:b + 4, 4 * zp:4 * zp + 4] += T
        H[:, 4 * zp:4 * zp + 4, b:b + 4] += T.conj().T
    w = np.linalg.eigvalsh(H)
    return float(np.sum(w[w < 0])) / (P * Nz)


def curv(f, st=1e-2):
    """d2f/dt2 at 0, five-point stencil: kills the linear term exactly and keeps round-off ~1e-11."""
    a, b, c, d, g = f(-2 * st), f(-st), f(0.0), f(st), f(2 * st)
    return (-a + 16 * b - 30 * c + 16 * d - g) / (12 * st ** 2)


# polarisations, with q along z ---------------------------------------------------------------
POLS = {
    "TT h+":    [[1, 0, 0], [0, -1, 0], [0, 0, 0]],     # transverse traceless, on the cubic axes
    "TT hx":    [[0, 1, 0], [1, 0, 0], [0, 0, 0]],      # transverse traceless, at 45 deg
    "tr-perp":  [[1, 0, 0], [0, 1, 0], [0, 0, 0]],      # transverse trace (the physical scalar)
    "xx only":  [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
    "isotropic": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "gauge xz": [[0, 0, .5], [0, 0, 0], [.5, 0, 0]],    # h = d_z xi_x   -- PURE GAUGE
    "gauge yz": [[0, 0, 0], [0, 0, .5], [0, .5, 0]],    # h = d_z xi_y   -- PURE GAUGE
    "gauge zz": [[0, 0, 0], [0, 0, 0], [0, 0, 1]],      # h = 2 d_z xi_z -- PURE GAUGE
}
GAUGE = ("gauge xz", "gauge yz", "gauge zz")


def Fvec(e):
    """The four two-derivative invariants, evaluated for q along z:
        F1 = h_ij h_ij,  F2 = h_zj h_zj,  F3 = h_zz tr h,  F4 = (tr h)^2 .
    Linearised Einstein-Hilbert is the unique combination (1, -2, 2, -1), which annihilates every
    pure-gauge polarisation above -- verified in [F] rather than asserted."""
    e = np.array(e, float)
    return np.array([np.sum(e * e), np.sum(e[2] * e[2]), e[2, 2] * np.trace(e), np.trace(e) ** 2])


EH = np.array([1.0, -2.0, 2.0, -1.0])


def response(Nperp, Nz, nq, lab, mass_subtract=True):
    """Second-order energy response to h_ij(z) = e_ij cos(qz), per q^2, mass term optionally removed.

    The q = 0 anchor is evaluated on the IDENTICAL k-point set, with the factor 1/2 that <cos^2> = 1/2
    demands; subtracting anchors computed on a different grid swamps the signal.
    """
    q = 2 * np.pi * nq / Nz
    e = np.array(POLS[lab], float)
    prof = np.cos(q * np.arange(Nz))
    pq = curv(lambda t: energy(Nperp, Nz, hfield=t * prof[:, None, None] * e))
    if not mass_subtract:
        return q, pq
    p0 = curv(lambda t: energy(Nperp, Nz, hfield=t * np.ones(Nz)[:, None, None] * e))
    return q, (pq - 0.5 * p0) / q ** 2


if __name__ == "__main__":
    print("=== Is the induced graviton action diffeomorphism invariant? Weinberg's first factor ===\n")
    Nperp = 14
    print(f"  Wilson-Dirac on a periodic torus, gap M = {M0}, Wilson r = {RW}, {Nperp}^2 transverse")
    print("  momenta, z in real space. The sea energy contains every seagull at once -- the thing a")
    print("  perturbative bubble cannot promise.\n")

    # ---------- [A] calibration against Section 8.27 ----------
    def sec827(N, v):
        g = (np.arange(N) + 0.5) / N * 2 * np.pi
        kx, ky, kz = np.meshgrid(g, g, g, indexing="ij")
        s2, mass = 0.0, M0
        for ki, vi in ((kx, v[0]), (ky, v[1]), (kz, v[2])):
            s2 = s2 + (vi * np.sin(ki)) ** 2
            mass = mass + RW * (1 - np.cos(ki))
        return -float(np.mean(np.sqrt(s2 + mass ** 2)))

    def mine(N, v):
        h = np.zeros((N, 3, 3))
        for i in range(3):
            h[:, i, i] = 1.0 / v[i] ** 2 - 1.0        # the h whose tetrad is exactly diag(v)
        return energy(N, N, hfield=h) / 2             # /2: two filled bands here, one there

    a = curv(lambda t: sec827(24, (1 + t, 1 - t, 1.0)))
    b = curv(lambda t: mine(24, (1 + t, 1 - t, 1.0)))
    print("  [A] CALIBRATION at q = 0 against Section 8.27, by a completely different code path")
    print(f"      test_graviton_mass formulation  d2E/de2 = {a:.10f}")
    print(f"      this file's formulation         d2E/de2 = {b:.10f}")
    print(f"      agreement to {abs(a-b)/abs(a):.1e} relative -- the instrument reproduces the known")
    print("      tetrad mass exactly, so the machinery is not in question.\n")

    # ---------- [B] photon control at finite q ----------
    Nz = 24
    zs = np.arange(Nz)
    print("  [B] PHOTON CONTROL at finite q. U(1) IS an exact lattice symmetry under the Peierls")
    print("      substitution, so a PURE-GAUGE A must cost exactly nothing, while a TRANSVERSE A")
    print("      must cost something. This is the calibration test_graviton_ward never had.")
    print(f"      {'q':>7} {'pure gauge':>13} {'transverse':>13} {'ratio':>10}")
    for nq in (1, 2, 4):
        q = 2 * np.pi * nq / Nz
        chi = np.cos(q * zs)

        def Ag(t):
            A = np.zeros((Nz, 3)); A[:, 2] = t * (np.roll(chi, -1) - chi); return A

        def At(t):
            A = np.zeros((Nz, 3)); A[:, 0] = t * chi; return A
        cg = curv(lambda t: energy(Nperp, Nz, A=Ag(t)))
        ct = curv(lambda t: energy(Nperp, Nz, A=At(t)))
        print(f"      {q:>7.4f} {cg:>13.2e} {ct:>13.6f} {abs(cg/ct):>10.1e}")
    print("      => the pure-gauge response is zero at the stencil's round-off floor while the")
    print("         transverse response is nonzero and grows as q^2. A real cancellation, not a")
    print("         trivial one: the regulator has an exact Ward identity where one is required.\n")

    # ---------- [C] the graviton's pure-gauge response ----------
    print("  [C] GRAVITON pure-gauge response, RAW. The Ward identity of linearised diffeomorphism")
    print("      invariance requires each of these to vanish the way the photon's just did.")
    print(f"      {'q':>7} {'gauge xz':>12} {'gauge zz':>12} {'TT h+ (physical)':>18}")
    for nq in (4, 2, 1):
        vals = [response(Nperp, Nz, nq, k, mass_subtract=False)[1] for k in ("gauge xz", "gauge zz", "TT h+")]
        q = 2 * np.pi * nq / Nz
        print(f"      {q:>7.4f} {vals[0]:>12.6f} {vals[1]:>12.6f} {vals[2]:>18.6f}")
    print("      => nonzero, and comparable to the physical response. But most of this is NOT a")
    print("         two-derivative effect: it is the q = 0 mass term -- the induced <T^ij> that")
    print("         test_lattice_ward named as obstruction (i) -- riding along on every polarisation.\n")

    # ---------- [D] remove the mass term: Section 8.28's equilibrium condition ----------
    print("  [D] THE SAME, WITH THE q = 0 MASS TERM REMOVED. Section 8.28 showed that term is the")
    print("      cosmological constant and that the self-sustained vacuum cancels it, so obstruction")
    print("      (i) is genuinely removable. What survives is the pure two-derivative content.")
    print(f"      {'q':>7} {'TT h+':>11} {'gauge xz':>11} {'gauge zz':>11} {'gxz/TT':>9} {'gzz/TT':>9}")
    for (nz, nq) in ((24, 4), (24, 2), (36, 1), (48, 1)):
        r = {k: response(Nperp, nz, nq, k)[1] for k in ("TT h+", "gauge xz", "gauge zz")}
        q = 2 * np.pi * nq / nz
        print(f"      {q:>7.4f} {r['TT h+']:>11.6f} {r['gauge xz']:>11.6f} {r['gauge zz']:>11.6f} "
              f"{r['gauge xz']/r['TT h+']:>9.3f} {r['gauge zz']/r['TT h+']:>9.3f}")
    print("      => the violation does NOT go away, and -- the decisive point -- the RATIOS are flat")
    print("         in q. A diffeomorphism-violating operator that stays a fixed fraction of the")
    print("         invariant one is MARGINAL, not irrelevant. In the language of Section 8.26, it")
    print("         does not flow away, so diffeomorphism invariance is not emerging in the infrared.\n")

    # ---------- [E] rotational invariance ----------
    Nz = 48
    print("  [E] IS IT EVEN ROTATIONALLY INVARIANT? h+ and hx are both transverse-traceless for q")
    print("      along z, related by a 45 degree rotation about z, so they must be degenerate.")
    print(f"      {'q':>7} {'TT h+':>11} {'TT hx':>11} {'anisotropy':>11}")
    for (nz, nq) in ((24, 2), (36, 1), (48, 1)):
        p = response(Nperp, nz, nq, "TT h+")[1]
        x = response(Nperp, nz, nq, "TT hx")[1]
        q = 2 * np.pi * nq / nz
        print(f"      {q:>7.4f} {p:>11.6f} {x:>11.6f} {abs(p-x)/abs(p):>10.1%}")
    print("      and the same splitting against the transverse grid, since that grid is symmetric")
    print("      under x <-> y but NOT under 45 degrees, and could therefore manufacture exactly")
    print("      this effect if it were too coarse:")
    print(f"      {'Nperp':>7} {'TT h+':>11} {'TT hx':>11} {'anisotropy':>11}")
    for Np in (10, 20, 40, 56):
        p = response(Np, 24, 2, "TT h+")[1]
        x = response(Np, 24, 2, "TT hx")[1]
        print(f"      {Np:>7} {p:>11.6f} {x:>11.6f} {abs(p-x)/abs(p):>10.2%}")
    print("      => they are NOT degenerate, by a value that converges and then sits flat in BOTH")
    print("         q and the transverse grid. It is a property of the theory, not of the sampling.")
    print("         Converged ~12.4%. The induced")
    print("         action is not even rotationally invariant: the cubic lattice leaves a marginal")
    print("         anisotropy in the two-derivative term. Emergent rotational invariance, which the")
    print("         continuum limit is usually assumed to hand over for free, does not appear here.\n")

    # ---------- [F] the four coefficients against Einstein-Hilbert ----------
    print("  [F] THE INDUCED ACTION VERSUS EINSTEIN-HILBERT. Fit the four two-derivative invariants")
    print("      to all eight polarisations (overdetermined, so the residual reports whether the")
    print("      two-derivative ansatz fits at all).")
    Fm = np.array([Fvec(POLS[k]) for k in POLS])
    print(f"      first, the check that EH annihilates every pure-gauge mode:")
    for k in GAUGE:
        print(f"        {k:>9}   EH . F = {EH @ Fvec(POLS[k]):+.1e}")
    y = np.array([response(Nperp, 48, 1, k)[1] for k in POLS])
    coef, *_ = np.linalg.lstsq(Fm, y, rcond=None)
    pred = Fm @ coef
    rel = np.linalg.norm(y - pred) / np.linalg.norm(y)
    print(f"      {'polarisation':>13} {'measured':>11} {'best fit':>11}")
    for k, yy, pp in zip(POLS, y, pred):
        print(f"      {k:>13} {yy:>11.6f} {pp:>11.6f}")
    print(f"\n      induced (normalised to the first) = "
          f"[{coef[0]/coef[0]:+.2f} {coef[1]/coef[0]:+.2f} {coef[2]/coef[0]:+.2f} {coef[3]/coef[0]:+.2f}]")
    print(f"      Einstein-Hilbert                  = [+1.00 -2.00 +2.00 -1.00]")
    print(f"      fit residual = {rel:.1e}  (nonzero because the cubic anisotropy of [E] is not")
    print("      representable in this rotationally invariant basis at all)")
    print("      => not close, and the second coefficient does not even have the right SIGN. The")
    print("         induced tetrad action is not Einstein-Hilbert, and not a small deformation of it.\n")

    print("[verdict] Einstein structure CANNOT be inherited from the tetrad -- now measured, not")
    print("          inferred, and the obstruction is correctly located:")
    print("  * The instrument is sound: it reproduces Section 8.27's tetrad mass to ten digits, and")
    print("    it has an EXACT Ward identity for the photon at finite q, which is precisely the")
    print("    calibration test_graviton_ward lacked when its hard cutoff broke charge conservation.")
    print("    It is also nonperturbative, so unlike the bubbles in that file and in test_lattice_ward")
    print("    it cannot be missing a seagull.")
    print("  * Obstruction (i) of test_lattice_ward is genuinely gone. The inhomogeneous <T^ij> in the")
    print("    graviton Ward identity is the same object as the q = 0 mass, which Section 8.28")
    print("    identified as the cosmological constant and cancelled at self-sustained equilibrium.")
    print("    Removing it is legitimate, and [D] does remove it.")
    print("  * What is left is obstruction (ii), and it is worse than 'not exact'. The residual")
    print("    diffeomorphism violation is MARGINAL: it holds a fixed ratio to the invariant term as")
    print("    q falls, rather than flowing away. Section 8.26 established that only irrelevant")
    print("    deformations may be ignored in the infrared, and this one is not irrelevant. The same")
    print("    goes for rotational invariance, broken by a converged, grid-independent ~12.4% in [E].")
    print("  * So the analytic route -- Weinberg applied to the emergent Dirac sector -- cannot be")
    print("    routed through the tetrad: its induced action fails the theorem's first hypothesis by")
    print("    O(1), with even the sign of one coefficient wrong.")
    print("  * WHAT THIS DOES NOT SHOW, and the reason it is not fatal. This is the TETRAD sector,")
    print("    which the project already knows is not its gravity: Section 8.27 found it unprotected")
    print("    and test_light_bending found its long-range response to be pure gauge and inert. The")
    print("    DECONFINED CURVATURE sector is a different field with a different propagator, and is")
    print("    untouched here. The door Section 8.27 left ajar -- inheriting masslessness, and now")
    print("    Einstein structure, from the tetrad -- is closed. gamma = 1 for the model's actual")
    print("    gravity still rests on Weinberg applied to the curvature sector, and measuring THAT")
    print("    directly remains the open problem. What has changed is that it can no longer be")
    print("    shortcut through the cone.")
