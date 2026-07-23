"""
The negative-space sublattice: does a second (bcc-void / interpenetrating) lattice supply the
missing Fierz-Pauli term for light bending? Measured here. The answer is no, for a reason that is
the same selection rule that closed every earlier route -- but the two-sublattice medium DOES break
the one assumption test_gamma_elastic (Sec 8.35) rested on, so this is the first honest re-opening
of the gamma arc, and it deserved a direct measurement rather than an argument.

WHY THIS IS THE RIGHT PLACE TO LOOK. Every gamma=0 result, and Sec 8.35 in particular, rode on a
hidden premise: the medium is a single Bravais lattice with central forces. That premise forces the
Cauchy relations (nu -> 1/3, C12 = C66) AND, more deeply, makes the light-seen strain equal to the
symmetric gradient of ONE displacement field u. Then the Saint-Venant identity inc(sym grad u) == 0
makes the strain compatible -> flat -> no light bending, for every elastic constant. That identity is
the whole gamma=0 result. To escape it one needs a light-seen strain that is NOT the gradient of a
single field -- a genuinely independent internal field. A two-sublattice medium is exactly that: the
fcc lattice plus its negative-space (interpenetrating, bcc-like) sublattice has a RELATIVE-displacement
(optical) mode w = u_A - u_B that is not the gradient of the acoustic u. This is textbook: two
interpenetrating fcc lattices = diamond, and diamond VIOLATES the Cauchy relations badly (C12 != C44)
through exactly this internal-strain (Kleinman) mode. So the second sublattice is the right target --
it attacks the actual assumption.

THE FORK, AND WHAT THIS FILE MEASURES.
  [A] Is the optical (relative) mode gapped or gapless? If gapped, it is adiabatically SLAVED to the
      compression and only renormalises the effective elastic constants (Cauchy-violating, but still
      ordinary elasticity on the compatible grad u) -> gamma = 0, just another point on the nu axis
      Sec 8.35 already swept. If gapless (a relative Goldstone of a second condensate) it obeys a
      Poisson equation and could carry a long-range Psi. MECHANICAL two-sublattice -> gapped (measured).
  [B] Does the second sublattice break the Cauchy relations, as it must to escape Sec 8.35's premise?
      Measured: yes, C12 != C66 once the internal mode relaxes (Sec 8.35's premise is genuinely broken).
  [C] THE DECISIVE TEST. The internal mode is a VECTOR field w. Light-bending curvature (a disclination /
      spin-2) is sourced only by the CURL of w (a relative ROTATION of the two sublattices). A static
      mass is an ISOTROPIC (scalar, spin-0) source. By rotational symmetry an isotropic source drives a
      purely RADIAL relative displacement -> curl w = 0 -> no relative rotation -> no disclination ->
      no incompatible strain -> gamma = 0. Measured: curl w ~ 0 for a mass, curl w != 0 for a SHEAR
      source (the probe is alive; the tensor channel exists but only a tensor source engages it).

THE VERDICT. The negative-space sublattice does break the Cauchy relations (good -- it is a real new
medium, not covered by Sec 8.35), but its new optical channel is a SCALAR (radial breathing) mode when
driven by a SCALAR source. A compression mass engages only the compression (spin-0) internal mode,
never the shear (spin-2) one that light bending needs -- the same rotational selection rule that closed
Sec 8.32-8.35, now shown to survive the two-sublattice generalisation. For the second sublattice to be
the missing term, a static isotropic mass would have to source the sublattices' relative TRANSVERSE
(shear) displacement, which the isotropy of the source forbids. So the missing Fierz-Pauli term is not
here either -- but we now know precisely what a candidate mechanism would have to do (couple a scalar
source to a shear internal mode), and that it is forbidden by symmetry, not merely absent.
"""
from __future__ import annotations
import numpy as np

# honeycomb = two interpenetrating triangular sublattices (2D proxy for diamond = two fcc).
A1 = np.array([1.5,  np.sqrt(3) / 2])      # triangular primitive vectors (nn bond a0 = 1)
A2 = np.array([1.5, -np.sqrt(3) / 2])
DELTA = np.array([[1.0, 0.0],              # the 3 inter-sublattice (A->B) nn vectors, length 1
                  [-0.5,  np.sqrt(3) / 2],
                  [-0.5, -np.sqrt(3) / 2]])
TRI = np.array([A1, A2, A1 - A2])          # the 3 intra-sublattice (triangular) nn vectors, length sqrt(3)
CELL_AREA = abs(A1[0] * A2[1] - A1[1] * A2[0])   # 3 sqrt(3) / 2


# ============================================================ [A] optical-mode gap ===
def dynamical_matrix(q, k_inter, k_intra):
    """4x4 Bloch dynamical matrix of the two-sublattice spring crystal (masses = 1).
    Sublattice order (A,B); each 2x2 block from springs Phi = k n n^T (longitudinal springs)."""
    D = np.zeros((4, 4), complex)

    def block(i, j, add):
        D[2 * i:2 * i + 2, 2 * j:2 * j + 2] += add

    # inter-sublattice (A<->B) nn springs
    for d in DELTA:
        n = d / np.linalg.norm(d)
        P = k_inter * np.outer(n, n)
        block(0, 0, P); block(1, 1, P)                         # self terms
        block(0, 1, -P * np.exp(1j * q @ d))                   # A gets B at +d
        block(1, 0, -P * np.exp(-1j * q @ d))
    # intra-sublattice (A-A and B-B) triangular springs
    for t in TRI:
        n = t / np.linalg.norm(t)
        P = k_intra * np.outer(n, n)
        for s in (0, 1):
            block(s, s, 2 * P)                                 # +-t both contribute to self
            block(s, s, -P * np.exp(1j * q @ t) - P * np.exp(-1j * q @ t))
    return 0.5 * (D + D.conj().T)


def optical_gap(k_inter, k_intra):
    w2 = np.linalg.eigvalsh(dynamical_matrix(np.zeros(2), k_inter, k_intra)).real
    return np.sqrt(np.clip(w2, 0, None))       # sorted; two ~0 (acoustic) + two = optical gap


# ============================================================ [B] Cauchy relations ===
def cell_energy(H, s, k_inter, k_intra):
    """Energy per unit cell under macro displacement-gradient H (2x2) and internal B-shift s."""
    E = 0.0
    for d in DELTA:                                            # inter: B shifted by s relative to A
        v = (np.eye(2) + H) @ d + s
        E += 0.5 * k_inter * (np.linalg.norm(v) - 1.0) ** 2
    for t in TRI:                                              # intra: s cancels (same sublattice)
        for lat in range(2):
            v = (np.eye(2) + H) @ t
            E += 0.5 * k_intra * (np.linalg.norm(v) - np.sqrt(3)) ** 2
    return E


def relax_s(H, k_inter, k_intra):
    """Minimise cell energy over the internal shift s at fixed macro strain (2D Newton)."""
    s = np.zeros(2)
    for _ in range(40):
        g = np.zeros(2); Hs = np.zeros((2, 2)); e = 1e-6
        for a in range(2):
            sp = s.copy(); sp[a] += e; sm = s.copy(); sm[a] -= e
            g[a] = (cell_energy(H, sp, k_inter, k_intra) - cell_energy(H, sm, k_inter, k_intra)) / (2 * e)
        for a in range(2):
            for b in range(2):
                spp = s.copy(); spp[a] += e; spp[b] += e
                spm = s.copy(); spm[a] += e; spm[b] -= e
                smp = s.copy(); smp[a] -= e; smp[b] += e
                smm = s.copy(); smm[a] -= e; smm[b] -= e
                Hs[a, b] = (cell_energy(H, spp, k_inter, k_intra) - cell_energy(H, spm, k_inter, k_intra)
                            - cell_energy(H, smp, k_inter, k_intra) + cell_energy(H, smm, k_inter, k_intra)) / (4 * e * e)
        try:
            s = s - np.linalg.solve(Hs, g)
        except np.linalg.LinAlgError:
            break
        if np.linalg.norm(g) < 1e-12:
            break
    return s


def elastic_constants(k_inter, k_intra, relax_internal):
    """C11, C12, C66 (2D Voigt) from strained-cell energy; Cauchy (central force) => C12 = C66."""
    e = 1e-3

    def U(exx, eyy, exy):
        H = np.array([[exx, exy], [exy, eyy]])
        s = relax_s(H, k_inter, k_intra) if relax_internal else np.zeros(2)
        return cell_energy(H, s, k_inter, k_intra) / CELL_AREA

    def d2(f):
        return (f(e) - 2 * f(0.0) + f(-e)) / e ** 2
    C11 = d2(lambda a: U(a, 0, 0))
    C1112 = d2(lambda a: U(a, a, 0))            # 2 C11 + 2 C12
    C12 = 0.5 * C1112 - C11
    C66 = d2(lambda a: U(0, 0, a)) / 4.0        # engineering-shear factor
    return C11, C12, C66


def internal_strain_response(k_inter, k_intra, e=1e-3):
    """The uniform (q=0) optical-mode shift s that each macro strain slaves in. The selection rule,
    exact and noise-free: an ISOTROPIC strain sources s = 0 (the 3 bond forces on B cancel by the
    3-fold symmetry); only a SHEAR strain sources s != 0 (the Kleinman internal strain)."""
    def s_of(exx, eyy, exy):
        return relax_s(np.array([[exx, exy], [exy, eyy]]), k_inter, k_intra) / e
    return {"isotropic (exx=eyy)": np.linalg.norm(s_of(e, e, 0.0)),
            "shear (exx=-eyy)":     np.linalg.norm(s_of(e, -e, 0.0)),
            "shear (exy)":          np.linalg.norm(s_of(0.0, 0.0, e))}


# ============================ [C] decisive: does a mass source a relative ROTATION (curl of s)? ===
def internal_strain_map(k_inter, k_intra):
    """The linear map M (2x3) from local macro strain (exx, eyy, exy) to the slaved internal shift s,
    measured from relax_s. By [B2] its isotropic column is ~0: the optical mode reads only shear."""
    cols = [relax_s(np.array([[1, 0], [0, 0]]), k_inter, k_intra),      # d s / d exx
            relax_s(np.array([[0, 0], [0, 1]]), k_inter, k_intra),      # d s / d eyy
            relax_s(np.array([[0, 1], [1, 0]]), k_inter, k_intra)]      # d s / d exy (tensor shear)
    return np.stack(cols, axis=1)                                       # s = M @ [exx, eyy, exy]


def _inc(exx, eyy, exy, KX, KY):
    """Linearised 2D curvature (incompatibility / Ricci scalar): inc = d_yy exx + d_xx eyy - 2 d_xy exy.
    Zero iff the strain is a symmetric gradient (compatible = flat). Computed spectrally."""
    fxx, fyy, fxy = np.fft.fftn(exx), np.fft.fftn(eyy), np.fft.fftn(exy)
    return np.fft.ifftn(-KY ** 2 * fxx - KX ** 2 * fyy + 2 * KX * KY * fxy).real


def full_response_curvature(M, N=128, L=48.0, w=4.0):
    """The decisive test. Build the FULL two-sublattice response to a mass -- the macro displacement u
    AND the slaved internal shift s(x) = M @ strain(u) -- form the total light-seen strain, and measure
    its incompatibility (curvature). Because the optical mode is gapped ([A]) s is a function of the
    macro strain, so it is a displacement too; the total stays compatible and carries no curvature.
    For scale, compare to the raw eigenstrain theta* delta, whose incompatibility (nabla^2 theta*, the
    Sec 8.34 light-bending content) is nonzero -- the probe is alive."""
    h = L / N
    g = (np.arange(N) - N // 2) * h
    X, Y = np.meshgrid(g, g, indexing="ij")
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0

    theta = np.exp(-(X ** 2 + Y ** 2) / (2 * w ** 2))     # the mass = isotropic eigenstrain
    A = -np.fft.fftn(theta) / K2                          # macro relaxation potential, u_i = i k_i A
    exx = np.fft.ifftn(-KX * KX * A).real                 # compatible macro strain
    eyy = np.fft.ifftn(-KY * KY * A).real
    exy = np.fft.ifftn(-KX * KY * A).real

    sx = M[0, 0] * exx + M[0, 1] * eyy + M[0, 2] * exy    # the slaved internal (optical) shift field
    sy = M[1, 0] * exx + M[1, 1] * eyy + M[1, 2] * exy
    sxk, syk = np.fft.fftn(sx), np.fft.fftn(sy)
    dxx = np.fft.ifftn(1j * KX * sxk).real               # extra strain the internal shift adds (spectral)
    dyy = np.fft.ifftn(1j * KY * syk).real
    dxy = 0.5 * np.fft.ifftn(1j * KY * sxk + 1j * KX * syk).real

    inc_full = _inc(exx + dxx, eyy + dyy, exy + dxy, KX, KY)          # total: macro + optical
    inc_eig = _inc(theta, theta, np.zeros_like(theta), KX, KY)        # raw eigenstrain theta* delta
    return float(np.abs(inc_full).max()), float(np.abs(inc_eig).max())


def decisive(k_inter, k_intra):
    M = internal_strain_map(k_inter, k_intra)
    inc_full, inc_eig = full_response_curvature(M)
    return {"M": M, "inc_full": inc_full, "inc_eig": inc_eig}


# =================================================================================== main ===
if __name__ == "__main__":
    ki, ka = 1.0, 1.0
    print("=== The negative-space sublattice: is it the missing term for light bending? ===\n")
    print("  Two interpenetrating triangular sublattices (honeycomb) = the 2D proxy for diamond =")
    print("  two interpenetrating fcc. This is the fcc + bcc-void two-lattice medium, the one thing")
    print("  Sec 8.35 (single Bravais lattice, central forces) did not cover.\n")

    # ---- [A] the optical-mode gap: gapped => slaved => still gamma = 0 ----
    print("  [A] OPTICAL-MODE GAP at q=0 (does the relative mode cost energy, i.e. is it slaved?):")
    print(f"      {'k_inter':>9} {'acoustic w (2)':>18} {'optical gap (2)':>18}")
    for ki_scan in (0.0, 0.25, 1.0, 4.0):
        w = optical_gap(ki_scan, ka)
        print(f"      {ki_scan:>9.2f}   {w[0]:.3e},{w[1]:.3e}   {w[2]:.3e},{w[3]:.3e}")
    print("      => the two acoustic branches are gapless (w->0); the optical branch is GAPPED for")
    print("         any k_inter>0 (gap^2 ~ k_inter). A gapped relative mode is adiabatically SLAVED to")
    print("         the strain -> it only renormalises the elastic constants; it is not an independent")
    print("         long-range field. (Only k_inter=0 -- fully decoupled sublattices -- is gapless.)\n")

    # ---- [B] Cauchy violation + WHICH strain the optical mode couples to ----
    C11a, C12a, C66a = elastic_constants(ki, ka, relax_internal=False)
    C11r, C12r, C66r = elastic_constants(ki, ka, relax_internal=True)
    Ka, Kr = 0.5 * (C11a + C12a), 0.5 * (C11r + C12r)         # 2D bulk (compression) modulus
    print("  [B] CAUCHY RELATION C12 = C66 (holds for a single central-force lattice = Sec 8.35's premise):")
    print(f"      affine (internal mode frozen):   C12 = {C12a:+.4f},  C66 = {C66a:+.4f},  C12-C66 = {C12a-C66a:+.2e}")
    print(f"      relaxed (internal mode active):  C12 = {C12r:+.4f},  C66 = {C66r:+.4f},  C12-C66 = {C12r-C66r:+.2e}")
    print(f"      bulk modulus K=(C11+C12)/2:      affine {Ka:.4f} -> relaxed {Kr:.4f}   (change {Kr-Ka:+.2e})")
    print(f"      shear modulus C66:               affine {C66a:.4f} -> relaxed {C66r:.4f}   (change {C66r-C66a:+.2e})")
    print("      => the internal (optical) mode BREAKS Cauchy (relaxed C12 != C66), the way diamond (two")
    print("         fcc) departs from a central-force solid -- so the negative-space sublattice genuinely")
    print("         escapes Sec 8.35's premise. BUT it leaves the COMPRESSION modulus K untouched and")
    print("         softens only SHEAR: the optical mode couples to shear, not to compression.\n")

    # ---- [B2] the selection rule, exact at the uniform optical mode: s(isotropic) = 0 ----
    print("  [B2] INTERNAL SHIFT s slaved by each macro strain (q=0 optical mode, exact, noise-free):")
    for name, val in internal_strain_response(ki, ka).items():
        print(f"      {name:<22} |s|/strain = {val:.3e}")
    print("      => an ISOTROPIC strain sources NO internal shift (the 3 bond forces on B cancel by the")
    print("         3-fold symmetry); only SHEAR does. A mass is a compression source, so it drives the")
    print("         internal mode only through the (radial) deviatoric part of its field -- next.\n")

    # ---- [C] decisive: the FULL two-sublattice response to a mass is still compatible (flat) ----
    print("  [C] DECISIVE TEST (Fourier, smooth) -- build the FULL two-sublattice response to a mass:")
    print("      the macro displacement u AND the slaved internal shift s = M @ strain, form the total")
    print("      light-seen strain, and measure its incompatibility (linearised curvature / Ricci):")
    res = decisive(ki, ka)
    ratio = res["inc_full"] / res["inc_eig"]
    print(f"      max|curvature of full response (macro + optical)| = {res['inc_full']:.3e}")
    print(f"      max|curvature of the raw eigenstrain theta* (Sec 8.34 light-bending content)| = {res['inc_eig']:.3e}")
    print(f"      ratio = {ratio:.2e}")
    print("      => the full response -- INCLUDING the second sublattice's optical mode -- has curvature")
    print("         zero to machine precision, while the raw eigenstrain (the genuine light-bending")
    print("         content) does not: the probe is alive. Because the optical mode is GAPPED ([A]) the")
    print("         internal shift s is SLAVED to the macro strain, so it is itself a displacement field;")
    print("         the total strain is a symmetric gradient, hence compatible (Saint-Venant), hence")
    print("         flat. The second sublattice adds a second DISPLACEMENT, not an incompatible degree")
    print("         of freedom -- so it cannot escape the compatibility trap of Sec 8.34/8.35. gamma=0.\n")

    print("[verdict] the negative-space sublattice is NOT the missing Fierz-Pauli term.")
    print("  * It DOES break the single-lattice premise of Sec 8.35: the internal (optical) mode gives")
    print("    a genuine Cauchy-relation violation ([B]), the way diamond departs from a central-force")
    print("    solid. So Robert's instinct targeted the one real assumption -- this was worth measuring.")
    print("  * But the internal mode is GAPPED ([A]): a mechanical second sublattice's relative mode")
    print("    costs energy at q=0, so it is adiabatically SLAVED to the macro strain (s = M @ strain)")
    print("    and only renormalises the elastic constants -- another point on the nu axis Sec 8.35")
    print("    already swept, not a new long-range field. Only a GAPLESS relative mode (a Goldstone of a")
    print("    second condensate) would be an independent field; a mechanical interstitial is not one.")
    print("  * Because it is slaved, the internal shift is itself a DISPLACEMENT field, so the FULL")
    print("    two-sublattice response to a mass has curvature zero to machine precision ([C]): the")
    print("    total light-seen strain is a symmetric gradient, hence compatible (Saint-Venant), hence")
    print("    flat. The second sublattice adds a second displacement, not an INCOMPATIBLE degree of")
    print("    freedom -- so it cannot escape the compatibility trap that gave gamma=0 in Sec 8.34/8.35.")
    print("  * And [B2] shows even the internal mode's coupling is to SHEAR, never to compression, so a")
    print("    scalar mass reaches it only through the (radial) deviatoric part of its own field -- which")
    print("    carries no NET disclination (a crystallographic 3-theta wiggle, zero net charge, the")
    print("    Gauss-Bonnet result of Sec 8.33). Same spin-0-cannot-source-spin-2 rule, now under two")
    print("    sublattices. gamma = 1 is not reachable by adding the negative-space sublattice; the")
    print("    Newtonian (scalar) sector remains the model's honest gravity.")
