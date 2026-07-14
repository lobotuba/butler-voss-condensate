"""
Can the medium SHIELD gravity?  (The unneutralizability principle.)

Observationally, gravity cannot be shielded: a slab of lead between you and the Sun changes
nothing, and gravitational waves cross the universe unattenuated. The usual explanation --
"gravity is lossless" -- is NOT the mechanism: screening is not dissipation. A superconductor
screens B over the London depth with exactly ZERO loss, and a Yukawa field e^{-r/lambda}
conserves energy exactly. Screening is EVANESCENT, not absorptive. So losslessness cannot be
what protects gravity.

The real reason is that gravity is UNNEUTRALIZABLE. Every screening mechanism in physics --
Debye, the Faraday cage, the Meissner effect -- works by the medium rearranging charge of the
OPPOSITE sign into a cancelling cloud. Mass is unipolar: there is no negative mass to build the
cloud out of. You cannot screen what you cannot cancel.

That converts into a sharp, box-independent test. Put a source at the origin, let the medium do
ANYTHING it likes in a surrounding shell, and ask: HOW MUCH CHARGE IS STILL VISIBLE OUTSIDE?
Two couplings, same question:

  A. DILATATION (gravity-by-density: the source is energy density). Non-topological.
  B. CURVATURE  (Route 1: the source is a topological defect / winding).  Topological.

The measurement is a Gauss-law integral on a contour OUTSIDE the shell -- no range fits, no
finite-box exponents, no IR fragility (the failure modes that have bitten this project before).
"""
from __future__ import annotations
import numpy as np

TWO_PI = 2 * np.pi


# ---------------------------------------------------------------- Part A: dilatation
def dilatation_outside(N=128, lam=1.0, mu=1.0, eps=1.0, src_sig=3.0, r_out=20.0):
    """Isotropic eigenstrain source s(x) in a 2D isotropic elastic medium. Solve elastic
    equilibrium in Fourier and measure the DENSITY perturbation (div u) outside the source.

    Equilibrium with eigenstrain e*_ij = eps*s*delta_ij:
        -[(lam+mu) k_i k_j + mu k^2 delta_ij] u_j = 2(lam+mu) eps i k_i s(k)
    which inverts (M^-1 = [I - a k^k^]/(mu k^2), a=(lam+mu)/(lam+2mu)) to the exact result
        div u (x) = [2(lam+mu)/(lam+2mu)] eps s(x)     -- STRICTLY LOCAL.
    The medium's response cancels the source's density perturbation EXACTLY outside its
    support: a density-coupled force sees NOTHING at range. That is Bitter-Crum screening.
    """
    x = np.arange(N) - N / 2
    X, Y = np.meshgrid(x, x, indexing="ij")
    R = np.hypot(X, Y)
    s = np.exp(-R ** 2 / (2 * src_sig ** 2))                     # the "mass": a smooth blob

    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0                                               # k=0 handled separately
    sk = np.fft.fft2(s)

    a = (lam + mu) / (lam + 2 * mu)
    # u = -M^-1 f, f_i = 2(lam+mu) eps i k_i s  ->  u_i = -2(lam+mu)eps/(lam+2mu) * i k_i s / k^2
    pref = -2 * (lam + mu) * eps / (lam + 2 * mu)
    ux = np.fft.ifft2(pref * 1j * KX * sk / K2).real
    uy = np.fft.ifft2(pref * 1j * KY * sk / K2).real

    # density perturbation delta_rho = -rho0 * div u ; compute div u SPECTRALLY (no FD bias)
    divu = np.fft.ifft2(1j * KX * np.fft.fft2(ux) + 1j * KY * np.fft.fft2(uy)).real

    # The k=0 mode of u is identically zero (u ~ i k s / k^2), so the periodic box drops the
    # uniform part of div u -- the usual neutralizing-background (jellium) zero mode. It is a
    # constant offset, NOT a long-range tail; remove it from both sides before comparing.
    C = 2 * (lam + mu) * eps / (lam + 2 * mu)
    divu -= divu.mean()
    pred = C * (s - s.mean())

    # correctness gate: div u must equal C*s(x) EXACTLY (strictly local, no long-range part)
    resid = np.abs(divu - pred).max() / np.abs(pred).max()

    # The gate above establishes div u == C*s(x) identically. So the PHYSICAL far field is just
    # C*s(x) -- i.e. the source's own tail, with ZERO long-range contribution from the medium.
    # (Measuring divu directly outside would report the mean-subtraction pedestal, a periodic-box
    # zero-mode artifact, not physics -- so read the far field off the verified local relation.)
    inside = R < 1.5 * src_sig
    outside = R > r_out
    peak = C * float(np.abs(s[inside]).max())
    leak = C * float(np.abs(s[outside]).max())
    return peak, leak, resid, C


# ---------------------------------------------------------------- Part B: topological
def wrap(d):
    return (d + np.pi) % TWO_PI - np.pi


def winding(theta, c, R):
    """EXACT lattice winding on a square contour of half-width R about site c.
    Sum of wrapped nearest-neighbour phase differences / 2pi -> an integer, to machine
    precision, for ANY single-valued phase field. This is the Gauss law for topology."""
    pts = []
    for x in range(c - R, c + R):       pts.append((x, c - R))
    for y in range(c - R, c + R):       pts.append((c + R, y))
    for x in range(c + R, c - R, -1):   pts.append((x, c + R))
    for y in range(c + R, c - R, -1):   pts.append((c - R, y))
    tot = 0.0
    for i in range(len(pts)):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % len(pts)]
        tot += wrap(theta[x1, y1] - theta[x0, y0])
    return tot / TWO_PI


def defect(N, cx, cy, q=1):
    x = np.arange(N)[:, None]
    y = np.arange(N)[None, :]
    return q * np.arctan2(y - cy, x - cx)


def smooth_field(N, rng, sigma=5.0):
    """A single-valued smooth random field: zero winding BY CONSTRUCTION. This is the most
    general thing the medium can do WITHOUT nucleating a defect."""
    F = np.fft.fft2(rng.standard_normal((N, N)))
    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    F *= np.exp(-0.5 * sigma ** 2 * (KX ** 2 + KY ** 2))
    g = np.fft.ifft2(F).real
    return g / (np.abs(g).max() + 1e-30)


def relax(theta, J, fixed, iters=400):
    """Relax the medium (stiffness field J) by local energy minimisation, holding the pinned
    sites. The medium is free to do whatever it wants everywhere else -- including inside a
    shell of wildly different stiffness. Can it cancel the charge seen outside?"""
    th = theta.copy()
    for _ in range(iters):
        s = np.zeros_like(th)
        c = np.zeros_like(th)
        for ax in (0, 1):
            for sh in (1, -1):
                w = 0.5 * (J + np.roll(J, sh, axis=ax))
                s += w * np.sin(np.roll(th, sh, axis=ax))
                c += w * np.cos(np.roll(th, sh, axis=ax))
        th = np.where(fixed, th, np.arctan2(s, c))
    return th


if __name__ == "__main__":
    print("=== Can the medium SHIELD gravity? The unneutralizability principle ===\n")
    print("  Screening is NOT loss (a superconductor screens B with zero dissipation).")
    print("  Screening is NEUTRALIZATION. So: what charge is still visible OUTSIDE a shell")
    print("  in which the medium has done everything it can to cancel the source?\n")

    # ---------------- A. dilatation (gravity-by-density) ----------------
    print("  [A] DILATATION source (gravity-by-density: couples to ENERGY DENSITY)")
    peak, leak, resid, C = dilatation_outside()
    print(f"      gate: div u == C*s(x) with C = {C:.4f}, residual {resid:.2e} (strictly LOCAL)")
    print(f"      medium's density response inside the source : {peak:.4e}")
    print(f"      density perturbation visible OUTSIDE (r>20) : {leak:.4e}")
    print(f"      => leakage / source = {leak/peak:.2e}  -- the medium cancels it to nothing.")
    print("      (the far field is only the source's own Gaussian tail; the medium's RESPONSE")
    print("      adds no long-range part at all. Bitter-Crum: div u is strictly local.)")
    print("      A density-coupled force sees NOTHING at range. This charge is NEUTRALIZABLE,")
    print("      hence SHIELDABLE, hence short-range.\n")

    # ---------------- B. topological curvature (Route 1) ----------------
    N, c = 128, 64
    rng = np.random.default_rng(0)
    xg = np.arange(N)[:, None]; yg = np.arange(N)[None, :]
    Rg = np.hypot(xg - c, yg - c)

    core = Rg < 2.5                      # the SOURCE ("the mass"), pinned
    edge = (Rg > 44)                     # fixed far-field boundary condition
    shell = (Rg > 8) & (Rg < 15)         # the SHIELD: a shell of intervening matter
    fixed = core | edge
    R_probe = 22                         # contour OUTSIDE the shell

    th0 = defect(N, c, c, q=1)
    base = winding(th0, c, R_probe)
    print("  [B] TOPOLOGICAL CURVATURE source (Route 1: couples to a defect / winding)")
    print(f"      bare source, charge on contour r={R_probe}: {base:+.9f}")
    print("      now let the medium try to shield it -- everything short of making a defect:\n")

    trials = []
    # (i) a violent smooth deformation of the medium in the shell
    th = th0 + 4.0 * shell * smooth_field(N, rng)
    trials.append(("violent smooth deformation in the shell", winding(th, c, R_probe)))

    # (ii) a shell of hugely different stiffness, medium fully relaxed
    J = np.ones((N, N)); J[shell] = 50.0
    th = relax(th0.copy(), J, fixed)
    trials.append(("shell 50x stiffer + full relaxation", winding(th, c, R_probe)))

    # (iii) a nearly-fluid shell (stiffness -> 0), medium fully relaxed
    J = np.ones((N, N)); J[shell] = 0.02
    th = relax(th0.copy(), J, fixed)
    trials.append(("shell 50x softer + full relaxation", winding(th, c, R_probe)))

    # (iv) everything at once
    J = np.ones((N, N)); J[shell] = 50.0
    th = relax(th0 + 4.0 * shell * smooth_field(N, rng), J, fixed)
    trials.append(("both: deformed AND stiffness-contrasted", winding(th, c, R_probe)))

    for name, w in trials:
        print(f"      {name:<42} -> {w:+.9f}")
    print("      => the charge visible outside is EXACTLY unchanged. The medium cannot")
    print("         neutralize it: any smooth response has ZERO winding, identically.\n")

    # ---------------- C. positive control: the ONLY way to shield ----------------
    print("  [C] CONTROL -- the only thing that WOULD shield it: a genuine ANTI-charge")
    th = th0 + defect(N, c + 11, c, q=-1)        # an anti-defect inside the shell
    w = winding(th, c, R_probe)
    print(f"      nucleate an anti-defect (q=-1) in the shell -> {w:+.9f}")
    print("      The probe IS sensitive -- so the invariance above is physics, not blindness.")
    print("      But an anti-defect is a QUANTIZED topological charge: you cannot make half of")
    print("      one, and in the ordered medium it costs a finite core energy to make one at all.")
    print("      There is no infinitesimal, continuously-polarizable version of this charge --")
    print("      hence no linear-response screening cloud, hence NO shielding.\n")

    print("[verdict] the two couplings sit on opposite sides of the shielding question:")
    print("  * DILATATION (energy density) is NEUTRALIZABLE: the medium cancels it exactly, the")
    print("    far field is zero, gravity-by-density is SHIELDABLE -- and therefore screened and")
    print("    short-range. This is the model's gravity problem, restated as its true cause.")
    print("  * CURVATURE (topological) is UNNEUTRALIZABLE: no smooth rearrangement of any amount")
    print("    of intervening matter can change the charge seen outside, because a single-valued")
    print("    medium response carries zero winding. Shielding would require a real anti-defect --")
    print("    a NEGATIVE MASS -- which is quantized and unavailable.")
    print("\n  => Topological quantization does in the model exactly what 'there is no negative mass'")
    print("     does in nature. Gravity is unshieldable not because it is lossless (screening costs")
    print("     nothing anyway) but because it is UNCANCELLABLE. That is the real argument for")
    print("     Route 1: couple gravity to curvature, and it inherits nature's unshieldability.")
