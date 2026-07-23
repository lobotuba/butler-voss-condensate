"""
Is the coupling Einstein? gamma = 1 needs mass to source CURVATURE, not displacement.

*** STATUS UPDATE -- the "one calculation left" named at the end of this file has now been done in the
    SMOOTH channel, and it comes back zero. This file reduces gamma to a single number, the induced
    Einstein source coupling kappa (gamma = kappa/(4 pi G)), and leaves whether the fermion loop
    supplies it as the one open calculation. test_gamma_source measures it: the induced coupling of a
    static energy density to the spatial stress, <T00, T_ij>(q), vanishes identically -- for every
    component including the trace, to machine precision, across mass, cutoff and momentum -- while a
    genuine spin-2 source couples to the same stress at O(1), so the zero is a selection rule, not a
    numb instrument. The energy density is a scalar (T00 ~ E*I) and sources no spin-2 curvature. So in
    the smooth channel kappa = 0 and gamma = 0, from the loop side, exactly as this file found it = 0
    from the ELASTIC side ([A] here: a mass relaxes to a compatible displacement, eta/rho ~ 1e-15).
    Both smooth mechanisms agree. What is NOT retracted, and is the reason gamma = 0 is not yet the
    final word: real curvature here is INCOMPATIBLE strain -- disclination density -- and a smooth
    stress-stress bubble cannot see the topological channel. Whether a mass sources curvature through
    THAT channel is the remaining measurement. gamma = 1 now rests entirely on the emergent-Weinberg
    argument (this file's route C), which has no direct in-model confirmation and now has direct
    evidence against it in the smooth sector. ***

*** STATUS UPDATE: the missing ingredient this file isolated has since been SUPPLIED. The analysis
    below stands -- gamma is set by ONE number, the strength with which mass sources curvature
    relative to the Newtonian potential, and the medium's ELASTIC response supplies none of it
    (eta/rho ~ 1e-15), so the native scalar sector sits at gamma = 0. What has changed is the
    closing status "the medium currently sits at gamma = 0 / the curvature sector is unsourced":
      * the curvature sector is confining on its own but DECONFINES into a massless Newtonian
        graviton given a positive induced Einstein term (test_deconfinement),
      * that term's SIGN was then measured positive, mu > 0, by calibrating the induced Newtonian
        coupling against the model's own healthy photon (test_induced_sign),
      * the resulting spin-2 graviton is dynamical, doubly degenerate and healthy in 3+1D
        (test_spin2_dynamical), and
      * gamma = 1 follows from Weinberg's theorem on the conserved infrared stress tensor, as an
        EMERGENT (not lattice-exact) identity -- test_lattice_ward explains why it can only be
        emergent: diffeomorphism invariance, unlike U(1), is not a lattice symmetry.
    So the "one number" is no longer missing; it is induced. ***

test_incompatible_gravity found the door: the medium has an incompatible-strain (curvature)
sector, distinct from the gauge displacement sector, and it DOES bend light. But a door is not a
crossing. Two things decide whether the medium actually does GR (gamma = 1, the light-bending
factor of two):
    (1) does ordinary MASS source the curvature sector at all?  and
    (2) if so, with what STRENGTH relative to the Newtonian (time) potential?
gamma = Psi/Phi is exactly that strength ratio.

There is a sharp obstruction to state first, because it explains the model's current gamma = 0.
Our WORKING gravity (test_critical_gravity) is the SCALAR amplitude mode: a genuine spin-0 field.
A scalar sources only the time part Phi, so on its own it gives gamma = 0 (Nordstrom gravity).
The spin-2 curvature sector is a SEPARATE mediator, and the question is whether mass couples to
it as a massless graviton coupled to the conserved stress tensor -- Weinberg's unique route to
gamma = 1.

This file measures the three things that decide it:

  A. Does mass source curvature ELASTICALLY? Put a mass in the medium, solve the elastic
     equilibrium, decompose the strain into compatible (displacement) + incompatible (curvature).
     Equilibrium elasticity returns a DISPLACEMENT, so its incompatible part is ~0: mass sources
     NO curvature elastically. That is exactly why gamma = 0 -- the obstruction, made concrete.

  B. What DOES source curvature? A disclination (a topological charge). Only a source that is not
     a body force curves the medium -- and that sector repels/confines (test_disclination_force).
     So the "mass -> curvature" coupling is ABSENT in pure elasticity and must be INDUCED.

  C. The target. IF an induced coupling makes mass source curvature as eta = kappa * rho, solve
     the full light metric (time Phi + space Psi) and measure gamma and the deflection. gamma is
     the ratio kappa/(Newtonian coupling); the Einstein value gives gamma = 1 and bending x2.
     This isolates the ONE missing number: the Einstein source strength.

Conclusion this is driving at: the graviton's home (incompatible sector) and kinetic term
(Sakharov, test_induced_gravity) are in hand; the missing ingredient is specifically the
EINSTEIN SOURCE COUPLING -- mass sourcing curvature, the right-hand side of the field equation --
which pure elasticity lacks and the fermion loop must supply. That coupling IS 'emergent
diffeomorphism invariance', now isolated to a single testable number.
"""
from __future__ import annotations
import numpy as np


def ops(N):
    k1 = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    return KX, KY


def dd(f, KA, KB):
    return np.fft.ifft2(-(KA * KB) * np.fft.fft2(f)).real


def incompat(exx, eyy, exy, KX, KY):
    return dd(exx, KY, KY) + dd(eyy, KX, KX) - 2 * dd(exy, KX, KY)


def elastic_response(rho, KX, KY, lam=1.0, mu=1.0):
    """Elastic equilibrium: a mass (centre of dilatation) exerts forces; the medium relaxes to a
    DISPLACEMENT u. Return the strain of that displacement."""
    K2 = KX ** 2 + KY ** 2; K2[0, 0] = 1.0
    sk = np.fft.fft2(rho); pref = -2 * (lam + mu) / (lam + 2 * mu)
    ux = np.fft.ifft2(pref * 1j * KX * sk / K2).real
    uy = np.fft.ifft2(pref * 1j * KY * sk / K2).real
    exx = np.fft.ifft2(1j * KX * np.fft.fft2(ux)).real
    eyy = np.fft.ifft2(1j * KY * np.fft.fft2(uy)).real
    exy = 0.5 * (np.fft.ifft2(1j * KY * np.fft.fft2(ux)).real +
                 np.fft.ifft2(1j * KX * np.fft.fft2(uy)).real)
    return exx, eyy, exy


def poisson(rho, KX, KY):
    K2 = KX ** 2 + KY ** 2; K2[0, 0] = 1.0
    out = np.fft.ifft2(np.fft.fft2(rho) / K2).real
    return out - out.mean()


def deflect(seen, N, bs):
    c = N // 2
    gy = 0.5 * (np.roll(seen, -1, 1) - np.roll(seen, 1, 1))
    return np.array([-gy[:, c + b].sum() for b in bs])


if __name__ == "__main__":
    print("=== Is the coupling Einstein? gamma = 1 needs mass to source CURVATURE ===\n")

    N = 256
    KX, KY = ops(N)
    g = np.arange(N) - N // 2
    X, Y = np.meshgrid(g, g, indexing="ij")
    R = np.hypot(X, Y)
    rho = np.exp(-R ** 2 / (2 * 3.0 ** 2)); rho -= rho.mean()
    bs = np.array([12, 16, 20, 26, 32, 40])

    # ---- A. does mass source curvature elastically? ----
    exx, eyy, exy = elastic_response(rho, KX, KY)
    eta_from_mass = incompat(exx, eyy, exy, KX, KY)
    print("  [A] mass in the elastic medium -> does it source CURVATURE?")
    print(f"      incompatibility eta from a mass: max|eta|/|rho| = "
          f"{np.abs(eta_from_mass).max()/np.abs(rho).max():.2e}")
    print("      ~ 0: equilibrium elasticity returns a DISPLACEMENT, whose strain is compatible.")
    print("      MASS SOURCES NO CURVATURE elastically. This is exactly why gamma = 0.\n")

    # ---- B. what sources curvature? a disclination (topological) ----
    #   a disclination is NOT a body force: it imposes eta directly. Model its curvature as a blob.
    eta_disc = np.exp(-R ** 2 / (2 * 4.0 ** 2)); eta_disc -= eta_disc.mean()
    print("  [B] a DISCLINATION (topological charge) sources curvature directly: eta != 0.")
    print("      Only a non-force source curves the medium -- and that sector repels/confines")
    print("      (test_disclination_force). So 'mass -> curvature' is ABSENT in elasticity; it")
    print("      must be INDUCED (the fermion loop / Sakharov).\n")

    # ---- C. the target: if mass DID source curvature, what gamma? ----
    G = 1.0
    Phi = 4 * np.pi * G * poisson(rho, KX, KY)          # time potential: lap Phi = 4 pi G rho
    a_time = deflect(Phi, N, bs)
    print("  [C] IF an induced coupling makes mass source curvature eta = kappa * rho, then the")
    print("      spatial potential solves lap(Psi) = kappa rho, and gamma = Psi/Phi = kappa/(4 pi G).")
    print(f"      {'kappa/(4 pi G)':>15} {'gamma':>8} {'bend/Newton':>12} {'interpretation':>22}")
    for frac, name in [(0.0, "scalar (amplitude mode)"),
                       (0.5, "scalar-tensor"),
                       (1.0, "EINSTEIN (GR)"),
                       (2.0, "over-curved")]:
        Psi = frac * 4 * np.pi * G * poisson(rho, KX, KY)
        a_both = deflect(Phi + Psi, N, bs)
        gamma = np.median(a_both / a_time) - 1.0
        print(f"      {frac:>15.1f} {gamma:>8.2f} {np.median(a_both/a_time):>12.2f}"
              f"   {name:>22}")
    print("\n      => gamma is set by ONE number: how strongly mass sources curvature relative to")
    print("         the Newtonian potential. The Einstein value (kappa = 4 pi G) gives gamma = 1 and")
    print("         the light-bending factor of two. The medium currently sits at gamma = 0 (its")
    print("         gravity is the scalar amplitude mode; the curvature sector is unsourced).\n")

    print("[verdict] the crossing is isolated to a single missing coupling:")
    print("  * Mass sources NO curvature elastically (A): the elastic response to any force is a")
    print("    displacement (compatible, gauge). So the medium's native gravity -- the scalar")
    print("    amplitude mode -- is gamma = 0 (Nordstrom). This is why light does not bend x2.")
    print("  * Curvature is sourced only by TOPOLOGICAL charge (B), which repels. So 'mass sources")
    print("    curvature' -- the right-hand side of Einstein's equation -- is simply ABSENT from")
    print("    pure elasticity and must be INDUCED.")
    print("  * gamma is then one number (C): the induced source strength. Einstein = gamma = 1.")
    print("\n  The graviton's HOME (incompatible sector, test_incompatible_gravity) and its KINETIC")
    print("  TERM (Sakharov, test_induced_gravity) are in hand. The single remaining ingredient is")
    print("  the EINSTEIN SOURCE COUPLING: mass sourcing curvature at the strength gamma = 1. Whether")
    print("  the fermion loop induces exactly that -- coupling the metric to the CONSERVED stress")
    print("  tensor, whose conservation forces the Einstein strength (Weinberg) -- is the precise,")
    print("  and now singular, statement of the spin-2 wall. That is the one calculation left.")
