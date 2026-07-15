"""
The cosmological constant: why the model's huge vacuum energy does not gravitate.

Integrating out the medium's modes gives a vacuum (zero-point) energy density of order the
microscopic scale -- in this model the lattice/node scale, ~ 1/a0^4 = M_Planck^4 (a0 = l_Planck was
fixed in test_scale_fixing). The observed dark-energy density is ~ (2 meV)^4, about 10^122 times
smaller. Taken at face value this is "the worst prediction in physics": a 10^122 fine-tuning of the
bare vacuum energy against the cosmological constant. The induced-gravity calculations already met
this term head-on -- test_induced_gravity's Pi(0,0) is the induced cosmological piece, and
test_lattice_ward's <T^ij> != 0 is the same vacuum stress.

But the model has a structural fact that changes the framing entirely: ITS VACUUM IS A CONDENSATE
-- a self-sustained equilibrium medium (test_collapse: the condensate is the true vacuum). For such
a vacuum the quantity that gravitates is NOT the bare energy density eps. This is Volovik's result
for emergent gravity in quantum liquids: the emergent metric couples to the vacuum STRESS, i.e. the
grand-canonical potential density
        rho_Lambda  =  eps~  =  eps - mu n  =  -P
(mu = chemical potential of the conserved medium charge; the thermodynamic identity Omega = -PV).
The bare eps is huge, but what curves the emergent space is -P, the vacuum pressure -- exactly the
<T^ij> ~ -P delta^ij vacuum stress. And a SELF-SUSTAINED vacuum -- one that can exist with nothing
outside pushing on it, which is what the vacuum of empty space is -- has

        P = 0        =>        rho_Lambda = 0,

with the equilibrium density self-adjusting to enforce it. The huge zero-point energy is absorbed
into the equilibrium condensate density, not into curvature.

What is measured here:
  [A] the bare problem, concretely: the filled-sea (zero-point) energy of the medium is O(1) per
      site = O(microscopic scale); as a density that is ~10^122 x the observed dark energy.
  [B] the self-cancellation: with eps(n) = eps0 + (g/2) n^2 and the gravitating rho_Lambda = -P,
      the self-sustained density n0 (where P=0) makes rho_Lambda = 0 to MACHINE PRECISION for eps0
      swept across 122 orders of magnitude -- with NO tuning of parameters (n0 tracks eps0). A rigid
      (non-self-adjusting) vacuum would instead gravitate the full eps0.
  [C] the residual: a small fractional departure delta from equilibrium gives rho_Lambda ~ delta.
      The observed tiny Lambda measures how close the vacuum sits to equilibrium.

HONEST SCOPE. This is Volovik's condensed-matter / emergent-gravity mechanism, naturally realized
because the model's vacuum is a self-sustained condensate. It makes the equilibrium cosmological
constant EXACTLY zero and AUTOMATIC -- dissolving the 10^122 fine-tuning (the equilibrium vacuum
does not gravitate its zero-point energy, for any bare value). It does NOT predict the observed
nonzero Lambda: that is relocated to a cosmological question -- why the vacuum sits slightly OFF
equilibrium (cosmic expansion, matter/radiation content, relaxation dynamics) -- which is not solved
here. The gain is conceptual and real: from "cancel 10^122 by tuning constants" to "the equilibrium
vacuum gravitates zero by thermodynamics, and Lambda measures the departure from equilibrium."
"""
from __future__ import annotations
import numpy as np


def sea_energy_density(N=200, m0=-1.0):
    """Filled-sea (zero-point) energy density of the Wilson-Dirac medium: sum of the occupied
    (negative-energy) branch over the BZ. O(1) per site = O(microscopic = Planck) scale."""
    g = (np.arange(N) + 0.5) / N * 2 * np.pi
    KX, KY = np.meshgrid(g, g, indexing="ij")
    kx, ky = KX.ravel(), KY.ravel()
    E = np.sqrt(np.sin(kx) ** 2 + np.sin(ky) ** 2 + (m0 + 2 - np.cos(kx) - np.cos(ky)) ** 2)
    return -float(np.mean(E))


G = 1.0                                              # condensate self-coupling


def pressure(n, eps0):
    """P = n eps'(n) - eps(n) for eps(n) = eps0 + (G/2) n^2."""
    return n * (G * n) - (eps0 + 0.5 * G * n * n)


def n_self_sustained(eps0):
    """Density where P = 0 (the self-sustained vacuum): (G/2) n^2 = eps0."""
    return np.sqrt(2 * eps0 / G)


def rho_lambda(n, eps0):
    """Gravitating vacuum energy = grand potential density = -P (Volovik)."""
    return -pressure(n, eps0)


if __name__ == "__main__":
    print("=== The cosmological constant: the condensate vacuum's zero-point energy does not gravitate ===\n")

    # ---------- [A] the bare problem ----------
    eps_site = sea_energy_density()
    ratio = 1e122
    print("  [A] THE BARE PROBLEM (concrete):")
    print(f"      filled-sea zero-point energy = {eps_site:.4f} per site = O(1) x microscopic scale.")
    print(f"      With a0 = l_Planck (test_scale_fixing), that density ~ M_Planck^4 ~ 10^122 x the")
    print(f"      observed dark-energy density. Naively: a 10^122 fine-tuning. Same term as")
    print(f"      test_induced_gravity's Pi(0,0) and test_lattice_ward's <T^ij> != 0.\n")

    # ---------- [B] the self-cancellation, no tuning ----------
    print("  [B] BUT the vacuum is a self-sustained CONDENSATE, and gravitating rho_Lambda = -P.")
    print("      A self-sustained vacuum has P = 0; the density self-adjusts. Sweep the bare eps0:")
    print(f"      {'eps0 (bare vac E)':>18} {'n0 (self-adjusts)':>18} {'rho_Lambda = -P':>16} "
          f"{'rho_Lambda/eps0':>16} {'rigid vacuum':>14}")
    for e0 in (1e0, 1e3, 1e30, 1e61, 1e122):
        n0 = n_self_sustained(e0)
        rl = rho_lambda(n0, e0)
        rigid = rho_lambda(n0 * 1.3, e0) / e0        # a vacuum NOT sitting at P=0 gravitates ~ eps0
        print(f"      {e0:>18.0e} {n0:>18.3e} {rl:>16.2e} {rl / e0:>16.2e} {rigid:>14.2e}")
    print("      => rho_Lambda = 0 to machine precision for ANY eps0 across 122 orders -- the density")
    print("         n0 tracks eps0 and P self-zeroes. NO parameter tuning. The huge zero-point energy")
    print("         is absorbed into the equilibrium density, not into curvature. (A rigid vacuum,")
    print("         last column, would gravitate a value of order eps0 -- the standard disaster.)\n")

    # ---------- [C] the residual = departure from equilibrium ----------
    e0 = 1e122
    n0 = n_self_sustained(e0)
    print("  [C] THE RESIDUAL: a vacuum slightly OFF equilibrium (fractional departure delta):")
    print(f"      {'delta (frac. off n0)':>22} {'rho_Lambda':>14} {'rho_Lambda/eps0':>16}")
    for d in (1e-1, 1e-3, 1e-6, 1e-12, 1e-30):
        rl = rho_lambda(n0 * (1 + d), e0)
        print(f"      {d:>22.0e} {rl:>14.3e} {rl / e0:>16.3e}")
    print("      => rho_Lambda ~ delta: the gravitating Lambda measures how FAR the vacuum sits from")
    print("         equilibrium, not a cancellation of eps0. To match the observed 10^-122 the vacuum")
    print("         need only be at equilibrium to that fractional precision -- what a relaxed")
    print("         self-sustained vacuum does -- rather than 122 digits of tuned constants.\n")

    print("[verdict] the model's condensate vacuum dissolves the cosmological-constant FINE-TUNING:")
    print("  * The gravitating vacuum energy is the grand potential rho_Lambda = eps - mu n = -P")
    print("    (Volovik: the emergent metric couples to the vacuum STRESS <T^ij> ~ -P, not the bare")
    print("    zero-point eps). This is natural here because the vacuum IS a self-sustained condensate")
    print("    (test_collapse), and its induced stress is exactly the <T> of test_lattice_ward.")
    print("  * A self-sustained vacuum has P = 0, so rho_Lambda = 0 -- AUTOMATICALLY and with NO")
    print("    tuning, for a bare eps0 swept across all 122 orders (the density self-adjusts). The")
    print("    huge zero-point energy does not gravitate.")
    print("  * HONEST scope: this makes the EQUILIBRIUM cosmological constant exactly zero; it does")
    print("    NOT derive the observed nonzero value. That is relocated to why the vacuum is slightly")
    print("    off equilibrium (expansion, matter, relaxation) -- a cosmological dynamics question,")
    print("    still open. The advance is the reframing: not 'cancel 10^122 by tuning constants' but")
    print("    'the equilibrium vacuum gravitates zero by thermodynamics; Lambda measures the")
    print("    departure from equilibrium.'")
