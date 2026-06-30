"""
H9 -- Sub-quantum hierarchy: do base particles bind into a higher layer?
========================================================================

Uses the complex-field prototype (prototype_complex.py).  Tests whether the
fundamental n=+/-1 vortices "confluence" into composites:

  A. STATIC interaction potential V(d): energy of a seeded pair vs separation d.
       opposite charges (+/-)  -> energy DROPS as they approach  = attraction = binding
       like charges    (+/+)  -> energy RISES  as they approach  = repulsion  = no bond
  B. DYNAMIC test: release a pair at rest and watch the separation d(t).
       attraction -> d shrinks; repulsion -> d grows.
  C. COMPOSITE vs FUNDAMENTAL: is n=2 one particle or two n=1?  Seed two like
       vortices close together and measure whether they fly apart (=> n=2 is a
       composite/unstable, so the true fundamentals are n=0, +/-1 + Noether).
  D. SCALE GAP: core size & energy of a single vortex vs the bound (+/-) molecule
       -> a first estimate of the base-layer / next-layer separation.

A +/- pair has zero net winding, so its far field cancels and total energy is
box-independent -- which makes its binding curve the clean, quantitative one.
"""
from __future__ import annotations
import numpy as np
from prototype_complex import ComplexFabric


# ----------------------------------------------------------- helpers ---------
def plaq_centers(f):
    return f.pos[f.plaq].mean(axis=1)                      # (P,2)


def vortex_spread(f, sign):
    """charge-weighted centroid and RMS spread of plaquettes of a given sign."""
    w = f.topological_charge(per_plaquette=True)
    sel = (np.sign(w) == sign) & (w != 0)
    if not sel.any():
        return None, 0.0, 0
    cen = plaq_centers(f)[sel]
    wt = np.abs(w[sel]).astype(float)
    centroid = (cen * wt[:, None]).sum(0) / wt.sum()
    spread = np.sqrt(((cen - centroid) ** 2).sum(1) @ wt / wt.sum())
    return centroid, float(spread), int(np.abs(w[sel]).sum())


def separation(f):
    cp, _, qp = vortex_spread(f, +1)
    cm, _, qm = vortex_spread(f, -1)
    if cp is None or cm is None:
        return None
    return float(np.linalg.norm(cp - cm))


def n_cores(f, amp_cut=0.4):
    """Robust count of real vortex cores: plaquettes that both carry a winding
    AND sit in a genuine amplitude dip (|psi| < amp_cut*v0), which rejects the
    spurious integer windings thrown off by radiation in the bulk."""
    w = f.topological_charge(per_plaquette=True)
    amp = (np.abs(f.psi)[f.plaq]).min(axis=1) / f.v0      # min |psi| on each plaquette
    return int(np.abs(w[(w != 0) & (amp < amp_cut)]).sum())


def core_radius(f, center_frac=(0.5, 0.5), recover=0.5):
    """radius at which |psi| recovers to `recover`*v0 around a core."""
    c = f._center(*center_frac)
    r = np.linalg.norm(f.pos - c, axis=1)
    amp = np.abs(f.psi) / f.v0
    order = np.argsort(r)
    for i in order:
        if amp[i] >= recover:
            return float(r[i])
    return float(r.max())


def seed_pair(charges, s, core=3.0, **kw):
    """+/- or +/+ pair centered, half-separation s (fraction of box)."""
    f = ComplexFabric(potential="mexicanhat", **kw)
    f.set_vortices([(0.5 - s, 0.5, charges[0]),
                    (0.5 + s, 0.5, charges[1])], core=core)
    return f


# ----------------------------------------------------------- A: static -------
def static_potential():
    print("A. Static interaction energy vs half-separation s")
    print(f"   {'s':>5} {'phys_sep':>9} {'E(+/-)':>9} {'E(+/+)':>9}")
    rows = []
    for s in [0.06, 0.10, 0.14, 0.18, 0.22, 0.28]:
        fpm = seed_pair((+1, -1), s); fpp = seed_pair((+1, +1), s)
        d = separation(fpm)
        epm, epp = fpm.energy(), fpp.energy()
        rows.append((s, d, epm, epp))
        print(f"   {s:>5.2f} {str(round(d,2) if d else None):>9} {epm:>9.3f} {epp:>9.3f}")
    e_close = rows[0][2]; e_far = rows[-1][2]
    print(f"   -> opposite charges: E rises with separation by "
          f"{e_far - e_close:+.3f} (E goes UP as they part => they ATTRACT/bind)")
    print(f"   -> binding energy (release from far->close) ~ {e_far - e_close:.3f}")
    de_pp = rows[0][3] - rows[-1][3]
    print(f"   -> like charges: E(close)-E(far) = {de_pp:+.3f} "
          f"(positive => costs energy to approach => REPEL)\n")


# ----------------------------------------------------------- B: dynamic ------
def dynamic_pair(charges, label, steps=1500, s=0.20, damping=0.999):
    """Release a pair at rest; track min|psi|/v0 (≈0 at a live core, ->1 once
    cores annihilate) and the spread of the like-sign cores."""
    f = seed_pair(charges, s, damping=damping)
    traj, every = [], max(1, steps // 8)
    for k in range(steps):
        f.step()
        if k % every == 0 or k == steps - 1:
            minamp = float(np.abs(f.psi).min() / f.v0)
            _, sp, _ = vortex_spread(f, int(np.sign(charges[0])))
            traj.append((round(f.time, 1), round(minamp, 2), round(sp, 1)))
    print(f"B. {label}: min|psi|/v0 (core depth) and like-core spread over time")
    print("   " + "  ".join(f"t{t}:amp={a},spr={s}" for t, a, s in traj))
    end_amp, end_spr, start_spr = traj[-1][1], traj[-1][2], traj[0][2]
    if end_spr > start_spr + 2.0:
        verdict = "cores SURVIVE and spread apart  =>  REPEL"
    elif end_amp > 0.6:
        verdict = "cores merge & ANNIHILATE (amp recovers to v0)  =>  ATTRACT / bind"
    else:
        verdict = "cores survive, ~static"
    print(f"   final core depth={end_amp}  spread {start_spr}->{end_spr}  =>  {verdict}\n")


# ----------------------------------------------- C: composite vs fundamental -
def composite_test(steps=1500):
    print("C. Is n=2 a fundamental or two n=1 stuck together?")
    f = seed_pair((+1, +1), s=0.05, damping=0.999)         # two like vortices, close
    _, sp0, q0 = vortex_spread(f, +1)
    spreads = []
    for k in range(steps):
        f.step()
        if k % 150 == 0 or k == steps - 1:
            _, sp, _ = vortex_spread(f, +1)
            spreads.append((round(f.time, 1), round(sp, 2)))
    _, spF, _ = vortex_spread(f, +1)
    print(f"   total +winding = {q0} (=2).  spread of the +charge over time:")
    print("   " + "  ".join(f"t{t}={s}" for t, s in spreads))
    grew = spF > sp0 + 0.5
    print(f"   spread {sp0:.2f} -> {spF:.2f}  => "
          f"{'SPLITS: n=2 is a COMPOSITE of two n=1' if grew else 'stays bound'}\n")


# ------------------------------------------------------------- D: scale gap --
def scale_gap():
    print("D. Scale gap (base particle vs bound composite)")
    single = ComplexFabric(potential="mexicanhat")
    single.set_vortices([(0.5, 0.5, +1)])
    rc = core_radius(single)
    e1 = single.energy()
    pair = seed_pair((+1, -1), s=0.10)                     # a bound +/- molecule
    dmol = separation(pair); emol = pair.energy()
    print(f"   single n=1 vortex : core radius ~ {rc:.2f},  energy ~ {e1:.2f}")
    print(f"   bound (+/-) pair  : size ~ {dmol:.2f},        energy ~ {emol:.2f}")
    if rc > 0:
        print(f"   size ratio (composite/base) ~ {dmol/rc:.1f}x   "
              f"energy ratio ~ {emol/e1:.2f}x")


if __name__ == "__main__":
    print("=== H9 :: do base particles bind into a higher layer? ===\n")
    static_potential()
    # opposite charges bind (dynamic); like-charge repulsion is shown by A + C
    dynamic_pair((+1, -1), "opposite charges (+/-)", steps=4000, s=0.10)
    composite_test()
    scale_gap()
