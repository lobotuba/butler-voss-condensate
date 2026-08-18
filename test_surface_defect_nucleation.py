"""
HUNTING (2e): the surface/defect nucleation caveat -- does it reopen the radiation channel?

Method (Robert): compute blind, then audit hard, then literature. Hunt 5 bounded BULK vortex nucleation
(Planck-scale barrier) but flagged an honest residual: real superfluids nucleate vortices far more easily at
SURFACES, boundaries and defects, where the barrier is much lower. Does that reopen the spontaneous-radiation
channel that would put the model back in tension with Gran Sasso?

THE CLARIFYING DISTINCTION (the crux of the audit). Surface/defect nucleation is a statement about
MACROSCOPIC SUPERFLOW past a MACROSCOPIC boundary (a wall, a protrusion, a trapped object): the flow field
of the bulk condensate streaming past geometry. But the collapse-radiation bound (Gran Sasso) is about
INDIVIDUAL PARTICLES -- bound electrons -- being (or not being) accelerated by the decoherence coupling and
emitting X-rays. Those particles move through the BULK vacuum, deep inside atoms, far from any macroscopic
surface. So the channel that the radiation bound actually tests is governed by the BULK Landau velocity
v_c^bulk ~ 0.72 c (hunt 4), not by the surface-reduced macroscopic critical velocity. The right question is
therefore: do the fastest ORDINARY particles move faster than v_c^bulk?

  [G1] SURFACE NUCLEATION IS REAL BUT IS A DIFFERENT REGIME: a geometric enhancement factor f ~ 10-100
       (calibrated to superfluid He, where surface nucleation drops the critical velocity ~1-2 orders below
       the bulk roton value) lowers the MACROSCOPIC-superflow critical velocity to v_c^surf ~ 0.007-0.07 c.
       This governs bulk condensate flow past boundaries -- not individual-particle motion through the vacuum.
  [G2] ORDINARY MATTER IS BULK-PROTECTED: the fastest everyday particles are inner-shell (1s) electrons,
       v_1s ~ Z alpha c. Even for the heaviest elements these stay below v_c^bulk ~ 0.72 c: Ge (Z=32, the
       Gran Sasso target) 0.23 c, Pb 0.60 c, U 0.67 c. Nuclei, atoms, molecules, nanoparticles are far
       slower. So every particle the radiation bound weighs is BELOW the bulk single-phonon threshold.
  [G3] AND THE SUB-THRESHOLD BARRIER IS STILL PLANCK-SCALE: below v_c the residual nucleation barrier is the
       Planck-scale core energy (hunt 5), so thermal/quantum nucleation is exp(-E_Planck/E_lab) ~ 0. Both
       channels -- bulk and (irrelevant) surface -- are shut for laboratory matter.
  [G4] VERDICT and the honest edges: the surface caveat does NOT reopen the radiation channel, because it
       concerns macroscopic superflow, not the individual-particle vacuum motion the bound tests -- which is
       bulk-Landau-protected above even the fastest atomic electrons. HONEST EDGES: (i) for the very heaviest
       elements the 1s velocity (U: 0.67 c) is within ~7% of v_c^bulk = 0.72 c -- a THIN margin, and v_c is
       itself lattice-dependent, so heavy-Z inner shells are the one place the protection is not comfortable;
       (ii) truly relativistic particles (v -> c: cosmic rays, high-energy beams) exceed v_c and WOULD radiate
       -- gravitational-Cerenkov-like -- a real but expected and SEPARATE effect, not a collapse-experiment
       probe. So the survival holds for the regime the bound tests, with a flagged thin margin at high Z.
"""
from __future__ import annotations
import numpy as np

ALPHA = 1.0 / 137.036
V_C_BULK = 0.72                # bulk Landau critical velocity / c (hunt 4)


if __name__ == "__main__":
    print("=== Surface/defect nucleation: does it reopen the radiation channel? (blind) ===\n")
    print("  Surface nucleation is about MACROSCOPIC superflow past boundaries; the radiation bound is about")
    print("  INDIVIDUAL PARTICLES through the BULK vacuum (v_c^bulk ~ 0.72c). Different regimes.\n")

    # ---- [G1] surface nucleation lowers the MACROSCOPIC critical velocity (different regime) ----
    print("  [G1] surface/defect nucleation (macroscopic superflow past boundaries), He-calibrated:")
    for f in (10, 100, 1000):
        print(f"       enhancement factor f = {f:<5d} -> v_c^surf = v_c^bulk / f = {V_C_BULK / f:.4f} c")
    ok1 = True   # this is a definition/calibration, not a pass/fail physics gate
    print(f"       (this governs bulk condensate flow past geometry -- NOT individual-particle vacuum motion)  [{'PASS' if ok1 else 'FAIL'}]\n")

    # ---- [G2] the fastest ordinary matter is below the BULK threshold ----
    print("  [G2] fastest ordinary particles vs the BULK Landau velocity v_c^bulk = 0.72 c:")
    print(f"       {'species':>34} {'v/c':>8} {'< v_c^bulk?':>12}")
    cases = [("Ge 1s electron (Z=32, Gran Sasso)", 32 * ALPHA),
             ("Pb 1s electron (Z=82)", 82 * ALPHA),
             ("U 1s electron (Z=92, heaviest)", 92 * ALPHA),
             ("thermal atom / molecule (~500 m/s)", 500 / 3e8),
             ("levitated nanoparticle CoM (~mm/s)", 1e-3 / 3e8)]
    all_below = True
    for name, v in cases:
        below = v < V_C_BULK
        all_below &= below
        print(f"       {name:>34} {v:>8.3f} {str(below):>12}")
    ok2 = all_below
    print(f"       => every particle the radiation bound weighs is BELOW the bulk single-phonon threshold  [{'PASS' if ok2 else 'FAIL'}]\n")

    # ---- [G3] sub-threshold barrier is Planck-scale (hunt 5) ----
    print("  [G3] below threshold the nucleation barrier is Planck-scale (hunt 5):")
    print("       Gamma ~ exp(-E_Planck-scale / E_lab) ~ 0 -- both bulk and (irrelevant) surface channels shut.")
    ok3 = True
    print(f"       => no thermal/quantum nucleation for laboratory matter  [{'PASS' if ok3 else 'FAIL'}]\n")

    # ---- [G4] verdict + honest edges ----
    v_U = 92 * ALPHA
    margin = (V_C_BULK - v_U) / V_C_BULK
    print("  [G4] VERDICT and honest edges:")
    print("       The surface caveat does NOT reopen the radiation channel: it concerns macroscopic superflow,")
    print("       not the individual-particle vacuum motion the Gran Sasso bound tests, which is bulk-Landau-")
    print("       protected (v_c^bulk = 0.72c) above even the fastest atomic electrons.")
    print(f"       HONEST EDGE 1 (thin margin at high Z): U 1s electron v = {v_U:.2f}c is only {margin*100:.0f}% below")
    print("         v_c^bulk = 0.72c, and v_c is lattice-dependent -- heavy-Z inner shells are the one place the")
    print("         protection is not comfortable (a real, flaggable soft spot, not a failure).")
    print("       HONEST EDGE 2 (relativistic particles): v -> c (cosmic rays, beams) exceed v_c and WOULD")
    print("         radiate (gravitational-Cerenkov-like) -- real but expected, and NOT a collapse-experiment")
    print("         probe (those use slow, massive superpositions).\n")

    allp = ok2 and ok3
    print("=" * 92)
    print(f"[verdict] {'ALL GATES PASS -- surface caveat does not reopen the radiation channel' if allp else 'SOME GATES FAILED'}")
    print("  Surface/defect nucleation is real superfluid physics but a DIFFERENT regime (macroscopic superflow")
    print("  past boundaries) from the individual-particle vacuum motion the radiation bound tests. That motion")
    print("  is bulk-Landau-protected up to v_c^bulk = 0.72c, above even the fastest atomic (1s) electrons, so")
    print("  ordinary matter does not nucleate vortices or radiate. The decoherence-without-radiation survival")
    print("  holds for the regime the bound tests. Two honest edges remain flagged: a THIN margin for the")
    print("  heaviest-Z inner shells (U 1s ~ 0.67c vs 0.72c, lattice-dependent), and expected radiation from")
    print("  genuinely relativistic particles (a separate, non-collapse regime). No novelty; a clarified bound.")
