# Butler-Voss Condensate

A lattice- and dimension-agnostic research engine for the **Butler-Voss Condensate**
(Robert Voss, formerly the "Grain Fabric Model").

Space is modeled as a network of **grains**, each carrying a scalar displacement `u`
and velocity `v`, joined to neighbors by tension-bearing "strings." Disturbances
propagate as ripples (a wave equation on the lattice graph); a nonlinear potential
lets strong overlaps self-trap into long-lived localized structures ("particles"),
and high-energy regions tighten their strings — the candidate mechanism for gravity.

> This is a physics-**inspired** toy model, not a claim about the real universe.
> The engine produces measurements so a hypothesis can be stated, then confirmed
> or refuted from the numbers.

## Install

```bash
pip install -r requirements.txt
```

Requires Python 3.11+.

## Run

```bash
# Live animation of a single breathing oscillon
python simulation.py --live --experiment lump

# Headless run that writes CSVs
python simulation.py --headless --experiment collide --lattice hex2d --steps 3000

# Energy-conservation correctness check (expect small drift)
python simulation.py --headless --experiment lump --damping 1 --gravity 0
```

Output is written to `condensate_runs/<experiment>_<lattice>_seed<seed>/`
(`summary.csv` + `particles.csv`), or to `--out DIR`.

See [CHEATSHEET.md](CHEATSHEET.md) for the full parameter, lattice, experiment,
and hypothesis reference.

## Layout

| File | Purpose |
|---|---|
| `simulation.py` | The research engine (lattices, fabric dynamics, particle tracking, CLI) |
| `CHEATSHEET.md` | Equations, parameters, experiments, and hypotheses (H1–H5) |
