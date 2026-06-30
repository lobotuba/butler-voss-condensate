"""
Butler-Voss Condensate  —  Research Engine (lattice- and dimension-agnostic)
============================================================================

Concept (Robert Voss "Butler-Voss Condensate", formerly the "Grain Fabric Model")
---------------------------------------------------------------------------------
- Space is a network of nodes. A node here is NOT a bare graph vertex but a
  physical locus: it carries a scalar displacement u, a velocity v, and a local
  tension, and is joined to its neighbors by tension-bearing "strings" -- the
  edges are real relations, not mere adjacency.
- Disturbances propagate as ripples (a wave equation on the lattice graph).
- A nonlinear term lets strong overlaps behave differently from small ripples,
  so localized structures ("particles") can in principle form.
- Where local energy is high, strings tighten. Excitations then drift up the
  tension gradient -> the candidate mechanism for gravity.

This is a physics-INSPIRED toy model, not a claim about the real universe.
The engine produces MEASUREMENTS so you can state a hypothesis and then confirm
or refute it from the numbers.

----------------------------------------------------------------------------
Why this version is lattice- AND dimension-agnostic
----------------------------------------------------------------------------
The lattice (hex, square, cubic, fcc) is only scaffolding for approximating
continuous space. A REAL effect should survive changing the lattice; an effect
that changes qualitatively with the lattice was a numerical artifact.

  * Run the SAME experiment on different --lattice values and compare.
  * 2D vs 3D genuinely differ:
      - gravity falloff is 1/r^2 only in 3D (in 2D it is ~1/r),
      - dispersion is faster in 3D, so lumps decay quicker,
      - Derrick's theorem: simple scalar solitons cannot be statically stable
        in 3D by energy scaling alone -> expect 3D particles to be HARDER.

Available lattices:
    hex2d     6 neighbors, isotropic 2D (default; matches earlier versions)
    square2d  4 neighbors, anisotropic 2D (waves prefer the axes)
    cubic3d   6 neighbors, 3D simple cubic
    fcc3d     12 neighbors, 3D close-packed (the true 3D analog of hex)

----------------------------------------------------------------------------
Hypotheses
----------------------------------------------------------------------------
H1 PARTICLE FORMATION   nonlinear overlap creates lumps that persist far longer
                        than linear ripples disperse.  (compare --nonlinear 0.0)
H2 GRAVITY              two compression regions drift together. (separation down)
H3 MOTION THROUGH MEDIUM a particle is a moving pattern; nodes only oscillate.
H4 LATTICE INDEPENDENCE the H1-H3 results are the SAME across --lattice choices.
H5 DIMENSIONALITY        2D vs 3D change the results (esp. particle stability).

----------------------------------------------------------------------------
Usage
----------------------------------------------------------------------------
    python simulation.py --live
    python simulation.py --live --lattice cubic3d --experiment gravity
    python simulation.py --headless --lattice square2d --experiment collide --steps 4000
    python simulation.py --headless --lattice fcc3d --experiment gravity --steps 6000
    # lattice-independence test: run the same experiment on each lattice and diff the CSVs
"""

from __future__ import annotations

import os
import csv
import math
import argparse
from dataclasses import dataclass

import numpy as np


# ============================================================
# Configuration
# ============================================================

@dataclass
class Config:
    lattice: str = "hex2d"
    size: int = 0                    # 0 = use the lattice's default size

    dt: float = 0.06
    base_wave_speed: float = 1.0
    base_tension: float = 1.0
    damping: float = 0.999           # 1.0 = energy-conserving (good for clean tests)

    # --- self-interaction potential V(u); dynamics use force = -V'(u) ---
    #   "oscillon": -m^2 u + focus*u^3 - sat*u^5   long-lived oscillons (DEFAULT)
    #   "soft":     -nonlinear_strength * u^3       old softening model
    #   "wave":     0                               massless linear control
    # H1 control: keep "oscillon" but set --focus 0  (a massive linear wave).
    potential: str = "oscillon"
    mass2: float = 1.0          # m^2: linear restoring -> oscillation/mass scale
    focus: float = 1.0          # +u^3 focusing term (self-trapping); 0.0 = control
    saturation: float = 0.5     # -u^5 term that prevents collapse / blow-up
    nonlinear_strength: float = 0.06   # only used by potential="soft"

    gravity_strength: float = 0.3    # 0.0 = pre-gravity model
    max_drift: float = 0.5           # CFL guard on drift speed

    energy_threshold: float = 0.20   # blob / tightening cutoff
    tightening_rate: float = 0.004
    relaxation_rate: float = 0.001
    min_tension: float = 0.6
    max_tension: float = 4.0

    match_radius: float = 4.0        # max centroid move to keep a particle id
    min_blob_cells: int = 3

    steps_per_measure: int = 5
    seed: int = 0


# ============================================================
# Lattice builders
# ------------------------------------------------------------
# Each returns (pos, neighbor_idx, valid):
#   pos          (N, D) float coordinates
#   neighbor_idx (N, K) int  flat index of each neighbor (self where invalid)
#   valid        (N, K) float 1.0 where the neighbor exists else 0.0
# ============================================================

def _assemble(coords, offsets):
    """Generic adjacency from an integer coordinate list and fixed offsets."""
    index_of = {tuple(c): i for i, c in enumerate(coords)}
    N, K, D = len(coords), len(offsets), len(coords[0])
    nidx = np.zeros((N, K), dtype=np.int64)
    valid = np.zeros((N, K), dtype=np.float64)
    for i, c in enumerate(coords):
        for k, off in enumerate(offsets):
            nc = tuple(c[d] + off[d] for d in range(D))
            j = index_of.get(nc)
            if j is not None:
                nidx[i, k] = j
                valid[i, k] = 1.0
            else:
                nidx[i, k] = i
    return nidx, valid


def make_hex2d(size):
    rows = size or 70
    cols = int(round(rows * 1.28))
    N = rows * cols
    pos = np.zeros((N, 2))
    coords = []
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            pos[i] = (c + 0.5 * (r % 2), r * math.sqrt(3) / 2.0)
            coords.append((r, c))
    # hex neighbor offsets depend on row parity, so build adjacency directly
    nidx = np.zeros((N, 6), dtype=np.int64)
    valid = np.zeros((N, 6), dtype=np.float64)
    for r in range(rows):
        if r % 2 == 0:
            offs = [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
        else:
            offs = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
        for c in range(cols):
            i = r * cols + c
            for k, (dr, dc) in enumerate(offs):
                rr, cc = r + dr, c + dc
                if 0 <= rr < rows and 0 <= cc < cols:
                    nidx[i, k] = rr * cols + cc
                    valid[i, k] = 1.0
                else:
                    nidx[i, k] = i
    return pos, nidx, valid


def make_square2d(size):
    rows = size or 72
    cols = int(round(rows * 1.28))
    coords = [(r, c) for r in range(rows) for c in range(cols)]
    pos = np.array([(c, r) for (r, c) in coords], dtype=float)
    offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]  # 4-neighbor (anisotropic)
    nidx, valid = _assemble(coords, offsets)
    return pos, nidx, valid


def make_cubic3d(size):
    n = size or 18
    coords = [(i, j, k) for i in range(n) for j in range(n) for k in range(n)]
    pos = np.array(coords, dtype=float)
    offsets = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    nidx, valid = _assemble(coords, offsets)
    return pos, nidx, valid


def make_fcc3d(size):
    # FCC = even-sum sites of a cubic grid; 12 nearest neighbors at (+/-1,+/-1,0)-type.
    L = size or 24
    coords = [(i, j, k)
              for i in range(L) for j in range(L) for k in range(L)
              if (i + j + k) % 2 == 0]
    pos = np.array(coords, dtype=float)
    offsets = []
    for a in (-1, 1):
        for b in (-1, 1):
            offsets += [(a, b, 0), (a, 0, b), (0, a, b)]  # 12 close-packed neighbors
    nidx, valid = _assemble(coords, offsets)
    return pos, nidx, valid


LATTICES = {
    "hex2d": make_hex2d,
    "square2d": make_square2d,
    "cubic3d": make_cubic3d,
    "fcc3d": make_fcc3d,
}


# ============================================================
# Fabric engine (works in any dimension on any neighbor graph)
# ============================================================

class Fabric:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed)

        pos, nidx, valid = LATTICES[cfg.lattice](cfg.size)
        self.pos = pos
        self.neighbor_idx = nidx
        self.valid = valid
        self.N, self.D = pos.shape

        self.u = np.zeros(self.N)
        self.v = np.zeros(self.N)
        self.tension = np.full(self.N, cfg.base_tension)
        self.step_count = 0
        self.time = 0.0

        self.lo = pos.min(axis=0)
        self.hi = pos.max(axis=0)

        self._build_operators()

    # ---- precompute graph operators ------------------------------------

    def _build_operators(self):
        self.neighbor_count = np.clip(self.valid.sum(axis=1), 1.0, None)
        # Normalize the lattice operators by a single CONSTANT (the typical
        # degree), not per-node counts. A constant keeps the graph Laplacian
        # SYMMETRIC at boundaries, so the spring force stays the exact gradient
        # of the spring energy -> energy is conserved (no boundary pumping),
        # while the wave speed stays comparable across lattices of different
        # coordination number.
        self.deg_norm = float(self.valid.sum(axis=1).max())

        # neighbor delta vectors, masked: (N, K, D)
        self.dvec = (self.pos[self.neighbor_idx] - self.pos[:, None, :]) \
            * self.valid[:, :, None]

        # per-node second-moment matrix M = sum_k dvec dvec^T  -> (N, D, D)
        M = np.einsum("nkd,nke->nde", self.dvec, self.dvec)
        M = M + np.eye(self.D)[None, :, :] * 1e-9   # regularize degenerate nodes
        self.Minv = np.linalg.inv(M)

    def _laplacian(self, field):
        # symmetric graph Laplacian  sum_j (u_j - u_i)  divided by a constant
        neigh = field[self.neighbor_idx] * self.valid
        return (neigh.sum(axis=1) - self.neighbor_count * field) / self.deg_norm

    def _gradient(self, field):
        """Least-squares spatial gradient -> (N, D)."""
        df = (field[self.neighbor_idx] - field[:, None]) * self.valid     # (N,K)
        b = np.einsum("nk,nkd->nd", df, self.dvec)                        # (N,D)
        return np.einsum("nde,ne->nd", self.Minv, b)                      # (N,D)

    # ---- excitations (placed by fractional position, any dimension) -----

    def _center(self, frac):
        f = np.full(self.D, 0.5)
        f[:len(frac)] = frac
        return self.lo + f * (self.hi - self.lo)

    def add_ripple(self, frac, amplitude=1.0, radius=5.0):
        c = self._center(frac)
        d2 = np.sum((self.pos - c) ** 2, axis=1)
        self.u += amplitude * np.exp(-d2 / (2 * radius ** 2))

    def add_velocity_pulse(self, frac, amplitude=0.2, radius=4.0):
        c = self._center(frac)
        d2 = np.sum((self.pos - c) ** 2, axis=1)
        self.v += amplitude * np.exp(-d2 / (2 * radius ** 2))

    # ---- measurements ---------------------------------------------------

    # ---- self-interaction potential (force and energy stay consistent) --

    def potential_force(self):
        """force = -dV/du for the active potential."""
        cfg, u = self.cfg, self.u
        if cfg.potential == "oscillon":
            return -cfg.mass2 * u + cfg.focus * u ** 3 - cfg.saturation * u ** 5
        if cfg.potential == "soft":
            return -cfg.nonlinear_strength * u ** 3
        return 0.0  # "wave": massless linear

    def potential_energy(self):
        """V(u) for the active potential (>= 0)."""
        cfg, u = self.cfg, self.u
        if cfg.potential == "oscillon":
            return (0.5 * cfg.mass2 * u ** 2
                    - 0.25 * cfg.focus * u ** 4
                    + (cfg.saturation / 6.0) * u ** 6)
        if cfg.potential == "soft":
            return 0.25 * cfg.nonlinear_strength * u ** 4
        return np.zeros_like(u)

    def energy_density(self):
        kinetic = 0.5 * self.v ** 2
        potential = self.potential_energy()
        diff = (self.u[self.neighbor_idx] - self.u[:, None]) * self.valid
        gradient = 0.5 * self.tension * np.sum(diff ** 2, axis=1) / self.deg_norm
        return kinetic + potential + gradient

    def _update_tension(self, e):
        cfg = self.cfg
        hi = e > cfg.energy_threshold
        self.tension[hi] += cfg.tightening_rate * e[hi]
        lo = ~hi
        self.tension[lo] -= cfg.relaxation_rate * (self.tension[lo] - cfg.base_tension)
        np.clip(self.tension, cfg.min_tension, cfg.max_tension, out=self.tension)

    # ---- time step ------------------------------------------------------

    def step(self):
        cfg = self.cfg
        e = self.energy_density()
        self._update_tension(e)

        lap = self._laplacian(self.u)
        tension_force = cfg.base_wave_speed ** 2 * self.tension * lap
        accel = tension_force + self.potential_force()

        self.v += accel * cfg.dt
        self.v *= cfg.damping
        self.u += self.v * cfg.dt

        # gravity: advect excitations up the tension gradient.
        # Drift w = G*grad(tension); transport eqn du/dt = -w.grad(u).
        # Proportional to grad(u), so empty fabric does not move: only mass
        # feels gravity. Also self-focuses single particles.
        if cfg.gravity_strength != 0.0:
            w = cfg.gravity_strength * self._gradient(self.tension)       # (N,D)
            speed = np.linalg.norm(w, axis=1)
            scale = np.where(speed > cfg.max_drift,
                             cfg.max_drift / (speed + 1e-12), 1.0)
            w *= scale[:, None]
            gu = self._gradient(self.u)
            gv = self._gradient(self.v)
            self.u -= cfg.dt * np.einsum("nd,nd->n", w, gu)
            self.v -= cfg.dt * np.einsum("nd,nd->n", w, gv)

        self.step_count += 1
        self.time += cfg.dt

    # ---- particle (blob) detection on the true graph --------------------

    def label_blobs(self, e):
        cfg = self.cfg
        mask = e > cfg.energy_threshold
        labels = -np.ones(self.N, dtype=np.int64)
        nidx, valid = self.neighbor_idx, self.valid
        K = nidx.shape[1]
        blobs = []
        for start in np.nonzero(mask)[0]:
            if labels[start] != -1:
                continue
            stack = [start]
            labels[start] = len(blobs)
            cells = []
            while stack:
                i = stack.pop()
                cells.append(i)
                for k in range(K):
                    if valid[i, k]:
                        j = nidx[i, k]
                        if mask[j] and labels[j] == -1:
                            labels[j] = len(blobs)
                            stack.append(j)
            if len(cells) < cfg.min_blob_cells:
                continue
            cells = np.array(cells)
            w = e[cells]
            mass = float(w.sum())
            centroid = (self.pos[cells] * w[:, None]).sum(axis=0) / mass
            blobs.append({"cells": cells, "mass": mass, "peak": float(w.max()),
                          "pos": centroid})
        return blobs


# ============================================================
# Particle tracker (persistent IDs across frames, any dimension)
# ============================================================

class Tracker:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.next_id = 0
        self.tracks = {}

    def update(self, blobs, frame_dt):
        cfg = self.cfg
        prev, used, new = self.tracks, set(), {}
        for b in blobs:
            best_id, best_d = None, cfg.match_radius
            for tid, t in prev.items():
                if tid in used:
                    continue
                d = float(np.linalg.norm(b["pos"] - t["pos"]))
                if d < best_d:
                    best_id, best_d = tid, d
            if best_id is None:
                tid = self.next_id
                self.next_id += 1
                new[tid] = {"pos": b["pos"], "vel": np.zeros_like(b["pos"]),
                            "lifetime": frame_dt, "mass": b["mass"], "peak": b["peak"]}
            else:
                used.add(best_id)
                t = prev[best_id]
                new[best_id] = {"pos": b["pos"], "vel": (b["pos"] - t["pos"]) / frame_dt,
                                "lifetime": t["lifetime"] + frame_dt,
                                "mass": b["mass"], "peak": b["peak"]}
        self.tracks = new
        return new


# ============================================================
# Experiments (positions are fractional, so they work in any lattice/dim)
# ============================================================

def setup_collide(f: Fabric):
    f.add_ripple((0.28, 0.40), amplitude=1.1, radius=5.5)
    f.add_ripple((0.62, 0.80), amplitude=-1.0, radius=6.5)
    f.add_ripple((0.45, 0.50), amplitude=0.8, radius=4.5)
    f.add_velocity_pulse((0.28, 0.40), amplitude=0.18, radius=5.0)
    f.add_velocity_pulse((0.62, 0.80), amplitude=-0.16, radius=5.0)
    f.add_velocity_pulse((0.45, 0.50), amplitude=0.12, radius=4.0)


def setup_gravity(f: Fabric):
    f.add_ripple((0.30, 0.50), amplitude=1.6, radius=3.5)
    f.add_ripple((0.70, 0.50), amplitude=1.6, radius=3.5)
    f.v += f.rng.normal(0, 0.01, f.N)


def setup_travel(f: Fabric):
    f.add_ripple((0.22, 0.50), amplitude=1.4, radius=4.0)
    f.add_velocity_pulse((0.30, 0.50), amplitude=0.25, radius=4.0)  # net +x momentum


def setup_lump(f: Fabric):
    # one clean centered lump to nucleate a single oscillon (best for H1)
    f.add_ripple((0.50, 0.50), amplitude=1.3, radius=4.0)


EXPERIMENTS = {"collide": setup_collide, "gravity": setup_gravity,
               "travel": setup_travel, "lump": setup_lump}


# ============================================================
# Logging
# ============================================================

class Logger:
    def __init__(self, out_dir, D):
        os.makedirs(out_dir, exist_ok=True)
        self.D = D
        self.summary_path = os.path.join(out_dir, "summary.csv")
        self.particles_path = os.path.join(out_dir, "particles.csv")
        self._sf = open(self.summary_path, "w", newline="")
        self._pf = open(self.particles_path, "w", newline="")
        self._sw = csv.writer(self._sf)
        self._pw = csv.writer(self._pf)
        self._sw.writerow(["step", "time", "total_energy", "max_energy",
                           "mean_tension", "max_tension", "n_particles",
                           "total_particle_mass", "max_lifetime", "separation"])
        self._pw.writerow(["step", "time", "id", "mass", "peak",
                           "x", "y", "z", "speed", "lifetime"])

    def log(self, fab, tracks, e):
        n = len(tracks)
        total_mass = sum(t["mass"] for t in tracks.values())
        max_life = max((t["lifetime"] for t in tracks.values()), default=0.0)
        sep = ""
        if n >= 2:
            two = sorted(tracks.values(), key=lambda t: t["mass"], reverse=True)[:2]
            sep = float(np.linalg.norm(two[0]["pos"] - two[1]["pos"]))
        self._sw.writerow([fab.step_count, round(fab.time, 4),
                           round(float(e.sum()), 5), round(float(e.max()), 5),
                           round(float(fab.tension.mean()), 5),
                           round(float(fab.tension.max()), 5),
                           n, round(total_mass, 5), round(max_life, 4),
                           round(sep, 4) if sep != "" else ""])
        for tid, t in tracks.items():
            p = t["pos"]
            x = round(float(p[0]), 4)
            y = round(float(p[1]), 4) if self.D >= 2 else ""
            z = round(float(p[2]), 4) if self.D >= 3 else ""
            speed = float(np.linalg.norm(t["vel"]))
            self._pw.writerow([fab.step_count, round(fab.time, 4), tid,
                               round(t["mass"], 5), round(t["peak"], 5),
                               x, y, z, round(speed, 5), round(t["lifetime"], 4)])

    def close(self):
        self._sf.close()
        self._pf.close()


# ============================================================
# Run loops
# ============================================================

def run_headless(fab: Fabric, steps, out_dir):
    cfg = fab.cfg
    tracker = Tracker(cfg)
    logger = Logger(out_dir, fab.D)
    frame_dt = cfg.dt * cfg.steps_per_measure
    print(f"[headless] lattice={cfg.lattice} D={fab.D} N={fab.N} | {steps} steps -> {out_dir}")
    for _ in range(0, steps, cfg.steps_per_measure):
        for _ in range(cfg.steps_per_measure):
            fab.step()
        e = fab.energy_density()
        tracks = tracker.update(fab.label_blobs(e), frame_dt)
        logger.log(fab, tracks, e)
    logger.close()
    print(f"[done] {logger.summary_path}")
    print(f"[done] {logger.particles_path}")


def _viz_indices(fab: Fabric):
    """For 3D, take a thin slice near the mid-plane; for <=2D, take everything."""
    if fab.D <= 2:
        return np.arange(fab.N)
    zmid = np.median(fab.pos[:, 2])
    return np.nonzero(np.abs(fab.pos[:, 2] - round(float(zmid))) < 0.5)[0]


def run_live(fab: Fabric, steps, out_dir):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    cfg = fab.cfg
    tracker = Tracker(cfg)
    logger = Logger(out_dir, fab.D)
    frame_dt = cfg.dt * cfg.steps_per_measure

    sl = _viz_indices(fab)
    px, py = fab.pos[sl, 0], fab.pos[sl, 1]
    slice_note = "" if fab.D <= 2 else "  (mid-plane slice)"

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    e0 = fab.energy_density()
    disp = axes[0].scatter(px, py, c=fab.u[sl], cmap="RdBu", s=18, vmin=-1, vmax=1)
    en = axes[1].scatter(px, py, c=e0[sl], cmap="inferno", s=18, vmin=0, vmax=0.6)
    ten = axes[2].scatter(px, py, c=fab.tension[sl], cmap="viridis", s=18,
                          vmin=cfg.min_tension, vmax=cfg.max_tension)
    for ax, title in zip(axes, ["Displacement", "Energy", "Tension"]):
        ax.set_title(title + slice_note)
        ax.set_aspect("equal")
        ax.axis("off")
    marks = axes[1].scatter([], [], s=120, facecolors="none", edgecolors="cyan", linewidths=1.5)
    for ax, m in zip(axes, [disp, en, ten]):
        fig.colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()

    n_frames = max(1, steps // cfg.steps_per_measure)

    def animate(_frame):
        for _ in range(cfg.steps_per_measure):
            fab.step()
        e = fab.energy_density()
        tracks = tracker.update(fab.label_blobs(e), frame_dt)
        logger.log(fab, tracks, e)

        disp.set_array(fab.u[sl])
        en.set_array(e[sl])
        ten.set_array(fab.tension[sl])
        if fab.D == 2 and tracks:
            marks.set_offsets(np.array([t["pos"][:2] for t in tracks.values()]))
        else:
            marks.set_offsets(np.empty((0, 2)))
        max_life = max((t["lifetime"] for t in tracks.values()), default=0.0)
        axes[1].set_title(f"Energy | particles: {len(tracks)} | max life: {max_life:.1f}{slice_note}")
        return disp, en, ten, marks

    ani = FuncAnimation(fig, animate, frames=n_frames, interval=30, blit=False)
    try:
        plt.show()
    finally:
        logger.close()
        print(f"[live] logged to {out_dir}")
    return ani


# ============================================================
# CLI
# ============================================================

def main():
    p = argparse.ArgumentParser(description="Butler-Voss Condensate research engine")
    p.add_argument("--lattice", choices=list(LATTICES), default="hex2d")
    p.add_argument("--size", type=int, default=0, help="lattice size (0 = default)")
    p.add_argument("--experiment", choices=list(EXPERIMENTS), default="collide")
    p.add_argument("--steps", type=int, default=3000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--potential", choices=["oscillon", "soft", "wave"], default=None)
    p.add_argument("--mass", type=float, default=None, help="m^2 for oscillon potential")
    p.add_argument("--focus", type=float, default=None, help="u^3 focusing (0 = control)")
    p.add_argument("--saturation", type=float, default=None, help="u^5 anti-collapse")
    p.add_argument("--nonlinear", type=float, default=None, help="strength for potential=soft")
    p.add_argument("--gravity", type=float, default=None)
    p.add_argument("--damping", type=float, default=None)
    p.add_argument("--out", type=str, default=None)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--live", action="store_true")
    mode.add_argument("--headless", action="store_true")
    args = p.parse_args()

    cfg = Config(lattice=args.lattice, size=args.size, seed=args.seed)
    if args.potential is not None:
        cfg.potential = args.potential
    if args.mass is not None:
        cfg.mass2 = args.mass
    if args.focus is not None:
        cfg.focus = args.focus
    if args.saturation is not None:
        cfg.saturation = args.saturation
    if args.nonlinear is not None:
        cfg.nonlinear_strength = args.nonlinear
    if args.gravity is not None:
        cfg.gravity_strength = args.gravity
    if args.damping is not None:
        cfg.damping = args.damping

    fab = Fabric(cfg)
    EXPERIMENTS[args.experiment](fab)

    out_dir = args.out or os.path.join(
        "condensate_runs", f"{args.experiment}_{args.lattice}_seed{args.seed}")

    if args.headless:
        run_headless(fab, args.steps, out_dir)
    else:
        run_live(fab, args.steps, out_dir)


if __name__ == "__main__":
    main()
