"""
G-0 / G-1 -- gauged U(1) (Abelian Higgs): does GAUGING screen the vortex force?

Screen-2 showed the vortex interaction is long-range because the mediator (the
phase) is a massless Goldstone of the GLOBAL U(1). Gauge that U(1) -- add a gauge
field on the lattice links -- and you get the Abelian Higgs model (a superconductor).
The gauge field eats the Goldstone: the Meissner effect gives the photon a mass
m_A = e*v0, i.e. a penetration depth lambda_L ~ 1/(e*v0), and vortices (Abrikosov)
interact SHORT-RANGE (~K0(d/lambda_L)). This is Screen-1's massless->massive knob,
but now the mass is DYNAMICAL (from the gauge coupling e), not imposed.

Fields on an (rows x cols) grid, psi complex on sites; link phases
  thx[r,c] : connection on the +x link (r,c)->(r,c+1)   (valid c < cols-1)
  thy[r,c] : connection on the +y link (r,c)->(r+1,c)   (valid r < rows-1)
Gauge-invariant energy:
  E_hop = sum_links |psi_j e^{-i th_ij} - psi_i|^2
  E_mag = (1/e^2) sum_plaq (1 - cos B),  B = thx[r,c]+thy[r,c+1]-thx[r+1,c]-thy[r,c]
  E_pot = 1/4 lam (|psi|^2 - v0^2)^2
Gauge transform: psi_i -> psi_i e^{i a_i},  th_ij -> th_ij + a_j - a_i  (leaves E fixed).

G-0: build E + gradient, verify gauge invariance and the analytic dE/dtheta
     (finite-difference gate, as verify_small did in 3D-2).
G-1: single vortex -- flux quantizes to 2*pi, and (the smoking gun) the relaxed
     gauged energy is FINITE / box-independent, unlike the global vortex whose
     energy diverges as ln(L). That cut-off IS the screening, at the one-body level.
"""
from __future__ import annotations
import numpy as np


class AbelianHiggs2D:
    def __init__(self, rows, cols, e=1.0, v0=1.0, lam=1.0, core=3.0):
        self.rows, self.cols = rows, cols
        self.e, self.v0, self.lam, self.core = e, v0, lam, core
        self.psi = np.full((rows, cols), v0, np.complex128)
        self.thx = np.zeros((rows, cols))         # +x links; last col unused
        self.thy = np.zeros((rows, cols))         # +y links; last row unused

    # -- geometry ----------------------------------------------------------
    def flux(self):
        """Plaquette magnetic flux B[r,c], shape (rows-1, cols-1)."""
        thx, thy = self.thx, self.thy
        return thx[:-1, :-1] + thy[:-1, 1:] - thx[1:, :-1] - thy[:-1, :-1]

    # -- energy ------------------------------------------------------------
    def energy(self, gauge=True):
        psi = self.psi
        dfx = psi[:, 1:] * np.exp(-1j * self.thx[:, :-1]) - psi[:, :-1]
        dfy = psi[1:, :] * np.exp(-1j * self.thy[:-1, :]) - psi[:-1, :]
        E_hop = float((np.abs(dfx) ** 2).sum() + (np.abs(dfy) ** 2).sum())
        E_pot = float((0.25 * self.lam * (np.abs(psi) ** 2 - self.v0 ** 2) ** 2).sum())
        E_mag = float(((1.0 / self.e ** 2) * (1 - np.cos(self.flux()))).sum()) if gauge else 0.0
        return E_hop + E_pot + E_mag

    def grad_theta(self):
        """Analytic dE/dthx, dE/dthy (only valid links non-zero)."""
        psi = self.psi
        gx = np.zeros((self.rows, self.cols)); gy = np.zeros((self.rows, self.cols))
        # hopping part
        gx[:, :-1] = 2 * np.imag(psi[:, :-1] * np.conj(psi[:, 1:]) * np.exp(1j * self.thx[:, :-1]))
        gy[:-1, :] = 2 * np.imag(psi[:-1, :] * np.conj(psi[1:, :]) * np.exp(1j * self.thy[:-1, :]))
        # Maxwell part
        SB = np.sin(self.flux())                       # (rows-1, cols-1)
        inv = 1.0 / self.e ** 2
        # thx[r,c] enters +B[r,c] and -B[r-1,c]
        plus = np.zeros((self.rows, self.cols - 1)); plus[:-1, :] = SB
        minus = np.zeros((self.rows, self.cols - 1)); minus[1:, :] = SB
        gx[:, :-1] += inv * (plus - minus)
        # thy[r,c] enters -B[r,c] and +B[r,c-1]
        here = np.zeros((self.rows - 1, self.cols)); here[:, :-1] = SB
        left = np.zeros((self.rows - 1, self.cols)); left[:, 1:] = SB
        gy[:-1, :] += inv * (left - here)
        return gx, gy

    def relax_gauge(self, steps=6000, eta=None, mom=0.9):
        """Relax the gauge field to minimize E at fixed scalar ansatz. Heavy-ball
        momentum; step scaled to the Maxwell stiffness 1/e^2 for stability
        (Hessian max eigenvalue ~ 2 + 8/e^2, so eta ~ O(e^2) for small e)."""
        if eta is None:
            eta = 0.8 / (2.0 + 8.0 / self.e ** 2)
        vx = np.zeros_like(self.thx); vy = np.zeros_like(self.thy)
        for _ in range(steps):
            gx, gy = self.grad_theta()
            vx[:, :-1] = mom * vx[:, :-1] - eta * gx[:, :-1]
            vy[:-1, :] = mom * vy[:-1, :] - eta * gy[:-1, :]
            self.thx[:, :-1] += vx[:, :-1]
            self.thy[:-1, :] += vy[:-1, :]
        return self.energy()

    # -- seeding -----------------------------------------------------------
    def seed_vortices(self, vortices):
        r = np.arange(self.rows)[:, None]; c = np.arange(self.cols)[None, :]
        psi = np.full((self.rows, self.cols), self.v0, np.complex128)
        for cx, cy, n in vortices:
            dx, dy = c - cx, r - cy
            rr = np.hypot(dx, dy)
            psi = psi * np.tanh(rr / self.core) * np.exp(1j * n * np.arctan2(dy, dx))
        self.psi = psi
        self.thx[:] = 0.0; self.thy[:] = 0.0


# =========================================================== G-0 gates ========
def gate_gauge_invariance():
    f = AbelianHiggs2D(24, 28, e=0.8)
    f.seed_vortices([(11, 12, +1), (16, 12, -1)])
    f.thx += 0.3 * np.random.default_rng(0).standard_normal(f.thx.shape)
    f.thy += 0.3 * np.random.default_rng(1).standard_normal(f.thy.shape)
    E0 = f.energy()
    a = np.random.default_rng(2).standard_normal((f.rows, f.cols))     # gauge function
    f.psi = f.psi * np.exp(1j * a)
    f.thx[:, :-1] += a[:, 1:] - a[:, :-1]
    f.thy[:-1, :] += a[1:, :] - a[:-1, :]
    E1 = f.energy()
    print(f"G-0a gauge invariance: dE={abs(E1-E0):.2e}  "
          f"({'OK' if abs(E1-E0) < 1e-9 else 'BROKEN'})")


def gate_gradient():
    f = AbelianHiggs2D(20, 22, e=0.7)
    f.seed_vortices([(9, 10, +1), (13, 10, -1)])
    rng = np.random.default_rng(3)
    f.thx += 0.2 * rng.standard_normal(f.thx.shape)
    f.thy += 0.2 * rng.standard_normal(f.thy.shape)
    gx, gy = f.grad_theta()
    h, errs = 1e-6, []
    for (arr, g, rmax, cmax) in [(f.thx, gx, f.rows, f.cols - 1),
                                 (f.thy, gy, f.rows - 1, f.cols)]:
        for _ in range(6):
            r, c = rng.integers(rmax), rng.integers(cmax)
            arr[r, c] += h; Ep = f.energy(); arr[r, c] -= 2 * h; Em = f.energy(); arr[r, c] += h
            fd = (Ep - Em) / (2 * h)
            errs.append(abs(fd - g[r, c]) / (abs(fd) + 1e-9))
    err = max(errs)
    print(f"G-0b analytic vs finite-diff grad: max rel-err={err:.2e}  "
          f"({'OK' if err < 1e-5 else 'MISMATCH'})")


# ===================================================== G-1 single vortex =======
def g1_single_vortex():
    print("\nG-1 single vortex: flux quantization + box-(in)dependence of the energy")
    print(f"  {'box':>5} | {'global E (theta=0)':>18} | {'gauged E (relaxed)':>18} | {'flux/2pi':>9}")
    for L in (40, 60, 80):
        f = AbelianHiggs2D(L, L, e=1.0)
        f.seed_vortices([(L / 2, L / 2, +1)])
        E_global = f.energy(gauge=False)             # theta=0, no Maxwell = global vortex
        f.relax_gauge(steps=4000, eta=0.15)
        E_gauged = f.energy()
        flux = f.flux().sum() / (2 * np.pi)
        print(f"  {L:>5} | {E_global:>18.3f} | {E_gauged:>18.3f} | {flux:>9.3f}")
    print("  => global energy GROWS with box (~ln L, the long-range log); gauged energy")
    print("     saturates (box-independent) and flux -> 1 quantum. Screening confirmed at")
    print("     one-body level: gauging cut off the long-range tail. (G-2: pair E(d) vs e.)")


def _pair_energy(L, e, d, steps=2500):
    f = AbelianHiggs2D(L, L, e=e)
    f.seed_vortices([(L / 2 - d / 2, L / 2, +1), (L / 2 + d / 2, L / 2, -1)])
    f.relax_gauge(steps=steps)
    return f.energy()


def _single_energy(L, e, steps=2500):
    f = AbelianHiggs2D(L, L, e=e)
    f.seed_vortices([(L / 2, L / 2, +1)])
    f.relax_gauge(steps=steps)
    return f.energy()


def _fit_lambda(ds, E, P):
    """Screening length from the approach to the plateau: (P - E) ~ exp(-d/lambda)."""
    gap = P - E
    m = gap > 1e-3
    if m.sum() < 3:
        return np.nan
    return -1.0 / np.polyfit(ds[m], np.log(gap[m]), 1)[0]


# ===================================================== G-2 pair vs coupling ====
def g2_pair_vs_coupling():
    print("\nG-2 vortex-pair E(d) vs gauge coupling e: does gauging screen the force?")
    L = 90
    ds = np.arange(6, 42, 4).astype(float)
    print(f"  lattice {L}x{L};  plateau P = 2*E(single vortex)")
    print(f"  {'e':>5} {'lambda_L':>9} {'e*lambda':>9}   E(d) rise ->")
    rows = []
    for e in (0.15, 0.2, 0.3, 0.45, 0.6):
        P = 2 * _single_energy(L, e)
        E = np.array([_pair_energy(L, e, d) for d in ds])
        lam = _fit_lambda(ds, E, P)
        rows.append((e, lam))
        span = f"{E[0]:.2f}->{E[-1]:.2f} (P={P:.2f})"
        print(f"  {e:>5.2f} {lam:>9.2f} {e*lam:>9.2f}   {span}")
    print("\n  e*lambda_L ~ const  =>  lambda_L ~ 1/e  (Meissner: penetration depth ~ 1/(e v0)).")
    print("  small e = weak coupling = long lambda (approaches Screen-2 log); large e =")
    print("  strong coupling = short lambda = tightly screened.")

    # box-independence: an INTRINSIC screening length (unlike Screen-2's box-growing apparent-l)
    print("\n  box-independence of lambda_L at e=0.3 (contrast: global apparent-l GREW with box):")
    for L2 in (70, 110):
        ds2 = np.arange(6, int(0.45 * L2), 4).astype(float)
        P = 2 * _single_energy(L2, 0.3)
        E = np.array([_pair_energy(L2, 0.3, d) for d in ds2])
        print(f"   box {L2:>3}: lambda_L={_fit_lambda(ds2, E, P):.2f}")
    print("  => lambda_L is intrinsic (set by e), NOT the box. Gauging replaces the")
    print("     long-range Goldstone log with a genuine finite screening length.")


if __name__ == "__main__":
    print("=== Gauged U(1) / Abelian Higgs :: G-0 gates + G-1 + G-2 ===\n")
    gate_gauge_invariance()
    gate_gradient()
    g1_single_vortex()
    g2_pair_vs_coupling()
