# certT: DP construction of the barrier (float) + rigorous cell-by-cell verification (arb).
import sys, math, time, json, hashlib
import numpy as np
from flint import arb, ctx
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from ctlib import ArcTangent, cell_Q_low, disc_low, sigma_prime_arb, sigma_prime_float, arc_float, PI, lo_f, hi_f

ctx.prec = 100


class Box:
    """parameter box: tau, A0 (= 2 m2), w (window bound), sigma'_low (arb lower bound for sigma'(m1)), Lam0"""
    def __init__(self, tau, m1, m2, wstar, Lam0='37.7', a='0.2595', Ns=4000, Nl=48, arct=None):
        self.tau_s, self.m1_s, self.m2_s, self.wstar_s, self.Lam0_s = tau, m1, m2, wstar, Lam0
        self.tau = arb(tau); self.Lam0 = arb(Lam0)
        self.A0 = 2 * arb(m2)
        wq = min(hi_f(arb(wstar)), hi_f(arb(m2) / arb(tau)))
        # window bound: any w_eff >= min(w*, m2/tau) is valid (larger region = stronger claim); wq is an exact float
        self.w = arb(wq)
        self.sig_low = lo_f(sigma_prime_arb(m1, a=a, lam=Lam0))
        self.Ns, self.Nl = Ns, Nl
        self.ds = self.A0 / Ns            # arb
        self.dl = self.w / Nl             # arb
        self.dsf, self.dlf = float(self.ds.mid()), float(self.dl.mid())
        self.sf = np.arange(Ns + 1) * self.dsf      # float grid (only for the DP; the check uses arb nodes)
        self.lf = np.arange(Nl + 1) * self.dlf
        self.arct = arct if arct is not None else ArcTangent(str(2 * float(tau)), hi_f(self.A0) / 2 + 0.01)
        # rigorous per-cell lower bounds
        self.Q = np.empty((Ns, Nl))       # lower bound of P~ over the l>0 part of cell (i,j)
        self.D = np.empty(Ns)             # lower bound of sqrt(4 pi s) on column i (l = 0 edge)
        for i in range(Ns):
            s_lo = arb(i) * self.ds
            self.D[i] = disc_low(s_lo) if i > 0 else 0.0
            for j in range(Nl):
                l_hi = arb(j + 1) * self.dl
                self.Q[i, j] = cell_Q_low(self.tau, s_lo, l_hi, self.arct)
        self.PL = self.Q.copy()           # P~ lower bound over the whole closed cell
        self.PL[:, 0] = np.minimum(self.Q[:, 0], self.D)

    def desc(self):
        return (f"tau={self.tau_s} m'=[{self.m1_s},{self.m2_s}] A0={float(self.A0.mid()):.6f} w={float(self.w.mid()):.6f} "
                f"sigma'_low={self.sig_low:.6f} Lam0={self.Lam0_s} Ns={self.Ns} Nl={self.Nl}")


def build_barrier(B, lam_c, sig_c, iters=6):
    """Float DP: minimal reachable Riccati value, pessimistic step (max vertex rho^2, min cell P).  Returns rho or None."""
    Ns, Nl, ds = B.Ns, B.Nl, B.dsf
    l = B.lf
    rho = np.empty((Ns + 1, Nl + 1))
    R = 2 * sig_c * l               # column 0: rho(0,0)=0, then jumps (r increases by 2 sig dl)
    R[0] = 0.0
    rho[0] = R
    for i in range(Ns):
        # P lower bound for the s-edge of row j: min over cells (i,j) [LR] and (i,j-1) [UL]
        PLe = np.full(Nl + 1, np.inf)
        PLc = B.PL[i] if i > 0 else B.Q[0]   # column 0: the l=0 edge of cell (0,0) is handled analytically (check)
        PLe[:Nl] = PLc
        PLe[1:] = np.minimum(PLe[1:], PLc)
        Rn = R + ds * (-lam_c - R ** 2 / PLe ** 2)
        for _ in range(iters):
            M = np.maximum(R ** 2, Rn ** 2)
            M[:Nl] = np.maximum(M[:Nl], Rn[1:] ** 2)      # LR(i,j) vertex (i+1,j+1)
            M[1:] = np.maximum(M[1:], R[:-1] ** 2)        # UL(i,j-1) vertex (i,j-1)
            Rn = R + ds * (-lam_c - M / PLe ** 2)
            Rn = np.minimum(Rn, np.minimum.accumulate(Rn - 2 * sig_c * l) + 2 * sig_c * l)
        if not np.all(np.isfinite(Rn)) or Rn.min() < -1e5:
            return None, i
        R = Rn
        rho[i + 1] = R
    return rho, Ns


def check_float(B, rho, Lam):
    """fast float screen of (b),(c): returns (worst slack b, worst slack c, #fail)"""
    Ns, Nl = B.Ns, B.Nl
    ds, dl = B.dsf, B.dlf
    a = rho[:-1, :-1]; b = rho[1:, :-1]; c = rho[1:, 1:]; d = rho[:-1, 1:]
    PL2 = B.PL ** 2
    # LR = {(i,j),(i+1,j),(i+1,j+1)}: rho_s=(b-a)/ds, rho_l=(c-b)/dl
    MLR = np.maximum(np.maximum(a ** 2, b ** 2), c ** 2)
    bLR = (b - a) / ds + Lam + MLR / PL2
    # UL = {(i,j),(i,j+1),(i+1,j+1)}: rho_l=(d-a)/dl, rho_s=(c-d)/ds
    MUL = np.maximum(np.maximum(a ** 2, d ** 2), c ** 2)
    bUL = (c - d) / ds + Lam + MUL / PL2
    # special first cell (i=0,j=0)
    Q00 = B.Q[0, 0]
    rs = (b[0, 0] - a[0, 0]) / ds
    bLR[0, 0] = rs + Lam + max(MLR[0, 0] / Q00 ** 2, rs ** 2 * ds / (4 * math.pi))
    bUL[0, 0] = (c[0, 0] - d[0, 0]) / ds + Lam + MUL[0, 0] / Q00 ** 2
    cLR = (c - b) / dl - 2 * B.sig_low
    cUL = (d - a) / dl - 2 * B.sig_low
    worst_b = max(bLR.max(), bUL.max()); worst_c = max(cLR.max(), cUL.max())
    nf = int((bLR > 0).sum() + (bUL > 0).sum() + (cLR > 0).sum() + (cUL > 0).sum())
    return worst_b, worst_c, nf


def check_arb(B, rho, Lam_s, log=None):
    """rigorous check of (a),(b),(c) on every triangle.  rho node values are the exact floats given.
    Returns (ok, n_fail, worst_b_upper, worst_c_upper)."""
    Lam = arb(Lam_s)
    Ns, Nl = B.Ns, B.Nl
    ds, dl = B.ds, B.dl
    two_sig = 2 * arb(B.sig_low)
    four_pi = 4 * PI
    nfail = 0; wb = -1e9; wc = -1e9
    if not (rho[0, 0] <= 0.0):
        return False, 1, float('inf'), float('inf')
    for i in range(Ns):
        ri = rho[i]; rn = rho[i + 1]
        for j in range(Nl):
            a, b, c, d = ri[j], rn[j], rn[j + 1], ri[j + 1]
            PL = arb(B.PL[i, j])
            # --- LR triangle
            rs = (arb(b) - arb(a)) / ds
            M = max(arb(a) ** 2, arb(b) ** 2, arb(c) ** 2, key=lambda z: float(z.upper()))
            if i == 0 and j == 0:
                # l>0 part: P~ >= Q00 (= 2 arc(2tau,0) = 4 tau);  l=0 edge: rho = rs*s (rho(0,0)=0), P^2 = 4 pi s
                Q = arb(B.Q[0, 0])
                t1 = M / Q ** 2
                t2 = rs ** 2 * ds / four_pi
                bound = arb(max(float(t1.upper()), float(t2.upper())))
                lhs = rs + Lam + bound
            else:
                lhs = rs + Lam + M / PL ** 2
            u = float(lhs.upper()); wb = max(wb, u)
            if not (u <= 0):
                nfail += 1
            rl = (arb(c) - arb(b)) / dl - two_sig
            u = float(rl.upper()); wc = max(wc, u)
            if not (u <= 0):
                nfail += 1
            # --- UL triangle
            rs = (arb(c) - arb(d)) / ds
            M = max(arb(a) ** 2, arb(d) ** 2, arb(c) ** 2, key=lambda z: float(z.upper()))
            if i == 0 and j == 0:
                lhs = rs + Lam + M / arb(B.Q[0, 0]) ** 2
            else:
                lhs = rs + Lam + M / PL ** 2
            u = float(lhs.upper()); wb = max(wb, u)
            if not (u <= 0):
                nfail += 1
            rl = (arb(d) - arb(a)) / dl - two_sig
            u = float(rl.upper()); wc = max(wc, u)
            if not (u <= 0):
                nfail += 1
    return nfail == 0, nfail, wb, wc
