# certT, Pruefer-angle form of the barrier.
# Barrier rho = K cot(Theta), Theta continuous P1 on the triangulated grid, 0 < Theta < pi.
# Conditions (equivalent to (a),(b),(c) for rho, since rho_s = -K Theta_s / sin^2, rho_l = -K Theta_l / sin^2):
#   (a') Theta(0,0) = pi/2 exactly (node (0,0) is treated as the exact real number pi/2)  => rho(0,0) = 0
#   (b') Theta_s >= sup_T [ Lam0 sin^2 Theta + K^2 cos^2 Theta / P~^2 ] / K
#   (c') Theta_l >= - 2 sigma' inf_T sin^2 Theta / K
#   (d') sup Theta < pi, inf Theta > 0
# First cell (0,0): on its l=0 edge, |cos Theta| <= |Theta - pi/2| = Theta_s * s and P~^2 = 4 pi s,
#   so the cos-term is <= K Theta_s^2 Delta_s/(4 pi); on its l>0 part P~ >= Q00 = 4 tau.
import sys, math, time
import numpy as np
from flint import arb, ctx
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from ctlib import PI, lo_f, hi_f
from ctcert import Box

ctx.prec = 100
HALFPI_F = math.pi / 2


def sin2_range(ta, tb):
    lo, hi = min(ta, tb), max(ta, tb)
    smax = 1.0 if lo <= HALFPI_F <= hi else max(math.sin(lo) ** 2, math.sin(hi) ** 2)
    smin = min(math.sin(lo) ** 2, math.sin(hi) ** 2)
    cmax = max(math.cos(lo) ** 2, math.cos(hi) ** 2)
    return smin, smax, cmax


def build_theta(B, lam_c, sig_c, K=20.0, iters=6, slack_s=0.0):
    """Float DP for the maximal reachable Pruefer angle (pessimistic over each triangle).  Returns Theta or (None, col)."""
    Ns, Nl, ds, dl = B.Ns, B.Nl, B.dsf, B.dlf
    Th = np.empty((Ns + 1, Nl + 1))
    T = np.empty(Nl + 1)
    T[0] = HALFPI_F
    jump = 2 * sig_c * dl / K
    for j in range(1, Nl + 1):
        T[j] = math.atan2(1.0, 1.0 / math.tan(T[j - 1]) + jump) if T[j - 1] != 0 else 0.0   # arccot(cot + jump)
    Th[0] = T
    for i in range(Ns):
        PLc = B.PL[i] if i > 0 else B.Q[0]
        PLe = np.full(Nl + 1, np.inf); PLe[:Nl] = PLc; PLe[1:] = np.minimum(PLe[1:], PLc)
        Tn = T + ds * (lam_c * np.sin(T) ** 2 + K * K * np.cos(T) ** 2 / PLe ** 2) / K
        for _ in range(iters):
            # vertex sets V_j = {(i,j),(i+1,j),(i+1,j+1),(i,j-1)}: theta range
            lo = np.minimum(T, Tn); hi = np.maximum(T, Tn)
            lo[:Nl] = np.minimum(lo[:Nl], Tn[1:]); hi[:Nl] = np.maximum(hi[:Nl], Tn[1:])
            lo[1:] = np.minimum(lo[1:], T[:-1]); hi[1:] = np.maximum(hi[1:], T[:-1])
            smax = np.where((lo <= HALFPI_F) & (hi >= HALFPI_F), 1.0, np.maximum(np.sin(lo) ** 2, np.sin(hi) ** 2))
            cmax = np.maximum(np.cos(lo) ** 2, np.cos(hi) ** 2)
            Tn = T + ds * ((lam_c * smax + K * K * cmax / PLe ** 2) / K + slack_s)
            # (c') along the vertical edge (i+1,j)->(i+1,j+1): Th_{j+1} >= Th_j - jump * inf sin^2 over the two triangles
            # using that edge: LR(i,j) [extra vertex (i,j)] and UL(i+1,j) [extra vertex (i+2,j+1), predicted by the
            # maximal one-step growth].  Raising theta is always admissible for the DP (max reachable angle).
            grow = ds * (lam_c + K * K / np.maximum(PLe, 1e-12) ** 2) / K
            for j in range(Nl):
                t0 = Tn[j]
                for _r in range(3):
                    t1 = Tn[j + 1]
                    pred = min(t1 + grow[j + 1], math.pi - 1e-12)
                    cands = (t0, t1, pred, T[j])
                    s2min = min(math.sin(x) ** 2 for x in cands)
                    need = t0 - jump * s2min
                    if Tn[j + 1] < need:
                        Tn[j + 1] = need
                    else:
                        break
        if not np.all(np.isfinite(Tn)) or Tn.max() >= math.pi - 1e-9:
            return None, i
        T = Tn
        Th[i + 1] = T
    return Th, Ns


def check_theta_float(B, Th, Lam, K=20.0):
    Ns, Nl, ds, dl = B.Ns, B.Nl, B.dsf, B.dlf
    a = Th[:-1, :-1]; b = Th[1:, :-1]; c = Th[1:, 1:]; d = Th[:-1, 1:]
    PL2 = B.PL ** 2
    out = {}
    for name, (v1, v2, v3, s_from, s_to, l_from, l_to) in {
            'LR': (a, b, c, a, b, b, c), 'UL': (a, d, c, d, c, a, d)}.items():
        lo = np.minimum(np.minimum(v1, v2), v3); hi = np.maximum(np.maximum(v1, v2), v3)
        smax = np.where((lo <= HALFPI_F) & (hi >= HALFPI_F), 1.0, np.maximum(np.sin(lo) ** 2, np.sin(hi) ** 2))
        smin = np.minimum(np.sin(lo) ** 2, np.sin(hi) ** 2)
        cmax = np.maximum(np.cos(lo) ** 2, np.cos(hi) ** 2)
        Ts = (s_to - s_from) / ds; Tl = (l_to - l_from) / dl
        with np.errstate(divide='ignore', invalid='ignore'):
            rhs = (Lam * smax + K * K * cmax / PL2) / K
        if name == 'LR':
            rhs[0, 0] = (Lam * smax[0, 0] + K * K * max(cmax[0, 0] / B.Q[0, 0] ** 2, Ts[0, 0] ** 2 * ds / (4 * math.pi))) / K
        else:
            rhs[0, 0] = (Lam * smax[0, 0] + K * K * cmax[0, 0] / B.Q[0, 0] ** 2) / K
        out[name + 'b'] = Ts - rhs                     # must be >= 0
        out[name + 'c'] = Tl + 2 * B.sig_low * smin / K  # must be >= 0
    wb = min(out['LRb'].min(), out['ULb'].min()); wc = min(out['LRc'].min(), out['ULc'].min())
    nf = int(sum((v < 0).sum() for v in out.values()))
    return wb, wc, nf, float(Th.max())


def check_theta_arb(B, Th, Lam_s, K_s='20'):
    """Rigorous check of (a')-(d') on every triangle.  Node (0,0) is the exact pi/2; other nodes are the given floats."""
    Ns, Nl = B.Ns, B.Nl
    Lam = arb(Lam_s); K = arb(K_s); K2 = K * K
    ds, dl = B.ds, B.dl
    two_sig_over_K = 2 * arb(B.sig_low) / K
    four_pi = 4 * PI
    # node enclosures: theta, sin^2, cos^2
    t0 = time.time()
    TA = [[None] * (Nl + 1) for _ in range(Ns + 1)]
    S2 = [[None] * (Nl + 1) for _ in range(Ns + 1)]
    C2 = [[None] * (Nl + 1) for _ in range(Ns + 1)]
    for i in range(Ns + 1):
        for j in range(Nl + 1):
            t = PI / 2 if (i == 0 and j == 0) else arb(float(Th[i, j]))
            TA[i][j] = t
            s, c = t.sin(), t.cos()
            S2[i][j] = s * s; C2[i][j] = c * c
    # (d'): range
    tmax = max(hi_f(TA[i][j]) for i in range(Ns + 1) for j in range(Nl + 1))
    tmin = min(lo_f(TA[i][j]) for i in range(Ns + 1) for j in range(Nl + 1))
    if not (tmax < math.pi * (1 - 1e-12) and tmin > 0):
        return False, -1, None, None, tmax
    halfpi = PI / 2
    hp_lo, hp_hi = lo_f(halfpi), hi_f(halfpi)
    nfail = 0; wb = 1e9; wc = 1e9
    ntri = 0
    for i in range(Ns):
        for j in range(Nl):
            PL = arb(B.PL[i, j]) if not (i == 0 and j == 0) else None
            tris = (('LR', (i, j), (i + 1, j), (i + 1, j + 1), ((i, j), (i + 1, j)), ((i + 1, j), (i + 1, j + 1))),
                    ('UL', (i, j), (i, j + 1), (i + 1, j + 1), ((i, j + 1), (i + 1, j + 1)), ((i, j), (i, j + 1))))
            for name, p1, p2, p3, sedge, ledge in tris:
                vs = (p1, p2, p3)
                ths = [TA[p[0]][p[1]] for p in vs]
                # sin^2 max: 1 if pi/2 may lie in [min, max] theta of the triangle (theta is linear on T)
                lo = min(lo_f(t) for t in ths); hi = max(hi_f(t) for t in ths)
                s2s = [S2[p[0]][p[1]] for p in vs]; c2s = [C2[p[0]][p[1]] for p in vs]
                if lo <= hp_hi and hi >= hp_lo:
                    smax = arb(1)
                else:
                    smax = arb(max(hi_f(z) for z in s2s))
                smin_low = min(lo_f(z) for z in s2s)   # sin^2 unimodal on (0,pi): min over T at a vertex
                cmax = arb(max(hi_f(z) for z in c2s))
                Ts = (TA[sedge[1][0]][sedge[1][1]] - TA[sedge[0][0]][sedge[0][1]]) / ds
                Tl = (TA[ledge[1][0]][ledge[1][1]] - TA[ledge[0][0]][ledge[0][1]]) / dl
                if i == 0 and j == 0:
                    Q = arb(B.Q[0, 0])
                    if name == 'LR':
                        ccos = cmax / (Q * Q)
                        cedge = Ts * Ts * ds / four_pi
                        cterm = arb(max(hi_f(ccos), hi_f(cedge)))
                    else:
                        cterm = cmax / (Q * Q)
                    rhs = (Lam * smax + K2 * cterm) / K
                else:
                    rhs = (Lam * smax + K2 * cmax / (PL * PL)) / K
                ntri += 1
                gb = Ts - rhs
                gc = Tl + two_sig_over_K * arb(smin_low)
                vb = lo_f(gb); vc = lo_f(gc)
                wb = min(wb, vb); wc = min(wc, vc)
                if not (gb.lower() >= 0):
                    nfail += 1
                if not (gc.lower() >= 0):
                    nfail += 1
    check_theta_arb.ntri = ntri
    return nfail == 0, nfail, wb, wc, tmax
