"""A1Dc cert_slice: interval (arb) certificate of the horizontal-slice constant by the LADDER (level-crossing) argument.

Claim certified: for every measurable c on [0,tau] with c >= w and int_0^tau c <= m_hi, and every m' in [m_lo, m_hi]
(sigma' = sigma'(m') >= sigma'(m_lo)), the half-interval slice eigenvalue
    mu_h[c] = inf { int_0^tau (F'^2 + pi^2 F^2/c^2) + sigma' F(tau)^2 } / int_0^tau F^2      (Neumann at 0)
satisfies mu_h[c] > L0 (= 37.7 > 12 pi).   No monotonicity of c is used.

Proof skeleton (author-level, see report):
  Pruefer angle th of F (F'' = (pi^2/c^2 - L0) F, F(0)=1, F'(0)=0): th' = f(th,c) = 1 + s (L0 - 1 - pi^2/c^2), s = sin^2 th, th(0)=pi/2.
  mu_h <= L0  <=>  th(tau) >= thR := pi/2 + arctan sigma'   (terminal condition; once th reaches pi it stays >= pi).
  A 'potential' Phi (piecewise linear, Phi' = -p on cells) with  g_nu(th,p) := min_{c>=w} [c - nu - p f(th,c)] >= 0  on an interval J
  gives, for any path piece inside J:  int c >= nu * (duration) + Phi(start) - Phi(end).
  g_nu(th,p) = G(p, s) := m(p, s) - nu - p (1 + (L0-1) s),  s = sin^2 th,
     m(p,s) = min_{c>=w}(c + p pi^2 s/c^2)  =  w + K/w^2 if (2K)^{1/3} <= w,  (3/2)(2K)^{1/3} otherwise (K = p pi^2 s),
     and (3/2)(2K)^{1/3} is a lower bound in all cases.  G is CONCAVE in s (min of affine functions), and s is monotone on
     each theta-cell, so  min over the cell = min(G(p, s(th_a)), G(p, s(th_b))):  ghat(p) := that minimum (exact, no slack).
     The one cell containing the true pi/2 uses s in {min(s_a,s_b), 1} instead (still exact by concavity).
  Levels thR = t_0 < ... < t_K, floor th_b.  Cases for a feasible path:
    dip  (min th < th_b):            mass >= w tau + int_{th_b}^{thR} pmax(.;w)
    top  (max th >= t_K):            mass >= A_K
    k    (max th in [t_k,t_{k+1}), min th >= th_b):   mass >= max(A_k, B_k),
         A_k = w tau + int_{pi/2}^{t_k} pmax(.;w)
         B_k = nu_k tau + int_{pi/2}^{t_k} pmax(.;nu_k) - int_{thR}^{t_k} pmin(.;nu_k)   (two potentials, before/after first hitting t_k)
       requires the feasible p-interval [pmin,pmax] of ghat(.;nu_k) to be nonempty on every cell of [th_b, t_{k+1}].
       (A_k and the top case do not need min th >= th_b: for nu = w, p = 0 is admissible everywhere.  Near th* the code may
        pick nu_k slightly below w; then B_k is merely weaker than A_k, still valid.)
  All verdicts are arb comparisons; floats are exact grid nodes / exact trial slopes; printed floats use directed rounding.
usage: cert_slice.py L0 w m_lo m_hi K [fine]       (decimal strings)
"""
import sys, math, time, json
import numpy as np
from flint import arb, ctx
from cert_common import PI, Q, sigp, up, lo
ctx.prec = 128

TAU = Q('45/100')

def G_f(p, s, nu, w, L0):
    K = p * math.pi ** 2 * s
    t = np.cbrt(2 * K)
    m = np.where(t <= w, w + K / w ** 2, 1.5 * t)
    return m - nu - p * (1 + (L0 - 1) * s)

def ghat_f(p, sa, sb, nu, w, L0):
    """float version (vectorised): min of G at the two cell ends"""
    return np.minimum(G_f(p, sa, nu, w, L0), G_f(p, sb, nu, w, L0))

def G_a(P, s, nu, w, L0):
    K = P * PI * PI * s
    t = (2 * K).root(3)
    if t < w:
        m = w + K / (w * w)
    else:
        m = arb(3) / 2 * t
    return m - nu - P * (1 + (L0 - 1) * s)

def ghat_a(p, sa, sb, nu, w, L0):
    """arb lower bound of min over the cell"""
    if p == 0:
        if nu is w:
            return arb(0)        # min_{c>=w} (c - w) = 0 exactly
        return w - nu            # min_{c>=w} (c - nu)
    P = arb(p)
    return G_a(P, sa, nu, w, L0).min(G_a(P, sb, nu, w, L0))

def feasible_interval(slo_f, shi_f, nu, w, L0):
    """float: for each cell the interval [pmin,pmax] where ghat >= 0 (nan if empty)"""
    n = len(slo_f)
    pg = np.exp(np.linspace(np.log(1e-7), np.log(50.0), 400))
    vals = ghat_f(pg[None, :], slo_f[:, None], shi_f[:, None], nu, w, L0)
    ib = np.argmax(vals, axis=1)
    vmax = vals[np.arange(n), ib]
    ppk = pg[ib]
    pmin = np.full(n, np.nan); pmax = np.full(n, np.nan)
    ok = vmax > 0
    # pmax: bisection on [ppk, 60]
    a = ppk.copy(); b = np.full(n, 60.0)
    for _ in range(70):
        mdl = 0.5 * (a + b)
        v = ghat_f(mdl, slo_f, shi_f, nu, w, L0)
        a = np.where(v >= 0, mdl, a); b = np.where(v >= 0, b, mdl)
    pmax[ok] = a[ok]
    # pmin: bisection on [0, ppk]
    a = np.zeros(n); b = ppk.copy()
    for _ in range(70):
        mdl = 0.5 * (a + b)
        v = ghat_f(mdl, slo_f, shi_f, nu, w, L0)
        a = np.where(v >= 0, a, mdl); b = np.where(v >= 0, mdl, b)
    pmin[ok] = b[ok]
    if nu <= w:            # p = 0 is always admissible when nu <= w
        pmax = np.where(ok, pmax, 0.0)
        pmin = np.zeros(n)
        ok = np.ones(n, dtype=bool)
    return pmin, pmax, ok

class Cells:
    def __init__(self, nodes):
        self.nodes = list(nodes)
        self.A = [arb(x) for x in self.nodes]
        S = [x.sin() ** 2 for x in self.A]
        # s = sin^2 is monotone on every cell not containing pi/2, so the two cell ends suffice (concavity of G in s).
        # A cell that arb cannot place on one side of the TRUE pi/2 (e.g. [fl(pi/2), .], since fl(pi/2) < pi/2 by 6e-17)
        # has s-range [min(s_a, s_b), 1]; by concavity the exact cell minimum is then min(G(min(s_a,s_b)), G(1)).
        # (revAC/S1D patch 1)
        self.slo, self.shi, self.dth = [], [], []      # the two s-values at which G is evaluated (arb balls)
        self.fixed = []
        half = PI / 2
        for i in range(len(nodes) - 1):
            a, b = self.A[i], self.A[i + 1]
            if (b <= half) or (a >= half):
                self.slo.append(S[i]); self.shi.append(S[i + 1])
            else:
                self.slo.append(S[i].min(S[i + 1])); self.shi.append(arb(1))
                self.fixed.append(i)
            self.dth.append(b - a)
        self.slo_f = np.array([float(x.mid()) for x in self.slo])
        self.shi_f = np.array([float(x.mid()) for x in self.shi])

def certified_p(cells, idx, nu, nu_a, w, w_a, L0, L0_a, which):
    """sum over cells idx of p*dth with p = pmax (which='max', lower-bounded) or pmin (which='min', upper-bounded);
    every chosen p is verified (arb) to satisfy ghat >= 0 on its cell.  Returns arb sum (lower bound for 'max',
    upper bound for 'min'), or None if some cell is infeasible."""
    if len(idx) == 0:
        return arb(0)
    slo = cells.slo_f[idx]; shi = cells.shi_f[idx]
    pmin, pmax, ok = feasible_interval(slo, shi, nu, w, L0)
    if not ok.all():
        return None
    tot = arb(0)
    for j, i in enumerate(idx):
        if which == 'max':
            p = pmax[j] * (1 - 1e-9)
        elif which == 'min':
            p = pmin[j] * (1 + 1e-9) + 1e-14 if pmin[j] > 0 else 0.0
        else:   # 'any': take the midpoint, only check feasibility
            p = 0.5 * (pmin[j] + pmax[j])
        tries = 0
        while not (ghat_a(p, cells.slo[i], cells.shi[i], nu_a, w_a, L0_a) >= 0):
            tries += 1
            if tries > 30:
                return None
            if which == 'max':
                p *= (1 - 1e-6)
            elif which == 'min':
                p = p * (1 + 1e-6) + 1e-12
            else:
                return None
        if which != 'any':
            tot += arb(p) * cells.dth[i]
    return tot

def run(L0s, ws, mlos, mhis, K, fine=1e-5, verbose=True, thb=0.25):
    t0 = time.time()
    L0_a, w_a, mhi_a = Q(L0s), Q(ws), Q(mhis)
    L0, w = float(L0s), float(ws)
    sig_lo = sigp(Q(mlos)).lower()
    thR_a = PI / 2 + sig_lo.atan()
    thR = float(thR_a.lower()); thR = np.nextafter(thR, 0)
    assert arb(thR) <= thR_a
    sstar = 1 / (math.pi ** 2 / w ** 2 - L0 + 1)
    ths = math.pi - math.asin(math.sqrt(sstar))       # c=w equilibrium (float, only used to place levels)
    ttop = ths - 2e-6
    # grid
    if not (ths > thR + 1e-4):
        raise ValueError('theta* <= thR: ladder not needed (use TF)')
    n1 = int(math.ceil((math.pi / 2 - thb) / 0.004))
    tm_ = thR - 0.15
    n2a = int(math.ceil((tm_ - math.pi / 2) / 0.0004))
    n2b = int(math.ceil((thR - tm_) / 0.00005))
    n2 = n2a + n2b
    n3 = int(math.ceil((ttop - thR) / fine))
    nodes = (list(np.linspace(thb, math.pi / 2, n1 + 1)) + list(np.linspace(math.pi / 2, tm_, n2a + 1)[1:])
             + list(np.linspace(tm_, thR, n2b + 1)[1:]) + list(np.linspace(thR, ttop, n3 + 1)[1:]))
    nodes[n1] = math.pi / 2
    cells = Cells(nodes)
    i_half, i_R = n1, n1 + n2
    i_top = len(nodes) - 1
    # levels: node indices in [i_R, i_top], denser near the top
    lev = sorted(set([i_R + int(round((i_top - i_R) * (1 - (1 - u) ** 2))) for u in np.linspace(0, 1, K + 1)]))
    if verbose:
        print('cells: %d (floor %.3f, pi/2, thR_lo=%.7f, top=%.7f), levels %d;  sigma\'_lo=%.6f  w=%s  m_hi=%s  L0=%s'
              % (len(nodes) - 1, thb, thR, ttop, len(lev), lo(sig_lo), ws, mhis, L0s), flush=True)
    rng = lambda a, b: np.arange(a, b)
    # A_k for all levels at once (nu = w): cumulative pmax(w) from pi/2
    base_w = certified_p(cells, rng(i_half, i_R), w, w_a, w, w_a, L0, L0_a, 'max')
    A = {}
    acc = base_w
    prev = i_R
    for li in lev:
        acc = acc + certified_p(cells, rng(prev, li), w, w_a, w, w_a, L0, L0_a, 'max'); prev = li
        A[li] = w_a * TAU + acc
    dip = w_a * TAU + certified_p(cells, rng(0, i_R), w, w_a, w, w_a, L0, L0_a, 'max')
    res = {'dip': lo(dip), 'top': lo(A[lev[-1]]), 'cases': [], 'pi2_cells': [(float(cells.nodes[i]), float(cells.nodes[i + 1])) for i in cells.fixed]}
    # all comparisons in arb; floats (directed rounding, revAC/S1D patch 2) only for printing / json
    worst = dip.lower().min(A[lev[-1]].lower())
    for k in range(len(lev) - 1):
        lk, lk1 = lev[k], lev[k + 1]
        # nu_k: below c0(t_{k+1}); search downward until all cells on [th_b, t_{k+1}] feasible
        s = math.sin(nodes[lk1]) ** 2
        nu = math.pi * math.sqrt(s / (1 + (L0 - 1) * s)) * (1 - 2e-4)
        Bk = None
        for attempt in range(12):
            nu_a = arb(nu)
            feas = certified_p(cells, rng(0, lk1), nu, nu_a, w, w_a, L0, L0_a, 'any')
            if feas is not None:
                pm = certified_p(cells, rng(i_half, lk), nu, nu_a, w, w_a, L0, L0_a, 'max')
                pn = certified_p(cells, rng(i_R, lk), nu, nu_a, w, w_a, L0, L0_a, 'min')
                if pm is not None and pn is not None:
                    Bk = nu_a * TAU + pm - pn
                    break
            nu = w + (nu - w) * 0.8
        Ak = A[lk]
        if Bk is None:
            val = Ak
        else:
            val = Ak.max(Bk)
        res['cases'].append(dict(k=k, t_k=nodes[lk], t_k1=nodes[lk1], A=lo(Ak), B=(lo(Bk) if Bk is not None else None), nu=nu))
        worst = worst.min(val.lower())
        if verbose and (k % 10 == 0 or k == len(lev) - 2):
            print('  case %d [%.6f,%.6f): A>=%.7f  B>=%s  nu=%.6f   (%.0fs)' % (k, nodes[lk], nodes[lk1], lo(Ak), ('%.7f' % lo(Bk)) if Bk is not None else 'n/a', nu, time.time() - t0), flush=True)
    margin = worst - mhi_a
    ok = bool(margin > 0)
    res.update(dict(L0=L0s, w=ws, m_lo=mlos, m_hi=mhis, K=K, fine=fine, thb=thb, bound=lo(worst), margin=lo(margin),
                    bound_arb=worst.str(20, radius=True), certified=ok))
    print('[SLICE] L0=%s w=%s sigma\'>=sigma\'(%s) mass<=%s : min-mass >= %.7f  dip=%.6f top=%.6f  margin=%.3e  %s  (pi/2-cells %s; %.0fs)'
          % (L0s, ws, mlos, mhis, lo(worst), res['dip'], res['top'], lo(margin), 'CERTIFIED' if ok else 'FAILED',
             res['pi2_cells'], time.time() - t0), flush=True)
    return res

if __name__ == '__main__':
    L0s, ws, mlos, mhis, K = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
    fine = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-5
    res = run(L0s, ws, mlos, mhis, K, fine)
    json.dump(res, open('cs_L%s_w%s_%s_%s_K%s.json' % (L0s, ws, mlos, mhis, K), 'w'), indent=1)
