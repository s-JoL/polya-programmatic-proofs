"""A1Dc cert_riccati: interval (arb) certificates for the three Talenti-type 1D constants at tau = 0.45.

(W)  Lemma W:  T(w; tau, A0) > lam_c.   1D problem  inf int I^2 g'^2 / int g^2,  g(A0)=0,  I = min(sqrt(4 pi s), P_lens(2tau,s) - 2w).
(N)  neck Talenti:  lam_1D(I_neck, 2a) > lam_c,  I_neck = max(sqrt(4 pi s), 2(sqrt(4 pi (s/2-m0)_+) - m0/tau)).
     Both by ground-state substitution with a Riccati SUPERsolution  zeta' >= lam + zeta^2/I^2,  zeta(0)=0, finite on [0,A0]:
     phi = exp(-int zeta/I^2) > 0 gives int I^2 g'^2 >= lam int g^2 (upward jumps of zeta allowed).
     On [0,s0] (where the profile used is sqrt(4 pi s) <= I) zeta = 2 pi X J1(X)/J0(X), X = sqrt(lam s/pi) (disc, exact);
     on each cell [s_i,s_{i+1}] with k = rigorous lower bound of I on the cell (I is increasing in s):
         zeta(s) = k sqrt(lam) tan( sqrt(lam)(s-s_i)/k + arctan(zeta_i/(k sqrt(lam))) )   (exact solution with I:=k).
(C)  Lemma C:  psi'' = -F(s,psi,psi'), F = psi f(s, psi'/psi),  f(s,r) = max{ lam/(4 pi s), max_{0<=l<=2w} (lam - sig r l)/(P(s)-l)^2 };
     need psi' > 0 on [0,A0].  Riccati r = psi'/psi:  r' = -r^2 - f(s,r).  On [0,s0] (4 pi s <= (P-2w)^2) psi is the Bessel solution,
     r = J0(X)/((2 pi/lam) X J1(X)).  On each cell, f(s,r) <= K := f-hat(s_i, P_lo(s_i), r_lo) for r >= r_lo
     (f nonincreasing in r and in P, P increasing in s, lam/(4 pi s) decreasing), and the SUBsolution
         r(s) = sqrt K tan( arctan(r_i/sqrt K) - sqrt K (s - s_i) )
     is valid while r >= r_lo.  r(A0) > 0  ==>  psi' > 0 on [0, A0]  ==>  lambda_red > lam (Lemma C, with sigma'_lo, A0 >= |D|, l_max = 2w).
usage:
  cert_riccati.py W  w  A0  lam_c  N
  cert_riccati.py N  m0 tau a lam_c N
  cert_riccati.py C  w  A0  m'_lo  lam_c  N          (sigma' = sigma'(m'_lo))
  (w, A0, m0, ... as exact decimal strings; lam_c decimal string)
"""
import sys, math, time
import numpy as np
from flint import arb, ctx
from cert_common import PI, LAM, KAPPA, SIG1, Q, sigp, up, lo
ctx.prec = 128

# ---------- lens perimeter lower bound ----------
def s_of_t(t, a):
    return 2 * a * a * (t - math.sin(t) * math.cos(t)) / math.sin(t) ** 2

def t_of_s_float(s, a):
    s = np.asarray(s, dtype=float)
    lo_, hi_ = np.full(s.shape, 1e-9), np.full(s.shape, math.pi / 2)
    for _ in range(80):
        mid = 0.5 * (lo_ + hi_)
        val = 2 * a * a * (mid - np.sin(mid) * np.cos(mid)) / np.sin(mid) ** 2 - s
        lo_ = np.where(val < 0, mid, lo_); hi_ = np.where(val < 0, hi_, mid)
    return lo_

def Plens_lower(s_list, a_arb, a_float):
    """for each exact float s, an arb lower bound of P_lens(2a, s) (lens case s < pi a^2), via t_lo with s(t_lo) < s."""
    ts = t_of_s_float(s_list, a_float)
    out = []
    for s, t in zip(s_list, ts):
        tl = t * (1 - 1e-12)
        for _ in range(20):
            T = arb(tl)
            sn = T.sin()
            st = 2 * a_arb * a_arb * (T - sn * T.cos()) / (sn * sn)
            if st < arb(s):
                break
            tl *= (1 - 1e-10)
        else:
            raise RuntimeError('P_lens bracket failed')
        out.append((4 * a_arb * T / sn).lower())
    return out

def grid(s0f, A0, N):
    """float grid s0f = s_0 < ... < s_N with exact arb(s_N) >= A0 (A0 arb exact rational)."""
    A0f = float(A0.upper())
    g = list(np.linspace(s0f, A0f, N + 1))
    g[-1] = np.nextafter(A0f, 2.0)
    while not (arb(g[-1]) >= A0):
        g[-1] = np.nextafter(g[-1], 2.0)
    return g

# ---------- (W) and (N): linear Riccati supersolution ----------
def bessel_zeta(lam, s0):
    X = (lam * arb(s0) / PI).sqrt()
    return 2 * PI * X * X.bessel_j(1) / X.bessel_j(0), X

def run_linear(Ilow_fun, s0f, A0, lam, N, label):
    lam = Q(lam) if isinstance(lam, str) else lam
    z, X = bessel_zeta(lam, s0f)
    assert X < arb('2.4048255576957727'), 'Bessel start beyond j01'
    assert z.is_finite()
    sg = grid(s0f, A0, N)
    Ilo = Ilow_fun(sg)
    sl = lam.sqrt()
    zmax = z
    for i in range(N):
        k = Ilo[i]
        if not (k > 0):
            return False, 'I_lo <= 0 at s=%.6f' % sg[i], None
        h = arb(sg[i + 1]) - arb(sg[i])
        arg = sl * h / k + (z / (k * sl)).atan()
        if not (arg < PI / 2):
            return False, 'blow-up at s=%.6f' % sg[i], None
        z = (k * sl * arg.tan()).upper()     # upward jump allowed
        z = arb(z)
    return True, 'zeta(A0) <= %.6g' % up(z), up(z)

def Ilow_W(w, tau):
    a_arb = tau; a_f = float(tau.mid())
    def f(sg):
        P = Plens_lower(sg, a_arb, a_f)
        out = []
        for s, p in zip(sg, P):
            c1 = (4 * PI * arb(s)).sqrt()
            c2 = p - 2 * w
            out.append(c1.min(c2).lower())
        return out
    return f

def Ilow_N(m0, tau):
    def f(sg):
        out = []
        for s in sg:
            S = arb(s)
            c1 = (4 * PI * S).sqrt()
            t = S / 2 - m0
            if t > 0:
                c2 = 2 * ((4 * PI * t).sqrt() - m0 / tau)
                out.append(c1.max(c2).lower())
            else:
                out.append(c1.lower())
        return out
    return f

def cert_W(w, A0, lam, N):
    w, A0 = Q(w), Q(A0)
    tau = Q('45/100')
    # on [0,s0] the profile must EQUAL sqrt(4 pi s): need sqrt(4 pi s) <= P(s) - 2w there.  P(s) >= 4 tau (t -> 0 limit,
    # P increasing in s), so s0 with 4 pi s0 <= (4 tau - 2w)^2 suffices.
    s0 = float(((4 * tau - 2 * w) ** 2 / (4 * PI)).lower()) * (1 - 1e-12)
    s0 = min(s0, float(A0.lower()) * 0.5)
    assert 4 * PI * arb(s0) < (4 * tau - 2 * w) ** 2
    ok, msg, zA = run_linear(Ilow_W(w, tau), s0, A0, lam, N, 'W')
    return ok, msg, s0

def cert_N(m0, tau, a, lam, N):
    m0, tau, a = Q(m0), Q(tau), Q(a)
    A0 = 2 * a
    s0 = float((2 * m0).lower())   # for s <= 2 m0 the second branch is negative; Bessel with sqrt(4 pi s) <= I is a supersolution anyway
    ok, msg, zA = run_linear(Ilow_N(m0, tau), s0, A0, lam, N, 'N')
    return ok, msg, s0

# ---------- (C): nonlinear Riccati subsolution ----------
def gstar_upper(lam, sig, P, r, w2):
    """upper bound of max_{0<=l<=w2} (lam - sig r l)/(P - l)^2 for r >= r (arb lower value), P = arb lower bound of P(s)."""
    a = sig * r
    g0 = lam / (P * P)
    g2 = (lam - a * w2) / ((P - w2) ** 2)
    if not (a > 0):
        # a could be ~0: max over l of lam/(P-l)^2 at l = w2
        return lam / ((P - w2) ** 2)
    ls = 2 * lam / a - P
    if ls < 0:
        return g0
    if ls > w2:
        return g2 if g2 > g0 else arb(max(up(g0), up(g2)))
    d = a * P - lam
    if d > 0:
        gs = a * a / (4 * d)
        return arb(max(up(gs), up(g0), up(g2)))
    return lam / ((P - w2) ** 2)          # crude but valid (numerator <= lam, denominator >= (P-w2)^2)

def cert_C(w, A0, mlo, lam, N, verbose=False):
    w, A0, mlo = Q(w), Q(A0), Q(mlo)
    lam = Q(lam) if isinstance(lam, str) else lam
    sig = sigp(mlo).lower()
    tau = Q('45/100')
    w2 = 2 * w
    # s0: 4 pi s0 <= (P_lo(s0) - 2w)^2, and P(s) >= P(s0)?? no: on [0,s0] need 4 pi s <= (P(s)-2w)^2 for all s<=s0:
    # P(s) >= 4 tau (t->0 limit, P increasing), so it suffices 4 pi s0 <= (4 tau - 2w)^2.
    s0 = float(((4 * tau - w2) ** 2 / (4 * PI)).lower()) * (1 - 1e-12)
    assert 4 * PI * arb(s0) <= (4 * tau - w2) ** 2
    s0 = min(s0, float(A0.lower()) * 0.5)
    X = (lam * arb(s0) / PI).sqrt()
    r = (X.bessel_j(0) / ((2 * PI / lam) * X * X.bessel_j(1))).lower()
    assert r > 0
    sg = grid(s0, A0, N)
    Plo = Plens_lower(sg, tau, 0.45)
    rmin_seen = None
    for i in range(N):
        S = arb(sg[i]); h = arb(sg[i + 1]) - S
        t0 = lam / (4 * PI * S)
        # r_lo: start optimistic (r), then just below the computed endpoint
        rl = r
        for it in range(30):
            if not (rl > 0):
                rl = arb(0)
            g = gstar_upper(lam, sig, Plo[i], rl, w2)
            K = t0.max(g)
            sK = K.sqrt()
            arg = (r / sK).atan() - sK * h
            if not (arg > 0):
                return False, 'psi\' may vanish in cell at s=%.6f (r=%.4g)' % (sg[i], lo(r)), sg[i]
            rn = (sK * arg.tan()).lower()
            if rn >= rl:
                break
            rl = (rn - (r - rn) * arb('0.02')).lower()
        else:
            return False, 'r_lo iteration failed at s=%.6f' % sg[i], sg[i]
        r = arb(rn)
        if verbose and i % (N // 10) == 0:
            print('   s=%.5f r>=%.6g K=%.5g' % (sg[i], lo(r), up(K)), flush=True)
    return (r > 0), 'r(A0) >= %.6g' % lo(r), lo(r)

if __name__ == '__main__':
    t0 = time.time()
    kind = sys.argv[1]
    if kind == 'W':
        w, A0, lam, N = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
        ok, msg, s0 = cert_W(w, A0, lam, N)
        print('[W] w=%s A0=%s lam_c=%s N=%d s0=%.5f : %s  %s  (%.1fs)' % (w, A0, lam, N, s0, 'CERTIFIED' if ok else 'FAILED', msg, time.time() - t0))
    elif kind == 'N':
        m0, tau, a, lam, N = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])
        ok, msg, s0 = cert_N(m0, tau, a, lam, N)
        print('[N] m0=%s tau=%s a=%s lam_c=%s N=%d : %s  %s  (%.1fs)' % (m0, tau, a, lam, N, 'CERTIFIED' if ok else 'FAILED', msg, time.time() - t0))
    elif kind == 'C':
        w, A0, mlo, lam, N = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])
        ok, msg, x = cert_C(w, A0, mlo, lam, N, verbose=('-v' in sys.argv))
        print('[C] w=%s A0=%s sigma\'=sigma\'(%s)=%.6f lam_c=%s N=%d : %s  %s  (%.1fs)' % (w, A0, mlo, lo(sigp(Q(mlo))), lam, N, 'CERTIFIED' if ok else 'FAILED', msg, time.time() - t0))
