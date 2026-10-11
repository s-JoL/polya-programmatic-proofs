# certT library: barrier certificate for constant (T) (Lemma T of certA-cp/report2.md).
#
# Claim certified for one parameter box (tau, A0, w, sigma'_low, Lam0):
#   for every nondecreasing path l: [0,A0] -> [0,w] and every q in H^1(0,A0) with q(A0) = 0,
#     F_l(q) = int P~(l,s)^2 q'^2 ds + 2 sigma' int q^2 dl - Lam0 int q^2 ds >= 0,
#   P~(0,s) = sqrt(4 pi s);  P~(l,s) = max( 2 arc(2tau, max(0,(s-2tau l)/2)), sqrt(4 pi s) - 2l )  (l > 0).
# Barrier criterion (Picone, measure form): a continuous P1 function rho on [0,A0]x[0,w] with
#   (a) rho(0,0) <= 0,  (b) rho_s <= -Lam0 - rho^2/P~^2,  (c) rho_l <= 2 sigma'   on every triangle.
# Proof sketch (report.md section 2): g(s) = rho(s,l(s)) satisfies dg <= -(Lam0 + g^2/P~^2) ds + 2 sigma' dl,
#   F >= int (P~ q' - g q/P~)^2 + int q^2 [ -dg - (Lam0 + g^2/P~^2)ds + 2 sigma' dl ] + [g q^2]_{0-}^{A0} >= -rho(0,0) q(0)^2 >= 0.
import math
import numpy as np
from flint import arb, ctx

ctx.prec = 100
PI = arb.pi()


def arbq(x):
    """exact arb from a python float or string"""
    return arb(x)


# ---------------------------------------------------------------- sigma' (arb lower bound), same formula as star-cert engine
def sigma_prime_arb(mprime, a='0.2595', lam='37.7'):
    j0 = arb('2.404825557695772768621631879326454643124244909145')
    j0 = arb(j0.mid(), arb('1e-45'))
    kappa = PI * j0 ** 2
    al = arb(a) - arb(mprime)
    s1 = arb(3) / 2 * (PI / arb(3).sqrt()).sqrt()
    nu = kappa / (2 * al)
    return s1 / al.sqrt() * (1 - arb(lam) / nu)


def sigma_prime_float(mprime, a=0.2595, lam=37.7):
    kappa = math.pi * 2.404825557695773 ** 2
    al = a - mprime
    s1 = 1.5 * math.sqrt(math.pi / math.sqrt(3))
    return s1 / math.sqrt(al) * (1 - lam / (kappa / (2 * al)))


# ---------------------------------------------------------------- arc(L, A): circular arc on chord L enclosing area A
def _A_of_phi(L, ph):
    h = L / 2
    return h * h * (ph - math.sin(ph) * math.cos(ph)) / math.sin(ph) ** 2


def arc_float(L, A):
    if A <= 0:
        return L
    lo, hi = 1e-12, math.pi / 2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _A_of_phi(L, mid) < A:
            lo = mid
        else:
            hi = mid
    ph = 0.5 * (lo + hi)
    return L * ph / math.sin(ph)


def lo_f(x):
    """float that is <= every point of the arb ball x (directed rounding)"""
    return float(np.nextafter(float(x.lower()), -np.inf))


def hi_f(x):
    """float that is >= every point of the arb ball x"""
    return float(np.nextafter(float(x.upper()), np.inf))


class ArcTangent:
    """Rigorous lower bound of arc(L, A) by tangent lines.  arc is convex in A on [0, pi L^2/8]
    (d arc/dA = 1/R, R = (L/2)/sin(phi) decreasing in A), hence for every table point phi_k:
        arc(A) >= arc(phi_k) + (sin(phi_k)/(L/2)) * (A - A(phi_k))      for all A in [0, pi L^2/8].
    Table points phi_k are floats; arc(phi_k), A(phi_k), slope_k are arb enclosures."""
    def __init__(self, L_str, Amax, n=4000):
        self.Lb = arb(L_str)
        L = float(self.Lb.mid())
        h = self.Lb / 2
        self.Lf_low = lo_f(self.Lb)
        self.Ak = np.linspace(0.0, Amax, n + 1)
        self.phi = []; self.Aphi = []; self.arcphi = []; self.slope = []; self.Aphi_up = np.empty(n + 1)
        for k in range(n + 1):
            A = self.Ak[k]
            if k == 0:
                self.phi.append(None); self.Aphi.append(arb(0)); self.arcphi.append(self.Lb); self.slope.append(arb(0))
                self.Aphi_up[0] = 0.0
                continue
            lo, hi = 1e-12, math.pi / 2
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if _A_of_phi(L, mid) < A:
                    lo = mid
                else:
                    hi = mid
            p = arb(lo)
            Ap = h * h * (p - p.sin() * p.cos()) / p.sin() ** 2
            self.phi.append(lo); self.Aphi.append(Ap); self.arcphi.append(self.Lb * p / p.sin())
            self.slope.append(p.sin() / h)
            self.Aphi_up[k] = hi_f(Ap)
        assert self.Aphi_up[-1] <= math.pi * L * L / 8

    def lower_arb(self, A):
        """rigorous lower bound (float) of arc(L, max(0, A)) for the arb A (A may straddle 0)."""
        if float(A.upper()) <= 0:
            return self.Lf_low
        k = int(np.searchsorted(self.Aphi_up, max(float(A.mid()), 0.0), side='right') - 1)
        k = max(k, 0)
        if float(A.lower()) < 0:
            return self.Lf_low                            # arc(max(0,A)) >= L always
        val = self.arcphi[k] + self.slope[k] * (A - self.Aphi[k])
        return max(lo_f(val), self.Lf_low)


class ArcTable:
    """Rigorous lower bounds of arc(L, A) on a grid A_k (arc is increasing in A).  lower(A) uses the largest A_k <= A.
    L is given as a decimal string (e.g. '0.9' = 2 tau); all enclosures are valid for the exact decimal value."""
    def __init__(self, L_str, Amax, n=4000):
        Lb = arb(L_str)
        L = float(Lb.mid())
        self.Lb = Lb
        self.L = float(Lb.lower())
        self.Ak = np.linspace(0.0, Amax, n + 1)
        self.low = np.empty(n + 1)
        h = Lb / 2
        self.low[0] = self.L                 # arc(L, 0) = L (degenerate arc = chord)
        for k in range(1, n + 1):
            A = self.Ak[k]
            # find phi_lo with A(phi_lo) <= A rigorously
            lo, hi = 1e-12, math.pi / 2
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if _A_of_phi(L, mid) < A:
                    lo = mid
                else:
                    hi = mid
            ph = lo * (1 - 1e-12)
            phb = arb(ph)
            Aph = h * h * (phb - phb.sin() * phb.cos()) / phb.sin() ** 2
            while not (Aph.upper() <= arb(A).lower()):
                ph *= (1 - 1e-9)
                phb = arb(ph)
                Aph = h * h * (phb - phb.sin() * phb.cos()) / phb.sin() ** 2
            arcl = Lb * phb / phb.sin()
            self.low[k] = float(arcl.lower())
        # sanity: arc is >= L and increasing
        assert np.all(np.diff(self.low) >= -1e-15)

    def lower(self, A):
        """float lower bound (rigorous) for arc(L, A') for all A' >= A (A float lower bound)."""
        if A <= 0:
            return self.L                    # float(arb(L).lower())
        k = np.searchsorted(self.Ak, A, side='right') - 1
        return self.low[max(k, 0)]


# ---------------------------------------------------------------- per-cell lower bound of P~ (float, rigorous via arb inside)
def cell_Q_low(tau, s_lo, l_hi, arct):
    """rigorous float lower bound of P~(l,s) over points of the cell with l > 0:
    P~ is increasing in s and decreasing in l, so evaluate at (s_lo, l_hi).  tau, s_lo, l_hi are arb.
    arct: ArcTangent (tangent-line lower bounds, directed rounding)."""
    A = (s_lo - 2 * tau * l_hi) / 2
    q1 = 2 * arct.lower_arb(A)                          # exact doubling of a rigorous float lower bound
    q2 = (4 * PI * s_lo).sqrt() - 2 * l_hi
    return max(q1, lo_f(q2))


def disc_low(s_lo):
    """rigorous float lower bound of sqrt(4 pi s) for s >= s_lo (arb)"""
    return lo_f((4 * PI * s_lo).sqrt())


# ---------------------------------------------------------------- float profile on the grid (for the DP)
def P_float(tau, s, l, arcL):
    if l <= 0:
        return math.sqrt(4 * math.pi * s)
    A = max(0.0, (s - 2 * tau * l) / 2)
    return max(2 * arcL(A), math.sqrt(4 * math.pi * s) - 2 * l)
