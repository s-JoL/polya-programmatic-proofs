"""A1Dc cert_common: rigorous constants (python-flint arb) for the tau=0.45 one-dimensional certificates.
All parameters are exact rationals; transcendental constants are arb balls.
"""
import numpy as np
from flint import arb, ctx, fmpq
ctx.prec = 128

PI = arb.pi()
LAM = 12 * PI                                   # Lambda = 12 pi
# j_{0,1}: sign change of J0 on a short interval (J0 > 0 on (0, j01), standard)
_jlo, _jhi = arb('2.40482555769577276'), arb('2.40482555769577277')
assert _jlo.bessel_j(0) > 0 and _jhi.bessel_j(0) < 0
J01 = _jlo.union(_jhi)
KAPPA = PI * J01 ** 2                           # kappa = pi j01^2
SIG1 = arb(3) / 2 * (PI / arb(3).sqrt()).sqrt() # sigma_1 = (3/2)(pi/sqrt 3)^{1/2}  (closed form, theorem)

def Q(x):
    """rational (string 'p/q' or exact decimal) -> arb ball enclosing it (exact when dyadic)"""
    if isinstance(x, str):
        if '/' in x:
            p, q = x.split('/')
            return arb(fmpq(int(p), int(q)))
        if 'e' not in x.lower():      # exact decimal
            neg = x.startswith('-'); x = x.lstrip('+-')
            ip, _, fp = x.partition('.')
            num = int(ip or '0') * 10 ** len(fp) + int(fp or '0')
            v = arb(fmpq(num, 10 ** len(fp)))
            return -v if neg else v
        return arb(x)
    return arb(x)

A_CAP = Q('2595/10000')                         # a = 0.2595
TAU = Q('45/100')

def sigp(mp, a=A_CAP):
    """sigma'(m') = sigma_1/sqrt(a-m') * (1 - 2 Lambda (a-m')/kappa)  (increasing in m')"""
    al = a - mp
    return SIG1 / al.sqrt() * (1 - 2 * LAM * al / KAPPA)

def up(x):
    """float certified >= the upper endpoint of the arb ball x (directed rounding)"""
    U = x.upper()
    f = float(U)
    while arb(f) < U:
        f = float(np.nextafter(f, np.inf))
    return f

def lo(x):
    """float certified <= the lower endpoint of the arb ball x (directed rounding)"""
    L = x.lower()
    f = float(L)
    while arb(f) > L:
        f = float(np.nextafter(f, -np.inf))
    return f

if __name__ == '__main__':
    print('12 pi  =', LAM)
    print('j01    =', J01)
    print('kappa  =', KAPPA)
    print('sigma1 =', SIG1)
    print('a - kappa/(2 Lambda) (sigma\' > 0 iff m\' above this) =', A_CAP - KAPPA / (2 * LAM))
    print('w* tau = 0.240 * 0.45 =', Q('0.240') * TAU)
    for m in ('0.03', '0.108', '0.15', '0.18', '0.20', '0.21', '0.22', '0.224', '0.2265', '0.2282'):
        print('sigma\'(%s) =' % m, sigp(Q(m)))
