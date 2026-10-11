"""Independent EXACT (pure Fraction) re-check of certE.py (same nodal values).
Uses: Psi_k(A; om) is increasing in k and A and decreasing in om; Psi is convex and every tangent
line of the cubic branch g(A)=2k om + A^3/om^2 at a point A_j in [k^{1/3}om, (k+1/2)^{1/3}om], as well as
T1(A)=3k^{2/3}A and T2(A)=3(k+1/2)^{2/3}A - om, is a supporting line of Psi.  With rational lower bounds
kap_i <= P_i^{1/3}, (k^{2/3})_lo, ((k+1/2)^{2/3})_lo and om <= om_hi, each line gives an affine-in-X lower
bound of F = Psi(A(X)) - 3 S(X); a sub-box passes when one such line is >= 0 at its worst corner.
usage: certE_exact.py file.npz delta eps_target"""
import sys, itertools, time
import numpy as np
from fractions import Fraction as Fr
from math import isqrt

fn, delta, eps_t = sys.argv[1], Fr(sys.argv[2]), Fr(sys.argv[3])
K = 44; SC = 2 ** K
d = np.load(fn)
wf, vid, Uf, Sf = d['w'], d['vid'], d['U'], d['S']
Wd = Fr(float(d['W'])); Lam0 = Fr(float(d['Lam0'])); n1 = int(d['n1'])
N, NS = len(Uf) - 1, len(Sf) - 1
U = [Fr(float(x)) for x in Uf]; S = [Fr(float(x)) for x in Sf]
PI_LO = Fr(314159265358979, 10 ** 14); PI_HI = Fr(314159265358980, 10 ** 14)


def sqrt_lo(x, k=64):
    if x <= 0: return Fr(0)
    return Fr(isqrt(x.numerator * x.denominator * 4 ** k), x.denominator * 2 ** k)


def sqrt_hi(x, k=64):
    if x <= 0: return Fr(0)
    return Fr(isqrt(x.numerator * x.denominator * 4 ** k) + 1, x.denominator * 2 ** k)


def icbrt(n):
    if n <= 0: return 0
    x = int(round(n ** (1 / 3))) if n.bit_length() < 1000 else 1 << ((n.bit_length() + 2) // 3)
    while x ** 3 > n: x -= 1
    while (x + 1) ** 3 <= n: x += 1
    return x


def cbrt_lo(x, k=50):
    if x <= 0: return Fr(0)
    return Fr(icbrt((x * 2 ** (3 * k)).__floor__()), 2 ** k)


S3_LO, S3_HI = sqrt_lo(Fr(3)), sqrt_hi(Fr(3))
OM_HI = sqrt_hi(Lam0) * Wd / (S3_LO * PI_LO)
ok = True
# sigma_top > sigma0 = Lam0/(3 sqrt3 pi) - 2
if not S[-1] > Lam0 / (3 * S3_LO * PI_LO) - 2: ok = False; print('FAIL sigma_top')
WI = {}
for i in range(N + 1):
    for j in range(N + 1):
        for l in range(NS + 1):
            v = Fr(float(wf[vid[i, j, l]])) - delta * (U[i] + U[j] + S[l])
            WI[i, j, l] = (v * SC).__floor__()
for j in range(N + 1):
    for l in range(NS + 1):
        WI[0, j, l] = min(WI[0, j, l], 0); WI[j, 0, l] = min(WI[j, 0, l], 0)
Wv = lambda i, j, l: Fr(WI[i, j, l], SC)
ws = Wv(n1, n1, NS)
eps_lo = 3 * S3_LO * PI_LO / Lam0 * ws
print('eps_cert >=', float(eps_lo), ' target', float(eps_t))
if not eps_lo > eps_t: ok = False; print('FAIL start')


def rng(a, b):
    F = lambda t: t * (2 - t)
    pts = [a, b] + ([Fr(1)] if a < 1 < b else [])
    return min(sqrt_lo(F(t)) for t in pts), max(sqrt_hi(F(t)) for t in pts)


XR = [rng(U[i], U[i + 1]) for i in range(N)]
nbox = [0]


def line_for(Ac, k, k23lo, kh23lo):
    """supporting line (alpha, beta) of Psi(.; OM_HI) near Ac (float)."""
    om = OM_HI
    a1f = float(k) ** (1 / 3) * float(om); a2f = (float(k) + 0.5) ** (1 / 3) * float(om)
    if Ac <= a1f:
        return 3 * k23lo, Fr(0)
    if Ac >= a2f:
        return 3 * kh23lo, -om
    Aj = Fr(Ac).limit_denominator(2 ** 40)
    if Aj ** 3 < k * om ** 3 or Aj ** 3 > (k + Fr(1, 2)) * om ** 3:
        return 3 * k23lo, Fr(0)
    return 3 * Aj * Aj / (om * om), 2 * k * om - 2 * Aj ** 3 / (om * om)


def check_box(P, kap, k, k23lo, kh23lo, b1, b2, depth):
    nbox[0] += 1
    Ac = float(kap[0]) * float(b1[0] + b1[1]) / 2 + float(kap[1]) * float(b2[0] + b2[1]) / 2
    cands = [line_for(Ac, k, k23lo, kh23lo)]
    Alo_f = float(kap[0] * b1[0] + kap[1] * b2[0]); Ahi_f = float(kap[0] * b1[1] + kap[1] * b2[1])
    cands += [line_for(Alo_f, k, k23lo, kh23lo), line_for(Ahi_f, k, k23lo, kh23lo)]
    for (al, be) in cands:
        # F >= al*(kap1 X1 + kap2 X2) + be - 3(P1 X1 + P2 X2): affine, coefficients g_i
        g1 = al * kap[0] - 3 * P[0]; g2 = al * kap[1] - 3 * P[1]
        val = be + g1 * (b1[0] if g1 >= 0 else b1[1]) + g2 * (b2[0] if g2 >= 0 else b2[1])
        if val >= 0:
            return True
    if depth >= 16:
        return False
    if b1[1] - b1[0] >= b2[1] - b2[0]:
        m = (b1[0] + b1[1]) / 2
        return check_box(P, kap, k, k23lo, kh23lo, (b1[0], m), b2, depth + 1) and check_box(P, kap, k, k23lo, kh23lo, (m, b1[1]), b2, depth + 1)
    m = (b2[0] + b2[1]) / 2
    return check_box(P, kap, k, k23lo, kh23lo, b1, (b2[0], m), depth + 1) and check_box(P, kap, k, k23lo, kh23lo, b1, (m, b2[1]), depth + 1)


t0 = time.time(); nt = 0; nfail = 0
for i in range(N):
    for j in range(N):
        for l in range(NS):
            hs = [U[i + 1] - U[i], U[j + 1] - U[j], S[l + 1] - S[l]]
            for sg in itertools.permutations(range(3)):
                v = [i, j, l]; verts = [tuple(v)]
                for ax in sg:
                    v[ax] += 1; verts.append(tuple(v))
                g = [None] * 3
                for s_, ax in enumerate(sg):
                    g[ax] = (Wv(*verts[s_ + 1]) - Wv(*verts[s_])) / hs[ax]
                k = -g[2]; P = [g[0] + k, g[1] + k]
                nt += 1
                if k < 0 or P[0] < 0 or P[1] < 0:
                    ok = False; nfail += 1; continue
                kap = [cbrt_lo(P[0]), cbrt_lo(P[1])]
                if not check_box(P, kap, k, cbrt_lo(k * k), cbrt_lo((k + Fr(1, 2)) ** 2), XR[i], XR[j], 0):
                    ok = False; nfail += 1
                    if nfail < 10: print('fail', (i, j, l), sg, float(P[0]), float(P[1]), float(k))
    if i % 8 == 7: print('  i', i, 'tets', nt, 'boxes', nbox[0], 'fails', nfail, f'{time.time()-t0:.0f}s', flush=True)
print('tets', nt, 'boxes', nbox[0], 'fails', nfail)
print('EXACT-CERTIFIED' if ok else 'FAILED', ' eps_cert >=', float(eps_lo))
