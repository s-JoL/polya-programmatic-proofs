"""A1Dc cert_slices2: constant (S) for the revised route (coordinator notice 2 + revGT): tau = 0.45, m = 0.2282, w* = 0.240.
For every m'-interval [m_lo, m_hi] of the partition of [w* tau, m] (for m' < w* tau = 0.108 the class {c >= w*, int c <= m'} is
empty): ladder certificate that  c >= w*, int_0^tau c <= m_hi, sigma' >= sigma'(m_lo)  ==>  mu_h[c] > L0 = 37.7.
usage: cert_slices2.py [w*] [L0]
"""
import sys, json, time
import cert_slice as CS

ws = sys.argv[1] if len(sys.argv) > 1 else '0.240'
L0 = sys.argv[2] if len(sys.argv) > 2 else '37.7'
part = ['0.108', '0.15', '0.18', '0.20', '0.21', '0.22', '0.224', '0.2265', '0.2282']
out = []
t0 = time.time()
for mlo, mhi in zip(part[:-1], part[1:]):
    try:
        r = CS.run(L0, ws, mlo, mhi, 60, fine=1e-5, verbose=False, thb=0.25)
    except Exception as e:
        r = dict(m_lo=mlo, m_hi=mhi, error=str(e), certified=False)
        print('[SLICE] [%s,%s] ERROR %s' % (mlo, mhi, e), flush=True)
    out.append(r)
    json.dump(out, open('cert_slices2_%s_%s.json' % (ws, L0), 'w'), indent=1)
print('ALL DONE (%.0fs)' % (time.time() - t0), flush=True)
