from fractions import Fraction as F
from pathlib import Path
from flint import arb,fmpq,ctx
import json
R=Path(__file__).resolve().parent.parent
(R/"build").mkdir(exist_ok=True)
# Alternating expansions; bound remainder by the first omitted term.
def atan_bounds(x,n):
    s=sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n));e=(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(s,s+e),max(s,s+e)
def j0_bounds(x,n=80):
    z=x*x/4;t=F(1);s=t
    for k in range(1,n):t=-t*z/F(k*k);s+=t
    e=-t*z/F(n*n)
    assert z/F(n*n)<1
    return min(s,s+e),max(s,s+e)
a0,a1=atan_bounds(F(1,5),150);b0,b1=atan_bounds(F(1,239),50)
pi0,pi1=16*a0-4*b1,16*a1-4*b0
assert F(314159265358979,10**14)<pi0<pi1<F(314159265358980,10**14)
x=F('2.404825557695772768621631879326454643124244909145');e=F(1,10**45)
lo,hi=x-e,x+e
assert j0_bounds(lo)[0]>0 and j0_bounds(hi)[1]<0
assert hi<F(241,100) and F(241,100)**2<8 # J1 >= x/2-x^3/16 >0 on (0,2.41]
ctx.prec=200
Q=lambda x:arb(fmpq(F(x).numerator,F(x).denominator))
PI=Q(pi0).union(Q(pi1));J=Q(lo).union(Q(hi));K=PI*J*J
Lambda=12*PI;m=Q(F(1141,5000));a=Q(F(519,2000));tau=Q(F(9,20))
sigma1=Q(F(3,2))*(PI/arb(3).sqrt()).sqrt()
checks={
 '37.7>12pi':Q(F(377,10))-Lambda,
 'core_half_mass<a':a-(Q(F(1,2))-K/(2*Lambda)),
 'gap_mass>1/2':2*K/Lambda-Q(F(1,2)),
 'FK_zero_window>37.7':K/(2*m)-Q(F(377,10)),
 'B_budget_slack':Q(F('0.0879293'))-(1-4*m),
 'nu_at_m0>37.7':K/(2*(a-Q(F(3,100))))-Q(F(377,10)),
 'sigma37.7_at_m0>0':sigma1/(a-Q(F(3,100))).sqrt()*(1-Q(F(377,10))*2*(a-Q(F(3,100)))/K),
 'arc_cap_valid':PI*tau*tau/2-m}
assert all(v>0 for v in checks.values())
rec={'status':'PASS','arithmetic':'Fraction alternating-series pi and first Bessel root; 200-bit arb interfaces', 'root_uniqueness':'J0 derivative=-J1, J1>=x/2-x^3/16>0 for 0<x<=2.41', 'pi_bounds_decimal':[str(Q(pi0)),str(Q(pi1))],'j01_bracket':[str(lo),str(hi)], 'strict_slacks':{k:str(v) for k,v in checks.items()}, 'not_formal_transcendental_proof':True}
(R/'build/constants-receipt.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2))
