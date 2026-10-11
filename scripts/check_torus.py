from fractions import Fraction as F
from pathlib import Path
from math import comb
import json
R=Path(__file__).resolve().parent.parent
# An independent sparse ordinary polynomial ring, not the author's Laurent code.
def add(a,b):
    d=a.copy()
    for k,v in b.items():d[k]=d.get(k,F(0))+v
    return {k:v for k,v in d.items() if v}
def scale(a,t):return {k:v*t for k,v in a.items() if v*t}
def mul(a,b):
    d={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():d[(i+k,j+l)]=d.get((i+k,j+l),F(0))+v*w
    return {k:v for k,v in d.items() if v}
def power(a,n):
    d={(0,0):F(1)}
    for _ in range(n):d=mul(d,a)
    return d
one={(0,0):F(1)};x={(1,0):F(1)};y={(0,1):F(1)}
class P:
    def __init__(self,a):self.a=a if isinstance(a,dict) else {(0,0):F(a)}
    def __add__(self,b):return P(add(self.a,(b if isinstance(b,P) else P(b)).a))
    __radd__=__add__
    def __neg__(self):return P(scale(self.a,-1))
    def __sub__(self,b):return self+-P(b) if not isinstance(b,P) else self+-b
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):return P(mul(self.a,(b if isinstance(b,P) else P(b)).a))
    __rmul__=__mul__
    def __pow__(self,n):return P(power(self.a,n))
x,y=P(x),P(y)
def contractions(t):return 5*(2+t),F(1,2)*(16-t-t*t),6*(1-t),t*t+24*t+39
px,qx,rx,vx=contractions(x);py,qy,ry,vy=contractions(y)
d1=5*(ry*qx+py*vx)-8*py*qx
d2=5*(vy*px+qy*rx)-8*qy*px
d3=5*(vy*qx+qy*vx)-8*qy*qx
u=F(25,2)*(3*(x+y)-10)
v=F(25,2)*vx+(15*y-35)*qx
w=F(25,2)*vy+(15*x-35)*qy
minors=[d1,d1*d2-4*(1+x)*(1+y)*u*u,
 d1*d2*d3+8*(1+x)*(1+y)*u*v*w-2*(1+x)*d1*w*w-2*(1+y)*d2*v*v-4*(1+x)*(1+y)*d3*u*u]
def substitute(a,X,Y):
    z=P(0)
    for (i,j),v in a.a.items():z+=v*(X**i)*(Y**j)
    return z
records=[]
for a in minors:
    n=max(i for i,j in a.a);m=max(j for i,j in a.a)
    A=substitute(a,2*x-1,2*y-1)
    b=[[sum(v*F(comb(k,i),comb(n,i))*F(comb(l,j),comb(m,j)) for (i,j),v in A.a.items() if i<=k and j<=l) for l in range(m+1)] for k in range(n+1)]
    B=P(0)
    for k in range(n+1):
        for l in range(m+1):B+=b[k][l]*comb(n,k)*comb(m,l)*(x**k)*((1-x)**(n-k))*(y**l)*((1-y)**(m-l))
    assert A.a==B.a and min(v for row in b for v in row)>0
    records.append({'bidegree':[n,m],'min':str(min(v for row in b for v in row)), 'bernstein':[[str(v) for v in row] for row in b], 'power':[[i,j,str(v)] for (i,j),v in sorted(a.a.items())]})
expected=json.loads((R/'certificates/torus.json').read_text())['records']
def validate(saved):
 assert records==saved, 'certificate differs from independent reconstruction'
validate(expected)
assert sum(len(row) for r in records for row in r['bernstein'])==58
bad=json.loads(json.dumps(expected));bad[0]['bernstein'][0][0]='0'
try:
 validate(bad)
except AssertionError:
 pass
else:
 raise AssertionError('altered certificate was accepted')
print('PASS: exact identities and 58 strictly positive coefficients; altered coefficient rejected')
