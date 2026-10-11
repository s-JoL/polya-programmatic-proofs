#!/usr/bin/env python3
"""Independent exact certificate, reviewed before execution; standard library only.

Start with the 3-node element, wrap in Laurent variables, form the full tensor
matrix and its three complement columns. No author's minor is an input.
Conjugation is inversion of the unit-circle Laurent variables. Only determinants
are converted to cosines, after checking both separate inversion symmetries.
"""
from fractions import Fraction as F
from itertools import permutations
from math import comb
from pathlib import Path
import hashlib
import json


class P:
    def __init__(self, terms=None):
        if isinstance(terms, (int, F)):
            terms = {(0, 0): F(terms)}
        self.a = {k: F(v) for k, v in (terms or {}).items() if v}

    def __add__(self, other):
        other = other if isinstance(other, P) else P(other)
        a = dict(self.a)
        for k, v in other.a.items():
            a[k] = a.get(k, F(0)) + v
        return P(a)

    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.a.items()})

    def __sub__(self, other):
        return self + -as_p(other)

    def __rsub__(self, other):
        return as_p(other) + -self

    def __mul__(self, other):
        other = as_p(other)
        a = {}
        for (i, j), v in self.a.items():
            for (k, l), w in other.a.items():
                key = (i+k, j+l)
                a[key] = a.get(key, F(0)) + v*w
        return P(a)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert n >= 0
        result = P(1)
        for _ in range(n):
            result = result*self
        return result

    def star(self):
        return P({(-i, -j): v for (i, j), v in self.a.items()})

    def __eq__(self, other):
        return self.a == as_p(other).a


def as_p(a):
    return a if isinstance(a, P) else P(a)


def matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def adjoint(a):
    return [[as_p(a[j][i]).star() for j in range(len(a))]
            for i in range(len(a[0]))]


def tensor(a, b):
    return [[a[i][j]*b[k][l] for j in range(len(a[0]))
             for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def det(a):
    n = len(a)
    out = P()
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        term = P((-1)**inversions)
        for i in range(n):
            term = term*a[i][p[i]]
        out = out+term
    return out


z, w = P({(1, 0): 1}), P({(0, 1): 1})
x, y = P({(1, 0): 1}), P({(0, 1): 1})  # ordinary polynomial ring


def derive(axis):
    c = (axis+axis.star())*F(1, 2)
    s = 1+axis
    wrap = [[P(1), P()], [P(), P(1)], [axis, P()]]
    # Raw unit-element integer matrices; scaling 3d and 15m is exact.
    raw_d = [[7, -8, 1], [-8, 16, -8], [1, -8, 7]]
    raw_m = [[4, 2, -1], [2, 16, 2], [-1, 2, 4]]
    D = matmul(matmul(adjoint(wrap), raw_d), wrap)
    M = [[entry*F(1, 2) for entry in row]
         for row in matmul(matmul(adjoint(wrap), raw_m), wrap)]
    assert D == [[14+2*c, -8*s.star()], [-8*s, P(16)]]
    assert M == [[4-c, s.star()], [s, P(8)]]
    j, g = [P(1), s*F(1, 2)], [-s.star()*F(1, 2), P(1)]

    def product(a, Q, b):
        return sum(a[i].star()*Q[i][k]*b[k] for i in range(2) for k in range(2))

    p, q = product(j, M, j), product(g, M, g)
    r, v = product(j, D, j), product(g, D, g)
    aa, bb = product(j, M, g), product(j, D, g)
    n = sum(t.star()*t for t in j)
    assert n == (3+c)*F(1, 2)
    assert sum(t.star()*t for t in g) == n
    assert sum(j[i].star()*g[i] for i in range(2)) == 0
    assert p == 5*(2+c) and q == (16-c-c*c)*F(1, 2)
    assert r == 6*(1-c) and v == c*c+24*c+39
    assert aa == s.star()*F(5, 2) and bb == -3*(1-c)*s.star()
    assert det(M) == 10*(3-c)
    assert det(D) == 96*(1-c)
    return D, M, j, g, (p, q, r, v, aa, bb)


def twoT(n, variable):
    # variable represents c: z^n+z^-n=2*T_n(c), including n=0.
    a, b = P(2), 2*variable
    if n == 0:
        return a
    for _ in range(1, n):
        a, b = b, 2*variable*b-a
    return b


def substitute(poly, X, Y):
    assert all(i >= 0 and j >= 0 for i, j in poly.a)
    return sum(v*(X**i)*(Y**j) for (i, j), v in poly.a.items())


def to_cos(poly):
    assert poly == P({(-i, j): v for (i, j), v in poly.a.items()})
    assert poly == P({(i, -j): v for (i, j), v in poly.a.items()})
    result = P()
    for (i, j), coefficient in poly.a.items():
        if i >= 0 and j >= 0:
            result += coefficient*(twoT(i, x) if i else P(1))*(twoT(j, y) if j else P(1))
    # Independent inverse substitution, identity in the entire Laurent ring.
    assert substitute(result, (z+z.star())*F(1, 2), (w+w.star())*F(1, 2)) == poly
    return result


def encode(poly):
    return [[i, j, str(v)] for (i, j), v in sorted(poly.a.items())]


def bernstein(poly):
    n = max(i for i, j in poly.a)
    m = max(j for i, j in poly.a)
    transformed = substitute(poly, 2*x-1, 2*y-1)
    coeff = [[sum(v*F(comb(k, i), comb(n, i))*F(comb(l, j), comb(m, j))
                  for (i, j), v in transformed.a.items() if i <= k and j <= l)
              for l in range(m+1)] for k in range(n+1)]
    expansion = P()
    for k in range(n+1):
        for l in range(m+1):
            expansion += (coeff[k][l]*comb(n, k)*comb(m, l)*
                          (x**k)*((1-x)**(n-k))*(y**l)*((1-y)**(m-l)))
    assert expansion == transformed
    minimum = min(v for row in coeff for v in row)
    assert minimum > 0
    return {"bidegree": [n, m], "minimum": str(minimum),
            "coefficients": [[str(v) for v in row] for row in coeff],
            "power_coefficients_cx_cy": encode(poly),
            "roundtrip_exact": True}


def main():
    Dx, Mx, jx, gx, cx = derive(z)
    Dy, My, jy, gy, cy = derive(w)
    dymx, mydx, mymx = tensor(Dy, Mx), tensor(My, Dx), tensor(My, Mx)
    A = [[5*(dymx[i][k]+mydx[i][k])-8*mymx[i][k]
          for k in range(4)] for i in range(4)]
    columns = [tensor([[t] for t in jy], [[t] for t in gx]),
               tensor([[t] for t in gy], [[t] for t in jx]),
               tensor([[t] for t in gy], [[t] for t in gx])]
    V = [[columns[k][i][0] for k in range(3)] for i in range(4)]
    C = matmul(matmul(adjoint(V), A), V)
    assert C == adjoint(C)
    px, qx, rx, vx, ax, bx = cx
    py, qy, ry, vy, ay, by = cy
    # Check every entry against a second contraction-level construction.
    diag = [5*(ry*qx+py*vx)-8*py*qx,
            5*(vy*px+qy*rx)-8*qy*px,
            5*(vy*qx+qy*vx)-8*qy*qx]
    e12 = 5*(by*ax.star()+ay*bx.star())-8*ay*ax.star()
    e13 = 5*(by*qx+ay*vx)-8*ay*qx
    e23 = 5*(vy*ax+qy*bx)-8*qy*ax
    for i in range(3):
        assert C[i][i] == diag[i]
    assert C[0][1] == e12 and C[0][2] == e13 and C[1][2] == e23
    hdet = (diag[0]*diag[1]*diag[2] +
            e12*e23*e13.star()+(e12*e23*e13.star()).star() -
            diag[0]*e23*e23.star()-diag[1]*e13*e13.star()-diag[2]*e12*e12.star())
    assert det(C) == hdet
    minors = [det([row[:k] for row in C[:k]]) for k in (1, 2, 3)]
    records = [bernstein(to_cos(p)) for p in minors]
    expected=json.loads((Path(__file__).resolve().parent.parent/"certificates/torus.json").read_text())["records"]
    for actual, saved in zip(records,expected):
        assert actual["power_coefficients_cx_cy"]==saved["power"]
        assert actual["coefficients"]==saved["bernstein"]
    result = {"status": "strict_positive_exact", "input_energy": 8,
              "input": "unit 3-node Q2 matrices and endpoint wrap; direct tensor restriction",
              "arithmetic": "Fraction Laurent polynomials; no floating point",
              "symbol_scaling": "225*A = 5*(Dy tensor Mx + My tensor Dx) - E*(My tensor Mx)",
              "basis": ["jy tensor gx", "gy tensor jx", "gy tensor gx"],
              "checks": ["original-element wrap", "all six 1D contractions",
                         "basis norms and orthogonality", "all complement entries",
                         "determinant by permutations and Hermitian formula",
                         "separate phase inversion symmetries",
                         "cosine Laurent roundtrip", "Bernstein power roundtrip"],
              "leading_principal_minors": records,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "original_problem_solved": False}
    path = Path(__file__).resolve().parent.parent / "build/torus-raw.json"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "minima": [r["minimum"] for r in records],
                      "bidegrees": [r["bidegree"] for r in records],
                      "checker_sha256": result["checker_sha256"]}))


if __name__ == "__main__":
    main()
