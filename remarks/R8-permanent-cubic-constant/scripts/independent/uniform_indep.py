"""Uniform constant: for every m >= m0, exhibit (k,q,d) admissible for the manuscript's chain and record
c^{-1}(m)=4m^3/((r-1)(d-1)).  Range [m0, M]: local search around k~32m/330, q~11k/32 (BOX) or q~k/3 (paper).
Tail m > M: explicit choice k=floor(32(m-2)/330), q=floor(11k/32), with an exact (sympy) monotonicity proof
of the resulting rational bound g(m).  Every reported value is an upper bound on the best c^{-1}(m)."""
import sys
from math import comb
import sympy as sp

def P(model, k, q):
    r = k-2*q
    if r < 2 or q < 0: return 0
    if model == 'paper':
        if 3*q > k-1: return 0
        D = comb(q+2, 2)-k
    else:
        D = (q+1)**2-k
    d = min(D, k*k-1)
    return (r-1)*(d-1) if d >= 2 else 0

def local(m, model):
    best = 0
    k0 = 32*m//330 if model == 'BOX' else 3*m//31
    for k in range(k0-6, k0+8):
        q0 = 11*k//32 if model == 'BOX' else (k-1)//3
        for q in range(q0-6, q0+7):
            if q < 0 or 11*k-2*q > m: continue
            best = max(best, P(model, k, q))
    return best

m0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1408
M = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
for model in ['paper', 'BOX']:
    worst, wm = 0.0, None
    for m in range(m0, M+1):
        v = local(m, model)
        assert v > 0, (model, m)
        c = 4*m**3/v
        if c > worst: worst, wm = c, m
    print(f"{model}: max over {m0}<=m<={M} of c^-1(m) = {worst:.1f} e  (attained m={wm})")

# tail for BOX
m = sp.symbols('m', positive=True)
K = sp.Rational(32, 330)*(m-2) - 1
A = sp.Rational(10, 32)*K - 1                    # <= r-1 = k-2q-1  (q<=11k/32)
B = (sp.Rational(11, 32)*K)**2 - K - 1           # <= d-1 = (q+1)^2-k-1  (q+1 > 11k/32), increasing in k
g = sp.together(4*m**3/(A*B))
dg = sp.together(sp.diff(g, m))
num = sp.Poly(sp.numer(dg), m)
roots = [complex(r) for r in num.nroots()]
real_roots = sorted(r.real for r in roots if abs(r.imag) < 1e-9)
print("tail: largest real root of numerator of g'(m):", max(real_roots) if real_roots else None,
      "; sign of g' at M:", sp.sign(dg.subs(m, M)))
print(f"tail: g(M) = {float(g.subs(m, M)):.1f} e ; lim g = {float(sp.limit(g, m, sp.oo)):.1f} e")
# also check A>0, B>0 at M and that 11k-2q<=m for the explicit choice (11k-2q <= 330k/32+62/32 <= m)
print("A(M), B(M) =", float(A.subs(m, M)), float(B.subs(m, M)))
