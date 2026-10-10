#!/usr/bin/env python3
"""
R2 numerical cross-check (FLOATING POINT; not part of any proof):
minimise P(n,n) subject to OpenAI's axioms (S),(B),(C),(T) on a finite grid
1 <= a,b <= N with N >= 4n-1, and compare with P*(n,n) = (3n-1)phi(n)/2.
The minimality P >= P* is proved in note.tex (Theorem 1(ii)); this LP only
confirms numerically that no smaller diagonal value is feasible on the grid.

Environment: /usr/bin/python3 with numpy + scipy (HiGHS). Single thread.
Run:  OMP_NUM_THREADS=1 /usr/bin/python3 r2_least_profile_lp.py
Expected output: expected_output/r2_least_profile_lp.txt (runtime ~1 min)
"""
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

def phi(m):
    r = Fr(1)
    for j in range(1, m):
        r *= Fr(3 * j + 1, 3 * j)
    return r

def lp_min_diag(n, N):
    idx = {}
    for a in range(1, N + 1):
        for b in range(a, N + 1):
            idx[(a, b)] = len(idx)
    V = lambda a, b: idx[(min(a, b), max(a, b))]
    rows = []                       # each row: dict var->coef, meaning sum <= 0
    for a in range(1, N + 1):
        for b in range(2, N):       # (C): P(a,b+1)+P(a,b-1)-2P(a,b) <= 0
            r = {}
            for (x, cf) in ((V(a, b + 1), 1), (V(a, b - 1), 1), (V(a, b), -2)):
                r[x] = r.get(x, 0) + cf
            rows.append(r)
        for h in range(1, N + 1):   # (T): 3P(a,h) - P(a,3h+a-1) <= 0
            if 3 * h + a - 1 <= N:
                r = {}
                for (x, cf) in ((V(a, h), 3), (V(a, 3 * h + a - 1), -1)):
                    r[x] = r.get(x, 0) + cf
                rows.append(r)
    A = lil_matrix((len(rows), len(idx)))
    for i, r in enumerate(rows):
        for x, cf in r.items():
            A[i, x] += cf
    bounds = [(0, None)] * len(idx)
    for b in range(1, N + 1):       # (B)
        bounds[V(1, b)] = (b, b)
    c = np.zeros(len(idx)); c[V(n, n)] = 1
    res = linprog(c, A_ub=A.tocsr(), b_ub=np.zeros(len(rows)), bounds=bounds,
                  method="highs", options={"presolve": True})
    return res.fun

print(" n   N   LP min P(n,n)        P*(n,n)          |diff|")
worst = 0
for n in (2, 3, 4, 5, 6, 8, 10, 12):
    N = 4 * n - 1
    v = lp_min_diag(n, N)
    ex = (3 * n - 1) * phi(n) / 2
    worst = max(worst, abs(v - float(ex)))
    print(f"{n:2d} {N:3d}  {v:16.10f}  {float(ex):16.10f}  {abs(v-float(ex)):.1e}  ({ex})")
print(f"max |LP - P*| = {worst:.1e}  ->", "AGREE (to 1e-8)" if worst < 1e-8 else "DISAGREE")
