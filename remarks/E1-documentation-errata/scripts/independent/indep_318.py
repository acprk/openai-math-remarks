"""Independent re-implementation (does not import or copy 318_ce_ranks.py).

Ranks of the Chevalley-Eilenberg differentials d_k : C^k(g;Q) -> C^{k+1}(g;Q),
k = 0..3, trivial coefficients, for g = sl_2 and sl_3 (family 318,
SUPPORT-COVERAGE.md: expected (0,3,0,0) and (0,8,20,35)).
Differences from the main script: works on COCHAINS directly with the formula
  (d w)(x_0..x_k) = sum_{i<j} (-1)^{i+j} w([x_i,x_j], x_0..^i..^j..x_k),
uses a different basis (E_ij for i!=j and D_i = E_ii - E_nn), expresses brackets
in that basis by reading off matrix entries (no linear solve), and computes
ranks with sympy's exact Matrix.rank.  Also prints the Betti numbers
dim H^k = dim C^k - rank d_k - rank d_{k-1}.  Environment: python3, sympy. ~2 s.
"""
from itertools import combinations
import sympy as sp

def basis(n):
    B = []
    for i in range(n):
        for j in range(n):
            if i != j:
                M = sp.zeros(n); M[i, j] = 1; B.append(M)
    for i in range(n - 1):
        M = sp.zeros(n); M[i, i] = 1; M[n - 1, n - 1] = -1; B.append(M)
    return B

def coords(X, n):
    # X traceless; off-diagonal coordinates are entries, diagonal part = sum_i X_ii D_i
    c = [X[i, j] for i in range(n) for j in range(n) if i != j]
    c += [X[i, i] for i in range(n - 1)]
    return c

def ce_ranks(n):
    B = basis(n); m = len(B)
    br = {(a, b): coords(B[a] * B[b] - B[b] * B[a], n) for a in range(m) for b in range(m)}
    ranks, dims = [], []
    for k in range(0, 4):
        Ck = list(combinations(range(m), k)); idx = {s: t for t, s in enumerate(Ck)}
        Ck1 = list(combinations(range(m), k + 1))
        dims.append(len(Ck))
        A = sp.zeros(len(Ck1), len(Ck))   # rows: (k+1)-tuples x, cols: basis k-cochains w_s
        for r, x in enumerate(Ck1):
            for i in range(k + 1):
                for j in range(i + 1, k + 1):
                    rest = [x[t] for t in range(k + 1) if t not in (i, j)]
                    for c, coef in enumerate(br[x[i], x[j]]):
                        if coef == 0 or c in rest: continue
                        args = [c] + rest
                        # sign of the permutation sorting args
                        sgn = 1
                        for p in range(len(args)):
                            for q in range(p + 1, len(args)):
                                if args[p] > args[q]: sgn = -sgn
                        A[r, idx[tuple(sorted(args))]] += (-1) ** (i + j) * sgn * coef
        ranks.append(A.rank())
    dims.append(len(list(combinations(range(m), 4))))
    betti = [dims[k] - ranks[k] - (ranks[k - 1] if k else 0) for k in range(4)]
    return ranks, betti

for n, exp in [(2, [0, 3, 0, 0]), (3, [0, 8, 20, 35])]:
    r, b = ce_ranks(n)
    print(f"sl{n}: ranks d0..d3 = {r}  expected {exp}  match={r == exp}  dim H^0..H^3 = {b}")
