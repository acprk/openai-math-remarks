#!/usr/bin/env python3
"""Item 2 (G1, Lemma R): rank W_D = hw(D) - zm(D), checked by direct linear algebra on wedge^4, NOT via the kernel.
rank W_D = dim( HW_{2Q-1-D}(wedge^4 U_Q) cap ran What ) where ran What cap weight = span{ v_p ^ e_j ^ e_k : p+1+j+k = D+1 }.
  = rank(U) - rank(J_+ U),  U = matrix of the vectors u_b in the 4-wedge basis of index sum D+1.
Two independent evaluations:
  (i) exact rational arithmetic in the rescaled basis f_i = sqrt(C(Q,i)) e_i: by paper eq. (pair)
      v_p is proportional to sum_{x+y=p+1, x<y} (x-y) f_x ^ f_y (integer, Q-free), and J_+ f_i = (Q-i+1) f_{i-1}
      (integer).  Rank is invariant under this invertible diagonal change of basis.  Exact at Q = 24, 300, 1000, 10^4, 10^8.
 (ii) float64 SVD in the ORTHONORMAL basis with the CG pair vectors (indep_fock.pair_vectors_cg is too slow at Q=10^4,
      so (ii) uses the paper formula eq. pair in log-gamma form) at Q = 300, 1000, 10^4, reporting the singular-value gap."""
import itertools, sys
from fractions import Fraction
from math import lgamma, log, exp, sqrt
import numpy as np
from flint import fmpq_mat, fmpq

def subsets(n):
    return sorted(c for c in itertools.combinations(range(n + 1), 4) if sum(c) == n)

def wsign(lst):
    srt = sorted(lst); perm = [lst.index(v) for v in srt]
    s = 1
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]: s = -s
    return tuple(srt), s

def build(D, pair, raise_coef):
    rows = subsets(D + 1); ri = {r: i for i, r in enumerate(rows)}
    lows = subsets(D); li = {r: i for i, r in enumerate(lows)}
    cols = [(p, j, k) for p in range(D + 1) for j in range(D + 1) for k in range(j + 1, D + 1) if p + j + k == D]
    U = [[0] * len(cols) for _ in rows]
    for c, (p, j, k) in enumerate(cols):
        for (x, y), a in pair(p).items():
            lst = [x, y, j, k]
            if len(set(lst)) < 4: continue
            t, s = wsign(lst)
            U[ri[t]][c] += s * a
    Jp = [[0] * len(rows) for _ in lows]
    for c, A in enumerate(rows):
        for pos, i in enumerate(A):
            if i >= 1 and (i - 1) not in A:
                t = list(A); t[pos] = i - 1
                Jp[li[tuple(t)]][c] += raise_coef(i)
    return U, Jp, len(rows)

def hw(D):
    c = lambda n: sum(1 for cc in itertools.combinations(range(n + 1), 4) if sum(cc) == n)
    return c(D + 1) - c(D)

def zm(D):
    def p4(n):
        if n < 0: return 0
        return sum(1 for a in range(n + 1) for b in range(a, n + 1) for c in range(b, n + 1) if n - a - b - c >= c)
    return p4(D - 17) - p4(D - 18)

def rank_q(M):
    if not M or not M[0]: return 0
    return fmpq_mat([[fmpq(v) for v in r] for r in M]).rank()

print("D  hw  zm  hw-zm | exact rank (same at Q=24,300,1e3,1e4,1e8; -1 if not) | float rank Q=300 / 1000 / 10000 (sigma_r, sigma_{r+1})")
bad = 0
for D in range(1, 24):
    mono_pair = lambda p: {(x, p + 1 - x): (x - (p + 1 - x)) for x in range(0, (p + 2) // 2) if x < p + 1 - x}
    rexs = []
    for Qe in (24, 300, 1000, 10**4, 10**8):
        U, Jp, nr = build(D, mono_pair, lambda i, Qe=Qe: Qe - i + 1)
        if nr == 0:
            rexs.append(0); continue
        JU = [[sum(Jp[a][b] * U[b][c] for b in range(len(U))) for c in range(len(U[0]))] for a in range(len(Jp))]
        rexs.append(rank_q(U) - rank_q(JU))
    rex = rexs[0] if len(set(rexs)) == 1 else -1
    fl = []
    for Q in (300, 1000, 10000):
        def onb_pair(p, Q=Q):
            out = {}
            for x in range(0, (p + 2) // 2):
                y = p + 1 - x
                if x >= y: continue
                lf = lambda s, n: sum(log(s - i) for i in range(n))
                lv = lf(Q, x) + lf(Q, y) + lgamma(p + 1) - log(Q) - lf(2 * Q - 2, p) - lgamma(x + 1) - lgamma(y + 1)
                out[(x, y)] = (x - y) * exp(0.5 * lv)
            return out
        Uf, Jf, nr = build(D, onb_pair, lambda i, Q=Q: sqrt(i * (Q - i + 1)) / sqrt(Q))
        if nr == 0:
            fl.append((0, float('nan'), float('nan'))); continue
        Uf = np.array(Uf, float); Jf = np.array(Jf, float)
        # orthonormal basis of ran U, then kernel of J_+ on it
        u, s, _ = np.linalg.svd(Uf, full_matrices=False)
        ru = int((s > 1e-9 * s[0]).sum())
        Ob = u[:, :ru]
        sv = np.linalg.svd(Jf @ Ob, compute_uv=False) if Jf.size else np.zeros(0)
        sv = np.concatenate([sv, np.zeros(max(0, ru - len(sv)))])
        nul = int((sv < 1e-9).sum())
        sorted_sv = np.sort(sv)
        lo = sorted_sv[nul - 1] if nul > 0 else float('nan')
        hi = sorted_sv[nul] if nul < len(sorted_sv) else float('nan')
        fl.append((nul, lo, hi))
    ok = rex == hw(D) - zm(D) and all(f[0] == rex for f in fl)
    bad += not ok
    print(f"{D:2d} {hw(D):3d} {zm(D):3d} {hw(D)-zm(D):4d} | {rex:4d} | " +
          " / ".join(f"{f[0]} ({f[1]:.1e},{f[2]:.1e})" for f in fl) + ("" if ok else "   MISMATCH"))
print("Lemma R mismatches:", bad)
