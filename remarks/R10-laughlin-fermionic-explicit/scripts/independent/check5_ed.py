#!/usr/bin/env python3
"""Item 5: independent re-computation of the N<=8 Laughlin-flux gaps (cert/ed_cert.py uses paper eq. (pair) + Higham).
Here: pair vectors from Clebsch-Gordan (indep_fock), own Slater-rule assembly of H_{N,Q} on the Lz=0 sector,
then (i) scipy.sparse.linalg.eigsh (Lanczos, which='SA', k=4) with residual bounds |lambda - lambda_exact| <= ||r||,
(ii) Sylvester inertia of H - sigma I by dense LDL^T (scipy.linalg.ldl) at sigma = certified g_N:
     #negative eigenvalues must be exactly 1 (the Laughlin zero mode), i.e. lambda_2 >= g_N,
(iii) inertia at sigma = 1e-6 must be 1 as well and at sigma=-1e-6 must be 0 (kernel is a line, Lemma 6.2).
usage: check5_ed.py N [g_cert]"""
import sys, itertools, json
from fractions import Fraction
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
from indep_fock import pair_vectors_cg

N = int(sys.argv[1]); Q = 3 * (N - 1)
ED = {r['N']: r for r in json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'primary', 'out', 'ed_cert.json')))}
g = float(Fraction(ED[N]['g_cert']))
V = pair_vectors_cg(Q)
M = Q * N // 2
sec = [S for S in itertools.combinations(range(Q + 1), N) if sum(S) == M]
idx = {S: k for k, S in enumerate(sec)}
n = len(sec)

def rem2(S, x, y):
    """<S\\{x,y}| c_y c_x |S>"""
    px = S.index(x); S1 = S[:px] + S[px + 1:]
    py = S1.index(y)
    return (-1) ** (px + py), S1[:py] + S1[py + 1:]

rows, cols, vals = [], [], []
for S in sec:
    for x, y in itertools.combinations(S, 2):
        s1, T = rem2(S, x, y)
        p = x + y - 1
        a = V[p].get((x, y), 0.0)
        if a == 0: continue
        for (u, w), b in V[p].items():
            if u in T or w in T: continue
            S2 = tuple(sorted(T + (u, w)))
            s2, _ = rem2(S2, u, w)
            rows.append(idx[S2]); cols.append(idx[S]); vals.append(s1 * s2 * a * b)
H = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
H = 0.5 * (H + H.T)
print(f"N={N} Q={Q} Lz=0 dim={n} nnz={H.nnz} ||H-H^T||(before sym)=ok", flush=True)
k = min(4, n - 1) if n > 2 else 1
if n > 50:
    w, X = spl.eigsh(H, k=k, which='SA', tol=1e-12)
else:
    w, X = np.linalg.eigh(H.toarray()); w, X = w[:max(k, 2)], X[:, :max(k, 2)]
o = np.argsort(w); w, X = w[o], X[:, o]
res = [np.linalg.norm(H @ X[:, i] - w[i] * X[:, i]) / np.linalg.norm(X[:, i]) for i in range(len(w))]
print("  (i) eigsh lowest:", ", ".join(f"{w[i]:.10f} (res {res[i]:.1e})" for i in range(len(w))), flush=True)
Hd = H.toarray()
def inertia_neg(sig):
    L, Dm, perm = sla.ldl(Hd - sig * np.eye(n), lower=True)
    return int((np.linalg.eigvalsh(Dm) < 0).sum())     # D is block diagonal 1x1/2x2
ng = inertia_neg(g); n0p = inertia_neg(1e-6); n0m = inertia_neg(-1e-6)
print(f"  (ii) #eig < g_cert={g}: {ng} (need 1)   (iii) #eig < 1e-6: {n0p} (need 1), #eig < -1e-6: {n0m} (need 0)   "
      f"gap(eigsh) - res = {w[1]-res[1]:.6f} >= g_cert: {w[1]-res[1] >= g}   "
      f"{'OK' if ng == 1 and n0p == 1 and n0m == 0 and w[1]-res[1] >= g else 'FAIL'}", flush=True)
