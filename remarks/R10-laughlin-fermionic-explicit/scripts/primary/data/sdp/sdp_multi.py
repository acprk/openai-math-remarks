#!/usr/bin/env python3
"""269: SOS/SDP with ONE comparison operator Gamma shared by several fluxes (interval robustness test).
Reduced 4-body LMI (scalar allowance rho_D * I_b, cost n_D rho_D):
     U^T S^T C S U + rho I >= 0   <=>   L^T C L + rho I >= 0,   S Pi_hatG S^T = L L^T  (size <= |R_D|).
maximise min_Q gamma(Q).
usage: python sdp_multi.py Q1,Q2,... [pmax M tmax] [tag]"""
import sys, time, numpy as np, scipy.sparse as sp
import os; sys.path.insert(0, os.environ.get('LAUGHLIN_STEP1_DIR', '.'))  # directory holding engine.py from the step-1 search (provenance only)
from engine import *
import cvxpy as cp

Qs = [int(float(x)) for x in sys.argv[1].split(',')]
pmax, M, tmax = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) >= 5 else (7, 8, 7)
tag = sys.argv[5] if len(sys.argv) >= 6 else ''
sup = Support(pmax, M, range(tmax + 1))
Gam = {t: cp.Variable((sup.dim(t), sup.dim(t)), symmetric=True) for t in sup.tset}
cons = [Gam[t] >> 0 for t in sup.tset]
eta = sum(Gam[t][0, 0] for t in sup.tset)
g = cp.Variable()
vecG = {t: cp.vec(Gam[t], order='C') for t in sup.tset}
T0 = time.time()
for Q in Qs:
    assert Q >= sup.Dmax + 1
    eng = Engine(Q, sup)
    da = eng.da
    cost = 0
    for z in eng.zlist:
        if z in (0, 1, 3):
            continue
        tau = cp.Variable(nonneg=True)
        ez = sum(cp.sum(cp.multiply(eng.C3[z][t], Gam[t])) for t in sup.tset)
        cons.append(eng.qz(z) - da / eng.d3(z) * ez + tau >= 0)
        cost = cost + tau * (z + 1) * eng.d3(z) / da
    # S arrays
    Sarr = {}
    for D in range(1, sup.Dmax + 1):
        R = eng.Rset[D]
        for T in range(D, sup.Dmax + 1):
            d = {}
            for ir, r in enumerate(R):
                for (p, j, k), val in eng.Scoef[(D, T, r)].items():
                    for key, sg in (((p, j, k), 1.0), ((p, k, j), -1.0)):
                        d.setdefault(key, np.zeros(len(R)))[ir] += sg * val
            Sarr[(D, T)] = d
    Lm = {}
    for D in range(1, sup.Dmax + 1):
        bl, Smat, vecs = eng.bspace(D)
        U = range_basis(gram(vecs))
        if U.shape[1] == 0:
            continue
        KK = Smat @ U @ U.T @ Smat.T
        w, V = np.linalg.eigh(KK)
        keep = w > 1e-10 * max(1.0, w.max())
        if keep.sum() == 0:
            continue
        Lm[D] = V[:, keep] * np.sqrt(w[keep])
    for D, L in Lm.items():
        R = eng.Rset[D]; nR = len(R); k = L.shape[1]
        P = np.kron(L, L).T   # vec_C(L^T E L) = kron(L,L)^T vec_C(E)
        expr = 0
        for t in sup.tset:
            S = sup.S[t]; n = sup.dim(t)
            rows, cols, vals = [], [], []
            for a, (p, j) in enumerate(S):
                for b, (q, l) in enumerate(S):
                    i = p + j - t; kk = q + l - t; T = p + j + kk
                    if T < D:
                        continue
                    x = Sarr[(D, T)].get((p, j, kk)); y = Sarr[(D, T)].get((q, l, i))
                    if x is None or y is None:
                        continue
                    o = -(P @ np.outer(x, y).ravel())
                    nz = np.nonzero(np.abs(o) > 1e-15)[0]
                    rows.extend(nz.tolist()); cols.extend([(a + 1) * n + (b + 1)] * len(nz)); vals.extend(o[nz].tolist())
            if rows:
                Mp = sp.csr_matrix((vals, (rows, cols)), shape=(k * k, n * n))
                expr = expr + Mp @ vecG[t]
        if isinstance(expr, int):
            expr = cp.Constant(np.zeros(k * k))
        chi = np.zeros((nR, nR)); chi[0, 0] = (4 * Q - 1 - 2 * D) / da
        rho = cp.Variable(nonneg=True)
        Lmat = L.T @ chi @ L - cp.reshape(expr, (k, k), order='C') + rho * np.eye(k)
        cons.append(0.5 * (Lmat + Lmat.T) >> 0)
        cost = cost + eng.nD(D) * rho
    cons.append(g <= 1 - eta - eng.delta_tail() - cost)
print(f"built {len(Qs)} fluxes in {time.time()-T0:.1f}s", flush=True)
prob = cp.Problem(cp.Maximize(g), cons)
T0 = time.time()
prob.solve(solver='CLARABEL')
print(f"solve {time.time()-T0:.1f}s status={prob.status}  min_Q gamma_SDP={g.value:.6f}  eta={eta.value:.6f}", flush=True)
fn = f'multi_{"_".join(map(str,Qs))}_p{pmax}M{M}t{tmax}{tag}.npz'
np.savez(fn, **{f'G{t}': Gam[t].value for t in sup.tset}, gamma=g.value)
print("saved", fn)
