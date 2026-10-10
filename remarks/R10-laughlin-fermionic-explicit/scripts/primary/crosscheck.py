#!/usr/bin/env python3
"""269 cert: INDEPENDENT cross-check (Jordan-Wigner sparse matrices, no spin/coupling reduction of the certificate).
For given Q and N: (a) smallest positive eigenvalue of H_{N,Q} (dense/sparse eig, Gamma-free) must be >= gamma_cert(Q);
(b) operator check  H^2 - gamma H - K_Q >= 0  on each spin-j highest-weight space with the RATIONAL Gamma of the
certificate (K_Q = d_a * Haar average of K, evaluated via descendant frames), gamma = gamma_cert(Q).
usage: crosscheck.py gamma.json Q N1,N2,.. gamma_cert
(original header of validate_sector.py follows)
269: brute-force check of an SDP comparison operator in the N-particle sectors (sparse, physical operators).
For each N and each total spin j, on the highest-weight space h (J_+ h = 0, J_z = j):
   X_j = h^T H^2 h - gamma h^T H h - d_a Tr_j(K)/(2j+1),    Tr_j = sum_mu f_mu^T K f_mu (descendant frames)
since K_Q = d_a * twirl(K) and all of H, K_Q are SU(2)-invariant.
Reports lambda_min(X_j) and the best gamma' in this sector  (min over j of generalised eig on ran H).
usage: python validate_sector.py file.npz Q N1,N2,.. [gamma] [pmax M tmax]"""
import sys, itertools, numpy as np, scipy.sparse as sp, scipy.linalg as sla
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from engine_float import pair_vec, Support
import aengine as ae, rationalize as rz
fn = sys.argv[1]; Q = int(sys.argv[2]); Ns = [int(x) for x in sys.argv[3].split(',')]
_G = ae.Gamma(rz.load(fn)).float()
dat = {f'G{t}': _G[t] for t in _G}
gamma = float(sys.argv[4])
sup = Support(*(int(x) for x in sys.argv[5:7]), range(int(sys.argv[7]) + 1)) if len(sys.argv) > 7 else Support()
nm = Q + 1
bases = {}; index = {}
def basis(N):
    if N not in bases:
        L = list(itertools.combinations(range(nm), N))
        bases[N] = L; index[N] = {s: i for i, s in enumerate(L)}
    return bases[N]
def ann(j, N):  # c_j : N -> N-1
    src = basis(N); basis(N - 1); idx = index[N - 1]
    r, c, v = [], [], []
    for col, s in enumerate(src):
        if j in s:
            pos = s.index(j)
            t = s[:pos] + s[pos + 1:]
            r.append(idx[t]); c.append(col); v.append((-1) ** pos)
    return sp.csr_matrix((v, (r, c)), shape=(len(bases[N - 1]), len(src)))
cache = {}
def C(j, N):
    if (j, N) not in cache:
        cache[(j, N)] = ann(j, N)
    return cache[(j, N)]
def Bp(p, N):  # N -> N-2 ; B(e_x^e_y) = c_y c_x
    M = None
    for (x, y), v in pair_vec(Q, p).items():
        term = v * (C(y, N - 1) @ C(x, N))
        M = term if M is None else M + term
    return M
for N in Ns:
    basis(N)
    dimN = len(bases[N])
    Bs = {p: Bp(p, N) for p in range(2 * Q - 1)}
    H = sum((Bs[p].T @ Bs[p]) for p in Bs).tocsr()
    K = sp.csr_matrix((dimN, dimN))
    for t in sup.tset:
        G = 0.5 * (dat[f'G{t}'] + dat[f'G{t}'].T)
        w, V = np.linalg.eigh(G)
        # F_k = sum_m g_k[m] O_m ; K_t = sum_{m,m'} G[m,m'] O_m^T O_m'
        ops = [Bs[t]]
        for (p, j) in sup.S[t]:
            i = p + j - t
            ops.append((C(i, N - 2).T @ C(j, N - 2) @ Bs[p]) if N >= 3 else None)
        for k in range(len(w)):
            if w[k] <= 0:
                continue
            g = V[:, k] * np.sqrt(w[k])
            F = g[0] * ops[0]
            for m in range(1, len(ops)):
                if ops[m] is not None and g[m] != 0:
                    F = F + g[m] * ops[m]
            K = K + (F.T @ F)
    K = K.tocsr()
    # spin structure
    msum = np.array([sum(s) for s in bases[N]])
    Jm = sp.csr_matrix((dimN, dimN))
    for j in range(Q):
        Jm = Jm + np.sqrt((j + 1) * (Q - j)) * (C(j, N).T @ C(j, N) * 0)  # placeholder (replaced below)
    # build J_- directly: e_j -> sqrt((j+1)(Q-j)) e_{j+1}
    r, c, v = [], [], []
    idx = index[N]
    for col, s in enumerate(bases[N]):
        for pos, j in enumerate(s):
            if j + 1 <= Q and (j + 1) not in s:
                t = list(s); t[pos] = j + 1
                # reorder: j+1 replaces j at same position; since j+1 not in s and j<j+1, order is preserved
                r.append(idx[tuple(t)]); c.append(col); v.append(np.sqrt((j + 1) * (Q - j)))
    Jm = sp.csr_matrix((v, (r, c)), shape=(dimN, dimN))
    Jp = Jm.T.tocsr()
    H2 = (H @ H).tocsr()
    worst = (1e9, None); bestg = 1e9; lamHW = float('inf')
    # exact spin-j subspaces level by level from the Casimir J^2 = J_- J_+ + Jz^2 + Jz (noise-free frames)
    levs = {}
    for ms in sorted(set(msum)):
        lev = np.where(msum == ms)[0]
        jz = (Q * N - 2 * ms) / 2.0
        A = Jp[:, lev].toarray()
        A = A[np.any(A != 0, axis=1)] if A.size else A
        J2 = (A.T @ A if A.shape[0] else np.zeros((len(lev), len(lev)))) + (jz * jz + jz) * np.eye(len(lev))
        w, V = np.linalg.eigh(J2)
        jj = np.round((-1 + np.sqrt(1 + 4 * np.maximum(w, 0))) / 2 * 2) / 2   # spin labels
        levs[ms] = (lev, jj, V)
    for ms in sorted(levs):
        lev, jjs, V = levs[ms]
        jz2 = Q * N - 2 * ms
        if jz2 < 0:
            continue
        sel = np.isclose(jjs, jz2 / 2)
        if not sel.any():
            continue
        hw = V[:, sel]
        h = np.zeros((dimN, hw.shape[1])); h[lev] = hw
        Tr = np.zeros((hw.shape[1], hw.shape[1])); f = h.copy()
        for mu in range(jz2 + 1):
            Tr += f.T @ (K @ f)
            if mu < jz2:
                g_ = Jm @ f
                lv, jl, Vl = levs[ms + mu + 1]
                E = Vl[:, np.isclose(jl, jz2 / 2)]
                f = np.zeros_like(f); f[lv] = E @ (E.T @ g_[lv])
                f /= np.linalg.norm(f, axis=0)
        KQ = (2 * Q - 1) * Tr / (jz2 + 1)
        Hh = h.T @ (H @ h); H2h = h.T @ (H2 @ h)
        X = H2h - gamma * Hh - KQ
        lm = np.linalg.eigvalsh(0.5 * (X + X.T)).min()
        if lm < worst[0]:
            worst = (lm, jz2 / 2)
        # best gamma' on ran Hh
        wh, Vh = np.linalg.eigh(0.5 * (Hh + Hh.T))
        if (wh > 1e-9).any():
            lamHW = min(lamHW, float(wh[wh > 1e-9].min()))
        P = Vh[:, wh > 1e-9]
        if P.shape[1]:
            Y = P.T @ (H2h - KQ) @ P; Z = P.T @ Hh @ P
            gp = sla.eigh(0.5 * (Y + Y.T), 0.5 * (Z + Z.T), eigvals_only=True).min()
            bestg = min(bestg, gp)
    lam1 = lamHW                      # min positive eigenvalue of H over all highest-weight spaces (= over the sector)
    if dimN <= 6000:                  # consistency: full dense spectrum
        ev = np.linalg.eigvalsh(H.toarray()); pos = ev[ev > 1e-8]
        assert abs(float(pos.min()) - lam1) < 1e-8, (pos.min(), lam1)
    sc = max(1.0, float(np.abs(H2.diagonal()).max()))
    okb = worst[0] >= -1e-9 * sc
    print(f"Q={Q} N={N} dim={dimN}: gamma_cert={gamma:.6f} | (a) min positive eig H={lam1:.6f} {'OK' if lam1 >= gamma else 'VIOLATION'}"
          f" | (b) min_j lambda_min(H^2-gH-K_Q)={worst[0]:.3e} (j={worst[1]}) {'OK' if okb else 'VIOLATION'}"
          f" best gamma' in sector={bestg:.6f}", flush=True)
