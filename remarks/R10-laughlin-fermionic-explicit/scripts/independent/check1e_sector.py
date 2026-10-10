#!/usr/bin/env python3
"""Item 1(e): independent sector test of the certified operator inequality  H_Q^2 - gamma H_Q - K_Q >= 0,
K_Q = d_a * Haar average of K (paper eq. (Adecomp), 04-comparison l.170-174), K = sum_t sum_{m,m'} Gamma_t[m,m'] O_m^dag O_m',
O_0 = B_t, O_(p,j) = c^dag_{p+j-t} c_j B_p (paper eq. (F), l.143-148), Gamma_t = the certificate's rational data.
Pair vectors from Clebsch-Gordan (indep_fock), own Jordan-Wigner, Haar average by Schur's lemma on highest-weight spaces
(04-comparison l.111-120): on the spin-J multiplicity space, K_Q = d_a/(2J+1) * sum_mu <f_mu, K f_mu>, f_mu = normalised
J_-^mu h  (normalisation prod sqrt((k+1)(2J-k)), no projections).  gamma = certified gamma_lo(Q) from cert/out/certQ_all.jsonl.
usage: check1e_sector.py Q N1,N2,..."""
import sys, json, itertools
from fractions import Fraction
import numpy as np, scipy.sparse as sp, scipy.linalg as sla
from indep_fock import pair_vectors_cg, Fock

import os
CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'primary') + '/'
Q = int(sys.argv[1]); Ns = [int(x) for x in sys.argv[2].split(',')]
FLIP = len(sys.argv) > 3 and sys.argv[3] == 'flip'   # negative control: v_p -> (-1)^p v_p (wrong phase convention)
rec = {json.loads(l)['Q']: json.loads(l) for l in open(CERT + 'out/certQ_all.jsonl')}[Q]
gamma = rec['gamma_lo']
d = json.load(open(CERT + rec['gamma_file']))
Gam = {}
for t, rows in d['F'].items():
    Fm = np.array([[float(Fraction(v)) for v in r] for r in rows])
    Gam[int(t)] = Fm @ Fm.T
S = {t: [(p, j) for p in range(8) for j in range(9) if 0 <= p + j - t <= 8] for t in range(8)}   # paper S_rho
V = pair_vectors_cg(Q)
if FLIP: V = {p: {k: (-1) ** p * a for k, a in dd.items()} for p, dd in V.items()}
F = Fock(Q)
da = 2 * Q - 1
for N in Ns:
    F.basis(N); F.basis(N - 1); F.basis(N - 2)
    n = len(F.bases[N])
    B = {p: F.B(V[p], N) for p in V}
    H = sum((B[p].T @ B[p] for p in B)).tocsr()
    K = sp.csr_matrix((n, n))
    for t in range(8):
        ops = [B[t]] + [(F.c(p + j - t, N - 2).T @ F.c(j, N - 2) @ B[p]) if N >= 3 else None for (p, j) in S[t]]
        G = Gam[t]
        assert G.shape[0] == len(ops)
        w, U = np.linalg.eigh(0.5 * (G + G.T))
        for k in range(len(w)):
            if w[k] <= 1e-300: continue
            g = U[:, k] * np.sqrt(w[k])
            Fk = None
            for m, o in enumerate(ops):
                if o is None or g[m] == 0: continue
                Fk = g[m] * o if Fk is None else Fk + g[m] * o
            if Fk is not None: K = K + Fk.T @ Fk
    K = K.tocsr()
    Jm = F.lowering(N); Jp = Jm.T.tocsr()
    lev = {}
    for i, s in enumerate(F.bases[N]): lev.setdefault(sum(s), []).append(i)
    worst, minpos, ndim, bestg = np.inf, np.inf, 0, np.inf
    # per weight level: eigen-decomposition of the Casimir J^2 = J_- J_+ + Jz^2 + Jz (dense, level-sized); used to
    # (i) get the spin-J highest-weight space and (ii) re-project every descendant onto spin J (removes the round-off
    # contamination by higher spins, which J_-^mu would otherwise amplify by up to sqrt(C(2J',2J'-2J)))
    cas = {}
    for s0 in sorted(lev):
        idx = lev[s0]; jz = (Q * N - 2 * s0) / 2.0
        A = Jp[lev[s0 - 1]][:, idx].toarray() if s0 - 1 in lev else np.zeros((0, len(idx)))
        C2 = A.T @ A + (jz * jz + jz) * np.eye(len(idx))
        w_, U_ = np.linalg.eigh(C2)
        cas[s0] = (w_, U_)
    def spinJ(s0, J2):
        w_, U_ = cas[s0]; JJ = J2 / 2.0
        return U_[:, np.abs(w_ - JJ * (JJ + 1)) < 1e-6]
    for s0 in sorted(lev):
        J2 = Q * N - 2 * s0          # 2J
        if J2 < 0: break
        idx = lev[s0]
        ns = spinJ(s0, J2)
        if ns.shape[1] == 0: continue
        h = np.zeros((n, ns.shape[1])); h[idx] = ns
        ndim += ns.shape[1]
        f = h.copy(); Tr = f.T @ (K @ f)
        for k in range(J2):
            g_ = (Jm @ f) / np.sqrt((k + 1) * (J2 - k))
            E = spinJ(s0 + k + 1, J2); li = lev[s0 + k + 1]
            f = np.zeros_like(f); f[li] = E @ (E.T @ g_[li])
            Tr += f.T @ (K @ f)
        KQ = da * Tr / (J2 + 1)
        Hh = h.T @ (H @ h); H2h = (H @ h).T @ (H @ h)
        X = H2h - gamma * Hh - KQ
        worst = min(worst, np.linalg.eigvalsh(0.5 * (X + X.T)).min())
        wh = np.linalg.eigvalsh(0.5 * (Hh + Hh.T))
        if (wh > 1e-9).any(): minpos = min(minpos, wh[wh > 1e-9].min())
        wv, Uv = np.linalg.eigh(0.5 * (Hh + Hh.T)); P = Uv[:, wv > 1e-9]
        if P.shape[1]:
            Y = P.T @ (H2h - KQ) @ P; Z = P.T @ Hh @ P
            bestg = min(bestg, sla.eigh(0.5 * (Y + Y.T), 0.5 * (Z + Z.T), eigvals_only=True).min())
    print(f"Q={Q} N={N} dim={n} HWtotal={ndim} gamma_cert={gamma:.6f} (file {rec['gamma_file']}) | "
          f"min_J lambda_min(H^2-gH-K_Q) = {worst:.3e} | min positive eig H = {minpos:.6f} "
          f"best gamma' in sector = {bestg:.6f} " f"{'OK' if worst > -1e-9 and minpos >= gamma else 'VIOLATION'}", flush=True)
