#!/usr/bin/env python3
"""269 finite-flux engine (float64).  Builds, for a given flux Q and auxiliary support,
  * pair vectors v_p (eq:pair),
  * 3-body block functionals  e_z^Q(Gamma)   (trace of the unaveraged 3-body matrix in spin block z),
  * 4-body copy matrices      E_D^Q(Gamma)   (partial trace in block D, r,s in R_D),
  * the T=D coupling matrix S^Q (rows r, cols b=(p,j,k), j<k) and the finite-Q Gram of the
    uncoupled 4-wedges  hatG^Q_D = < v_p^ e_j ^ e_k , v_q ^ e_l ^ e_m >.
Gamma is a dict t -> PSD matrix indexed by [B_t] + S_t, S_t = {(p,j): p<=pmax, j<=M, 0<=p+j-t<=M}.
K = sum_t (row vectors of Gamma_t factorisation) reproduces the paper's F_rho when Gamma_t = sum g g^T.
"""
import numpy as np
from math import lgamma, sqrt, log, exp
from functools import lru_cache

import os
ROWS = os.environ.get('OPENAI_MATH', '../../../../../../../openai-math') + '/preprints/A-Fock-space-inequality-and-the-Laughlin-spectral-gap-September-24-2026/verification/rows.txt'


def lff(s, p):  # log of falling factorial (s)_p, s>=p
    if p < 200:
        return sum(log(s - i) for i in range(p))
    return lgamma(s + 1) - lgamma(s - p + 1)


def pair_vec(Q, p):
    """v_p as dict {(x,y): coef}, x<y, x+y=p+1, x,y<=Q (eq:pair)."""
    out = {}
    for x in range(0, (p + 1) // 2 + 1):
        y = p + 1 - x
        if x >= y or y > Q:
            continue
        lv = lff(Q, x) + lff(Q, y) + lgamma(p + 1) - log(Q) - lff(2 * Q - 2, p) - lgamma(x + 1) - lgamma(y + 1)
        out[(x, y)] = (x - y) * exp(0.5 * lv)
    return out


def pair_vec_star(p):
    out = {}
    for x in range(0, (p + 1) // 2 + 1):
        y = p + 1 - x
        if x >= y:
            continue
        out[(x, y)] = (x - y) * exp(0.5 * (lgamma(p + 1) - p * log(2) - lgamma(x + 1) - lgamma(y + 1)))
    return out


def couple(n1, n2, z, Tmax):
    """Coupled vectors of spin (n1+n2)/2 - z in V_{n1} (x) V_{n2}, levels T=z..Tmax.
    returns list over T of dict {(p,m): c} with p+m=T (normalized, sign (-1)^z at p=0)."""
    c = {}
    val = (-1.0) ** z
    for p in range(0, z + 1):
        m = z - p
        if p <= n1 and m <= n2:
            c[(p, m)] = val
        if p < z:
            num = (z - p) * (n2 - z + p + 1)
            den = (p + 1) * (n1 - p)
            if den <= 0 or num < 0:
                val = 0.0
            else:
                val = -val * sqrt(num / den)
    nrm = sqrt(sum(v * v for v in c.values()))
    c = {k: v / nrm for k, v in c.items()}
    out = {z: c}
    for T in range(z, Tmax):
        new = {}
        for (p, m), v in c.items():
            if p + 1 <= n1:
                new[(p + 1, m)] = new.get((p + 1, m), 0.0) + v * sqrt((p + 1) * (n1 - p))
            if m + 1 <= n2:
                new[(p, m + 1)] = new.get((p, m + 1), 0.0) + v * sqrt((m + 1) * (n2 - m))
        nrm = sqrt(sum(v * v for v in new.values()))
        c = {k: v / nrm for k, v in new.items()}
        out[T + 1] = c
    return out


def parse_rows():
    rows = []
    cur = None
    for line in open(ROWS):
        if not line.strip():
            continue
        toks = line.split()
        if not line.startswith(' '):
            cur = {'t': int(toks[0]), 'l': int(toks[1]), 'a': {}}
            rows.append(cur)
            toks = toks[2:]
        for tk in toks:
            k, v = tk.split(':')
            cur['a'][(int(k[0]), int(k[1]))] = int(v)
    return rows


class Support:
    def __init__(self, pmax=7, M=8, tset=range(8)):
        self.pmax, self.M, self.tset = pmax, M, list(tset)
        self.S = {}
        for t in self.tset:
            self.S[t] = [(p, j) for p in range(pmax + 1) for j in range(M + 1) if 0 <= p + j - t <= M]
        self.T3max = pmax + M
        self.Dmax = pmax + 2 * M

    def dim(self, t):
        return 1 + len(self.S[t])


def paper_gamma(sup):
    P = 1e7
    G = {t: np.zeros((sup.dim(t), sup.dim(t))) for t in sup.tset}
    for r in parse_rows():
        t = r['t']
        g = np.zeros(sup.dim(t))
        g[0] = r['l'] / P
        idx = {m: k + 1 for k, m in enumerate(sup.S[t])}
        for (p, j), a in r['a'].items():
            al = a / P * sqrt(2.0 ** (p - t) * exp(lgamma(t + 1) + lgamma(p + j - t + 1) - lgamma(p + 1) - lgamma(j + 1)))
            g[idx[(p, j)]] = al
        G[t] += np.outer(g, g)
    return G


class Engine:
    def __init__(self, Q, sup):
        self.Q, self.sup = Q, sup
        self.da = 2 * Q - 1
        self.v = {p: pair_vec(Q, p) for p in range(0, min(2 * Q - 1, sup.Dmax + 1))}
        self._three()
        self._four_coupling()

    # ---------------- three body ----------------
    def _three(self):
        Q, sup = self.Q, self.sup
        self.zlist = [z for z in range(0, sup.T3max + 1) if z <= Q]
        self.C3 = {}   # z -> t -> matrix (dim x dim):  e_z = sum_t <C3[z][t], Gamma_t>
        for z in self.zlist:
            vecs = couple(2 * Q - 2, Q, z, sup.T3max)
            self.C3[z] = {}
            for t in sup.tset:
                S = sup.S[t]
                Cm = np.zeros((sup.dim(t), sup.dim(t)))
                for T in range(max(z, t), sup.T3max + 1):
                    xi = vecs[T]
                    u = np.zeros(sup.dim(t))
                    u[0] = xi.get((t, T - t), 0.0)
                    for k, (p, j) in enumerate(S):
                        if p + j == T:
                            u[k + 1] = xi.get((p, j), 0.0)
                    Cm += np.outer(u, u)
                    Cm[0, 0] -= u[0] ** 2
                self.C3[z][t] = Cm

    def qz(self, z):
        Q = self.Q
        r = 1.0
        for i in range(z):
            r *= (Q - i) / (2 * Q - 2 - i)
        return (-1) ** z * r * (3 * z - 1 - z * (z + 1) / Q)

    def d3(self, z):
        return 3 * self.Q - 1 - 2 * z

    def delta_tail(self, z0=None):
        """three-body tail cost for blocks z>z0 (default T3max): sum over odd z with q_z<0."""
        Q = self.Q
        z0 = self.sup.T3max if z0 is None else z0
        s = 0.0
        for z in range(z0 + 1, min(Q, 600) + 1):   # terms <= 6z(z+1)(2/3)^z: beyond 600 < 1e-95
            q = self.qz(z)
            if q < 0:
                s += self.d3(z) / self.da * (-q) * (z + 1)
        return s

    # ---------------- four body ----------------
    def _four_coupling(self):
        Q, sup = self.Q, self.sup
        Dmax = sup.Dmax
        # spin Q-r summands of wedge^2 U_Q, as wedge-coefficient dicts per level n=j+k
        self.f = {}
        for r in range(1, Dmax + 1, 2):
            if r > Q:
                break
            vecs = couple(Q, Q, r, Dmax)
            fr = {}
            for n, c in vecs.items():
                d = {}
                for (x, y), val in c.items():
                    if x < y:
                        d[(x, y)] = sqrt(2.0) * val
                fr[n - r] = d   # index m = n - r
            self.f[r] = fr
        self.Rset = {}
        self.Scoef = {}   # (D,T,r) -> dict {(p,j,k): val}, j<k
        for D in range(1, Dmax + 1):
            R = [r for r in range(1, D + 1, 2) if r in self.f]
            self.Rset[D] = R
            for r in R:
                outer = couple(2 * Q - 2, 2 * Q - 2 * r, D - r, Dmax - r)
                for T in range(D, Dmax + 1):
                    d = {}
                    for (p, m), sig in outer[T - r].items():
                        for (j, k), fv in self.f[r][m].items():
                            d[(p, j, k)] = d.get((p, j, k), 0.0) + sig * fv
                    self.Scoef[(D, T, r)] = d

    def Sval(self, D, T, r, p, j, k):
        if j == k:
            return 0.0
        if j < k:
            return self.Scoef[(D, T, r)].get((p, j, k), 0.0)
        return -self.Scoef[(D, T, r)].get((p, k, j), 0.0)

    def E_of(self, D, G):
        """E_D^Q for a given Gamma (dict t->matrix)."""
        R = self.Rset[D]
        nR = len(R)
        E = np.zeros((nR, nR))
        for t in self.sup.tset:
            S = self.sup.S[t]
            Gt = G[t]
            for a, (p, j) in enumerate(S):
                for b, (q, l) in enumerate(S):
                    g = Gt[a + 1, b + 1]
                    if g == 0.0:
                        continue
                    i = p + j - t
                    k = q + l - t
                    T = p + j + k
                    if T < D:
                        continue
                    for ir, r in enumerate(R):
                        x = self.Sval(D, T, r, p, j, k)
                        if x == 0.0:
                            continue
                        for js, s in enumerate(R):
                            y = self.Sval(D, T, s, q, l, i)
                            if y:
                                E[ir, js] -= g * x * y
        return E

    def bspace(self, D):
        """b=(p,j,k), j<k, p+j+k=D; S^Q at T=D (nR x nb); 4-wedge vectors; Gram."""
        Q = self.Q
        bl = [(p, j, k) for p in range(0, D + 1) for j in range(0, D + 1) for k in range(j + 1, D + 1)
              if p + j + k == D and k <= Q and p <= 2 * Q - 2]
        R = self.Rset[D]
        Smat = np.array([[self.Sval(D, D, r, *b) for b in bl] for r in R]).reshape(len(R), len(bl))
        vecs = [wedge4(self.v[p], j, k) for (p, j, k) in bl]
        return bl, Smat, vecs

    def nD(self, D):
        return (D + 1) ** 2 // 4


def wedge4(vp, j, k):
    """v_p ^ e_j ^ e_k as dict over sorted 4-tuples."""
    out = {}
    for (x, y), c in vp.items():
        lst = [x, y, j, k]
        if len(set(lst)) < 4:
            continue
        # sign of sorting permutation
        sgn = 1
        a = lst[:]
        for i1 in range(4):
            for i2 in range(3 - i1):
                if a[i2] > a[i2 + 1]:
                    a[i2], a[i2 + 1] = a[i2 + 1], a[i2]
                    sgn = -sgn
        key = tuple(a)
        out[key] = out.get(key, 0.0) + sgn * c
    return out


def gram(vecs):
    keys = sorted(set().union(*[set(v) for v in vecs])) if vecs else []
    ki = {k: i for i, k in enumerate(keys)}
    A = np.zeros((len(keys), len(vecs)))
    for c, v in enumerate(vecs):
        for k, x in v.items():
            A[ki[k], c] = x
    return A


def range_basis(A, tol=1e-9):
    """orthonormal basis (in column space of A^T, i.e. b-space) of ran(A^T A)."""
    if A.size == 0:
        return np.zeros((A.shape[1], 0))
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    rk = int((s > tol * max(1.0, s.max())).sum())
    return Vt[:rk].T
