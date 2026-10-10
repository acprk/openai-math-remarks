#!/usr/bin/env python3
"""269 certificate engine in arb ball arithmetic (python-flint).

Everything is written in the variable x = 1/Q.  For an integer flux Q >= 24 pass x = arb(fmpq(1,Q)) (exact);
for Q >= 300 one may pass a ball x = [x1,x2] subset of [0,1/300]: every formula below is a rational/algebraic
function of x that is analytic on [0,1/24] and agrees with the finite-flux object at x = 1/Q, so ball
evaluation encloses all Q with 1/Q in the ball (x = 0 is the Q -> infinity limit).

Index conventions follow engine.py (float engine used by the SDP) and the paper (sections 3-5):
  Support: pmax=7, M=8, t in 0..7;  S_t = {(p,j): p<=pmax, j<=M, 0<=p+j-t<=M};  Gamma_t on [B_t] + S_t.
All index guards of the float engine (p<=n1, m<=n2, k<=Q, ...) are automatically satisfied for Q>=24
(asserted in comments; see README section "finite-Q side conditions").
"""
from math import factorial
from fractions import Fraction
import itertools
import numpy as np
from flint import arb, arb_mat, fmpq, fmpz, fmpq_mat, ctx

ctx.prec = 256
PMAX, MM, TSET = 7, 8, list(range(8))
T3MAX = PMAX + MM          # 15
DMAX = PMAX + 2 * MM       # 23


def S_t(t):
    return [(p, j) for p in range(PMAX + 1) for j in range(MM + 1) if 0 <= p + j - t <= MM]


def dim_t(t):
    return 1 + len(S_t(t))


def A(v):
    """exact arb from int / Fraction / fmpq"""
    if isinstance(v, Fraction):
        return arb(fmpq(v.numerator, v.denominator))
    return arb(v)


# ----------------------------------------------------------------- coupling (Lemma 3.3, finite flux)
def couple_x(a1, b1, a2, b2, z, Tmax, x):
    """Normalised coupled vectors of spin J=(n1+n2)/2 - z in V_{n1} (x) V_{n2}, n_i = a_i Q - b_i, levels z..Tmax.
    Same vectors as engine.couple (sign (-1)^z at p=0, normalised lowering), but computed without per-step
    renormalisation: the highest vector is a closed product, descendants are J_-^{T-z} applied unnormalised and
    divided by the exact norm prod_{mu<T-z} sqrt((mu+1)(2J-mu)).  All quantities divided by the appropriate power
    of Q, so the formulas are analytic at x=0 (Q=infinity).  Valid when all indices stay inside the
    representations, which holds for every use below once Q>=24."""
    c = {}
    acc = arb(1)
    for p in range(0, z + 1):
        m = z - p
        c[(p, m)] = (arb(-1) ** (z + p)) * acc.sqrt()
        if p < z:
            acc *= (z - p) * (a2 - (b2 + z - p - 1) * x) / ((p + 1) * (a1 - (b1 + p) * x))
    nrm = sum((v * v for v in c.values()), arb(0)).sqrt()
    c = {k: v / nrm for k, v in c.items()}
    out = {z: c}
    u = c
    scale = arb(1)
    for T in range(z, Tmax):
        mu = T - z
        new = {}
        for (p, m), v in u.items():
            w1 = v * ((p + 1) * (a1 - (b1 + p) * x)).sqrt()
            new[(p + 1, m)] = new.get((p + 1, m), arb(0)) + w1
            w2 = v * ((m + 1) * (a2 - (b2 + m) * x)).sqrt()
            new[(p, m + 1)] = new.get((p, m + 1), arb(0)) + w2
        scale *= ((mu + 1) * ((a1 + a2) - (b1 + b2 + 2 * z + mu) * x)).sqrt()
        u = new
        out[T + 1] = {k: v / scale for k, v in u.items()}
    return out


def falling_x(n, a, b, x):
    """prod_{i<n} (a - (b+i) x)"""
    r = arb(1)
    for i in range(n):
        r *= (a - (b + i) * x)
    return r


# ----------------------------------------------------------------- Gamma
class Gamma:
    """Gamma_t = F_t F_t^T with F_t rational (so Gamma_t >= 0 and K = sum_t sum_cols (row F)^dag (row F) is an exact SOS)."""

    def __init__(self, F):   # F: dict t -> list of rows (list of Fraction), shape dim_t x rank
        self.F = F
        self.G = {}
        for t in TSET:
            Fm = fmpq_mat([[fmpq(v.numerator, v.denominator) for v in row] for row in F[t]])
            Gq = Fm * Fm.transpose()
            self.G[t] = [[arb(Gq[i, j]) for j in range(Gq.ncols())] for i in range(Gq.nrows())]
            self.Gq = getattr(self, 'Gq', {})
            self.Gq[t] = Gq
        self.eta_q = sum((self.Gq[t][0, 0] for t in TSET), fmpq(0))

    def float(self):
        return {t: np.array([[float(self.Gq[t][i, j]) for j in range(self.Gq[t].ncols())]
                             for i in range(self.Gq[t].nrows())]) for t in TSET}


# ----------------------------------------------------------------- the engine
class AEngine:
    def __init__(self, x):
        self.x = x
        self._three()
        self._four()

    # ---- scalars
    def da_over(self, d_num_a, d_num_b):  # (d_num_a Q - d_num_b)/(2Q-1)
        x = self.x
        return (d_num_a - d_num_b * x) / (2 - x)

    def qz(self, z):
        x = self.x
        r = arb(1)
        for i in range(z):
            r *= (1 - i * x) / (2 - (2 + i) * x)
        return arb(-1) ** z * r * (3 * z - 1 - z * (z + 1) * x)

    def d3_over_da(self, z):   # (3Q-1-2z)/(2Q-1)
        return self.da_over(3, 1 + 2 * z)

    def dD_over_da(self, D):   # (4Q-1-2D)/(2Q-1)
        return self.da_over(4, 1 + 2 * D)

    def tail(self, Zc=250):
        """upper bound for delta_Q = sum_{odd 17<=z<=Q} (d3z/da)|q_z|(z+1), valid for every Q>=max(Zc,4) with 1/Q in x.
        terms z<=Zc evaluated exactly (z<=Zc<=Q so the block exists); z>Zc bounded by 6z(z+1)(2/3)^z (paper, Lemma 6.1)."""
        s = arb(0)
        for z in range(17, Zc + 1, 2):
            s += self.d3_over_da(z) * (-self.qz(z)) * (z + 1)   # z odd: q_z < 0 (bracket >= 2z-2 > 0)
        # sum_{z>Zc} 6 z (z+1) (2/3)^z  (all z, upper bound): closed form via ratio bound
        z0 = Zc + 1
        r = arb(2) / 3
        # term_z <= term_{z0} * rho^{z-z0}, rho = (2/3)*((z0+2)/z0)  (ratio of successive terms <= rho for z>=z0)
        rho = r * arb(z0 + 2) / z0
        term0 = 6 * z0 * (z0 + 1) * r ** z0
        s += term0 / (1 - rho)
        return s

    # ---- three body
    def _three(self):
        x = self.x
        self.zlist = list(range(0, T3MAX + 1))
        self.xi3 = {z: couple_x(2, 2, 1, 0, z, T3MAX, x) for z in self.zlist}   # V_Q (2Q-2) x U_Q (Q)

    def ez(self, z, Gam):
        """trace of the unaveraged 3-body matrix in spin block z (paper eq. (A3)); exact finite-flux value."""
        vecs = self.xi3[z]
        tot = arb(0)
        for t in TSET:
            S = S_t(t)
            G = Gam.G[t]
            for T in range(max(z, t), T3MAX + 1):
                xi = vecs[T]
                u = [xi.get((t, T - t), arb(0))]
                idx = [0]
                for k, (p, j) in enumerate(S):
                    if p + j == T:
                        u.append(xi.get((p, j), arb(0)))
                        idx.append(k + 1)
                n = len(u)
                for a in range(n):
                    for b in range(n):
                        if a == 0 and b == 0:
                            continue
                        g = G[idx[a]][idx[b]]
                        tot += g * u[a] * u[b]
        return tot

    # ---- four body
    def _four(self):
        x = self.x
        self.f = {}
        for r in range(1, DMAX + 1, 2):
            vecs = couple_x(1, 0, 1, 0, r, DMAX, x)       # U_Q x U_Q, spin Q-r summand
            self.f[r] = {n - r: {(a, b): arb(2).sqrt() * v for (a, b), v in c.items() if a < b}
                         for n, c in vecs.items()}
        self.Rset = {D: list(range(1, D + 1, 2)) for D in range(1, DMAX + 1)}
        self.Scoef = {}
        for D in range(1, DMAX + 1):
            for r in self.Rset[D]:
                outer = couple_x(2, 2, 2, 2 * r, D - r, DMAX - r, x)   # V_Q x (spin Q-r), deficit D-r
                for T in range(D, DMAX + 1):
                    d = {}
                    for (p, m), sig in outer[T - r].items():
                        for (j, k), fv in self.f[r].get(m, {}).items():
                            d[(p, j, k)] = d.get((p, j, k), arb(0)) + sig * fv
                    self.Scoef[(D, T, r)] = d

    def Sval(self, D, T, r, p, j, k):
        if j == k:
            return None
        if j < k:
            return self.Scoef[(D, T, r)].get((p, j, k))
        v = self.Scoef[(D, T, r)].get((p, k, j))
        return None if v is None else -v

    def E_of(self, D, Gam):
        R = self.Rset[D]
        nR = len(R)
        E = [[arb(0)] * nR for _ in range(nR)]
        for t in TSET:
            S = S_t(t)
            G = Gam.G[t]
            # group slots by level n = p+j
            lev = {}
            for a, (p, j) in enumerate(S):
                lev.setdefault(p + j, []).append((a, p, j))
            for na, LA in lev.items():
                for nb, LB in lev.items():
                    T = na + nb - t
                    if T < D:
                        continue
                    i = na - t
                    k = nb - t
                    X = []   # rows a: vector over r
                    for (a, p, j) in LA:
                        X.append([self.Sval(D, T, r, p, j, k) for r in R])
                    Y = []
                    for (b, q, l) in LB:
                        Y.append([self.Sval(D, T, s, q, l, i) for s in R])
                    # E_rs -= sum_{a,b} G_ab X[a][r] Y[b][s]
                    # first Z[a][s] = sum_b G_ab Y[b][s]
                    Zs = []
                    for ia, (a, p, j) in enumerate(LA):
                        if all(v is None for v in X[ia]):
                            Zs.append(None)
                            continue
                        row = [arb(0)] * nR
                        for ib, (b, q, l) in enumerate(LB):
                            g = G[a + 1][b + 1]
                            if g == 0:
                                continue
                            yb = Y[ib]
                            for s in range(nR):
                                if yb[s] is not None:
                                    row[s] += g * yb[s]
                        Zs.append(row)
                    for ia in range(len(LA)):
                        if Zs[ia] is None:
                            continue
                        xa = X[ia]
                        for rr in range(nR):
                            if xa[rr] is None:
                                continue
                            for s in range(nR):
                                E[rr][s] -= xa[rr] * Zs[ia][s]
        return E


# ----------------------------------------------------------------- Q-independent integer data for the wedge map
def wedge_data(D):
    """b-list, wedge list A (increasing 4-tuples, sum D+1), integer matrix N[A,b] = sigma*(x-y) with
    v_p^*-pattern (Lemma 3.1): u_b = v_p ^ e_j ^ e_k = Delta'_A * N * Lambda''_b exactly at every Q>=24."""
    bl = [(p, j, k) for p in range(D + 1) for j in range(D + 1) for k in range(j + 1, D + 1) if p + j + k == D]
    As = sorted({tuple(sorted(c)) for c in itertools.combinations(range(D + 2), 4) if sum(c) == D + 1})
    Ai = {a: i for i, a in enumerate(As)}
    N = [[0] * len(bl) for _ in As]
    for ib, (p, j, k) in enumerate(bl):
        for xx in range(0, (p + 1) // 2 + 1):
            yy = p + 1 - xx
            if xx >= yy:
                continue
            lst = [xx, yy, j, k]
            if len(set(lst)) < 4:
                continue
            sgn = 1
            a = lst[:]
            for i1 in range(4):
                for i2 in range(3 - i1):
                    if a[i2] > a[i2 + 1]:
                        a[i2], a[i2 + 1] = a[i2 + 1], a[i2]
                        sgn = -sgn
            N[Ai[tuple(a)]][ib] += sgn * (xx - yy)
    return bl, As, N


def rowspace_basis(N):
    """rational basis (columns) of ran N^T = row space of N."""
    if not N or not N[0]:
        return []
    M = fmpq_mat([[fmpq(v) for v in row] for row in N])
    R, rk = M.rref()
    return [[R[i, c] for i in range(rk)] for c in range(M.ncols())]   # nb x rk


def hw_mult(D):
    """dim of highest-weight space of spin 2Q-1-D in wedge^4 U_Q (Q>=24): c(D+1)-c(D), c(n)=#4-subsets sum n."""
    def c(n):
        return sum(1 for cc in itertools.combinations(range(n + 1), 4) if sum(cc) == n)
    return c(D + 1) - c(D)


def zm_mult(D):
    """dim of the highest-weight space of spin 2Q-1-D inside ker H_{4,Q} (Q>=24).  By the divisibility argument of
    Lemma 6.2 (paper) ker H_{4,Q} = prod_{i<j} b_ij^3 * Sym^4(U_{Q-9}); the Laughlin factor is SU(2)-invariant, so this is
    the HW count of spin 2(Q-9)-(D-17) in Sym^4 U_{Q-9}: p4(d)-p4(d-1), d=D-17, p4 = partitions into <=4 parts."""
    def p4(n):
        if n < 0:
            return 0
        return sum(1 for a in range(n + 1) for b in range(a, n + 1) for c in range(b, n + 1)
                   if n - a - b - c >= c)
    d = D - 17
    return p4(d) - p4(d - 1)


def rank_W(D):
    """rank of the wedged copy vectors {w_{Dr}}: = dim(HW_{2Q-1-D} cap ran What) = hw_mult - zm_mult, since
    ran What = (ker H_4)^perp in wedge^4 U_Q (B(v_p)psi=0 for all p <=> psi orthogonal to v_p ^ wedge^2)."""
    return hw_mult(D) - zm_mult(D)


def nD(D):
    return (D + 1) ** 2 // 4


def lam_tilde(b, x):
    """Q^{(D+1)/2} Lambda''_b = sqrt(p! j! k! / (prod_{i<p}(2-(2+i)x) prod_{i<j}(1-ix) prod_{i<k}(1-ix)))"""
    p, j, k = b
    return (arb(factorial(p) * factorial(j) * factorial(k)) /
            (falling_x(p, 2, 2, x) * falling_x(j, 1, 0, x) * falling_x(k, 1, 0, x))).sqrt()


# ----------------------------------------------------------------- small arb linear algebra
def mat(rows):
    return arb_mat(rows)


def to_rows(M):
    return [[M[i, j] for j in range(M.ncols())] for i in range(M.nrows())]


def ldl_pd(M):
    """certify a (ball) symmetric matrix positive definite by LDL^T in ball arithmetic.
    Returns (ok, list_of_pivots).  ok=True => every symmetric matrix in the ball is positive definite."""
    n = len(M)
    Aa = [row[:] for row in M]
    piv = []
    for k in range(n):
        d = Aa[k][k]
        piv.append(d)
        if not d > 0:
            return False, piv
        for i in range(k + 1, n):
            l = Aa[i][k] / d
            for j in range(k + 1, n):
                Aa[i][j] -= l * Aa[k][j]
    return True, piv
