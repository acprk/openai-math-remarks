#!/usr/bin/env python3
"""269: rigorous evaluation of the finite-flux comparison for one x-ball (x = 1/Q) and one rational Gamma.

Certified inequality (README, Prop. C):  on the whole Fock space F_Q,
    H_Q^2 - gamma H_Q  >=  d_a * int K(g) dg  >= 0,
    gamma = 1 - eta - delta_Q - cost3 - sum_D n_D rho_D,
where every term is an arb upper bound and each rho_D is certified by a strict LDL^T in ball arithmetic of
    T_D = rho_D G_J + G_r[J,:] C_D G_r[:,J]       (m_D x m_D, m_D = dim of the highest-weight space)
(the r-space form of  W_D C_D W_D^T + rho_D What What^T >= 0 on wedge^4, see README "four-body reduction").
"""
import numpy as np
from fractions import Fraction
from math import ceil
from flint import arb, arb_mat, fmpq
import aengine as ae

_WD = {}


def wedge_cache(D):
    if D not in _WD:
        bl, As, N = ae.wedge_data(D)
        Z = ae.rowspace_basis(N)            # nb x k, fmpq
        _WD[D] = (bl, As, N, Z, ae.rank_W(D))
    return _WD[D]


def amid(a):
    return float(a.mid())


def fl(M):
    return np.array([[amid(v) for v in row] for row in M])


def arb_chol_lower(M):
    """plain Cholesky in arb (used only to build an exact preconditioner from midpoints)."""
    n = len(M)
    L = [[arb(0)] * n for _ in range(n)]
    for j in range(n):
        s = M[j][j] - sum((L[j][k] * L[j][k] for k in range(j)), arb(0))
        if not s > 0:
            raise ValueError("preconditioner Gram not PD at %d" % j)
        L[j][j] = s.sqrt()
        for i in range(j + 1, n):
            L[i][j] = (M[i][j] - sum((L[i][k] * L[j][k] for k in range(j)), arb(0))) / L[j][j]
    return L


def falling_ratio(n, a, b, x0, x):
    """prod_{i<n} (a-(b+i)x0)/(a-(b+i)x), tight for a ball x around x0"""
    r = arb(1)
    for i in range(n):
        r *= (a - (b + i) * x0) / (a - (b + i) * x)
    return r


def lam_ratio(bv, x0, x):
    p, j, k = bv
    return (falling_ratio(p, 2, 2, x0, x) * falling_ratio(j, 1, 0, x0, x) * falling_ratio(k, 1, 0, x0, x)).sqrt()


def evaluate(x, Gam, Qint=None, x0=None, rho_given=None, rho_margin=1e-12, rho_rel=1e-9, verbose=False):
    """x: arb ball (or exact point).  x0: exact center used for the preconditioner (default = x.mid()).
    rho_given: optional dict D -> Fraction (use these allowances instead of choosing them).
    Returns dict with gamma lower bound and all pieces."""
    if x0 is None:
        x0 = arb(x.mid())
    eng = ae.AEngine(x)
    out = {'fail': []}
    # ---------------- two body
    eta = arb(Gam.eta_q)
    # ---------------- three body
    cost3 = arb(0)
    slack = {}
    for z in eng.zlist:
        if z in (0, 1, 3):
            continue   # W_3 R_{3,z} = 0 (q_z = -1 exactly)
        r3 = eng.d3_over_da(z)
        ez = eng.ez(z, Gam)
        s = eng.qz(z) - ez / r3
        slack[z] = s
        if s > 0:
            continue
        neg = -s
        up = max(0.0, float(neg.upper()))
        cost3 += arb(up) * (z + 1) * r3
        # upper() is rounded outward by arb; use upper of product below
    cost3_up = float(cost3.upper()) * (1 + 1e-15) + 1e-300
    # ---------------- tail
    if Qint is not None:
        tl = arb(0)
        for z in range(17, Qint + 1, 2):
            tl += eng.d3_over_da(z) * abs(eng.qz(z)) * (z + 1)
    else:
        tl = eng.tail()
    tail_up = float(tl.upper()) * (1 + 1e-15)
    # ---------------- four body
    rho = {}
    cost4 = arb(0)
    info4 = {}
    for D in range(1, ae.DMAX + 1):
        bl, As, N, Z, m = wedge_cache(D)
        R = eng.Rset[D]
        nR = len(R)
        if m == 0:
            rho[D] = Fraction(0)
            info4[D] = (0, None)
            continue
        E = eng.E_of(D, Gam)
        C = [[-E[i][j] for j in range(nR)] for i in range(nR)]
        C[0][0] += eng.dD_over_da(D)                    # chi = e_{r=1} e_{r=1}^T
        S0 = [[eng.Sval(D, D, r, *b) or arb(0) for b in bl] for r in R]
        k = len(Z[0])
        # preconditioned basis of ran What^T = Lambda'' ran N^T :  Bn(x) = diag(lam(x)/lam(x0)) * B0 R0
        lam0 = [ae.lam_tilde(b, x0) for b in bl]
        B0 = [[lam0[ib] * arb(Z[ib][c]) for c in range(k)] for ib in range(len(bl))]
        G0 = (arb_mat(B0).transpose() * arb_mat(B0))
        L0 = arb_chol_lower(ae.to_rows(G0))
        Linv = arb_mat(L0).inv()
        R0 = arb_mat([[arb(Linv[j, i].mid()) for j in range(k)] for i in range(k)])  # exact ~ L0^{-T}
        U0 = arb_mat(B0) * R0
        rat = [lam_ratio(b, x0, x) for b in bl]
        Bn = arb_mat([[rat[ib] * U0[ib, c] for c in range(k)] for ib in range(len(bl))])
        Gn = Bn.transpose() * Bn
        Am = arb_mat(S0) * Bn                           # nR x k
        Gr = Am * Gn.solve(Am.transpose())              # = S Pi S^T  (nR x nR)
        Grf = fl(ae.to_rows(Gr))
        # choose J (m indices) by greedy pivoted Cholesky on the midpoint
        J = []
        Wk = Grf.copy()
        for _ in range(m):
            dg = np.array([Wk[i, i] if i not in J else -1 for i in range(nR)])
            i = int(np.argmax(dg))
            J.append(i)
            if Wk[i, i] <= 0:
                break
            v = Wk[:, i] / np.sqrt(Wk[i, i])
            Wk = Wk - np.outer(v, v)
        J = sorted(J)
        GJ = [[Gr[a, b] for b in J] for a in J]
        okG, _ = ae.ldl_pd(GJ)
        if not okG:
            out['fail'].append(('GJ', D))
            info4[D] = ('GJ-fail', None)
            rho[D] = None
            continue
        GrJ = arb_mat([[Gr[a, b] for b in range(nR)] for a in J])     # m x nR
        P = GrJ * arb_mat(C) * GrJ.transpose()
        Pf = fl(ae.to_rows(P))
        GJf = fl(GJ)
        Lg = np.linalg.cholesky(0.5 * (GJf + GJf.T))
        Li = np.linalg.inv(Lg)
        lam = np.linalg.eigvalsh(Li @ (0.5 * (Pf + Pf.T)) @ Li.T).min()
        if rho_given is not None:
            rD = rho_given[D]
        else:
            need = max(0.0, -lam)
            if need == 0.0 and lam > 1e3 * rho_margin:
                rD = Fraction(0)
            else:
                val = need * (1 + rho_rel) + rho_margin
                rD = Fraction(ceil(val * 1e15), 10 ** 15)
        T = [[arb(fmpq(rD.numerator, rD.denominator)) * GJ[a][b] + P[a, b] for b in range(m)] for a in range(m)]
        ok, piv = ae.ldl_pd(T)
        info4[D] = (m, lam)
        if not ok:
            out['fail'].append(('T', D, float(rD), lam))
            rho[D] = None
            continue
        rho[D] = rD
        cost4 += ae.nD(D) * arb(fmpq(rD.numerator, rD.denominator))
    if out['fail']:
        out['gamma_lo'] = None
    else:
        g = 1 - eta - arb(tail_up) - arb(cost3_up) - cost4
        out['gamma_lo'] = float(g.lower()) - 1e-14
    out.update(eta=float(eta.mid()), tail_up=tail_up, cost3_up=cost3_up, cost4=float(cost4.upper()), rho=rho,
               info4=info4, minslack=min((float(s.mid()), z) for z, s in slack.items()))
    if verbose:
        print(out)
    return out


def bspace_check(Q, Gam, rho_r, tiny=Fraction(1, 10 ** 12)):
    """Second, rank-free certificate of the four-body step at integer Q (b-space form):
         Am^T C Am + rho Gn > 0   (k x k, k = rank N),   Am = S_D Bn,  Gn = Bn^T Bn,  Bn = Lambda'' Z R0,
    which is exactly  What(S^T C S + rho I)What^T >= 0 restricted to ran What^T (no use of rank_W).
    rho = rho_r[D] (r-space allowance) + tiny.  Returns (ok_all, cost4_upper, failures)."""
    x = arb(fmpq(1, Q))
    eng = ae.AEngine(x)
    cost = arb(0)
    fails = []
    for D in range(1, ae.DMAX + 1):
        bl, As, N, Z, m = wedge_cache(D)
        if not Z or not Z[0]:
            continue
        R = eng.Rset[D]
        nR = len(R)
        E = eng.E_of(D, Gam)
        C = [[-E[i][j] for j in range(nR)] for i in range(nR)]
        C[0][0] += eng.dD_over_da(D)
        S0 = arb_mat([[eng.Sval(D, D, r, *b) or arb(0) for b in bl] for r in R])
        k = len(Z[0])
        lam0 = [ae.lam_tilde(b, x) for b in bl]
        B0 = arb_mat([[lam0[ib] * arb(Z[ib][c]) for c in range(k)] for ib in range(len(bl))])
        L0 = arb_chol_lower(ae.to_rows(B0.transpose() * B0))
        Linv = arb_mat(L0).inv()
        R0 = arb_mat([[arb(Linv[j, i].mid()) for j in range(k)] for i in range(k)])
        Bn = B0 * R0
        Am = S0 * Bn
        rD = (rho_r.get(D) or Fraction(0)) + tiny
        rq = arb(fmpq(rD.numerator, rD.denominator))
        T = Am.transpose() * arb_mat(C) * Am + (Bn.transpose() * Bn) * rq
        ok, _ = ae.ldl_pd(ae.to_rows(T))
        if not ok:
            fails.append(D)
        cost += ae.nD(D) * rq
    return (not fails), float(cost.upper()), fails
