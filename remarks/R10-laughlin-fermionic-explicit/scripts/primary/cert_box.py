#!/usr/bin/env python3
"""269: rigorous certificate on an x-box [x0-h, x0+h] (x = 1/Q) by second-order Taylor models.

Every quantity f(x) used by the certificate is analytic on [0,1/24].  We evaluate the whole pipeline in truncated
power-series arithmetic (python-flint arb_series, 3 coefficients) twice:
  (A) at the exact centre x0              -> f0 = f(x0), f1 = f'(x0)               (tight balls)
  (B) at the ball base point X=[x0-h,x0+h] -> F2 encloses f''(xi)/2 for every xi in X  (Taylor coefficient of a
      series computed at a ball base point encloses the coefficient at every point of the ball).
Lagrange:  f(x0+t) = f0 + f1 t + (f''(xi)/2) t^2,  |t|<=h.
 * scalars: upper bounds via f0 + |f1| h + |F2| h^2.
 * matrices: T(t) >= T0 + T1 t - h^2 ||F2||_F I; the right side is affine in t, so strict positivity at t=+-h
   (certified by ball LDL^T) gives T(x) > 0 on the whole box.
Four-body reduction per D (see README):  K = {r in R_D : r <= 16, not (r=1 and D even)};  copies outside K have
w_{Dr}=0 (r=1, D even) or an identically zero row/column of C_D (r>=17), so W C W^T = W_K C_KK W_K^T.
rank W_K <= mK := min(rank_W(D), |K|);  J subset K with |J| = mK and G_J > 0 on the box gives rank = mK, and then
   W_K C_KK W_K^T + rho What What^T >= 0   <=>   rho G_J + G[J,K] C_KK G[K,J] >= 0.
"""
import numpy as np
from fractions import Fraction
from math import ceil
from flint import arb, arb_series, arb_mat, fmpq, ctx
import aengine as ae
import cert_core as cc

NC = 7


def co(s, i):
    if isinstance(s, arb_series):
        c = s.coeffs()
        return c[i] if i < len(c) else arb(0)
    return s if i == 0 else arb(0)


# ---------------- polynomial (in t) matrices: list of arb_mat coefficients, truncated at t^2
def pm_from(rows):   # rows: list of list of scalars/series
    n, m = len(rows), len(rows[0])
    return [arb_mat([[co(rows[i][j], c) for j in range(m)] for i in range(n)]) for c in range(NC)]


def pm_mul(A, B):
    return [sum((A[a] * B[c - a] for a in range(c + 1)), arb_mat(A[0].nrows(), B[0].ncols())) for c in range(NC)]


def pm_T(A):
    return [M.transpose() for M in A]


def pm_add(A, B):
    return [A[c] + B[c] for c in range(NC)]


def pm_scale(s, A):   # scalar series s times matrix
    return [sum((A[c - a] * co(s, a) for a in range(c + 1)), arb_mat(A[0].nrows(), A[0].ncols())) for c in range(NC)]


def pm_inv(A):
    X = [A[0].inv()]
    for c in range(1, NC):
        acc = sum((A[a] * X[c - a] for a in range(1, c + 1)), arb_mat(A[0].nrows(), A[0].ncols()))
        X.append(-(X[0] * acc))
    return X


def pm_sub(A, rows, cols):
    return [arb_mat([[M[i, j] for j in cols] for i in rows]) for M in A]


def series_x(xbase):
    return arb_series([xbase, arb(1)])


def pipeline(xbase, x0, Gam, J_choice=None):
    """Evaluate all certificate ingredients as 3-term Taylor series in t at base point xbase (x = xbase + t)."""
    ctx.cap = NC
    x = series_x(xbase)
    x0 = arb(fmpq(x0.numerator, x0.denominator))      # exact centre (Fraction -> arb)
    eng = ae.AEngine(x)
    res = {}
    # three body slacks, as cost-relevant series:  c_z = -(q_z - e_z / r3) * r3 * (z+1)
    res['c3'] = {}
    for z in eng.zlist:
        if z in (0, 1, 3):
            continue
        r3 = eng.d3_over_da(z)
        s = eng.qz(z) - eng.ez(z, Gam) / r3
        res['c3'][z] = -s * r3 * (z + 1)
    res['tail'] = eng.tail()
    res['four'] = {}
    for D in range(1, ae.DMAX + 1):
        bl, As, N, Z, m = cc.wedge_cache(D)
        R = eng.Rset[D]
        K = [i for i, r in enumerate(R) if r <= 16 and not (r == 1 and D % 2 == 0)]
        mK = min(m, len(K))
        if mK == 0:
            continue
        E = eng.E_of(D, Gam)
        Crows = [[-E[i][j] for j in K] for i in K]
        if R[K[0]] == 1:
            Crows[0][0] = Crows[0][0] + eng.dD_over_da(D)
        C = pm_from(Crows)
        S0 = pm_from([[eng.Sval(D, D, R[i], *b) or arb(0) for b in bl] for i in K])
        k = len(Z[0])
        # preconditioner at the exact centre x0 (Q-independent rational Z)
        lam0 = [ae.lam_tilde(b, x0) for b in bl]
        B0 = arb_mat([[lam0[ib] * arb(Z[ib][c]) for c in range(k)] for ib in range(len(bl))])
        key = (D, str(x0))
        if key not in _PRE:
            G0 = B0.transpose() * B0
            L0 = cc.arb_chol_lower(ae.to_rows(G0))
            Linv = arb_mat(L0).inv()
            _PRE[key] = arb_mat([[arb(Linv[j, i].mid()) for j in range(k)] for i in range(k)])
        R0 = _PRE[key]
        U0 = B0 * R0
        rat = [cc.lam_ratio(b, x0, x) for b in bl]          # series ratio lam(x)/lam(x0)
        Bn = [arb_mat([[co(rat[ib], c) * U0[ib, cc_] for cc_ in range(k)] for ib in range(len(bl))]) for c in range(NC)]
        Gn = pm_mul(pm_T(Bn), Bn)
        Am = pm_mul(S0, Bn)
        Gr = pm_mul(pm_mul(Am, pm_inv(Gn)), pm_T(Am))         # |K| x |K|
        if J_choice is None:
            Grf = np.array([[float(Gr[0][i, j].mid()) for j in range(len(K))] for i in range(len(K))])
            J = []
            Wk = Grf.copy()
            for _ in range(mK):
                dg = np.array([Wk[i, i] if i not in J else -1 for i in range(len(K))])
                i = int(np.argmax(dg)); J.append(i)
                if Wk[i, i] <= 0:
                    break
                v = Wk[:, i] / np.sqrt(Wk[i, i]); Wk = Wk - np.outer(v, v)
            J = sorted(J)
        else:
            J = J_choice[D]
        GJ = pm_sub(Gr, J, J)
        GrJ = pm_sub(Gr, J, range(len(K)))
        P = pm_mul(pm_mul(GrJ, C), pm_T(GrJ))
        res['four'][D] = dict(J=J, GJ=GJ, P=P, mK=mK)
    return res


_PRE = {}


def EYE(m):
    return arb_mat([[arb(1) if i == j else arb(0) for j in range(m)] for i in range(m)])


def fnorm_up(M):
    s = arb(0)
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            v = M[i, j]
            u = max(abs(float(v.upper())), abs(float(v.lower())))
            s += arb(u) * arb(u)
    return float(s.sqrt().upper()) * (1 + 1e-12)


def hi_terms(Apt, Bball, hf):
    """sum_{i=2}^{K-1} ||A_i|| h^i + ||B_K|| h^K  (matrix Taylor remainder bound), K = NC-1"""
    K = NC - 1
    return (sum(fnorm_up(Apt[i]) * hf ** i for i in range(2, K)) + fnorm_up(Bball[K]) * hf ** K) * (1 + 1e-12)


def upper_scalar(f0s, F2s, h):
    """upper bound of f over the box: f0 + |f1| h + sum_{2<=i<K}|f_i(x0)| h^i + |F_K(ball)| h^K."""
    K = NC - 1
    hh = arb(fmpq(h.numerator, h.denominator))
    up = co(f0s, 0) + abs(co(f0s, 1)) * hh
    for i in range(2, K):
        up += abs(co(f0s, i)) * hh ** i
    up += abs(co(F2s, K)) * hh ** K
    return float(up.upper()) * (1 + 1e-15) + 1e-300


def certify_box(xa, xb, Gam, rho_margin=1e-10):
    """xa<xb Fractions. Returns dict (gamma_lo or None, rho, details)."""
    xm = (xa + xb) / 2
    h = (xb - xa) / 2
    x0 = arb(fmpq(xm.numerator, xm.denominator))
    X = arb(fmpq(xa.numerator, xa.denominator)).union(arb(fmpq(xb.numerator, xb.denominator)))
    A = pipeline(x0, xm, Gam)
    Jc = {D: v['J'] for D, v in A['four'].items()}
    B = pipeline(X, xm, Gam, J_choice=Jc)
    hh = arb(fmpq(h.numerator, h.denominator))
    out = {'fail': []}
    # scalars
    eta = arb(Gam.eta_q)
    tail_up = upper_scalar(A['tail'], B['tail'], h)
    cost3 = 0.0
    for z in A['c3']:
        u = upper_scalar(A['c3'][z], B['c3'][z], h)
        cost3 += max(0.0, u)
    cost3 *= (1 + 1e-14)
    rho = {}
    cost4 = arb(0)
    hf = float(h)
    for D, a in A['four'].items():
        b = B['four'][D]
        m = a['mK']
        # (i) G_J > 0 on the box
        g2 = hi_terms(a['GJ'], b['GJ'], hf)
        okG = True
        for sg in (1, -1):
            M = a['GJ'][0] + a['GJ'][1] * (sg * hh) - EYE(m) * arb(g2)
            ok, _ = ae.ldl_pd(ae.to_rows(M))
            okG = okG and ok
        if not okG:
            out['fail'].append(('GJ', D)); continue
        # (ii) choose rho from float data at the two endpoints, then certify
        def fl(M):
            return np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])
        need = 0.0
        p2 = hi_terms(a['P'], b['P'], hf)
        for sg in (1, -1):
            Pf = fl(a['P'][0]) + sg * hf * fl(a['P'][1]) - p2 * np.eye(m)
            Gf = fl(a['GJ'][0]) + sg * hf * fl(a['GJ'][1])
            Li = np.linalg.inv(np.linalg.cholesky(0.5 * (Gf + Gf.T)))
            lam = np.linalg.eigvalsh(Li @ (0.5 * (Pf + Pf.T)) @ Li.T).min()
            need = max(need, -lam)
        rq = None
        cands = ([Fraction(0)] if need < -1e-9 else []) + \
            [Fraction(ceil((max(need, 0.0) * (1 + 1e-6) + dm) * 1e15), 10 ** 15) for dm in (rho_margin, 1e-8, 1e-7, 1e-6, 1e-5)]
        okT = False
        for rD in cands:
            rq = arb(fmpq(rD.numerator, rD.denominator))
            okT = True
            for sg in (1, -1):
                g2r = float(rD) * g2 + p2
                M = (a['GJ'][0] + a['GJ'][1] * (sg * hh)) * rq + a['P'][0] + a['P'][1] * (sg * hh) - EYE(m) * arb(g2r)
                ok, _ = ae.ldl_pd(ae.to_rows(M))
                okT = okT and ok
            if okT:
                break
        if not okT:
            out['fail'].append(('T', D, float(rD))); continue
        rho[D] = rD
        cost4 += ae.nD(D) * rq
    out.update(eta=float(eta.mid()), tail_up=tail_up, cost3_up=cost3, cost4=float(cost4.upper()), rho=rho)
    if out['fail']:
        out['gamma_lo'] = None
    else:
        g = 1 - eta - arb(tail_up) - arb(cost3) - cost4
        out['gamma_lo'] = float(g.lower()) - 1e-14
    return out
