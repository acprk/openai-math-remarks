# Exact certificate for a small instance of the angular cylinder cover of family 100
# (Finite-angular-cylinder-covers-below-the-half-area-bound, construction.tex:52-123) with coarse
# parameters: rational eps, n sectors per family (2n cylinders), uniform partition q_j=-1+2j/n,
# alpha_j,beta_j as in eq. (axis-coefficients), cutoffs T_j = 1/2 + eps^2 (eta + max(d(q_j),d(q_{j+1}))).
#
# Coverage (exact): the first-crossing selection (coverage.tex:13-45, valid for every eps) puts every
# point of K in an angular wedge W_j of family 1 and W_k of family 2.  For every pair (j,k) we certify
#   max { z : p in K cap W_j cap W_k, eps^2 z <= t0 - 1/2 - eps^2 d_j, eps^2 z <= s0 - 1/2 - eps^2 d_k } <= U_jk
# by an exact rational dual vector y >= 0 (y^T A_x = 0, y^T a_z = 1  =>  z <= y^T b).
# If eta > max U_jk, no point fails both selected cylinders, so the 2n cylinders cover K.
# Area (exact): S/A_min = sum_j (2/n) T_j^2 / sqrt(1+eps^2(alpha_j^2+2 beta_j^2)) bounded above with
# rational upper bounds r_j >= 1/sqrt(.) (checked by r_j^2 * (.) >= 1).
import numpy as np, sys, itertools, time
from fractions import Fraction as Fr
from scipy.optimize import linprog

def build(eps, n):
    qs = [Fr(-1) + Fr(2 * j, n) for j in range(n + 1)]
    a = [(1 + qs[j] * qs[j + 1]) / 4 for j in range(n)]
    b = [(qs[j] + qs[j + 1]) / 4 for j in range(n)]
    d = lambda q: q * q * (1 + q * q) / 16
    dj = [max(d(qs[j]), d(qs[j + 1])) for j in range(n)]
    return qs, a, b, dj

def system(eps, qs, a, b, dj, j, k):
    e2 = eps * eps
    A = [[0, 0, -1, 0], [0, 0, 1, 0], [1, 0, 1, 0], [-1, 0, 1, 0], [0, 1, -1, 0], [0, -1, -1, 0]]
    B = [0, 1, 1, 1, 0, 0]
    for qq, sg in ((qs[j], 1), (qs[j + 1], -1)):
        A.append([sg * (-qq * eps * b[j] + eps * a[j]), -sg, sg * qq, 0]); B.append(0)
    for qq, sg in ((qs[k], 1), (qs[k + 1], -1)):
        A.append([-sg, sg * (qq * eps * b[k] - eps * a[k]), -sg * qq, 0]); B.append(-sg * qq)
    A.append([eps * b[j], 0, -1, e2]); B.append(Fr(-1, 2) - e2 * dj[j])
    A.append([0, -eps * b[k], 1, e2]); B.append(Fr(1, 2) - e2 * dj[k])
    return A, B

def pair_bound(eps, qs, a, b, dj, j, k):
    A, B = system(eps, qs, a, b, dj, j, k)
    Af = np.array([[float(v) for v in r] for r in A]); Bf = np.array([float(v) for v in B])
    r = linprog([0, 0, 0, -1], A_ub=Af, b_ub=Bf, bounds=[(None, None)] * 4, method='highs')
    if r.status == 2:  # infeasible region (wedges disjoint inside K): certify by Farkas
        A2 = [r_[:3] + [-1] for r_ in A]; Af2 = np.array([[float(v) for v in r_] for r_ in A2])
        r2 = linprog([0, 0, 0, -1], A_ub=Af2, b_ub=Bf, bounds=[(None, None)] * 4, method='highs')
        du = -r2.ineqlin.marginals; idx = list(np.argsort(-du)[:6])
        for I in itertools.combinations(sorted(idx), 4):
            y = solve_exact([[A2[i][c] for i in I] for c in range(4)], [Fr(0), Fr(0), Fr(0), Fr(1)])
            if y is None or any(v < 0 for v in y): continue
            if sum(y[t] * B[I[t]] for t in range(4)) < 0: return Fr(-10**9), 'infeasible-certified'
        return None, 'infeasible-unc'
    duals = -r.ineqlin.marginals  # >= 0
    idx = list(np.argsort(-duals)[:6])
    # exact: find 4 constraints among top candidates whose exact dual solves A_I^T y = (0,0,0,1), y>=0
    best = None
    for I in itertools.combinations(sorted(idx), 4):
        M = [[A[i][c] for i in I] for c in range(4)]
        rhs = [Fr(0), Fr(0), Fr(0), Fr(1)]
        y = solve_exact(M, rhs)
        if y is None or any(v < 0 for v in y): continue
        ub = sum(y[t] * B[I[t]] for t in range(4))
        if best is None or ub < best: best = ub
    return best, -r.fun

def solve_exact(M, rhs):
    n = len(M); A = [row[:] + [rhs[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return None
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]

def area_ub(eps, n, eta, qs, a, b, dj):
    tot = Fr(0)
    for j in range(n):
        T = Fr(1, 2) + eps * eps * (eta + dj[j])
        u = 1 + eps * eps * (a[j] ** 2 + 2 * b[j] ** 2)
        r = Fr(1 / float(u) ** 0.5).limit_denominator(10 ** 12) + Fr(1, 10 ** 11)
        assert r * r * u >= 1
        tot += Fr(2, n) * T * T * r
    return tot

if __name__ == "__main__":
    eps = Fr(sys.argv[1]); n = int(sys.argv[2])
    qs, a, b, dj = build(eps, n)
    per = len(sys.argv) > 3
    if per:  # per-sector margins e_j (cyl_persector.py output), rounded UP to rationals; then need max U_jk < 0 with eta = 0
        ev = np.load(sys.argv[3])[0]
        eR = [Fr(float(v)).limit_denominator(10 ** 7) + Fr(1, 10 ** 7) for v in ev]
        dj = [dj[j] + eR[j] for j in range(n)]
    t0 = time.time(); U = None; worst = None; bad = 0; uncert = 0
    for j in range(n):
        for k in range(n):
            ub, fl = pair_bound(eps, qs, a, b, dj, j, k)
            if fl == 'infeasible-unc': uncert += 1; continue
            if fl == 'infeasible-certified': continue
            if ub is None: bad += 1; continue
            if U is None or ub > U: U = ub; worst = (j, k, fl)
    print(f"eps={eps} n={n} cylinders={2*n} pairs={n*n} exact max U_jk = {float(U):.8f} (float LP {worst[2]:.8f} at {worst[:2]}) "
          f"uncertified={bad} infeasible-pairs-uncertified={uncert} time={time.time()-t0:.0f}s", flush=True)
    if per:
        print(f"per-sector mode: max U_jk = {float(U):.3e} < 0 : {U < 0}  (coverage certified iff True)", flush=True)
        eta = Fr(0)
    else:
        eta = Fr(float(U)).limit_denominator(10 ** 6) + Fr(1, 10 ** 6)
        assert eta > U
    S = area_ub(eps, n, eta, qs, a, b, dj)
    print(f"eta = {eta} > max U_jk ; exact area upper bound S/A_min <= {float(S):.10f} ; 1/2 - S >= {float(Fr(1,2)-S):.3e} ; BELOW HALF: {S < Fr(1,2)}", flush=True)
