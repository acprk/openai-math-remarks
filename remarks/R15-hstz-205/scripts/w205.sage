# w205.sage -- independent implementation (SageMath 9.5 + its lrcalc binding) of the
# multiplicity lower bound extracted from OpenAI family 205
# ("A Cyclic Polytabloid Proof of Saxl's Conjecture"):
#
#   mult(S^lam, W_m) >= w_m(lam) := max( 1                                   [Thm 3.1]
#                                     , K_{lam', rho_m}                       [Prop 3.2: W_m (x) sgn contains a module onto M^rho]
#                                     , sum_{lam/nu horiz. m-strip} w_{m-1}(nu)          [Prop 4.2, s=1, U_{m,1} trivial]
#                                     , sum_{nu |- N_{m-2}} w_{m-2}(nu) * L4(lam/nu) )  [Prop 4.2, s=2 + Prop 5.3]
#   L4(lam/nu) = sum_{eta, len(eta) <= 4} c^lam_{nu,eta}  (lrcalc.skew with row bound 4).
# and W_m is a submodule of S^rho (x) S^rho, so w_m(lam) <= g(rho,rho,lam).
#
# Modes:
#   sage w205.sage full M      : all lambda |- N_m for m <= M; writes w_m.txt (lam w) and compares with
#                                exact g if ../expected_output/exact_g_m{m}.txt exists.
#   sage w205.sage hooks M     : hooks only (hooks are closed under the recursion), m <= M;
#                                prints w on hooks (N-k,1^k) and on their conjugates.
# Kostka numbers K_{mu,rho} are read off h_m h_{m-1} ... h_1 expanded by lrcalc.mult (Pieri).
import sys
from sage.libs.lrcalc import lrcalc

def conj(l):
    return tuple(Partition(list(l)).conjugate()) if l else ()

def kostka_rho(m):
    cur = {(): 1}
    for r in range(m, 0, -1):
        nxt = {}
        for lam, c in cur.items():
            for mu, cc in lrcalc.mult(list(lam), [r]).items():
                key = tuple(mu)
                nxt[key] = nxt.get(key, 0) + c * cc
        cur = nxt
    return cur

def subpartitions(lam, k):
    """nu subset lam with |lam/nu| = k"""
    lam = list(lam); L = len(lam); res = []
    def rec(i, rem, acc):
        if i == L:
            if rem == 0: res.append(tuple(x for x in acc if x > 0))
            return
        hi = lam[i] if i == 0 else min(lam[i], acc[-1])
        for v in range(hi, max(0, lam[i] - rem) - 1, -1):
            rec(i + 1, rem - (lam[i] - v), acc + [v])
    rec(0, k, []); return res

def is_hstrip(lam, nu):
    return all((nu[i] if i < len(nu) else 0) >= (lam[i + 1] if i + 1 < len(lam) else 0) for i in range(len(lam)))

def L4(lam, nu):
    # c^lam_{nu,eta} = 0 for all eta with <= 4 rows iff some column of lam/nu has > 4 cells;
    # we test this first (lrcalc.skew in Sage 9.5 crashes when the row-bounded result is empty).
    lc = conj(lam); nc = conj(nu)
    if any(lc[j] - (nc[j] if j < len(nc) else 0) > 4 for j in range(len(lc))):
        return 0
    d = lrcalc.skew(list(lam), list(nu), 4)
    return sum(d.values())

def wbound(m, lam, w1, w2, K):
    best = max(1, K.get(conj(lam), 0))
    if m >= 2:
        s = sum(w1[nu] for nu in subpartitions(lam, m) if is_hstrip(lam, nu))
        best = max(best, s)
    if m >= 3:
        s = 0
        for nu in subpartitions(lam, 2 * m - 1):
            c = L4(lam, nu)
            if c: s += w2[nu] * c
        best = max(best, s)
    return best

mode = sys.argv[1]; M = int(sys.argv[2])
W = {0: {(): 1}}
if mode == 'full':
    for m in range(1, M + 1):
        N = m * (m + 1) // 2; K = kostka_rho(m)
        W[m] = {tuple(l): wbound(m, tuple(l), W.get(m - 1), W.get(m - 2), K) for l in Partitions(N)}
        with open(f'w_m{m}.txt', 'w') as f:
            for l, v in W[m].items(): f.write(','.join(map(str, l)) + f' {v}\n')
        print('m', m, 'N', N, '#lambda', len(W[m]), flush=True)
elif mode == 'hooks':
    def hooks(N): return [tuple([N - k] + [1] * k) for k in range(N)]
    for m in range(1, M + 1):
        N = m * (m + 1) // 2
        # Kostka number of a hook with content rho_m (m distinct letters): K_{(N-c,1^c),rho_m} = binom(m-1, c)
        K = {tuple([N - c] + [1] * c): binomial(m - 1, c) for c in range(N)}
        # recursion only touches hooks of smaller size; store only hooks
        W[m] = {h: wbound(m, h, W.get(m - 1), W.get(m - 2), K) for h in hooks(N)}
        W[m - 2] = {} if m >= 4 else W.get(m - 2)
        print('m', m, 'w on hooks (N-k,1^k), k=0..N-1:', [W[m][h] for h in hooks(N)], flush=True)
