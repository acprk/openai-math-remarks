#!/usr/bin/env sage -python
"""Independent re-implementation (adversarial review, 2026-10-09) for OpenAI preprint
'A Torsion-Free Group Algebra with Zero Divisors' (family 196), Lemma 2.1 (squared words).
Written from the paper (types.tex), NOT from kz196_spectral.py / kz196_exact.py.

For q = 2^a we
 (1) build PG(2,q) from GF(q) in Sage and the letter alphabet T = points + 3 extras with an
     explicit fixed-point-free involution (extras paired with 3 ordinary letters);
 (2) build the vertex types of Section 2.1 at the smallest admissible n = 2 v^2, with exact
     multiplicities per line class, and verify: integrality, letter frequency p on both sides,
     odd cross-intersections |S_x cap S_y| for all type pairs;
 (3) compute the turn weights w(t,u) EMPIRICALLY from the B-side type counts (frequency of u
     among B-vertices containing bar t), check A-side frequencies are <= w, and form
     M_{t,u} = w(t,u)^2 (0 for u = bar t) as exact rationals;
 (4) q=16: certificate max_t (M f)(t)/f(t) < 1 with the paper's f (5 if bar t extra, else 1);
     q=8 : Collatz--Wielandt lower bound min_t (M g)(t)/g(t) > 1 with a rational positive g;
 (5) for every q, the exact Perron root via the 3-class equitable partition, cross-checked
     against numpy eigenvalues of the full matrix (q <= 16).
Usage: sage -python kz196_indep.py
"""
from fractions import Fraction as Fr
from sage.all import GF, QQ, PolynomialRing, matrix, AA
import itertools, numpy as np

def plane(q):
    F = GF(q, 'z'); els = list(F)
    pts = []
    for x in els:
        for y in els: pts.append((x, y, F(1)))
    for x in els: pts.append((x, F(1), F(0)))
    pts.append((F(1), F(0), F(0)))
    lines = [frozenset(i for i, P in enumerate(pts) if a*P[0] + b*P[1] + c*P[2] == 0) for (a, b, c) in pts]
    return pts, lines

def build(q):
    pts, lines = plane(q); v = len(pts); assert v == q*q + q + 1
    assert all(len(l) == q + 1 for l in lines)
    inter = {len(l1 & l2) for l1, l2 in itertools.combinations(lines, 2)}
    assert inter == {1}, inter                               # two distinct lines meet in one point
    p = Fr(q + 1, v)
    nT = v + 3; X = [v, v + 1, v + 2]                       # extras
    # involution: extras <-> ordinary letters 0,1,2 ; remaining ordinary letters paired (3,4),(5,6),...
    bar = list(range(nT))
    for k in range(3): bar[X[k]] = k; bar[k] = X[k]
    assert (v - 3) % 2 == 0
    for i in range(3, v, 2): bar[i] = i + 1; bar[i + 1] = i
    assert all(bar[bar[t]] == t and bar[t] != t for t in range(nT))
    # vertex types at n = 2 v^2 : per line class 2v vertices
    n = 2 * v * v; per = n // v
    A_types, B_types = [], []                                # (frozenset S, multiplicity)
    for l in lines:
        for pair in itertools.combinations(X, 2):
            m = p / 2 * per; assert m.denominator == 1; A_types.append((l | frozenset(pair), int(m)))
        m0 = per - 3 * int(p / 2 * per); assert m0 >= 0; A_types.append((l, m0))
        m = p * per; assert m.denominator == 1; B_types.append((l | frozenset(X), int(m)))
        B_types.append((l, per - int(m)))
    for side in (A_types, B_types):
        assert sum(m for _, m in side) == n
        cnt = [0] * nT
        for S, m in side:
            for t in S: cnt[t] += m
        assert all(Fr(c, n) == p for c in cnt), 'letter frequency'
    # odd cross-intersections (all type pairs)
    Al = sorted({S for S, _ in A_types}, key=sorted); Bl = sorted({S for S, _ in B_types}, key=sorted)
    if q <= 16:                                              # exhaustive over all type pairs
        assert all(len(SA & SB) % 2 == 1 for SA in Al for SB in Bl)
    else:
        import random; rng = random.Random(0)
        assert all(len(rng.choice(Al) & rng.choice(Bl)) % 2 == 1 for _ in range(200000))
    # empirical joint counts J_side[a][u] = #vertices containing a and u
    def joint(side):
        J = np.zeros((nT, nT), dtype=np.int64)
        for S, m in side:
            idx = np.fromiter(S, dtype=np.int64)
            J[np.ix_(idx, idx)] += m
        return J
    JA, JB = joint(A_types), joint(B_types)
    w = [[Fr(0)] * nT for _ in range(nT)]
    for t in range(nT):
        a = bar[t]
        for u in range(nT):
            if u == a: continue                                 # cancelling turn
            w[t][u] = Fr(int(JB[a][u]), int(JB[a][a]))
            assert Fr(int(JA[a][u]), int(JA[a][a])) <= w[t][u]
    vals = sorted({x for r in w for x in r if x > 0})
    M = [[x * x for x in r] for r in w]
    return dict(q=q, v=v, p=p, bar=bar, X=X, M=M, wvals=vals, nT=nT)

def quotient_rho(q):
    """Exact Perron root via equitable partition O0 (t,bar t ordinary), O1 (bar t extra), E (t extra)."""
    v = q*q + q + 1; p = QQ(q + 1) / v; o = QQ(1) / (q + 1)**2; e = p * p
    Q = matrix(QQ, [[(v - 4) * o, 3 * o, 3 * e],
                    [(v - 3) * e, 3 * e, 2],
                    [(v - 3) * o, 2 * o, 3 * e]])
    return max(r for r in Q.charpoly().roots(AA, multiplicities=False)), Q

def row_ratio(M, f):
    return max(sum(M[t][u] * f[u] for u in range(len(f))) / f[t] for t in range(len(f)))

def cw_lower(M, g):
    return min(sum(M[t][u] * g[u] for u in range(len(g))) / g[t] for t in range(len(g)))

if __name__ == '__main__':
    print('positive control q=128 (class rows, exact):')
    q = 128; v = q*q + q + 1; p = Fr(q + 1, v)
    rO0 = Fr(v - 4, (q + 1)**2) + Fr(3 * 5, (q + 1)**2) + 3 * p * p
    rO1 = ((v - 3) * p * p + 3 * 5 * p * p + 2) / 5
    print('  row O0 =', rO0, rO0 == Fr(500731261911, 504183783481))
    print('  row O1 =', rO1, rO1 == Fr(820350863, 1363395845))
    for q in [2, 4, 8, 16, 32]:
        d = build(q); M = d['M']; nT = d['nT']; bar = d['bar']; X = d['X']
        f = [Fr(5) if bar[t] in X else Fr(1) for t in range(nT)]
        lamf = row_ratio(M, f)
        rho, Qm = quotient_rho(q)
        # check the 3-class partition is equitable on the actual matrix and its quotient equals Qm
        cls = lambda t: 2 if t in X else (1 if bar[t] in X else 0)
        for t in range(nT):
            rs = [Fr(0)] * 3
            for u in range(nT): rs[cls(u)] += M[t][u]
            assert all(QQ(rs[j].numerator) / rs[j].denominator == Qm[cls(t)][j] for j in range(3)), (q, t)
        A = np.array([[float(x) for x in r] for r in M])
        rnum = max(abs(np.linalg.eigvals(A))) if nT <= 400 else float('nan')
        line = f"q={q:3d} v={d['v']:5d} |T|={nT:5d}  paper-f max(Mf)/f={float(lamf):.6f} ({'<1' if lamf < 1 else '>=1'})  rho(exact)={float(rho):.6f}  rho(numpy)={rnum:.6f}  min/max weight={float(d['wvals'][0]):.5f}/{float(d['wvals'][-1])}"
        print(line, flush=True)
        if q == 16:
            print('  q=16 certificate exact value:', lamf, ' sum f =', sum(f))
        if q in (8, 4, 2):
            # rational positive test vector from the quotient Perron vector (class-constant)
            ev = (Qm - rho).right_kernel().basis()[0]
            ev = [QQ(round(float(abs(c)) * 10**9) + 1) / 10**9 for c in ev]
            g = [Fr(int(ev[cls(t)].numerator()), int(ev[cls(t)].denominator())) for t in range(nT)]
            mu = cw_lower(M, g)
            print(f'  q={q} Collatz-Wielandt: min (Mg)/g = {float(mu):.6f} ({mu}) > 1: {mu > 1}')
    print('Perron root of the type-construction matrix (formula, any projective plane of order q):')
    for q in [2, 4, 8, 10, 12, 14, 16, 32, 64, 128]:
        print(f'  q={q:3d}: rho={float(quotient_rho(q)[0]):.6f}')
