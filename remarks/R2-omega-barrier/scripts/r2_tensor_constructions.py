#!/usr/bin/env python3
"""
R2 tensor-level check: the degenerations behind the inequalities are verified at
the level of coefficients for small sizes (independent re-implementation).

Part 1 (monomial weight degenerations; support level):
  * OpenAI Lemma 4.2 (shifted tripling) of C(a, 3h+a-1)          paper.tex l.533-602
  * scalar peeling C(n,n) ~> 1 (+) C(n-1,n-1)                     CODEX-ROUND3 (5)
  * truncated convolution ~> big CW_{n-2}                         CODEX-ROUND2 (14)
  * zero/rest partition support S_0                               CODEX-ROUND3 (7)
Part 2 (D-adapted bases; coefficient level, exact rationals):
  bivariate multiplication W_{a-1,c-1} x W_{b-1,d-1} -> W_{a+b-2,c+d-2}, with
  W_{e,f} = V_e(u,v) (x) V_f(s,w) and D = uw - vs.  Checks that in D-adapted bases
  the product of grades (i,j) has no output component of grade < i+j and that its
  grade-(i+j) component is exactly the convolution C(A_i,B_j) in quotient indices.
  c=1, d=m gives C(a,b) (x) B_X(m) (the paper's Lemma 4.1 is m=2; family (B));
  general (a,b,c,d) gives families (C) and (D).  For (D) the second-stage weights
  2i^2, 2j^2, -k^2 are checked on the computed support.

Environment: Python >= 3.8, standard library only.
Run:  /usr/bin/python3 r2_tensor_constructions.py
Expected output: expected_output/r2_tensor_constructions.txt  (runtime < 1 min)
"""
from fractions import Fraction as Fr
from itertools import product
import sys

fails = 0
def report(name, n, bad):
    global fails
    fails += bad
    print(f"{name:62s} cases={n:6d} failures={bad}")

# ------------------------------------------------------------------ Part 1
def C_support(a, b):
    return [(i, j, i + j) for i in range(a) for j in range(b)]

def degenerate(support, wx, wy, wz):
    """Return (min weight, weight-min sub-support) for a monomial weighting."""
    ws = [wx(i) + wy(j) + wz(k) for (i, j, k) in support]
    m = min(ws)
    return m, sorted(t for t, w in zip(support, ws) if w == m)

# OpenAI tripling, all 1<=a,h<=12
n = bad = 0
for a in range(1, 13):
    for h in range(1, 13):
        Bv = 3 * h + a - 1
        midY = lambda y: h <= y <= 2 * h + a - 2
        midZ = lambda z: h + a - 1 <= z <= 2 * h + a - 2
        m, S = degenerate(C_support(a, Bv), lambda i: 0,
                          lambda j: 1 if midY(j) else 0, lambda k: -1 if midZ(k) else 0)
        left = sorted((x, y, x + y) for x in range(a) for y in range(h))
        right = sorted((x, y + 2 * h + a - 1, x + y + 2 * h + a - 1) for x in range(a) for y in range(h))
        # middle: x_{a-1-u} y_{h+u+r} z_{h+a-1+r}, u<a, r<h  (C with legs 2,3 swapped, leg 1 reversed)
        middle = sorted((a - 1 - u, h + u + r, h + a - 1 + r) for u in range(a) for r in range(h))
        good = (m == 0 and S == sorted(left + middle + right))
        # sector matching: Y-sector determines Z-sector on surviving terms
        sec = lambda y: 0 if y < h else (1 if y <= 2 * h + a - 2 else 2)
        secz = lambda z: 0 if z <= h + a - 2 else (1 if z <= 2 * h + a - 2 else 2)
        good &= all(sec(y) == secz(z) for (_, y, z) in S)
        n += 1; bad += (not good)
report("OpenAI Lemma 4.2 tripling degeneration (support)", n, bad)

# alternating k-sector chain (k odd), shift (k-1)(a-1)/2: Y-sectors of lengths
# h1, h2+a-1, h3, h4+a-1, ...; Z-sectors h1+a-1, h2, h3+a-1, ...; odd-position
# sectors get Y weight +1, Z weight -1 (k=3 is OpenAI's Lemma 4.2).
import random
random.seed(1)
n = bad = 0
for k in (3, 5, 7):
    for a in range(1, 7):
        for _ in range(40):
            hs = [random.randint(1, 5) for _ in range(k)]
            ylen = [hs[s] + (a - 1 if s % 2 else 0) for s in range(k)]
            zlen = [hs[s] + (0 if s % 2 else a - 1) for s in range(k)]
            ystart = [sum(ylen[:s]) for s in range(k)]
            zstart = [sum(zlen[:s]) for s in range(k)]
            Bv = sum(ylen)
            assert Bv == sum(hs) + (k - 1) * (a - 1) // 2 and sum(zlen) == Bv + a - 1
            ysec = lambda y: max(s for s in range(k) if ystart[s] <= y)
            zsec = lambda z: max(s for s in range(k) if zstart[s] <= z)
            m, S = degenerate(C_support(a, Bv), lambda i: 0,
                              lambda j: ysec(j) % 2, lambda z: -(zsec(z) % 2))
            good = m == 0 and all(ysec(y) == zsec(z) for (_, y, z) in S)
            # each sector s carries exactly a*h_s terms, i.e. a full C(a,h_s) (or its sigma-twin)
            for s in range(k):
                good &= sum(1 for (_, y, z) in S if ysec(y) == s) == a * hs[s]
            n += 1; bad += not good
report("alternating k-sector chains k=3,5,7 (support)", n, bad)

# scalar peeling: x0,y0 weight 1; z0 weight -2
n = bad = 0
for N_ in range(2, 30):
    m, S = degenerate(C_support(N_, N_), lambda i: 1 if i == 0 else 0,
                      lambda j: 1 if j == 0 else 0, lambda k: -2 if k == 0 else 0)
    exp = sorted([(0, 0, 0)] + [(i, j, i + j) for i in range(1, N_) for j in range(1, N_)])
    n += 1; bad += not (m == 0 and S == exp)
report("scalar peeling C(n,n) ~> 1 (+) C(n-1,n-1)", n, bad)

# truncated convolution t_n ~> CW_{n-2}
n = bad = 0
for N_ in range(3, 40):
    w = lambda i: 0 if i == 0 else (2 if i == N_ - 1 else 1)
    supp = [(i, j, i + j) for i in range(N_) for j in range(N_) if i + j < N_]
    m, S = degenerate(supp, w, w, lambda k: -w(k))
    exp = sorted({(0, j, j) for j in range(N_)} | {(i, 0, i) for i in range(N_)}
                 | {(i, N_ - 1 - i, N_ - 1) for i in range(1, N_ - 1)})
    n += 1; bad += not (m == 0 and S == exp)
report("truncated convolution ~> big CW_{n-2} (support)", n, bad)

# zero/rest partition of C(n,n): X,Y {0}|[1,n-1] weights (0,1); Z {0}|[1,n-1]|[n,2n-2] weights (0,-1,-2)
n = bad = 0
for N_ in range(2, 30):
    sx = lambda i: 0 if i == 0 else 1
    sz = lambda k: 0 if k == 0 else (1 if k <= N_ - 1 else 2)
    m, S = degenerate(C_support(N_, N_), sx, sx, lambda k: -sz(k))
    sectors = sorted({(sx(i), sx(j), sz(k)) for (i, j, k) in S})
    exp_sec = [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 2)]
    ok = (m == 0 and sectors == exp_sec)
    # block 112 is the window i+j >= n with i,j>=1 (truncated convolution after reversal)
    blk = sorted((i, j) for (i, j, k) in S if sz(k) == 2)
    ok &= blk == sorted((i, j) for i in range(1, N_) for j in range(1, N_) if i + j >= N_)
    n += 1; bad += not ok
report("zero/rest partition support S_0 = {000,011,101,112}", n, bad)

# ------------------------------------------------------------------ Part 2
# polynomials in (u,v,s,w): dict {(pu,pv,ps,pw): Fraction}
def pmul(p, q):
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            k = tuple(x + y for x, y in zip(m1, m2))
            out[k] = out.get(k, 0) + c1 * c2
    return {k: c for k, c in out.items() if c}

def ppow(p, e):
    r = {(0, 0, 0, 0): Fr(1)}
    for _ in range(e):
        r = pmul(r, p)
    return r

Dpoly = {(1, 0, 0, 1): Fr(1), (0, 1, 1, 0): Fr(-1)}      # uw - vs

def adapted_basis(e, f):
    """D-adapted basis of W_{e,f}: list of (grade k, quotient index l, polynomial)."""
    B = []
    for k in range(min(e, f) + 1):
        Dk = ppow(Dpoly, k)
        for l in range(e + f - 2 * k + 1):
            al = min(l, e - k)
            be = l - al
            mono = {(e - k - al, al, f - k - be, be): Fr(1)}
            B.append((k, l, pmul(Dk, mono)))
    return B

def monomials(e, f):
    return [(e - x, x, f - y, y) for x in range(e + 1) for y in range(f + 1)]

def solve_coords(basis, e, f):
    """Return function mapping a polynomial in W_{e,f} to its coordinates in basis."""
    mons = monomials(e, f)
    idx = {m: r for r, m in enumerate(mons)}
    nB = len(basis)
    assert nB == len(mons)
    # matrix columns = basis vectors; invert by Gauss-Jordan
    A = [[Fr(0)] * nB for _ in range(nB)]
    for c, (_, _, p) in enumerate(basis):
        for m, val in p.items():
            A[idx[m]][c] = val
    # augment with identity
    M = [row[:] + [Fr(int(r == c)) for c in range(nB)] for r, row in enumerate(A)]
    for col in range(nB):
        piv = next(r for r in range(col, nB) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(nB):
            if r != col and M[r][col] != 0:
                fct = M[r][col]
                M[r] = [x - fct * y for x, y in zip(M[r], M[col])]
    inv = [row[nB:] for row in M]
    def coords(p):
        vec = [Fr(0)] * nB
        for m, val in p.items():
            vec[idx[m]] += val
        return [sum(inv[r][c] * vec[c] for c in range(nB)) for r in range(nB)]
    return coords

n = bad = 0
supports = {}
for a, b, c, d in product(range(1, 5), repeat=4):
    if (a + b + c + d) > 13:
        continue
    BX = adapted_basis(a - 1, c - 1)
    BY = adapted_basis(b - 1, d - 1)
    eo, fo = a + b - 2, c + d - 2
    BZ = adapted_basis(eo, fo)
    coords = solve_coords(BZ, eo, fo)
    supp = set()
    for (i, l, px), (j, m, py) in product(BX, BY):
        co = coords(pmul(px, py))
        ok = True
        for (k, deg, _), val in zip(BZ, co):
            if val != 0:
                supp.add((i, j, k))
            if k < i + j and val != 0:
                ok = False
            if k == i + j and val != int(deg == l + m):
                ok = False
        n += 1; bad += not ok
    supports[(a, b, c, d)] = supp
report("D-adapted bivariate product: grade<i+j vanish, grade=i+j is C", n, bad)

# (D) second stage: weights K(-i-j+k) + (2i^2 + 2j^2 - k^2) on the computed support
n = bad = 0
for key, supp in supports.items():
    K = 1 + max(abs(2 * i * i + 2 * j * j - k * k) for (i, j, k) in supp)
    ws = {t: K * (t[2] - t[0] - t[1]) + 2 * t[0] ** 2 + 2 * t[1] ** 2 - t[2] ** 2 for t in supp}
    surv = sorted(t for t, wv in ws.items() if wv == 0)
    a, b, c, d = key
    R = min(a, b, c, d)
    ok = min(ws.values()) == 0 and surv == [(i, i, 2 * i) for i in range(R)]
    n += 1; bad += not ok
report("(D) two-stage weights keep exactly grades (i,i,2i), i<min", n, bad)

print("ALL PASS" if fails == 0 else f"FAILURES: {fails}")
sys.exit(1 if fails else 0)
