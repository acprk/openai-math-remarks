#!/usr/bin/env python3
"""Independent re-implementation (pure Python + numpy, no Sage) of the numerical check of
Theorem G:  d_r(C_d(X)) >= mu_r(a, s-d-1) = d_r(C_d(X_Car))
for reduced complete intersections X of type a in P^n(F_p), p prime.

Written from the statement only; it shares no code with the primary Sage scripts.
Differences from the primary computation:
  * GHW computed from the column matroid: d_r = N - max{|S| : rank(G_S) <= k - r},
    ranks of all 2^N column subsets by depth-first search with an incremental echelon basis;
  * mu_r computed by filtering all 2^|box| subsets of the box for down-closedness (bitmasks);
  * CI instances certified by a Hilbert-function test (no Groebner bases):
      P^2:  f1, f2 have no common factor  <=>  dim_k (S/(f1,f2))_T = a1*a2 at T = a1*a2+a1+a2;
      P^3 split case f3 = l_1...l_s: X finite <=> on every plane l_j the restrictions of f1, f2
      have no common factor (same P^2 test after a linear change of coordinates);
    then |X(F_p)| = deg X with distinct points  =>  X reduced with all points rational;
  * new instance types not in the primary run: (2,5), (3,5), (4,4) in P^2 and the split
    family G3 in P^3: (2,2,2), (2,2,3), (2,3,3) with f3 a product of planes.
Also a "teeth" control: the off-by-one shift m = s-d (instead of s-d-1) must produce violations.

Usage:  python3 ghw_ci_indep.py [seed]      (single thread; about 5-15 min)
"""
import sys, itertools, random
import numpy as np

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 2026
rng = random.Random(SEED)

# ---------------------------------------------------------------- linear algebra mod p
def rank_mod_p(M, p):
    M = np.array(M, dtype=np.int64) % p
    if M.size == 0:
        return 0
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c] != 0:
                piv = i; break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r

def monomials(nv, deg):
    return [e for e in itertools.product(range(deg + 1), repeat=nv) if sum(e) == deg]

def ev_mono(e, pt, p):
    v = 1
    for ei, xi in zip(e, pt):
        v = v * pow(xi, ei, p) % p
    return v

# polynomials: dict exponent-tuple -> coeff mod p
def poly_eval(f, pt, p):
    return sum(c * ev_mono(e, pt, p) for e, c in f.items()) % p

def poly_mul(f, g, p):
    h = {}
    for e1, c1 in f.items():
        for e2, c2 in g.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            h[e] = (h.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in h.items() if c}

def random_form(nv, deg, p):
    f = {e: rng.randrange(p) for e in monomials(nv, deg)}
    return {e: c for e, c in f.items() if c}

def proj_points(nv, p):
    pts = []
    for v in itertools.product(range(p), repeat=nv):
        if any(v):
            i = next(j for j in range(nv) if v[j])
            if v[i] == 1:
                pts.append(v)
    return pts

def hilb_quot(forms, nv, T, p):
    """dim_k (S/(forms))_T, S = k[x_1..x_nv]."""
    mons = monomials(nv, T); idx = {e: i for i, e in enumerate(mons)}
    rows = []
    for f in forms:
        dg = sum(next(iter(f)))
        if dg > T:
            continue
        for m in monomials(nv, T - dg):
            row = [0] * len(mons)
            for e, c in f.items():
                row[idx[tuple(a + b for a, b in zip(e, m))]] = c
            rows.append(row)
    return len(mons) - (rank_mod_p(rows, p) if rows else 0)

# ---------------------------------------------------------------- GHW via column matroid
def ghw_matroid(G, p):
    """G: list of k rows (k = rank), N columns.  d_r = N - max{|S| : rank(G_S) <= k - r}."""
    k, N = len(G), len(G[0])
    cols = [[int(G[i][j]) % p for i in range(k)] for j in range(N)]
    best = [-1] * (k + 1)          # best[t] = max |S| with rank(S) == t

    def reduce(v, basis):
        v = v[:]
        for piv, row in basis:
            c = v[piv]
            if c:
                v = [(a - c * b) % p for a, b in zip(v, row)]
        return v

    def dfs(j, basis, size):
        rk = len(basis)
        if size > best[rk]:
            best[rk] = size
        if j == N:
            return
        # include column j
        v = reduce(cols[j], basis)
        piv = next((i for i, x in enumerate(v) if x), None)
        if piv is None:
            dfs(j + 1, basis, size + 1)
        else:
            inv = pow(v[piv], p - 2, p)
            v = [(x * inv) % p for x in v]
            nb = [(pv, [(a - a_piv * b) % p for a, b in zip(row, v)]) if (a_piv := row[piv]) else (pv, row)
                  for pv, row in basis]
            nb.append((piv, v))
            dfs(j + 1, nb, size + 1)
        # exclude column j
        dfs(j + 1, basis, size)

    dfs(0, [], 0)
    # max |S| with rank <= t
    mx = [max(best[: t + 1]) for t in range(k + 1)]
    return [N - mx[k - r] for r in range(1, k + 1)]

def eval_code(pts, d, nv, p):
    mons = monomials(nv, d)
    M = [[ev_mono(e, pt, p) for pt in pts] for e in mons]
    # row-reduce to a basis
    A = np.array(M, dtype=np.int64) % p
    basis = []
    for row in A:
        cand = basis + [list(row)]
        if rank_mod_p(cand, p) > len(basis):
            basis.append([int(x) for x in row])
    return basis

# ---------------------------------------------------------------- mu_r by brute force over subsets
DS_CACHE = {}
def downsets_bruteforce(a):
    box = list(itertools.product(*[range(x) for x in a]))
    idx = {b: i for i, b in enumerate(box)}
    preds = []
    for b in box:
        m = 0
        for j in range(len(b)):
            if b[j] > 0:
                m |= 1 << idx[b[:j] + (b[j] - 1,) + b[j + 1:]]
        preds.append(m)
    res = []
    for mask in range(1 << len(box)):
        ok = True
        mm = mask
        while mm:
            low = mm & -mm; i = low.bit_length() - 1
            if preds[i] & ~mask:
                ok = False; break
            mm ^= low
        if ok:
            res.append([box[i] for i in range(len(box)) if mask >> i & 1])
    return res

def mu(a, m, r):
    a = tuple(a)
    if a not in DS_CACHE:
        DS_CACHE[a] = downsets_bruteforce(a)
    vals = [len(D) for D in DS_CACHE[a] if sum(1 for b in D if sum(b) > m) >= r]
    return min(vals) if vals else None

# ---------------------------------------------------------------- instance generators
def ci_P2(a1, a2, p, PTS):
    """random reduced CI of type (a1,a2) in P^2(F_p) with all points rational, or None."""
    f1 = random_form(3, a1, p)
    on = [pt for pt in PTS if poly_eval(f1, pt, p) == 0]
    if len(on) < a1 * a2:
        return None
    mons = monomials(3, a2)
    for _ in range(20):
        c = rng.randint(max(1, a1 * a2 - 4), a1 * a2)
        seed = rng.sample(on, c)
        M = [[ev_mono(e, pt, p) for e in mons] for pt in seed]
        ker = nullspace_mod_p(M, len(mons), p)
        if not ker:
            continue
        lam = [rng.randrange(p) for _ in ker]
        coeff = [sum(l * v[i] for l, v in zip(lam, ker)) % p for i in range(len(mons))]
        f2 = {e: c_ for e, c_ in zip(mons, coeff) if c_}
        if not f2:
            continue
        Z = [pt for pt in on if poly_eval(f2, pt, p) == 0]
        if len(Z) != a1 * a2:
            continue
        T = a1 * a2 + a1 + a2
        if hilb_quot([f1, f2], 3, T, p) != a1 * a2:
            continue                      # common factor
        return [f1, f2], Z
    return None

def nullspace_mod_p(M, ncols, p):
    A = [list(map(lambda x: x % p, row)) for row in M]
    rows = len(A); piv_cols = []; r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, rows) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv_cols.append(c); r += 1
        if r == rows:
            break
    free = [c for c in range(ncols) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * ncols; v[fc] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = (-A[i][fc]) % p
        basis.append(v)
    return basis

def plane_restriction_ok(f1, f2, plane, p, a1, a2):
    """restrict f1,f2 to the plane {plane . x = 0} in P^3 via a parametrisation P^2 -> plane;
    return True iff the restrictions have no common factor (P^2 Hilbert test)."""
    # basis of the plane (kernel of 1x4 matrix)
    B = nullspace_mod_p([plane], 4, p)            # 3 vectors in F_p^4
    def restrict(f):
        g = {}
        for e, c in f.items():
            # prod_i (sum_j B[j][i] y_j)^{e_i}
            term = {(0, 0, 0): c}
            for i, ei in enumerate(e):
                lin = {}
                for j in range(3):
                    if B[j][i] % p:
                        ex = [0, 0, 0]; ex[j] = 1
                        lin[tuple(ex)] = B[j][i] % p
                for _ in range(ei):
                    term = poly_mul(term, lin, p) if lin else {}
            for ee, cc in term.items():
                g[ee] = (g.get(ee, 0) + cc) % p
        return {e: c for e, c in g.items() if c}
    g1, g2 = restrict(f1), restrict(f2)
    if not g1 or not g2:
        return False
    T = a1 * a2 + a1 + a2
    return hilb_quot([g1, g2], 3, T, p) == a1 * a2

def ci_split_P3(a1, a2, s, p, PTS):
    """CI of type (a1,a2,s) in P^3: f1,f2 random through some points, f3 = product of s planes."""
    f1 = random_form(4, a1, p)
    mons2 = monomials(4, a2)
    on1 = [pt for pt in PTS if poly_eval(f1, pt, p) == 0]
    seed = rng.sample(on1, min(len(on1), len(mons2) - 2))
    M = [[ev_mono(e, pt, p) for e in mons2] for pt in seed]
    ker = nullspace_mod_p(M, len(mons2), p)
    if not ker:
        return None
    lam = [rng.randrange(p) for _ in ker]
    coeff = [sum(l * v[i] for l, v in zip(lam, ker)) % p for i in range(len(mons2))]
    f2 = {e: c for e, c in zip(mons2, coeff) if c}
    if not f2:
        return None
    curve = [pt for pt in on1 if poly_eval(f2, pt, p) == 0]
    planes, Z = [], []
    tries = 0
    while len(planes) < s and tries < 400:
        tries += 1
        if len(curve) < 3:
            return None
        three = rng.sample(curve, 3)
        pl = nullspace_mod_p(three, 4, p)
        if len(pl) != 1:
            continue
        pl = pl[0]
        on = [pt for pt in curve if sum(x * y for x, y in zip(pl, pt)) % p == 0]
        if len(on) != a1 * a2 or any(pt in Z for pt in on):
            continue
        if not plane_restriction_ok(f1, f2, pl, p, a1, a2):
            continue
        planes.append(pl); Z += on
    if len(planes) < s:
        return None
    f3 = {(0, 0, 0, 0): 1}
    for pl in planes:
        lin = {tuple(1 if i == j else 0 for i in range(4)): c % p for j, c in enumerate(pl) if c % p}
        f3 = poly_mul(f3, lin, p)
    assert all(poly_eval(f3, pt, p) == 0 and poly_eval(f1, pt, p) == 0 and poly_eval(f2, pt, p) == 0 for pt in Z)
    assert len(set(Z)) == a1 * a2 * s
    return [f1, f2, f3], Z

def cartesian(a, p):
    return [tuple([1] + list(t)) for t in itertools.product(*[list(range(ai)) for ai in a])]

def cartesian_ok(a, p, nv):
    cart = cartesian(a, p)
    s = sum(x - 1 for x in a)
    for d in range(1, s + 1):
        ghw = ghw_matroid(eval_code(cart, d, nv, p), p)
        if any(dr != mu(a, s - d - 1, r) for r, dr in enumerate(ghw, start=1)):
            return False
    return True

# ---------------------------------------------------------------- main loop
def check_instance(a, Z, nv, p, stats, label):
    s = sum(x - 1 for x in a)
    for d in range(1, s + 1):
        G = eval_code(Z, d, nv, p)
        ghw = ghw_matroid(G, p)
        for r, dr in enumerate(ghw, start=1):
            m_ = mu(a, s - d - 1, r)
            stats['tests'] += 1
            if m_ is None or dr < m_:
                stats['viol'] += 1
                print('VIOLATION', label, a, 'd', d, 'r', r, 'd_r', dr, 'mu', m_)
            elif dr == m_:
                stats['eq'] += 1
            else:
                stats['strict'] += 1
            # teeth control: wrong shift m = s - d
            mw = mu(a, s - d, r)
            if mw is not None and dr < mw:
                stats['wrong_shift_viol'] += 1

def main():
    stats = dict(tests=0, viol=0, eq=0, strict=0, wrong_shift_viol=0)
    plan_P2 = [((2, 2), 11, 6), ((2, 3), 11, 6), ((3, 3), 11, 6), ((2, 4), 13, 6), ((3, 4), 13, 6),
               ((2, 5), 13, 4), ((3, 5), 17, 2), ((4, 4), 17, 2)]
    plan_P3 = [((2, 2, 2), 7, 6), ((2, 2, 3), 11, 4), ((2, 3, 3), 13, 2)]
    for a, p, want in plan_P2:
        PTS = proj_points(3, p)
        found = 0; tries = 0
        while found < want and tries < 5000:
            tries += 1
            res = ci_P2(a[0], a[1], p, PTS)
            if res is None:
                continue
            found += 1
            check_instance(a, res[1], 3, p, stats, 'P2')
        cart_ok = cartesian_ok(a, p, 3)
        print(f'G1  P^2 type {a} p={p}: instances {found}; Cartesian GHW == mu_r for all (d,r): {cart_ok}', flush=True)
    for a, p, want in plan_P3:
        PTS = proj_points(4, p)
        found = 0; tries = 0
        while found < want and tries < 3000:
            tries += 1
            res = ci_split_P3(a[0], a[1], a[2], p, PTS)
            if res is None:
                continue
            found += 1
            check_instance(a, res[1], 4, p, stats, 'P3split')
        cart_ok = cartesian_ok(a, p, 4)
        print(f'G3  P^3 type {a} (f3 = product of planes) p={p}: instances {found}; Cartesian GHW == mu_r: {cart_ok}', flush=True)
    print('total (d,r) tests', stats['tests'], '| violations', stats['viol'],
          '| equality', stats['eq'], '| strict', stats['strict'])
    print('teeth control (wrong shift m=s-d): violations', stats['wrong_shift_viol'], '(must be > 0)')

if __name__ == '__main__':
    main()
