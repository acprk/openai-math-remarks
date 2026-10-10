# Theorem G check: for a reduced zero-dim complete intersection X of type a in P^n(F_q),
#   d_r(C_d(X)) >= mu_r(a, s-d-1),  s = sum(a_i - 1)  (CB7 shift: m = s-d-1),
#   mu_r(a,m) = min{ |D| : D down-set of box prod[0,a_i-1], #{beta in D : |beta| > m} >= r },
# and mu_r equals the GHW of the Cartesian code of the same type.
# GHW by brute force: d_r = min{|T| : dim{c in C : supp c subset T} >= r}.
# Usage: sage ghw_check.sage SEED
import sys, itertools, random
SEED = int(sys.argv[1]); set_random_seed(SEED); random.seed(SEED)

def downsets(a):
    box = list(itertools.product(*[range(x) for x in a]))
    box.sort(key=lambda b: (sum(b), b))
    res = []
    def rec(i, cur):
        if i == len(box):
            res.append(frozenset(cur)); return
        b = box[i]
        rec(i+1, cur)
        # can add b only if all predecessors in cur
        if all(tuple(b[:j]+(b[j]-1,)+b[j+1:]) in cur for j in range(len(b)) if b[j] > 0):
            cur.add(b); rec(i+1, cur); cur.remove(b)
    rec(0, set())
    return res

DS = {}
def mu(a, m, r):
    a = tuple(a)
    if a not in DS: DS[a] = downsets(a)
    best = None
    for D in DS[a]:
        if sum(1 for b in D if sum(b) > m) >= r:
            if best is None or len(D) < best: best = len(D)
    return best

def ghw_all(G):
    # G: k x N generator matrix over F_q; returns list d_1..d_k
    k, N = G.nrows(), G.ncols()
    best = [None]*(k+1)
    cols = list(range(N))
    for size in range(N+1):
        for T in itertools.combinations(cols, size):
            comp = [c for c in cols if c not in T]
            rk = G.matrix_from_columns(comp).rank() if comp else 0
            dimT = k - rk
            for r in range(1, dimT+1):
                if best[r] is None: best[r] = size
        if all(best[r] is not None for r in range(1, k+1)): break
    return best[1:]

def eval_code(K, pts, d, nv):
    R = PolynomialRing(K, nv, 'x')
    mons = [R.monomial(*e) for e in itertools.product(range(d+1), repeat=nv) if sum(e) == d]
    G = matrix(K, [[mm(*p) for p in pts] for mm in mons])
    return G.matrix_from_rows(G.pivot_rows())

def proj_points(K, nv):
    pts = []
    for v in itertools.product(K, repeat=nv):
        v = vector(K, v)
        if v.is_zero(): continue
        i = min(j for j in range(nv) if v[j] != 0)
        if v[i] == 1: pts.append(tuple(v))
    return pts

def ci_from_points(K, nv, a, seedpts):
    # forms of degrees a through seedpts; return (forms, full zero set) if zero set is a reduced CI of size prod(a)
    R = PolynomialRing(K, nv, 'x')
    fs = []
    for deg in a:
        mons = [R.monomial(*e) for e in itertools.product(range(deg+1), repeat=nv) if sum(e) == deg]
        M = matrix(K, [[mm(*p) for mm in mons] for p in seedpts])
        ker = M.right_kernel()
        if ker.dimension() == 0: return None
        f = sum(c*mm for c, mm in zip(ker.random_element(), mons))
        if f == 0: return None
        fs.append(f)
    I = R.ideal(fs)
    if I.dimension() != 1: return None
    Z = [p for p in ALLPTS[(K.order(), nv)] if all(f(*p) == 0 for f in fs)]
    if len(Z) != prod(a): return None
    return fs, Z

ALLPTS = {}
cases = [(11, 3, (2,2), 4), (11, 3, (2,3), None), (11, 3, (3,3), 8), (7, 3, (2,4), None), (7, 4, (2,2,2), 7)]
nviol = 0; ntests = 0
for (q, nv, a, nseed) in cases:
    K = GF(q); ALLPTS[(q, nv)] = proj_points(K, nv)
    s = sum(x-1 for x in a)
    found = 0
    for trial in range(400):
        if nseed is None:
            # points on a random smooth conic (P^2), then a curve of degree a[1] through prod(a) of them
            R = PolynomialRing(K, nv, 'x')
            C = R.random_element(degree=2, terms=Infinity).homogeneous_components().get(2, R(0))
            if C == 0 or len(C.factor()) != 1 or C.factor()[0][1] != 1: continue
            onC = [p for p in ALLPTS[(q, nv)] if C(*p) == 0]
            if len(onC) < prod(a): continue
            seed = random.sample(onC, prod(a))
            mons = [R.monomial(*e) for e in itertools.product(range(a[1]+1), repeat=nv) if sum(e) == a[1]]
            M = matrix(K, [[mm(*p) for mm in mons] for p in seed])
            g = sum(c*mm for c, mm in zip(M.right_kernel().random_element(), mons))
            if g == 0 or C.divides(g): continue
            Z = [p for p in ALLPTS[(q, nv)] if C(*p) == 0 and g(*p) == 0]
            if len(Z) != prod(a) or R.ideal([C, g]).dimension() != 1: continue
            res = ([C, g], Z)
        else:
            seed = random.sample(ALLPTS[(q, nv)], nseed)
            res = ci_from_points(K, nv, a, seed)
            if res is None: continue
        fs, Z = res
        found += 1
        for d in range(1, s+1):
            G = eval_code(K, Z, d, nv)
            ghw = ghw_all(G)
            for r, dr in enumerate(ghw, start=1):
                m_ = mu(a, s-d-1, r); ntests += 1
                if m_ is None or dr < m_:
                    nviol += 1; print('VIOLATION', q, a, d, r, dr, m_, fs)
        if found >= 6: break
    # Cartesian code of same type: affine grid A_1 x ... x A_n with |A_i| = a_i, embedded with x_0 = 1
    grid = [tuple([K(1)] + list(t)) for t in itertools.product(*[list(K)[:ai] for ai in a])]
    cart = []
    for d in range(1, s+1):
        ghw = ghw_all(eval_code(K, grid, d, nv))
        cart.append([(r, dr, mu(a, s-d-1, r)) for r, dr in enumerate(ghw, start=1)])
    eq = all(dr == m_ for row in cart for (_, dr, m_) in row)
    print('type', a, 'q', q, ': CI instances', found, '| Cartesian GHW == mu_r for all (d,r):', eq)
    if not eq: print(cart)
    sys.stdout.flush()
print('tests', ntests, 'violations', nviol)
