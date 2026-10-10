# R3 independent check. Run from scripts/: /usr/bin/python3 independent/check197.py (needs sympy; ~1 s)
# Expected output: expected_output/check197.txt
# Independent adversarial check of the T-OPT parameters for openai/math #197.
# Graph: k-subsets of [N], adjacent iff |I cap J| < s. Fit polynomial
# Q(y) = prod_{c=s}^{k-1} (y-c)/(k-c), degree k-s; t = sum_{j<=k-s} C(N,j).
from fractions import Fraction as Fr
from math import comb, factorial
import random, sys

def stirling2(n, r):
    S = [[0]*(n+2) for _ in range(n+2)]
    S[0][0] = 1
    for a in range(1, n+1):
        for b in range(1, a+1):
            S[a][b] = b*S[a-1][b] + S[a-1][b-1]
    return S

def check(N, k, s, p=None, label=""):
    out = {}
    m = comb(N, k)
    deg = k - s
    t = sum(comb(N, j) for j in range(deg+1))
    out['m'] = m; out['t'] = t; out['deg'] = deg
    out['t/m'] = t/m
    out['96t<m'] = 96*t < m
    out['100t<m'] = 100*t < m
    # triangle-free: three pairwise-adjacent k-sets need union >= 3k-3(s-1) > N
    out['triangle_free_bound'] = 3*k - 3*(s-1) > N
    # also exact: max over pairwise intersections <= s-1 of minimal union: inclusion-exclusion
    # |A u B u C| >= 3k - (|AB|+|BC|+|CA|) ; this is attainable bound check only.
    # nonadjacent distinct: intersection in [s, k-1]
    # coefficients of Q in y (exact rationals)
    poly = [Fr(1)]
    for c in range(s, k):
        den = k - c
        new = [Fr(0)]*(len(poly)+1)
        for i, a in enumerate(poly):
            new[i+1] += a/den
            new[i] -= a*c/den
        poly = new
    S = stirling2(deg, deg)
    # multilinear coefficient of a fixed size-r monomial in Q(sum X_x): c_r = sum_n a_n r! S(n,r)
    cr = [sum(poly[n]*factorial(r)*S[n][r] for n in range(r, deg+1)) for r in range(deg+1)]
    # F(x) = sum_{T subset of I_i cap I_j, |T|<=deg} c_|T| = sum_r C(x,r) c_r
    ok_diag = sum(comb(k, r)*cr[r] for r in range(deg+1)) == 1
    ok_non = all(sum(comb(x, r)*cr[r] for r in range(deg+1)) == 0 for x in range(s, k))
    ok_eval = all(sum(comb(x, r)*cr[r] for r in range(deg+1)) ==
                  eval_Q(x, s, k) for x in range(0, k+1))
    out['fit_diag'] = ok_diag; out['fit_nonedge'] = ok_non; out['fit_eval_all_x'] = ok_eval
    dens = sorted({k-c for c in range(s, k)})
    out['denominators'] = (min(dens), max(dens))
    if p is not None:
        out['p>m'] = p > m
        out['p%4'] = p % 4
        out['p>max_den'] = p > max(dens)
        # cr mod p well-defined (denominators of cr coprime to p)
        out['cr_den_coprime_p'] = all(c.denominator % p != 0 for c in cr)
    print(label, out)
    return out, cr

def eval_Q(x, s, k):
    v = Fr(1)
    for c in range(s, k):
        v *= Fr(x - c, k - c)
    return v

def brute_set_check(N, k, s, cr, trials=200, seed=1):
    """Pick actual k-subsets, compute F_ij = sum_T alpha_i(T) beta_j(T) with
    beta_j(T)=c_|T| if T subset I_j else 0, by enumerating T in I_i cap I_j
    only by size-count (enumeration of 2^x impossible for x up to 41) --
    so instead verify the support claim on small random sub-instances."""
    rnd = random.Random(seed)
    deg = k - s
    bad = 0
    for _ in range(trials):
        I = set(rnd.sample(range(N), k)); J = set(rnd.sample(range(N), k))
        x = len(I & J)
        F = sum(comb(x, r)*cr[r] for r in range(deg+1))
        if x >= s and I != J and F != 0: bad += 1
        if I == J and F != 1: bad += 1
    return bad

if __name__ == '__main__':
    # original parameters (sanity)
    check(1200, 600, 200, label="orig(1200,600,200)") if '--orig' in sys.argv else None
    p = 6054252965430758779982953
    out, cr = check(86, 41, 13, p, label="new(86,41,13)")
    print("random-set check bad =", brute_set_check(86, 41, 13, cr))
    # Tiny honest brute force of the multilinearisation on an actual set system
    # (N,k,s)=(10,4,2): enumerate all T explicitly.
    from itertools import combinations
    N, k, s = 10, 4, 2
    _, crs = check(N, k, s, label="tiny(10,4,2)")
    # build beta_j(T) by explicitly expanding Q(sum X) with multilinear reduction
    import sympy as sp
    X = sp.symbols('x0:%d' % N)
    verts = [frozenset(c) for c in combinations(range(N), k)]
    def beta_poly(J):
        y = sum(X[a] for a in J)
        Q = sp.Integer(1)
        for c in range(s, k):
            Q = sp.expand(Q*(y - c)/sp.Integer(k - c))
        P = sp.Poly(Q, *X)
        d = {}
        for mon, coef in P.terms():
            T = frozenset(i for i, e in enumerate(mon) if e > 0)
            d[T] = d.get(T, 0) + coef
        return d
    J0 = verts[0]
    bJ = beta_poly(J0)
    assert max(len(T) for T in bJ) <= k - s
    errs = 0
    for I in verts:
        F = sum(c for T, c in bJ.items() if T <= I)
        x = len(I & J0)
        want = eval_Q(x, s, k)
        if F != want: errs += 1
    print("tiny explicit multilinear fit errors:", errs, "over", len(verts), "vertices")
