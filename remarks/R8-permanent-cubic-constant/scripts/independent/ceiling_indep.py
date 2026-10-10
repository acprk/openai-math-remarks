"""Ceiling for the manuscript's pipeline (selector m>=11k-2q, linear section of G_D, thm:smooth-border):
a smooth d-dimensional linear section exists only if d <= codim Gamma_D (projective dimension theorem), and
codim Gamma_D <= codim of any family inside Gamma_D.  Families (Prop. F, verified char-free):
coordinate block-upper-triangular X with m blocks, diagonal-block ranks rho_i, sum rho_i + m <= r-1.
We use one ranked block (size n1, rank r-1-m) and m-1 unranked blocks of near-equal sizes; m=1 is the
rank<=r-2 family (codim (2q+2)^2 = (k-r+2)^2).
Prints the best achievable c^{-1}(m) = 4m^3/max_{k,q}(r-1)(codimUB-1)  -- no proof inside the pipeline can beat it --
and the asymptotic value min_beta 4(10+beta)^3/(beta*min((1-beta)^2, 1/2-beta^2)), beta=r/k."""
import sys
def famdim(k, r):
    best = k*k - (k - r + 2)**2                       # m = 1
    for m in range(2, r):
        rho = r - 1 - m; t = m - 1
        for n1 in range(max(rho, 1), k - t + 1):
            a, b = divmod(k - n1, t)
            if a == 0: break
            sq = n1*n1 + b*(a+1)**2 + (t-b)*a*a
            d = (k*k - sq)//2 + rho*(2*n1 - rho)
            if d > best: best = d
    return best
def ceiling(m, kwin=None):
    bv, arg = 0, None
    kc = int(m/10.45)
    for k in range(kc - 25, kc + 26):
        for q in range(1, (k-2)//2 + 1):
            if 11*k - 2*q > m: continue
            r = k - 2*q
            if not (0.30*k <= r <= 0.45*k): continue     # optimum is near r/k ~ 0.39 (checked below)
            cod = k*k - famdim(k, r)
            v = (r-1)*(cod-1)
            if v > bv: bv, arg = v, (k, q, r, cod)
    return 4*m**3/bv, arg
def asym():
    best = (1e18, None)
    for i in range(1, 100000):
        b = i/100000
        g = min((1-b)**2, 0.5 - b*b)
        if g <= 0: continue
        v = 4*(10+b)**3/(b*g)
        if v < best[0]: best = (v, b)
    return best
if __name__ == '__main__':
    v, b = asym()
    print(f"asymptotic ceiling: c^-1 >= {v:.1f} e at beta=r/k={b:.4f} ; stated/ceiling = {5529600/v:.1f}x ; BOX(118800)/ceiling = {118800/v:.3f}")
    vr, br = min((4*(10+b)**3/(b*(1-b)**2), b) for b in [i/100000 for i in range(1, 100000)])
    print(f"(rank-family-only ceiling, i.e. old 'Lemma A': {vr:.1f} e at beta={br:.4f})")
    for m in [int(a) for a in sys.argv[1:]] or [1408, 5000]:
        c, a = ceiling(m)
        print(f"m={m}: ceiling c^-1 >= {c:.1f} e at (k,q,r,codimUB)={a} ; stated/ceiling={5529600/c:.1f}x")
