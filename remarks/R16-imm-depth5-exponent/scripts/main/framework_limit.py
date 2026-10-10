#!/usr/bin/env python3
"""Family 135: framework limit of  'product-rank (Lemma 3.2) + expansion of any set of
factor occurrences (cost S^{|X|})'.

(a) continuum: tau*(sigma) = min(sigma/4, (sqrt(sigma^2+4)-sigma)/8), max at sigma=1/sqrt2,
    value 1/(4 sqrt2).
(b) exact finite n: H_e computed exactly (q = n^{-1/2} exactly as in the limit, lambda=k/n),
    adversary = uniform-degree factor lists (+ remainder), we take the best number of
    expanded copies.  Upper bound on what the framework can prove (no IMM-side loss K).
(c) dropping the constraint sum i_j = k (Lemma 3.2 line 111) is exponentially harmless:
    exact constrained sum [x^k] prod_j P_e(x) vs prod H_e for all-degree-2 lists.
(d) side conditions of Lemma 2.1 used with the new s at n = e^14.
"""
import math
import numpy as np
from mpmath import mp, mpf, sqrt as msqrt

def tau_star(sg):
    return min(sg / 4, (math.sqrt(sg * sg + 4) - sg) / 8)

def H(e, lam, q):
    # exact H_e = sum_{i=0}^e q^{|i - lam e|}; only i near lam e matter numerically
    c = lam * e
    i0 = int(math.floor(c))
    tot = 0.0
    for i in range(max(0, i0 - 60), min(e, i0 + 61) + 1):
        tot += q ** abs(i - c)
    return tot

def framework_value(n, sigma):
    s = int(math.floor(sigma * math.sqrt(n)))
    if (s - n) % 2:
        s -= 1
    k = (n - s) // 2
    lam = k / n
    q = n ** -0.5
    lH1 = -math.log(H(1, lam, q))
    best = None
    worst_e = None
    for e in range(2, int(4 * math.sqrt(n) / max(sigma, 0.2)) + 3):
        r, rem = divmod(n, e)
        lHe = -math.log(H(e, lam, q))
        lrem = -math.log(H(rem, lam, q)) if rem else 0.0
        # value = max_x [ (r-x) lHe + x e lH1 + lrem ] / (x+1)
        v = max(((r - x) * lHe + x * e * lH1 + lrem) / (x + 1) for x in range(0, r + 1))
        if best is None or v < best:
            best, worst_e = v, e
    return best / (math.sqrt(n) * math.log(n)), worst_e, s

def constrained_vs_H(n, sigma=1 / math.sqrt(2)):
    s = int(math.floor(sigma * math.sqrt(n)))
    if (s - n) % 2:
        s -= 1
    k = (n - s) // 2
    lam = k / n
    q = n ** -0.5
    r = n // 2
    # P(x) = sum_{i=0}^2 q^{|i-2lam|} x^i ; [x^k] P^r computed in log domain by tilting
    w = np.array([q ** abs(i - 2 * lam) for i in range(3)])
    # tilt x0 so that mean of tilted law = k/r
    from scipy.optimize import brentq
    def mean(t):
        p = w * np.exp(t * np.arange(3)); return (p * np.arange(3)).sum() / p.sum()
    t = brentq(lambda t: mean(t) - k / r, -50, 50)
    p = w * np.exp(t * np.arange(3)); Z = p.sum(); p = p / Z
    # exact coefficient of tilted multinomial at total k via log-binomial sum over (n0,n1,n2)
    from math import lgamma
    lp = np.log(p)
    best = -1e300; terms = []
    for n2 in range(0, r + 1):
        n0 = r - k + n2  # n1 = r - n0 - n2 ; n1 + 2 n2 = k
        n1 = k - 2 * n2
        if n0 < 0 or n1 < 0:
            continue
        terms.append(lgamma(r + 1) - lgamma(n0 + 1) - lgamma(n1 + 1) - lgamma(n2 + 1)
                     + n0 * lp[0] + n1 * lp[1] + n2 * lp[2])
    terms = np.array(terms); m = terms.max()
    logcoef_tilted = m + math.log(np.exp(terms - m).sum())
    log_constrained = logcoef_tilted + r * math.log(Z) - t * k
    log_H = r * math.log(w.sum())
    return log_constrained, log_H, s

if __name__ == "__main__":
    print("(a) continuum tau*(sigma):")
    for sg in [0.4, 0.5, 0.6, 1 / math.sqrt(2), 0.8, 1.0, 1.5]:
        print(f"   sigma={sg:.4f}  tau*={tau_star(sg):.6f}")
    grid = np.linspace(0.05, 3, 300001)
    vals = [tau_star(g) for g in grid]
    i = int(np.argmax(vals))
    print(f"   argmax sigma={grid[i]:.5f} (1/sqrt2={1/math.sqrt(2):.5f}) max={vals[i]:.6f} (1/(4sqrt2)={1/(4*math.sqrt(2)):.6f})")
    print("   paper's sigma=1:", tau_star(1.0), "  single-expansion 1/8 = 0.125")
    print("\n(b) exact finite-n framework value / (sqrt n log n)  [upper bound on provable c, no IMM loss]")
    for n in [10**4, 10**5, 10**6, 10**7]:
        row = []
        for sg in [0.5, 1 / math.sqrt(2), 1.0]:
            v, e, s = framework_value(n, sg)
            row.append(f"sigma={sg:.3f}: {v:.4f} (worst e={e}, s={s})")
        print(f"   n=1e{int(math.log10(n))}: " + " | ".join(row))
    print("   (finite-n values carry O(1/log n) corrections; limit tau*(sigma))")
    print("\n(c) constrained assignment sum vs H-product, all-degree-2 factors:")
    for n in [10**4, 10**5, 10**6]:
        lc, lh, s = constrained_vs_H(n)
        print(f"   n=1e{int(math.log10(n))}: log constrained={lc:.2f} log prodH={lh:.2f} "
              f"diff/sqrt n={(lh-lc)/math.sqrt(n):.4f}  saving prodH/(sqrt n log n)={-lh/(math.sqrt(n)*math.log(n)):.4f}")
    print("\n(d) Lemma 2.1 side conditions with s<=sqrt(n/2), n=e^14:")
    n = int(math.exp(14)); n += (n % 2)
    s = int(math.floor(math.sqrt(n / 2))); s -= (s - n) % 2
    k, m = (n - s) // 2, (n + s) // 2
    v, u = k * n * n, m * n * n
    a = math.ceil(v / (math.sqrt(n) - 1)); al = a / (a + v); rho = k / m
    b = math.ceil(u * al ** rho / (1 - al ** rho)); be = b / (b + u)
    q = al ** ((1 + rho) / 2)
    print(f"   n={n} s={s} k={k} m={m} lambda={k/n:.4f} s/k={s/k:.4f} a={a:.3e} b={b:.3e}")
    print(f"   2n<=min(a,b)/2: {2*n <= min(a,b)/2};  a>=k: {a>=k};  q={q:.3e} (<=0.1: {q<=0.1});"
          f" beta*sqrt n={be*math.sqrt(n):.4f};  t1=n/(2s)={n/(2*s):.1f}>2;  (1-1/(a+v))^-n={(1-1/(a+v))**-n:.6f}<=2")
