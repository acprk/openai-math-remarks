#!/usr/bin/env python3
"""Unit test: exact H_e vs the bounds used in the new proof, real parameters (eq:parameters) with s<=sqrt(n/2).
even e<t1: H_e <= q^x (1+(q^{1-2x}+q)/(1-q));  odd e<t1: H_e <= q^x (q^{1/2-2x}+q^{1/2})/(1-q), x=se/(2n);
H_1 <= 2 q^lambda.  Exact sum over all i in [0,e] with mpmath."""
from mpmath import mp, mpf, ceil, sqrt, floor
mp.dps = 40
for n in [10**6, 4*10**6, 1202604]:
    s = int(floor(sqrt(mpf(n)/2))); s -= (s-n)%2
    k, m = (n-s)//2, (n+s)//2; v = k*n*n
    a = int(ceil(mpf(v)/(sqrt(n)-1))); al = mpf(a)/(a+v); rho = mpf(k)/m
    q = al**((1+rho)/2); lam = mpf(k)/n; t1 = mpf(n)/(2*s)
    worst = 0
    for e in range(1, int(t1)+1):
        if e >= t1: break
        x = mpf(s)*e/(2*n)
        H = sum(q**abs(i-lam*e) for i in range(e+1))
        if e == 1:
            assert H <= 2*q**lam
        if e % 2 == 0:
            bd = q**x*(1+(q**(1-2*x)+q)/(1-q))
        else:
            bd = q**x*(q**(mpf(1)/2-2*x)+q**(mpf(1)/2))/(1-q)
        worst = max(worst, H/bd)
        assert H <= bd*(1+mpf(10)**-30), (n, e)  # bound tight up to tails q^{e/2}
    print(f"n={n}: s={s} t1={float(t1):.1f} q*sqrt n={float(q*sqrt(n)):.5f}  max H/bound={float(worst):.6f}  OK")
