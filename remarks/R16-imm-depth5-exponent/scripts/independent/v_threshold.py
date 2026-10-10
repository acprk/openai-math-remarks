#!/usr/bin/env python3
"""Independent threshold recomputation.  At integer n we use EXACT parameters (eq:parameters) with the
modified s, exact C_I(n) from v_imm.CI_parts (paper's own Sec.4-5 chain, all terms), the claimed
majorants Ce, Co (verified in v_params_H), and log S >= min(A,B):
   A = Lq s/2 - (C_I + Ce + Co + log 2),   B = t1 (Lq (n-2s)/(2n) - log 2).
Worst case over the floor of s: we also evaluate with s replaced by s-2 (same parity),
which is the smallest value s can take in a neighbourhood of n."""
import sys, math
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, sqrt, log, ceil, floor, exp
import v_params_H as VP
from v_imm import CI_parts
mp.dps = 80

def params_s(n, s):
    n = int(n); k, m = (n - s) // 2, (n + s) // 2
    v, u = k * n * n, m * n * n
    a = int(ceil(mpf(v) / (sqrt(mpf(n)) - 1))); al = mpf(a) / (a + v); rho = mpf(k) / m
    ar = al ** rho; b = int(ceil(u * ar / (1 - ar))); be = mpf(b) / (b + u)
    q = al ** ((1 + rho) / 2)
    return dict(n=n, s=s, k=k, m=m, v=v, u=u, a=a, b=b, al=al, be=be, rho=rho, q=q,
                lam=mpf(k) / n, t1=mpf(n) / (2 * s), Lq=-log(q))

def value(P):
    Ce, Co = VP.costs(P)
    CI = CI_parts(P)['CI']
    n, s, Lq = P['n'], P['s'], P['Lq']
    A = Lq * s / 2 - (CI + Ce + Co + log(2))
    B = P['t1'] * (Lq * (n - 2 * s) / (2 * n) - log(2))
    sc = sqrt(mpf(n)) * log(mpf(n))
    return A / sc, B / sc

def c_at(L):
    n = int(exp(mpf(L)))
    out = []
    for nn in (n, n + 1):
        P = VP.params(nn)
        for s in (P['s'], P['s'] - 2):
            a, b = value(params_s(nn, s))
            out.append(min(a, b))
    return min(out)

def thresh(c, lo=20, hi=2000):
    assert c_at(hi) >= c
    for _ in range(60):
        mid = (lo + hi) / 2
        if c_at(mid) >= c: hi = mid
        else: lo = mid
    return hi

if __name__ == "__main__":
    for L in [20, 30, 40, 50, 57.94, 60, 100, 200, 296.7, 300, 1000, 5000]:
        print(f"log n={L}: c(n) >= {float(c_at(L)):.6f}")
    # monotonicity scan
    prev, bad = None, 0
    for i in range(0, 3000):
        L = 20 + i * 0.5
        c = c_at(L)
        if prev is not None and c < prev - mpf(10)**-12: bad += 1
        prev = c
    print("non-monotone steps on grid logn in [20,1520] step .5:", bad)
    for c in [mpf(1)/400, mpf(3)/320, mpf(1)/8, mpf(1)/6]:
        print(f"c={float(c):.6f}: threshold log n = {float(thresh(c)):.3f}")
    print("1/(4 sqrt2) =", float(1/(4*sqrt(2))))
