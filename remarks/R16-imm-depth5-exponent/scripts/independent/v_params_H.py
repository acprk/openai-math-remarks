#!/usr/bin/env python3
"""Independent re-implementation, family 135.  Written from the tex, not
from main/.  Part 1: exact parameters (eq:parameters) with modified s<=sqrt(n/2); exact H_e;
check every inequality claimed in M2; check that the per-product inequality
   -sum_{j notin X} log H_{e_j} - E_X log H_1  >=  Lq s/2 + E_X c' - Ce - Co        (*)
for many adversarial degree lists (X = {e_j >= t1}), and that the resulting per-S value
is >= min(A,B) (with C_I+log2 set to the same K on both sides)."""
import random, sys
from mpmath import mp, mpf, sqrt, log, ceil, floor, exp
mp.dps = 80

def params(n, sfrac=mpf(1)/2):
    n = int(n)
    s = int(floor(sqrt(mpf(n) * sfrac)))
    if (s - n) % 2: s -= 1
    k, m = (n - s) // 2, (n + s) // 2
    v, u = k * n * n, m * n * n
    a = int(ceil(mpf(v) / (sqrt(mpf(n)) - 1)))
    al = mpf(a) / (a + v); rho = mpf(k) / m
    ar = al ** rho
    b = int(ceil(u * ar / (1 - ar)))
    be = mpf(b) / (b + u)
    q = al ** ((1 + rho) / 2)
    return dict(n=n, s=s, k=k, m=m, v=v, u=u, a=a, b=b, al=al, be=be, rho=rho, q=q,
                lam=mpf(k) / n, t1=mpf(n) / (2 * s), Lq=-log(q))

_Hcache = {}
def H(P, e):
    key = (P['n'], e)
    if key in _Hcache: return _Hcache[key]
    q, lam = P['q'], P['lam']
    c = lam * e
    i0 = int(floor(c))
    lo, hi = max(0, i0 - 40), min(e, i0 + 41)   # q<=1e-2 => q^40 negligible but add tail bound
    t = sum(q ** abs(i - c) for i in range(lo, hi + 1))
    # rigorous tail: remaining terms <= 2 q^{40}/(1-q)
    t += 2 * q ** 40 / (1 - q) if (lo > 0 or hi < e) else 0
    _Hcache[key] = t
    return t

def check_H(P):
    q, s, n, t1, lam = P['q'], P['s'], P['n'], P['t1'], P['lam']
    worst = 0
    assert H(P, 1) <= 2 * q ** lam
    e = 1
    while e < t1:
        x = mpf(s) * e / (2 * n)
        assert x < mpf(1) / 4
        # nearest-integer claims
        c = lam * e
        if e % 2 == 0:
            assert abs(abs(c - e // 2) - x) < mpf(10)**-60
            bd = q ** x * (1 + (q ** (1 - 2 * x) + q) / (1 - q))
        else:
            assert abs(abs(c - (e - 1) // 2) - (mpf(1) / 2 - x)) < mpf(10)**-60
            bd = q ** x * (q ** (mpf(1) / 2 - 2 * x) + q ** (mpf(1) / 2)) / (1 - q)
        r = H(P, e) / bd
        assert r <= 1 + mpf(10)**-60, (n, e, r)
        worst = max(worst, r)
        e += 1
    return worst

def costs(P):
    """Claimed majorants Ce (even corrections) and Co (odd corrections), M2."""
    q, s, n, t1 = P['q'], P['s'], P['n'], P['t1']
    f = lambda z: q ** (1 - mpf(s) * z / n) / z
    Ce = (n * max(f(2), f(t1)) + mpf(n) / 2 * q) / (1 - q)
    assert (q ** mpf('0.25') + q ** mpf('0.5')) / (1 - q) <= 1
    Co = 4 * s * (sqrt(q) + q / (1 - q))
    return Ce, Co

def lhs(P, degs):
    t1 = P['t1']
    from collections import Counter
    C = Counter(degs)
    nX = sum(c for e, c in C.items() if e >= t1)
    EX = sum(e * c for e, c in C.items() if e >= t1)
    val = -sum(c * log(H(P, e)) for e, c in C.items() if e < t1) - EX * log(H(P, 1))
    return val, nX, EX

def rhs(P, EX):
    Ce, Co = costs(P)
    n, s, Lq = P['n'], P['s'], P['Lq']
    cp = Lq * (n - 2 * s) / (2 * n) - log(2)
    return Lq * s / 2 + EX * cp - Ce - Co, cp

def lists(P, rng):
    n, t1 = P['n'], int(P['t1'])
    out = []
    for e in list(range(1, 3 * t1 + 3)) + [t1 - 1, t1, t1 + 1, 2 * t1, 2 * t1 + 1, 4 * t1]:
        if e < 1: continue
        r, rem = divmod(n, e)
        out.append([e] * r + ([rem] if rem else []))
    out.append([n])
    for _ in range(300):
        L, tot = [], 0
        pool = [1, 2, 2, 2, 4, rng.randint(1, t1), rng.randint(1, 3 * t1), t1, t1 + 1]
        while tot < n:
            e = min(rng.choice(pool), n - tot)
            L.append(e); tot += e
        out.append(L)
    return out

if __name__ == "__main__":
    rng = random.Random(1)
    for n in [10**4, 10**5 + 1, 10**6, 2 * 10**6 + 1, 10**7]:
        P = params(n)
        w = check_H(P)
        Ce, Co = costs(P)
        mins = None
        for degs in lists(P, rng):
            assert sum(degs) == n
            lv, nX, EX = lhs(P, degs)
            rv, cp = rhs(P, EX)
            assert lv >= rv - mpf(10) ** -20, (n, degs[:5], lv, rv)
            sl = lv - rv
            mins = sl if mins is None else min(mins, sl)
            # per-S value with common K'=0 (C_I+log2 shift is identical on both sides)
            A = P['Lq'] * P['s'] / 2 - Ce - Co
            B = P['t1'] * cp
            if A > 0 and B > 0:
                assert rv / (1 + nX) >= min(A, B) - mpf(10) ** -20
        # all-degree-2 list: how tight is Ce?
        r = n // 2
        exact_even = r * (log(H(P, 2)) + P['Lq'] * P['s'] * 2 / (2 * n))
        print(f"n={n}: s={P['s']} t1={float(P['t1']):.2f} q*sqrt(n)={float(P['q']*sqrt(n)):.5f} "
              f"maxH/bound={float(w):.6f}  Ce/sqrt n={float(Ce/sqrt(n)):.4f} "
              f"(exact all-2 even loss/sqrt n={float(exact_even/sqrt(n)):.4f})  Co={float(Co):.2f}  "
              f"min slack in (*)={float(mins):.3f}  [OK]")
