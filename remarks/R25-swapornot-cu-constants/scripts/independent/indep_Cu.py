#!/usr/bin/env python3
"""
Independent re-verification ("route B") of the values of
    C_u = sup_{m>=1} Z_u(m),   Z_u(m) = sum_{lambda |- m} f_lambda^{-u},
for u in {1, 1/2, 3/10, 1/5, 3/20}, and of the lower bound C_{1/10} >= Z_{1/10}(50).

Written from scratch; it does not import or reuse the primary scripts. Differences from the primary
route (paperH-style lemma_Cu: hook formula + MPFI, ballot-word peeling lemma, Catalan tail):
  * f_lambda from the branching rule f_lambda = sum_{corners} f_{lambda - corner} (exact integers);
  * exact region: every f^{-p/q} is enclosed by exact integer q-th roots at 2^-80 resolution
    (no floating point, no MPFR) and the enclosures are summed as exact integers;
  * peeling inequality: the exact first-row hook identity
        f_lambda = C(m,k) f_mu / Pi,  Pi = prod_{j<=mu_1} (1 + mu'_j/(L-j+1)),
    bounded by Chebyshev's sum inequality + AM-GM (Lemma B in the note) instead of ballot words;
  * a different analytic tail (Lemma T in the note: Robbins' binomial bound, three ranges).
Arithmetic in the recursive regions is IEEE double with an explicit outward inflation (see INF);
the tail inequality is evaluated in mpmath interval arithmetic.

Run: python3 indep_Cu.py            (single thread; ~1.5-3 min; ~0.8 GB RAM at peak)
Optional: ONLY=1/5 python3 indep_Cu.py   runs a single u (debugging).
Environment tested: CPython 3.10.12, numpy 2.2.6, mpmath 1.3.0.
"""
import math, sys, time, random
from fractions import Fraction
import numpy as np
import mpmath
from mpmath import iv

P = 80                      # exact region resolution 2^-P per term
INF = 1.0 + 2.0 ** -30      # outward inflation per float stage (>> accumulated IEEE/libm error)
TINY = 1e-280               # floor for positive entries (keeps upper bounds, avoids subnormals)

# (u, M0 exact, ME dp, M1 one-dimensional / tail start)
CASES_ALL = [(Fraction(1), 50, 120, 600), (Fraction(1, 2), 50, 120, 600), (Fraction(3, 10), 50, 150, 600),
         (Fraction(1, 5), 50, 250, 1000), (Fraction(3, 20), 50, 300, 1000)]
M0 = 50
import os
CASES = [c for c in CASES_ALL if os.environ.get('ONLY') in (None, str(c[0]))]

# ---------------------------------------------------------------- partitions + branching rule
def accel_asc(n):
    a = [0] * (n + 1); k = 1; y = n - 1
    while k != 0:
        x = a[k - 1] + 1; k -= 1
        while 2 * x <= y:
            a[k] = x; y -= x; k += 1
        l = k + 1
        while x <= y:
            a[k] = x; a[l] = y; yield a[:k + 2]; x += 1; y -= 1
        a[k] = x + y; y = x + y - 1
        yield a[:k + 1]

PNUM = [1]  # Euler pentagonal recurrence for p(n), used as a count check
def pcount(n):
    while len(PNUM) <= n:
        m = len(PNUM); s = 0; j = 1
        while True:
            g1 = j * (3 * j - 1) // 2
            if g1 > m: break
            sg = 1 if j % 2 else -1
            s += sg * PNUM[m - g1]
            g2 = j * (3 * j + 1) // 2
            if g2 <= m: s += sg * PNUM[m - g2]
            j += 1
        PNUM.append(s)
    return PNUM[n]

def iroot_floor(n, q):     # floor(n^(1/q)) for integers n >= 0, q >= 1 (integer Newton from above)
    if n < 2: return n
    bl = n.bit_length(); sh = max(bl - 64, 0)
    l2 = sh + math.log2(n >> sh)                       # float estimate of log2 n (start value only)
    r = l2 / q; ip = int(r)
    x = (int(2.0 ** (r - ip + 60) * (1 + 2.0 ** -30)) + 2) << ip >> 60 if ip >= 60 else \
        int(2.0 ** r * (1 + 2.0 ** -30)) + 2
    if x ** q <= n: x = 1 << -(-bl // q)               # start must be >= the root
    while True:
        y = ((q - 1) * x + n // x ** (q - 1)) // q
        if y >= x: return x
        x = y
def iroot_up(N, D, q):      # least integer Y with Y^q * D >= N   (N, D > 0)
    t = -(-N // D); y = iroot_floor(t, q)
    return y if y ** q >= t else y + 1
def iroot_dn(N, D, q):      # largest integer Y with Y^q * D <= N
    return iroot_floor(N // D, q)

t0 = time.time()
US = [c[0] for c in CASES] + [Fraction(1, 10)]
# Sup[u][m] : dict (lambda_1, ell) -> integer sum of upper enclosures (scaled 2^P);  Slo[u][m] : lower total
Sup = {u: {} for u in US}; Slo = {u: {} for u in US}
exact1 = {}                 # exact Z_1(m) as Fractions for m <= 12
prev = {(): 1}
fchk_ok = True
for m in range(1, M0 + 1):
    cur = {}
    cache = {u: {} for u in US}
    for u in US: Sup[u][m] = {}; Slo[u][m] = 0
    ex = Fraction(0)
    for asc in accel_asc(m):
        lam = tuple(reversed(asc)); r = len(lam); s = 0
        for i in range(r):
            if i == r - 1 or lam[i] > lam[i + 1]:
                mu = lam[:i] + ((lam[i] - 1,) if lam[i] > 1 else ()) + lam[i + 1:]
                s += prev[mu]
        cur[lam] = s
        key = (lam[0], r)
        for u in US:
            c = cache[u].get(s)
            if c is None:
                p, q = u.numerator, u.denominator
                N = 1 << (P * q); D = s ** p
                c = (iroot_up(N, D, q), iroot_dn(N, D, q) if (m <= 35 or u == Fraction(1, 10)) else 0)
                cache[u][s] = c
            d = Sup[u][m]; d[key] = d.get(key, 0) + c[0]; Slo[u][m] += c[1]
        if m <= 12: ex += Fraction(1, s)
    assert len(cur) == pcount(m)
    assert sum(v * v for v in cur.values()) == math.factorial(m)   # sum f^2 = m!
    if m <= 12: exact1[m] = ex
    prev = cur
print(f"exact region: all partitions of m <= {M0} enumerated ({sum(pcount(m) for m in range(1, M0+1))} shapes), "
      f"counts = p(m) and sum f^2 = m! checked for every m")
print(f"exact rational check: Z_1(6) = {exact1[6]}  (== 473/144: {exact1[6] == Fraction(473, 144)})")

# ---------------------------------------------------------------- unit tests of the peeling identity/lemma
def conj(lam):
    return [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
def hookdim(lam):
    if not lam: return 1
    cj = conj(lam); p = 1
    for i, r in enumerate(lam):
        for j in range(r): p *= (r - j - 1) + (cj[j] - i - 1) + 1
    return math.factorial(sum(lam)) // p
def lemmaB_test(Mt=24):
    n = viol = 0
    for m in range(1, Mt + 1):
        for asc in accel_asc(m):
            lam = tuple(reversed(asc)); L = lam[0]; mu = lam[1:]; k = m - L
            fl, fm = hookdim(lam), hookdim(mu)
            Pi = Fraction(1)
            for j, cj in enumerate(conj(mu), start=1): Pi *= 1 + Fraction(cj, L - j + 1)
            assert Fraction(math.comb(m, k) * fm) / Pi == fl            # exact identity
            if k:
                b = mu[0]
                Hd = sum(Fraction(1, L - j + 1) for j in range(1, b + 1))
                bndB = (1 + k * Hd / (b * b)) ** b                         # Chebyshev + AM-GM
                n += 1
                if Pi > bndB: viol += 1
                if Pi > Fraction(m, L) ** L: viol += 1                    # Pi <= (m/L)^L
                if math.log(Pi) > k / (L - min(k, L) + 1) + 1e-12: viol += 1  # Pi <= e^{k/(L-min(k,L)+1)}
    return n, viol
n, viol = lemmaB_test()
print(f"Lemma B unit test (all lambda |- m <= 24, exact rationals): identity holds; {n} shapes, {viol} violations")

# ---------------------------------------------------------------- helpers
def up_float(fr):            # float upper bound of a nonnegative Fraction
    x = float(fr)
    return math.nextafter(x, math.inf) if x > 0 else 0.0
def fmt_dn(fr, d=10): return f"{math.floor(fr * 10**d) / 10**d:.{d}f}"
def fmt_up(fr, d=10): return f"{math.ceil(fr * 10**d) / 10**d:.{d}f}"
SCALE = Fraction(1, 1 << P)

NMAX = 1001
mpmath.mp.dps = 40            # ln n! and H_n tabulated at 40 digits, then rounded to double (error <= 1 ulp)
lnfact = np.array([float(mpmath.loggamma(n + 1)) for n in range(NMAX + 1)])
Hn = np.array([float(mpmath.harmonic(n)) for n in range(NMAX + 1)])

def spot_check(u, samples):
    """recompute a sample of float weights in mpmath interval arithmetic; return #violations"""
    bad = 0; uiv = iv.mpf(u.numerator) / u.denominator
    for (kind, a1, a2, a3, wf) in samples:
        if kind == 'dp':
            k, a, b = a1, a2, a3; kp = k - a
            D = sum(iv.mpf(1) / (a - j + 1) for j in range(1, b + 1))
            lw = uiv * (b * iv.log(1 + kp * D / (b * b)) - iv.log(iv.mpf(math.comb(k, a))))
        else:
            m, L = a1, a2; k = m - L
            l1 = L * iv.log(iv.mpf(m) / L); l2 = iv.mpf(k) / (L - min(k, L) + 1)
            lp = iv.mpf([min(l1.a, l2.a), min(l1.b, l2.b)])
            lw = uiv * (lp - iv.log(iv.mpf(math.comb(m, k))))
        if iv.exp(lw).b > wf: bad += 1
    return bad

iv.dps = 30
results = []
for (u, M0u, ME, M1) in CASES:
    uf = float(u); tu = time.time()
    # ---------- exact region
    Zup = {}; Zlo = {}
    E = np.zeros((ME + 1, ME + 1, ME + 1)); E[0, 0, 0] = 1.0
    for m in range(1, M0u + 1):
        tot = 0
        for (a, c), S in Sup[u][m].items():
            E[m, a, c] = up_float(S * SCALE); tot += S
        Zup[m] = tot * SCALE; Zlo[m] = Slo[u][m] * SCALE
    mstar = max(Zup, key=lambda m: Zup[m])
    others = max(Zup[m] for m in Zup if m != mstar)
    assert others < Zlo[mstar]
    Zb = {m: up_float(Zup[m]) for m in Zup}           # float upper bounds of Z_u(m)
    # ---------- DP region  M0 < k <= ME
    rng = random.Random(20261009); samples = []
    for k in range(M0u + 1, ME + 1):
        E[k, k, 1] = E[k, 1, k] = 1.0
        for a in range(2, k):
            kp = k - a; bmax = min(a, kp); cmax = min(a, kp + 1)
            if cmax < 2: continue
            b = np.arange(1, bmax + 1)
            D = Hn[a] - Hn[a - b]
            lw = uf * (b * np.log1p(kp * D / (b * b)) - (lnfact[k] - lnfact[a] - lnfact[kp]))
            w = np.exp(lw) * INF
            if rng.random() < 0.02:
                j = rng.randrange(bmax); samples.append(('dp', k, a, j + 1, float(w[j])))
            vals = (w @ E[kp, 1:bmax + 1, 1:cmax]) * INF
            vals[(vals > 0) & (vals < TINY)] = TINY
            E[k, a, 2:cmax + 1] = vals
            E[k, 2:cmax + 1, a] = vals
        Zb[k] = math.nextafter(math.fsum(E[k].ravel().tolist()) * INF, math.inf)
    dpmax_m = max(range(M0u + 1, ME + 1), key=lambda k: Zb[k])
    # restricted cumulative sums R(k, L) = sum_{b<=min(L,k), c<=min(L-1,k)} E[k,b,c]
    CS = np.cumsum(np.cumsum(E, axis=1), axis=2) * (INF * INF)
    # ---------- one-dimensional region  ME < m <= M1
    for m in range(ME + 1, M1 + 1):
        s = math.isqrt(m - 1) + 1
        L = np.arange(s, m); k = m - L
        lp = np.minimum(L * np.log(m / L), k / (L - np.minimum(k, L) + 1))
        lw = uf * (lp - (lnfact[m] - lnfact[L] - lnfact[k]))
        w = np.exp(lw) * INF
        R = np.array([CS[kk, min(LL, kk), min(LL - 1, kk)] if kk <= ME else Zb[kk]
                      for LL, kk in zip(L.tolist(), k.tolist())])
        if rng.random() < 0.05:
            j = rng.randrange(len(L)); samples.append(('1d', m, int(L[j]), 0, float(w[j])))
        Zb[m] = math.nextafter(2.0 * math.fsum([1.0] + (w * R).tolist()) * INF * INF, math.inf)
    odmax_m = max(range(ME + 1, M1 + 1), key=lambda k: Zb[k])
    bad = spot_check(u, samples)
    # ---------- tail (Lemma T) at M = M1, interval arithmetic
    Cup = Zup[mstar]
    C = iv.mpf(math.ceil(Cup * 10**4)) / 10**4   # rounded up at 4th decimal
    viol_C = [m for m in Zb if Zb[m] > float(C.a)]
    if viol_C: print(f'  WARNING: bounds exceed C at m={viol_C[:10]}: {[Zb[m] for m in viol_C[:5]]}')
    uiv = iv.mpf(u.numerator) / u.denominator
    M = M1; s = math.isqrt(M - 1) + 1
    def lnbin(n, r): return iv.log(iv.mpf(math.comb(n, r)))
    best = None
    for tau in [Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(3, 10)]:
        tv = iv.mpf(tau.numerator) / tau.denominator
        kap = 1 - tv / 2 - tv * tv / 2
        ctau = min((1 - tv) * (-iv.log(1 - tv)), iv.log(3) / 3, key=lambda z: z.a)
        h = (8 * iv.mpf(s)) ** (uiv / 2) * iv.exp(-uiv * kap * s)
        rho = (1 + iv.mpf(1) / s) ** (uiv / 2) * iv.exp(-uiv * kap)
        assert rho.b < 1
        piece1 = h / (1 - rho)
        piece2 = M * (8 * iv.mpf(M)) ** (uiv / 2) * iv.exp(-uiv * ctau * M)
        mono2 = M >= float(((1 + uiv / 2) / (uiv * ctau)).b)
        for K in range(1, M0u + 1):
            if not (u * (K + 1) >= 1 + Fraction(1, M) and 3 * K <= M and mono2): continue
            A = sum(iv.exp(-uiv * lnbin(M, kk) + uiv * kk / (M - 2 * kk + 1)) * iv.mpf(up_float(Zup[kk]))
                    for kk in range(1, K + 1))
            gb = C * iv.exp(uiv) * (iv.mpf(M) / 3) * iv.exp(-uiv * lnbin(M, K + 1))
            T = 2 * (1 + A + gb + C * (piece1 + piece2))
            if best is None or T.b < best[0].b: best = (T, K, tau, A, gb, piece1, piece2)
    T, K, tau, A, gb, p1, p2 = best
    tail_ok = T.b <= C.a
    print(f"\n== u = {u}  (M0={M0u}, ME={ME}, M1={M1})")
    print(f"  exact region: sup over m<={M0u} at m*={mstar}: Z_u(m*) in [{fmt_dn(Zlo[mstar])}, {fmt_up(Zup[mstar])}]"
          f"; next largest m<={M0u} <= {fmt_up(others, 6)}")
    print(f"  DP region {M0u}<m<={ME}: max bound {Zb[dpmax_m]:.6f} at m={dpmax_m}")
    print(f"  1-D region {ME}<m<={M1}: max bound {Zb[odmax_m]:.6f} at m={odmax_m}; bound at M1: {Zb[M1]:.6f}")
    print(f"  weight spot-check vs mpmath.iv: {len(samples)} samples, {bad} violations")
    print(f"  tail (Lemma T, K={K}, tau={tau}): T(M1) <= {float(T.b):.6f} vs C = {float(C.a):.4f}: "
          f"{'OK' if tail_ok else 'FAIL'}  [A={float(A.b):.4f}, (b)={float(gb.b):.2e}, (c1)={float(p1.b):.2e}, (c2)={float(p2.b):.2e}]")
    sup_ok = tail_ok and float(T.b) < float(Zlo[mstar]) and bad == 0 and not viol_C and max(Zb[m] for m in Zb if m > M0u) < float(Zlo[mstar])
    print(f"  RESULT u={u}: C_u = Z_u({mstar}) in [{fmt_dn(Zlo[mstar])}, {fmt_up(Zup[mstar])}]"
          f"  ({'VERIFIED' if sup_ok else 'NOT VERIFIED'})   [{time.time()-tu:.0f}s]")
    results.append((u, mstar, sup_ok))
    del E, CS

u = Fraction(1, 10)
mm = max(Slo[u], key=lambda m: Slo[u][m])
print(f"\n== u = 1/10: lower bound only: max_m<={M0} Z(m) at m={mm}: Z in "
      f"[{fmt_dn(Slo[u][mm]*SCALE)}, {fmt_up(sum(Sup[u][mm].values())*SCALE)}]  => C_(1/10) >= {fmt_dn(Slo[u][mm]*SCALE, 4)}")
print("\nSUMMARY: " + ", ".join(f"u={u}: m*={m} {'ok' if ok else 'FAIL'}" for u, m, ok in results))
print(f"ALL {'VERIFIED' if all(ok for _, _, ok in results) else 'NOT ALL VERIFIED'}")
