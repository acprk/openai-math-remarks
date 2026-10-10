#!/usr/bin/env python3
# Usage (from scripts/): /usr/bin/python3 independent/indep_scan.py   -- ~7 min single thread, stdlib only.
# Expected output: expected_output/indep_scan.txt
"""Independent re-derivation of the prime conditions in operator.tex (Hahn-Wilson height 2, paper 319).
Does NOT import hw319_*.  Exact integer/Fraction arithmetic.  Centers are computed exactly
(d_i = floor((2i+1)p/(2M))) instead of the paper's p/M +- 1 approximations, and every band
bound is computed from the actual center differences.  All depth quantities are in units of X
(X = p^r), lower-portion quantities in units of X/p^2.

Conditions (operator.tex line refs):
  G   p >= 11 (comparison.tex l.19/41 p>8; geometry l.14 p>5)
  P0  M odd, 3<=M<p, centers distinct, d_{s0}=delta (l.497-499)
  CE  1 <= s0-R0+1, s0+R0 <= M-2 (l.507; needed for I-1 to be center indices / no wrap l.613)
  CB  (2R0)^3 > p (l.503, used obstruction.tex l.1082/1123)
  BASE d=floor(p/M): p <= q d < p^2, d <= delta, 2Md >= p (l.547-551)
  LOW  Lam = max(diam(d(I u (I-1)))+1, d): |lower portion| <= Lam*X/p^2 (l.590-596; r=2 gives d)
  CS   Lam*q < p^2  (same-cluster diffs < X/q, l.597-598)
  BAND band_k = [gmin_k/p - Lam/p^2, gmax_k/p + Lam/p^2]; need band_1 low >= 1/q (inside stage),
       band_1 low > 0 (no wrap l.613), band_1 high < min_{k>=2} band_k low (l.601-603)
  HYP  paper: band_1 low >= 1/(2M) (lemma op:coefficient-control hypothesis b >= X/(2M))
  C1   paper: q/(2M) - kappa(M+1) > 1, kappa = p/(p-1)  (l.322-326, l.91)
  C1v  refined: q*bmin - kappa(M+1) >= 1 with bmin = min(band_1 low, d/p)   (requires restating the lemma)
  J    every cluster index i reachable in l.680-694 lies in I (exact); literal variant: s0-3..s0+3 in I
  S2   R0 >= 2 (l.705-708: d(I) has a digit < delta, d(I-1) has a digit > delta)
"""
from fractions import Fraction as F
import sys

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True

def check(p, M, R0, c1mode='paper', jmode='exact', Rmin=2, why=False):
    fail = []
    if p < 11: fail.append('G')
    if M % 2 == 0 or M < 3 or M >= p: return ['P0']
    delta, q, s0 = (p - 1) // 2, p + 1, (M - 1) // 2
    c = [((2 * i + 1) * p) // (2 * M) for i in range(M)]
    if len(set(c)) < M or c[s0] != delta: fail.append('P0')
    lo, hi = s0 - R0 + 1, s0 + R0
    if R0 < 1 or lo < 1 or hi > M - 2: return fail + ['CE']
    if R0 < Rmin: fail.append('S2/Rmin')
    if (2 * R0) ** 3 <= p: fail.append('CB')
    d = p // M
    if not (p <= q * d < p * p and d <= delta and 2 * M * d >= p): fail.append('BASE')
    idx = range(lo - 1, hi + 1)                      # I u (I-1)
    Lam = max(max(c[i] for i in idx) - min(c[i] for i in idx) + 1, d)
    if not Lam * q < p * p: fail.append('CS')
    def band(k):
        g = [(c[i] - c[(i - k) % M]) % p for i in range(M)]
        return F(min(g), p) - F(Lam, p * p), F(max(g), p) + F(Lam, p * p)
    b1lo, b1hi = band(1)
    if not (b1lo >= F(1, q) and b1lo > 0): fail.append('BAND-in')
    if M > 2 and not all(b1hi < band(k)[0] for k in range(2, M)): fail.append('BAND-sep')
    kappa = F(p, p - 1)
    if c1mode == 'paper':
        if b1lo < F(1, 2 * M): fail.append('HYP')
        if not F(q, 2 * M) - kappa * (M + 1) > 1: fail.append('C1')
    else:
        bmin = min(b1lo, F(d, p))
        if not q * bmin - kappa * (M + 1) >= 1: fail.append('C1v')
    if jmode == 'literal':
        if not (lo <= s0 - 3 and s0 + 3 <= hi): fail.append('J')
    else:
        U = F(1, 2) + b1hi
        reach = [i for i in range(M) if F(c[i] + 1, p) > F(1, 2) and F(c[i], p) <= U]
        if not all(lo <= i <= hi for i in reach): fail.append('J')
    if R0 < 2: fail.append('S2')
    return fail

def ok(*a, **k): return not check(*a, **k)

def pmin(c1mode, jmode, Rmin, R0rule='free', limit=100000):
    for p in range(11, limit):
        if not is_prime(p): continue
        for M in range(3, p, 2):
            # cheap necessary prefilter: CB needs (2*Rmax)^3 > p; C1/C1v need p > ~2M^2 (so M^2 < p)
            Rmax = (M - 3) // 2 if R0rule == 'free' else M // 8
            if (2 * Rmax) ** 3 <= p: continue
            if M * M >= p: break
            Rs = range(1, (M - 3) // 2 + 1) if R0rule == 'free' else [M // 8]
            for R0 in Rs:
                if ok(p, M, R0, c1mode, jmode, Rmin):
                    return p, M, R0
    return None

if __name__ == '__main__':
    print('== minimal primes (independent) ==')
    for tag, args in [('A paper-C1, literal J (R0>=4)', ('paper', 'literal', 4)),
                      ('B paper-C1, exact J, R0>=2', ('paper', 'exact', 2)),
                      ('C C1v, literal J (R0>=4)', ('C1v', 'literal', 4)),
                      ('D C1v, exact J, R0>=2', ('C1v', 'exact', 2))]:
        print(tag, '->', pmin(*args))
    print('R0=floor(M/8) kept, drop M>=81 & p>16M^2, paper C1, literal J ->',
          pmin('paper', 'literal', 4, R0rule='floor8', limit=200000))
    # failure reasons for all primes below 293 in mode A at best M
    print('== mode A: why each prime 11..293 fails (union over M,R0 of minimal failure sets) ==')
    for p in range(11, 294):
        if not is_prime(p): continue
        best = None
        for M in range(3, min(p, 41), 2):
            for R0 in range(1, max(2, (M - 3) // 2 + 1)):
                f = check(p, M, R0, 'paper', 'literal', 4)
                if best is None or len(f) < len(best[0]): best = (f, M, R0)
        print(p, best)
