#!/usr/bin/env python3
"""R7 adversarial check (independent of hw319_pmin.py / indep_scan.py).

Instead of re-encoding the inequalities of operator.tex, this script builds the
actual surviving lower-label sets of the induction table (operator.tex l.566-574)
for a given (p, M, R0) and stage X = p^r, enumerates EVERY (column, row) pair,
computes its exact depth b = (c - rho) mod X, and checks each combinatorial
claim the proof of Proposition op:operator makes about these pairs:

  ALLOW   depths of eligible pairs inside the stage satisfy X <= q b < pX (l.588-589)
  SAME    same-cluster pairs never have an allowed depth (l.597-598)
  BAND1   every k1=1 pair (column cluster i, row cluster i-1 cyclically) has an
          allowed depth, and b >= X/(2M) (hypothesis of Lemma op:coefficient-control)
  BANDS   max depth of band k1=1  <  min depth of all allowed pairs with k1>=2 (l.601-603)
  PAPER   band-1 depths lie in [X(1/M-2/p), X(1/M+2/p)] as printed in (op:first-band)
  NOWRAP  for column cluster i in I, every band-1 pair has c > rho as integers (l.613-615)
  EQUAL   both sides of each adjacent cluster pair have B0^(r-2) labels (l.628-629)
  JCUT    every surviving column label within one band-1 width above (X-1)/2+J
          has cluster index in I (l.680-695), for several J
  LIT3    the same cluster indices are within 3 of s0 (literal wording of l.693)
  S2      d(I) has a digit < delta, d(I-1) has a digit > delta (l.705-708)
  INIT    the X=p stage: d=floor(p/M) satisfies p <= q d < p^2, d >= p/(2M), d <= delta
  C1      (q/(2M) - kappa(M+1)) > 1 with kappa = p/(p-1)   (l.322-329, exact rational)
  CB      B0^3 > p  i.e. beta = log_p B0 > 1/3             (l.503)

Usage:  python3 hw319_enumerate.py            (default: p=293, M=11, R0=4, r=2..5,
                                               plus negative controls)
Runtime: ~1 min, needs numpy.  Exit code 0 iff the main case passes all checks.
"""
import sys
from fractions import Fraction as F
import numpy as np


def centers(p, M):
    return [((2 * i + 1) * p) // (2 * M) for i in range(M)]


def labels(p, M, r, first, mid_set, last_set):
    """All lower labels mod p^r with digit0=first, digits 1..r-2 in mid_set,
    digit r-1 in last_set (with the cluster index of the last digit)."""
    vals = np.array([first], dtype=np.int64)
    for pos in range(1, r - 1):
        vals = (vals[:, None] + np.array(mid_set, dtype=np.int64)[None, :] * p ** pos).ravel()
    out, clus = [], []
    for i, dl in enumerate(last_set):
        out.append(vals + dl * p ** (r - 1))
        clus.append(np.full(vals.shape, i, dtype=np.int64))
    return np.concatenate(out), np.concatenate(clus), len(vals)


def check_stage(p, M, R0, r, Js=(-10**6, -1000, 0, 1000, 10**6)):
    fails = []
    q, delta, s0 = p + 1, (p - 1) // 2, (M - 1) // 2
    X = p ** r
    c = centers(p, M)
    I = list(range(s0 - R0 + 1, s0 + R0 + 1))
    Im1 = [i - 1 for i in I]
    if min(I) < 1 or max(I) > M - 2:
        return ['CE(I meets endpoint)']
    d = p // M
    col, ccl, ncol = labels(p, M, r, delta, [c[i] for i in I], c)
    row, rcl, nrow = labels(p, M, r, delta - d, [c[i] for i in Im1], c)
    if ncol != nrow or ncol != (2 * R0) ** (r - 2):
        fails.append('EQUAL')
    b = (col[:, None] - row[None, :]) % X
    k1 = (ccl[:, None] - rcl[None, :]) % M
    allowed = (q * b >= X) & (q * b < p * X)
    if np.any(allowed & (k1 == 0)):
        fails.append('SAME')
    band1 = (k1 == 1)
    if not np.all(allowed[band1]):
        fails.append('BAND1-allowed')
    b1 = b[band1]
    b1min, b1max = int(b1.min()), int(b1.max())
    if F(b1min) < F(X, 2 * M):
        fails.append('BAND1-hyp b>=X/2M')
    later = allowed & (k1 >= 2)
    if later.any() and not b1max < int(b[later].min()):
        fails.append('BANDS')
    if not (F(b1min) >= X * (F(1, M) - F(2, p)) and F(b1max) <= X * (F(1, M) + F(2, p))):
        fails.append('PAPER-band')
    Iset = set(I)
    inI = np.isin(ccl, I)
    sel = band1 & inI[:, None]
    if not np.all((col[:, None] - row[None, :])[sel] > 0):
        fails.append('NOWRAP')
    # J-cut: surviving column labels in ((X-1)/2+J, (X-1)/2+J+b1max]
    reach = set()
    for J in Js:
        lo, hi = (X - 1) // 2 + J, (X - 1) // 2 + J + b1max
        m = (col > lo) & (col <= hi)
        reach |= set(int(t) for t in np.unique(ccl[m]))
    if not reach <= Iset:
        fails.append('JCUT reach=%s' % sorted(reach))
    if not all(abs(i - s0) <= 3 for i in reach):
        fails.append('LIT3')
    if not (any(c[i] < delta for i in I) and any(c[i] > delta for i in Im1)):
        fails.append('S2')
    info = dict(r=r, X=X, ncol=len(col), b1=(b1min, b1max),
                b1_over_X=(float(F(b1min, X)), float(F(b1max, X))),
                k2min=(int(b[later].min()) if later.any() else None), reach=sorted(reach))
    return fails, info


def check_global(p, M, R0):
    fails = []
    q, delta = p + 1, (p - 1) // 2
    if M % 2 == 0 or M >= p:
        fails.append('M odd, M<p')
    c = centers(p, M)
    if len(set(c)) != M or c[(M - 1) // 2] != delta:
        fails.append('centers')
    d = p // M
    if not (p <= q * d < p * p and F(d) >= F(p, 2 * M) and d <= delta):
        fails.append('INIT')
    kappa = F(p, p - 1)
    if not F(q, 2 * M) - kappa * (M + 1) > 1:
        fails.append('C1')
    if not (2 * R0) ** 3 > p:
        fails.append('CB')
    if R0 < 4:
        fails.append('R0>=4 (literal l.694)')
    if p <= 8:
        fails.append('p>8')
    # analytic forms (valid for every stage r) of the band claims of l.588-603:
    L = F(2 * R0 * p, M) + 3                      # lower-portion bound, replacing p/4+3 (l.592)
    if not L * q <= p * p:
        fails.append('CS (same-cluster outside stage)')
    if not L <= p:
        fails.append('DEV (|lower|/X <= 1/p)')
    if not 4 * M < p:
        fails.append('SEP (band 1 below band 2; b >= X/2M)')
    # limiting J-cut (l.689-694): clusters whose window [d_i/p,(d_i+1)/p) meets [1/2, 1/2+1/M+2/p]
    s0 = (M - 1) // 2
    I = set(range(s0 - R0 + 1, s0 + R0 + 1))
    reach = {i for i in range(M) if F(c[i] + 1, p) > F(1, 2) and F(c[i], p) <= F(1, 2) + F(1, M) + F(2, p)}
    if not (reach <= I and all(abs(i - s0) <= 3 for i in reach)):
        fails.append('JLIM')
    return fails


def run(p, M, R0, rs=(2, 3, 4, 5), verbose=True):
    allf = check_global(p, M, R0)
    if verbose:
        print('(p,M,R0)=(%d,%d,%d) global:' % (p, M, R0), allf or 'OK',
              ' C1 margin =', F(p + 1, 2 * M) - F(p, p - 1) * (M + 1) - 1)
    for r in rs:
        res = check_stage(p, M, R0, r)
        if isinstance(res, list):
            allf += res
            if verbose:
                print('   r=%d' % r, res)
            continue
        f, info = res
        allf += ['r=%d:%s' % (r, x) for x in f]
        if verbose:
            print('   r=%d' % r, f or 'OK', info)
    return allf


COVER = 2 * 10 ** 6

if __name__ == '__main__':
    print('== main case: p=293, M=11, R0=(M-3)/2=4 ==')
    main = run(293, 11, 4)
    print('== negative controls (each must FAIL something) ==')
    neg = {
        '(283,11,4) prime just below 293': (283, 11, 4),
        '(293,11,5) R0=(M-1)/2: I meets endpoint': (293, 11, 5),
        '(37,11,4) p<4M: bands collide': (37, 11, 4),
        '(293,9,3) R0=3 <4 (literal three-places)': (293, 9, 3),
    }
    negok = True
    for k, args in neg.items():
        f = run(*args, rs=(2, 3, 4), verbose=False)
        print('  ', k, '->', f[:6] if f else 'PASSES (control broken!)')
        negok &= bool(f)
    print('== exhaustive p_min under the modified argument (literal J: R0>=4, R0<=(M-3)/2) ==')

    def isprime(n):
        if n < 2:
            return False
        i = 2
        while i * i <= n:
            if n % i == 0:
                return False
            i += 1
        return True
    found = None
    for p in range(11, 400):
        if not isprime(p):
            continue
        for M in range(11, p, 2):
            R0 = (M - 3) // 2
            if R0 >= 4 and not check_global(p, M, R0) and not run(p, M, R0, rs=(2, 3), verbose=False):
                found = (p, M, R0)
                break
        if found:
            break
    print('   p_min =', found)
    # also: no prime < 293 passes with ANY admissible R0 in [4,(M-3)/2]
    bad = [(p, M, R0) for p in range(11, 293) if isprime(p) for M in range(11, p, 2)
           if F(p + 1, 2 * M) - F(p, p - 1) * (M + 1) > 1      # cheap necessary condition (C1) first
           for R0 in range(4, (M - 3) // 2 + 1) if not check_global(p, M, R0)]
    print('   primes <293 passing the global conditions with any R0 in [4,(M-3)/2]:', bad)
    print('== original parameters (op:parameters): smallest prime with M>=81 odd, p>16M^2, 2*floor(M/8)>p^(1/3) ==')
    best = None
    for M in range(81, 2001, 2):
        a = M // 8
        lo, hi = 16 * M * M, (2 * a) ** 3     # need lo < p < hi
        if hi <= lo:
            continue
        n = lo + 1
        while n < hi and not isprime(n):
            n += 1
        if n < hi and (best is None or n < best[0]):
            best = (n, M)
    print('   original p_min =', best)
    print('== coverage: every prime 293 <= p <= %d has an admissible M (M_p = least odd M>=11 with (M-3)^3 > p, R0=(M-3)/2) ==' % COVER)
    sieve = bytearray([1]) * (COVER + 1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(COVER ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    uncovered, Ms = [], set()
    for p in range(293, COVER + 1):
        if not sieve[p]:
            continue
        M = 11
        while (M - 3) ** 3 <= p:
            M += 2
        if check_global(p, M, (M - 3) // 2):
            uncovered.append(p)
        Ms.add(M)
    print('   uncovered primes:', uncovered[:10], '(count %d); M values used: %d..%d' % (len(uncovered), min(Ms), max(Ms)))
    print('   analytic tail: for M>=13, p >= (M-5)^3 >= 2M(M+3) implies C1; checked (M-5)^3 >= 2M(M+3) for 13<=M<=10^5:',
          all((M - 5) ** 3 >= 2 * M * (M + 3) for M in range(13, 10 ** 5, 2)), '(and it is increasing in M beyond)')
    print('== stage enumeration at sample (p,M_p) from different windows ==')
    samp = [(509, 11, 4), (521, 13, 5), (1009, 15, 6), (10007, 25, 11)]
    sampok = True
    for (p, M, R0) in samp:
        f = run(p, M, R0, rs=(2, 3) if M < 20 else (2,), verbose=False)
        print('  ', (p, M, R0), f or 'OK')
        sampok &= not f
    ok = (not main) and negok and not uncovered and sampok and found == (293, 11, 4) and not bad and best == (17073437, 1033)
    print('VERDICT:', 'ALL CHECKS PASS' if ok else 'FAILURE')
    sys.exit(0 if ok else 1)
