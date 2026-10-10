#!/usr/bin/env python3
"""Independent re-verification of a product relation for OpenAI preprint
'Ambiently homeomorphic isolated hypersurfaces of multiplicities two and three'
(Lemma seed:product, Proposition seed:seeds, eq. topo:dimension).

Written from the paper text only; does NOT import or reuse explicit_relation.py.
Checks:
  (1) support: d = 3^k, 0<b<d/2, 3∤b, e_b nonzero integers, E = sum e_b > 0;
  (2) prod A_b^{e_b} = prod R_b^{e_b} = 1, by (a) direct big-integer products
      (no factorisation at all) when feasible, and (b) prime valuations via
      a smallest-prime-factor sieve;
  (3) the seed construction of Prop. seed:seeds built as explicit weight
      multisets (chains via the decreasing-remainder algorithm, tails on the
      opposite side, 4E squares, 2 cubes), giving s0 and N = 4 s0 + 16;
  (4) sampled Fourier test of 'spectra equal mod 2' (eq. seed:spectrum-product)
      on the full weight multisets (optional, --spectral).
Standard library only (plus cmath).  Usage: indep_check_059.py relation.json [--spectral]
"""
import json, sys, math, cmath, random
from collections import Counter
from fractions import Fraction

def load(path):
    J = json.load(open(path))
    d = int(J["d"])
    rel = [(int(r["b"]), int(r["e"])) for r in J["relation"]]
    return d, rel

def check_support(d, rel):
    k = 0; x = d
    while x % 3 == 0: x //= 3; k += 1
    assert x == 1, "d not a power of 3"
    bs = [b for b, _ in rel]
    assert len(set(bs)) == len(bs)
    for b, e in rel:
        assert 0 < b and 2 * b < d and b % 3 != 0 and e != 0
    E = sum(e for _, e in rel)
    assert E > 0, "E must be positive"
    return k, E

def products_direct(d, rel):
    # A_b = (d-b)/b ; R_b = (2d+b)(2d-b)/((3d-b)(d+b))
    numA = []; denA = []; numR = []; denR = []
    for b, e in rel:
        a_num, a_den = d - b, b
        r_num, r_den = (2*d + b) * (2*d - b), (3*d - b) * (d + b)
        if e < 0:
            a_num, a_den, r_num, r_den = a_den, a_num, r_den, r_num
        f = abs(e)
        numA.append(pow(a_num, f)); denA.append(pow(a_den, f))
        numR.append(pow(r_num, f)); denR.append(pow(r_den, f))
    PA = (math.prod(numA), math.prod(denA)); PR = (math.prod(numR), math.prod(denR))
    return PA[0] == PA[1], PR[0] == PR[1], PA[0].bit_length(), PR[0].bit_length()

def spf_sieve(n):
    s = bytearray(0)  # use list of ints
    spf = list(range(n + 1))
    i = 2
    while i * i <= n:
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j: spf[j] = i
        i += 1
    return spf

def products_valuation(d, rel):
    spf = spf_sieve(3 * d)
    def fac(n):
        out = Counter()
        while n > 1:
            p = spf[n]; out[p] += 1; n //= p
        return out
    vA = Counter(); vR = Counter()
    for b, e in rel:
        for n, s in ((d - b, 1), (b, -1)):
            for p, a in fac(n).items(): vA[p] += s * a * e
        for n, s in ((2*d + b, 1), (2*d - b, 1), (3*d - b, -1), (d + b, -1)):
            for p, a in fac(n).items(): vR[p] += s * a * e
    return all(v == 0 for v in vA.values()), all(v == 0 for v in vR.values())

def chain_weights(D, c):
    """Numerators (over D) of the weights of the chain starting at c/D."""
    out = []
    while c:
        a = D // c; r = D - a * c
        assert a >= 2 and 0 <= r < c
        # ordinary degree of the monomial: a+1 if r else a; must be >=3
        assert (a + 1 if r else a) >= 3
        out.append(c); c = r
    return out

def seed_multisets(d, rel, want_multisets=False):
    D = 3 * d
    side = [Counter(), Counter()] if want_multisets else None
    nvars = [0, 0]
    for b, e in rel:
        for c, cnt in ((3*b, 2*e), (b, -2*e), (d - b, 2*e), (d + b, 2*e)):
            assert 0 < 2 * c < D
            ch = chain_weights(D, c)
            j = 1 if cnt > 0 else 0          # positive count -> side 2
            m = abs(cnt)
            nvars[j] += m * len(ch)          # full chains on side j
            nvars[1 - j] += m * (len(ch) - 1)  # tails on the opposite side
            if want_multisets:
                for w in ch: side[j][w] += m
                for w in ch[1:]: side[1 - j][w] += m
    E = sum(e for _, e in rel)
    assert nvars[1] - nvars[0] == 4 * E, (nvars, E)
    nsq = 4 * E
    s_before = nvars[0] + nsq
    assert s_before == nvars[1]
    assert s_before % 2 == 0
    s0 = s_before + 2                         # two common pure cubes
    N = 4 * (s0 + 2) + 6 + 2
    assert N == 4 * s0 + 16 and N % 8 == 0
    return dict(D=D, E=E, nsq=nsq, s0=s0, N=N, side=side)

def spectral_test(d, rel, info, nsamples=400, seed=1):
    """Characters alpha -> exp(pi i k alpha), k mod 4D.  Factor of weight w
    (as a Fraction) in P_p(t)=prod (t^w - t)/(1 - t^w)."""
    D = info["D"]; side = info["side"]; E = info["E"]
    # full multisets: chains/tails (numerators over D), squares (w=1/2), cubes (1/3 x2)
    def factor_log(k, w):
        # returns complex log of factor value; None if zero (only odd k, 3|k for cubes)
        if k % 2 == 0:
            m = k // 2
            if (m * w).denominator == 1:
                return cmath.log(complex((1 - w) / w))
            return cmath.log(-1)
        # exact reduction of k*w modulo 2 before the floating exponential
        r = (k * w.numerator) % (2 * w.denominator)
        z = cmath.exp(1j * math.pi * r / w.denominator)
        if (k * w).denominator == 1 and (k * w).numerator % 2 == 0:
            raise ValueError("pole")
        val = (z + 1) / (1 - z)   # t = (-1)^k = -1
        return cmath.log(val)
    net = Counter(side[0]); net.subtract(side[1])
    net = {c: m for c, m in net.items() if m}
    rng = random.Random(seed)
    ks = list(range(0, 60)) + [rng.randrange(4 * D) for _ in range(nsamples)]
    tested = 0; maxerr = 0.0
    for k in ks:
        if k % 2 == 1 and k % 3 == 0:
            continue  # both sides contain a cube: character value 0 (paper Case 3)
        # common weights have finite nonzero factors here (paper Cases 1-2), so
        # compare the net signed multiset side1 - side2 (incl. 4E squares on side 1;
        # the two common cubes cancel).  The net is computed from the FULL multisets.
        diff = 0j
        for c in net:
            diff += net[c] * factor_log(k, Fraction(c, D))
        diff += info["nsq"] * factor_log(k, Fraction(1, 2))
        # equal values <=> real parts equal and imaginary parts equal mod 2pi
        re = abs(diff.real); im = abs((diff.imag + math.pi) % (2 * math.pi) - math.pi)
        maxerr = max(maxerr, re, im); tested += 1
    return tested, maxerr

if __name__ == "__main__":
    path = sys.argv[1]
    d, rel = load(path)
    k, E = check_support(d, rel)
    print("d = 3^%d, support %d, E = %d, sum|e| = %d, max|e| = %d" %
          (k, len(rel), E, sum(abs(e) for _, e in rel), max(abs(e) for _, e in rel)))
    maxdig = max(len(str(abs(e))) for _, e in rel)
    if maxdig <= 6:
        okA, okR, bA, bR = products_direct(d, rel)
        print("direct big-integer products: A ok=%s (%d bits), R ok=%s (%d bits)" % (okA, bA, okR, bR))
        assert okA and okR
    else:
        print("direct products skipped (|e| has %d digits)" % maxdig)
    vA, vR = products_valuation(d, rel)
    print("sieve valuations: A ok=%s, R ok=%s" % (vA, vR)); assert vA and vR
    spec = "--spectral" in sys.argv
    info = seed_multisets(d, rel, want_multisets=spec)
    print("squares 4E = %d, s0 = %d, N = 4 s0 + 16 = %d (N/8 = %d)" % (info["nsq"], info["s0"], info["N"], info["N"] // 8))
    if spec:
        t, err = spectral_test(d, rel, info)
        print("spectral mod-2 Fourier test: %d characters, max |log diff| = %.2e" % (t, err))
        assert err < 1e-6
    print("ALL CHECKS PASSED")
