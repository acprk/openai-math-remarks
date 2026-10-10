#!/usr/bin/env python3
"""Interval-arithmetic certificate for the numerical constants of remark R5 (Thorp shuffle).

Environment: Python >= 3.8 with mpmath (tested: CPython 3.10.12, mpmath 1.3.0).
Run:         python3 constants_certificate.py > expected_output/constants_certificate.txt
Runtime:     about 2 minutes, single thread.

Every inequality that the note relies on is checked here either with exact rationals
(fractions.Fraction) or with mpmath interval arithmetic (mpmath.iv, outward rounding);
the script raises AssertionError if any check fails.  Nothing in the note's proof depends
on floating point beyond what is certified below.

Notation (note, Section 2):
  L(b)     = log2(2b/(2b-1))             asymptotic per-step energy decay exponent for a
                                         factor whose lists leave >= n/b free positions
  kappa(b) = 2/L(b)                      per-factor cost (in units of d) of the first stage
  e(b)     = 1 - 2/b                     per-factor Fourier exponent gain (sigma -> 0)
Admissibility of a multiset {b_1..b_r} of block numbers (all b_i >= 3):
  route O (paper's orbit-span operator bound + ||B||_HS <= sqrt(D)||B||_op):  sum e(b_i) > 2
  route H (companion's isotypic Hilbert-Schmidt bound for one factor):        sum e(b_i) > 1
Constant of a route = min over admissible multisets of sum kappa(b_i).
Factors with b <= 2 have e(b) <= 0 and positive cost, so they never help and are excluded.
"""
from fractions import Fraction as F
from mpmath import iv, mpf, mp, nstr


def S(x, k):
    # fixed-point string with k decimals (x is an mpf endpoint, already outward-rounded)
    x = mpf(x)
    return nstr(x, k + len(str(int(abs(x)))) , min_fixed=-1e9, max_fixed=1e9)
import math

iv.dps = 50
LN2 = iv.log(2)


def L(b):
    b = iv.mpf(b)
    return iv.log(2 * b / (2 * b - 1)) / LN2


def kappa(b):
    return 2 / L(b)


def kappa_lower_elementary(b):
    # L(b) = log2(1 + 1/(2b-1)) <= 1/((2b-1) ln 2)  =>  kappa(b) >= 2 (2b-1) ln 2.
    return 2 * (2 * b - 1) * math.log(2) * (1 - 1e-12)


def e(b):
    return 1 - F(2, b)


def enumerate_route(threshold, X, label):
    """All multisets (b_i >= 3) with sum e(b_i) > threshold and certified cost lower end <= X."""
    k3 = kappa(3)
    rmin = int(threshold) + 1          # each e(b) < 1, so need more than `threshold` factors
    # any factor in a multiset of cost <= X satisfies kappa(b) <= X - (rmin-1) kappa(3)
    budget = X - (rmin - 1) * float(k3.a)
    Bmax = 3
    while kappa_lower_elementary(Bmax + 1) <= budget:
        Bmax += 1
    # check elementary lower bound is consistent with the interval values on the range
    for b in range(3, Bmax + 2):
        assert kappa(b).a >= kappa_lower_elementary(b)
        assert kappa(b + 1).a > kappa(b).b            # kappa strictly increasing
    rmax = int(X // float(k3.a))
    K = {b: kappa(b) for b in range(3, Bmax + 1)}
    out = []

    def rec(start, chosen, cost_lo, ex):
        if ex > threshold:
            cost = sum((K[b] for b in chosen), iv.mpf(0))
            out.append((tuple(chosen), cost, ex))
        if len(chosen) == rmax:
            return
        for b in range(start, Bmax + 1):
            c = cost_lo + float(K[b].a)
            if c > X:
                break
            rec(b, chosen + [b], c, ex + e(b))

    rec(3, [], 0.0, F(0))
    out.sort(key=lambda t: float(t[1].a))
    print(f"--- route {label}: sum(1-2/b_i) > {threshold};  search: b in [3,{Bmax}], r <= {rmax}, cost <= {X}")
    print(f"    completeness: kappa(b) >= 2(2b-1)ln2 > {budget:.3f} for b > {Bmax}; r > {rmax} costs > {X}")
    for comb, cost, ex in out[:8]:
        print(f"    {str(comb):22s} cost in [{S(cost.a,10)}, {S(cost.b,10)}]   sum e = {ex} = {float(ex):.5f}")
    best, second = out[0], out[1]
    assert best[1].b < second[1].a, "optimum not separated from runner-up"
    print(f"    certified optimum {best[0]}: {S(best[1].a,12)} <= c <= {S(best[1].b,12)};"
          f" runner-up {second[0]} >= {S(second[1].a,12)}")
    return best


def sigma_window(comb, threshold):
    # exponent with sigma>0: sum(1-(2+sigma)/b) - threshold > 0  <=>  sigma < (sum e - threshold)/sum(1/b)
    ex = sum(e(b) for b in comb) - threshold
    inv = sum(F(1, b) for b in comb)
    s = ex / inv
    print(f"    {comb}: slack sum e - {threshold} = {ex}; admissible for every sigma < {s} = {float(s):.6f}")
    return s


def fractional_relaxation(threshold):
    """Fractional relaxation: inf over real b >= 3 of kappa(b)/e(b), certified two-sided,
    plus the minimum over integer b >= 3 (attained at b = 4).
    Any integer multiset with sum e > threshold has cost = sum (kappa/e) e > threshold * min ratio.
    This is a strict LOWER BOUND for the routes' constants, not a limit that is attained."""
    # On [x, x+h] both kappa(b) and e(b) = 1-2/b are increasing, hence
    #   kappa(b)/e(b) >= kappa(x)/e(x+h)   (endpoints evaluated in interval arithmetic).
    # Coarse cells of width 1/1024 cover [3, 200] outside [3.75, 4]; fine cells of width 2^-17
    # cover [3.75, 4].  Tail b > 200: ratio >= kappa(b) >= 2(2b-1)ln2 > 550.
    def cell_min(a, ncell, h):
        best = None
        for i in range(ncell):
            x = a + i * h
            v = (kappa(x) / (1 - 2 / iv.mpf(x + h))).a
            best = v if best is None or v < best else best
        return best
    hc = mpf(1) / 1024
    coarse = min(cell_min(mpf(3), 768, hc), cell_min(mpf(4), 196 * 1024, hc))
    fine = cell_min(mpf(3.75), 2**15, mpf(1) / 2**17)
    assert coarse > fine
    lo_best = fine
    assert kappa_lower_elementary(200) > 550
    # upper end of the infimum: evaluate at a point near the minimiser (interval upper endpoint)
    xm = mpf("3.8646511")
    up = (kappa(xm) / (1 - 2 / iv.mpf(xm))).b
    assert lo_best <= up and up - lo_best < 1e-4
    # integer points: kappa(b)/e(b) for b = 3..200 (tail as above); minimum at b = 4
    ratios = {b: kappa(b) / (1 - 2 / iv.mpf(b)) for b in range(3, 201)}
    bint = min(ratios, key=lambda b: ratios[b].a)
    assert bint == 4 and all(ratios[4].b < ratios[b].a for b in ratios if b != 4)
    print(f"    fractional relaxation: inf_(real b>=3) kappa(b)/(1-2/b) in [{S(lo_best,6)}, {S(up,6)}]"
          f"  (minimiser b* ~ {S(xm,4)})")
    print(f"    integer points: min_(b>=3) kappa(b)/(1-2/b) = kappa(4)/(1/2) in [{S(ratios[4].a,6)}, {S(ratios[4].b,6)}]")
    print(f"    => every integer multiset admissible for threshold {threshold} costs > {threshold} x {S(ratios[4].a,6)}"
          f" (strict lower bound; not attained)")
    return lo_best, up, ratios[4]


print("=== per-factor costs kappa(b) = 2/log2(2b/(2b-1)), certified intervals ===")
for b in range(3, 9):
    k = kappa(b)
    print(f"  b={b}: kappa in [{S(k.a,12)}, {S(k.b,12)}]")

print("\n=== original constants for comparison ===")
g = iv.mpf(31) / 32
T0 = 2 + 3 / (-iv.log(g) / LN2)
print(f"  paper (1600d): 8 factors x 200d.  Paper's own inequalities need T > (2 + 3/log2(32/31)) d = "
      f"[{S(T0.a,6)}, {S(T0.b,6)}] d per factor")
print("  family's best explicit constant before this note: 512d (Random-coordinate-frames, Cor. 11.4,"
      " label shear:cor:512; 2 factors x 256d, already uses the companion's two-factor HS transfer)")

print("\n=== route O: orbit-span operator transfer of the main paper (four factors) ===")
bestO = enumerate_route(2, 47.0, "O")
sO = sigma_window(bestO[0], 2)
sigma_window((4, 4, 4, 5), 2)
print("\n=== route H: companion's isotypic Hilbert--Schmidt transfer (two factors) ===")
bestH = enumerate_route(1, 27.0, "H")
sH = sigma_window(bestH[0], 1)
assert bestO[0] == (3, 4, 5, 5) and bestH[0] == (4, 5)
assert 44.30 < bestO[1].a and bestO[1].b < 44.301
assert 23.539 < bestH[1].a and bestH[1].b < 23.540

print("\n=== fractional relaxation (strict lower bound for each route with these stage-1 rates) ===")
rlo, rup, rint = fractional_relaxation(1)
assert 20.738 < rlo and rup < 20.7383 and 20.7635 < rint.a and rint.b < 20.7637
print(f"    route H >= {S(rint.a,6)} (> {S(rlo,6)});  route O >= {S(2*rint.a,6)} (> {S(2*rlo,6)})")

print("\n=== ratios (512d = family's best, which already has the x4 two-factor step) ===")
for name, c in (("route H optimum", bestH[1]), ("route O optimum", bestO[1])):
    print(f"  1600 / {name} in [{S((1600 / c).a,4)}, {S((1600 / c).b,4)}];"
          f"  512 / {name} in [{S((512 / c).a,4)}, {S((512 / c).b,4)}]")

print("\n=== loss ledger, route H (cumulative; each line keeps all previous changes) ===")
def Lbeta(beta):
    return -iv.log(1 - iv.mpf(beta) / 2) / LN2
rows = [
    ("paper as written: 8 factors x 200d", iv.mpf(1600)),
    ("paper's own inequalities, minimal T per factor (8 x (2 + 3/log2(32/31)))", 8 * (2 + 3 / Lbeta(iv.mpf(1) / 16))),
    ("companion two-factor HS transfer, b=8 [already in 512d corollary; not ours]", 2 * (2 + 3 / Lbeta(iv.mpf(1) / 16))),
    ("+ concentration theta->0: beta 1/16 -> 1/8", 2 * (2 + 3 / Lbeta(iv.mpf(1) / 8))),
    ("+ KL chain rule: need n^2 E_T -> 0 instead of n^3 E_T -> 0", 2 * (2 + 2 / Lbeta(iv.mpf(1) / 8))),
    ("+ local look-back: offset 2d -> O(log d)", 2 * (2 / Lbeta(iv.mpf(1) / 8))),
    ("+ mixed block numbers {4,5} instead of {8,8}", kappa(4) + kappa(5)),
]
prev = None
for name, v in rows:
    step = "" if prev is None else f"   step factor {S((prev / v).a, 3)}"
    print(f"  {name:76s} {S(v.a, 3):>9s} d{step}")
    prev = v
print("\n=== loss ledger, route O (paper-internal transfer only; cumulative) ===")
rowsO = [
    ("paper as written: 8 factors x 200d", iv.mpf(1600)),
    ("paper's own inequalities; 5 factors suffice for its Lemma 5.1", 5 * (2 + 3 / Lbeta(iv.mpf(1) / 16))),
    ("+ concentration theta->0: beta 1/16 -> 1/8", 5 * (2 + 3 / Lbeta(iv.mpf(1) / 8))),
    ("+ KL chain rule", 5 * (2 + 2 / Lbeta(iv.mpf(1) / 8))),
    ("+ local look-back", 5 * (2 / Lbeta(iv.mpf(1) / 8))),
    ("+ any sigma>0 in the degree lemma: b=8 needs only 3 factors", 3 * (2 / Lbeta(iv.mpf(1) / 8))),
    ("+ mixed block numbers {3,4,5,5}", kappa(3) + kappa(4) + 2 * kappa(5)),
]
prev = None
for name, v in rowsO:
    step = "" if prev is None else f"   step factor {S((prev / v).a, 3)}"
    print(f"  {name:76s} {S(v.a, 3):>9s} d{step}")
    prev = v
print("\nALL CHECKS PASSED")
