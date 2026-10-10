"""Independent re-implementation (does not import or copy 238_network_bounds_recheck.py).

Family 238, 'Routing densities and representation contraction for Thorp sweeps',
build/sections/04-certificate.tex.  Works only from the paper's DEFINITIONS:
  h_j = max{ E(1/L) : L positive integer, P(L=l) <= l/2^{j-1} }          (l.134-136)
  H(rho) = max{ (1/2) sum_j h_j x_j : 0<=x_j<=1, sum_j 2^{-j} x_j <= rho } (l.138-140)
  bracket(rho) = H(rho) + (-rho-(1-rho)ln(1-rho))/(rho ln 2)              (l.148)
and certifies the five claimed inequalities (l.272-282)
  H(1/8)<.791402, H(13/50)<1.098069, H(1)<1.958069,
  bracket(7/25)<.904054, bracket(1)<.515374.
Method: h_j exactly by the greedy solution of its LP (mass l/M on l=1,2,...; brute force for j<=12);
H by a generic fractional-knapsack over j<=J=100 (items sorted by value/weight),
plus a rigorous tail: E(1/L) <= K/M + 1/(K+1) <= 2/sqrt(M) with K=isqrt(M), so
(1/2) sum_{j>J} h_j <= 2^{-J/2}/(1-2^{-1/2}) < 4*2^{-J/2}.  Logs by mpmath interval
arithmetic (mpmath.iv, 60 digits).  Environment: python3, mpmath>=1.2.  < 1 s.
"""
from fractions import Fraction as Fr
from math import isqrt
from mpmath import iv, mpf
iv.dps = 60

def h(j):
    # greedy solution of the LP: P(L=l)=l/M for l=1..q (largest q with sum <= 1),
    # remaining mass on l=q+1.  q is found by bisection (no closed formula used).
    M = 2 ** (j - 1)
    lo, hi = 0, 2 * M + 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid * (mid + 1) <= 2 * M: lo = mid
        else: hi = mid
    q = lo
    val = Fr(q, M) + (1 - Fr(q * (q + 1), 2 * M)) / (q + 1)
    K = isqrt(M); assert val <= Fr(K, M) + Fr(1, K + 1)   # tail lemma, checked where computable
    if j <= 12:  # brute-force cross-check of the greedy against the definition
        rem, bf, l = Fr(1), Fr(0), 1
        while rem > 0:
            p = min(Fr(l, M), rem); bf += p / l; rem -= p; l += 1
        assert bf == val
    return val

J = 100
items = [(j, h(j) / 2, Fr(1, 2 ** j)) for j in range(1, J + 1)]
TAIL = 4 * Fr(1, 2 ** (J // 2))          # rigorous bound on (1/2) sum_{j>J} h_j

def H_up(rho):
    rho = Fr(rho); val = Fr(0); budget = rho
    for j, v, w in sorted(items, key=lambda t: t[1] / t[2], reverse=True):
        x = min(Fr(1), budget / w); val += v * x; budget -= x * w
        if budget <= 0: break
    return val + TAIL                     # upper bound (tail items at most fully used)

def H_low(rho):                           # lower bound: truncated LP alone
    return H_up(rho) - TAIL

def bracket_iv(rho):
    r = iv.mpf(Fr(rho).numerator) / Fr(rho).denominator
    ent = (-r) if rho == 1 else (-r - (1 - r) * iv.log(1 - r))
    Hu = iv.mpf(H_up(rho).numerator) / H_up(rho).denominator
    return Hu + ent / (r * iv.log(2))

ok_all = True
for name, rho, claim in [("H(1/8)", Fr(1, 8), "0.791402"), ("H(13/50)", Fr(13, 50), "1.098069"),
                         ("H(1)", Fr(1), "1.958069")]:
    up = H_up(rho); ok = up < Fr(claim); ok_all &= ok
    print(f"{name:14s} in [{float(H_low(rho)):.10f}, {float(up):.10f}]  claim <{claim}  OK={ok}")
for name, rho, claim in [("bracket(7/25)", Fr(7, 25), "0.904054"), ("bracket(1)", Fr(1), "0.515374")]:
    b = bracket_iv(rho); ok = b.b < mpf(claim); ok_all &= bool(ok)
    print(f"{name:14s} upper={float(b.b):.10f}  claim <{claim}  OK={ok}")
print("tail bound used:", float(TAIL))
print("ALL FIVE OK" if ok_all else "FAILURE")
