"""Exact re-check of the five rational inequalities of family 238,
"Routing densities and representation contraction for Thorp sweeps"
(build/sections/04-certificate.tex, l.259-294), whose accompanying script
support/network_bounds.py is not in the repository.  Exact rational
arithmetic (fractions.Fraction); floats only for the informational lines."""
from fractions import Fraction as F
import math

def h(j):
    M = 2**(j-1)
    from math import isqrt
    q = max(0, (isqrt(8*M+1)-1)//2)   # largest q with q(q+1) <= 2M
    while (q+1)*(q+2) <= 2*M: q += 1
    while q*(q+1) > 2*M: q -= 1
    assert q*(q+1) <= 2*M < (q+1)*(q+2)
    v = F(q, M) + (1 - F(q*(q+1), 2*M)) / (q+1)
    assert v == F(q, 2*M) + F(1, q+1)          # simplified form
    assert v*v <= F(2, M)           # h_j <= sqrt(2/M)
    return v

# monotonicity of allocation ratios 2^{j-1} h_j
r = [2**(j-1)*h(j) for j in range(1, 120)]
assert all(r[i] <= r[i+1] for i in range(len(r)-1)), "ratio not monotone"

TAIL_BOUND = F(169, 70*1024)                    # paper's bound on (1/2) sum_{j>20} h_j
assert 1 + math.sqrt(2) < 169/70
# exact-ish tail value for comparison
tail_true = sum(float(h(j)) for j in range(21, 160)) / 2

def H_upper(rho):
    """Upper bound used in the paper: spend 2^-20 on all j>20 (value <= TAIL_BOUND),
    then greedy from j=20 downward."""
    budget = F(rho) - F(1, 2**20)
    val = TAIL_BOUND
    for j in range(20, 0, -1):
        w = F(1, 2**j)
        x = min(F(1), budget / w) if budget > 0 else F(0)
        val += F(1, 2) * h(j) * x
        budget -= x * w
    return val

def H_exact_float(rho):
    budget = rho - 2.0**-20; val = tail_true
    for j in range(20, 0, -1):
        w = 2.0**-j; x = min(1.0, budget/w) if budget > 0 else 0.0
        val += 0.5*float(h(j))*x; budget -= x*w
    return val

def ln1m_upper(x):   # upper bound for -ln(1-x)
    x = F(x); return sum(x**i / i for i in range(1, 81)) + x**81 / (81*(1-x))
def ln2_upper():
    t = F(1, 3); s = 2*sum(t**(2*j+1)/(2*j+1) for j in range(40))
    return s + 2*t**81/(81*(1-F(1, 9)))

def bracket_upper(rho):
    rho = F(rho)
    if rho == 1:
        num = -rho              # (1-rho)ln(1-rho) -> 0
    else:
        num = -rho + (1-rho)*ln1m_upper(rho)   # -rho-(1-rho)ln(1-rho) <= this (negative)
    assert num < 0
    corr = num / (rho * ln2_upper())           # bigger denominator => upper bound
    return H_upper(rho) + corr

claims = [("H(1/8)",   H_upper(F(1, 8)),    F("0.791402")),
          ("H(13/50)", H_upper(F(13, 50)),  F("1.098069")),
          ("H(1)",     H_upper(F(1)),       F("1.958069")),
          ("bracket(7/25)", bracket_upper(F(7, 25)), F("0.904054")),
          ("bracket(1)",    bracket_upper(F(1)),     F("0.515374"))]
for name, val, bound in claims:
    print(f"{name:14s} upper={float(val):.9f}  claimed<{float(bound)}  OK={val < bound}")
print("tail true/2 = %.3e  vs bound %.3e" % (tail_true, float(TAIL_BOUND)))
for rho in (0.125, 0.26, 1.0):
    print("H_float(%g) = %.9f" % (rho, H_exact_float(rho)))
print("slopes: 4*h3 =", 4*h(3), " 2*h2 =", 2*h(2))
print("25/(36 ln2) =", 25/(36*math.log(2)), "< 25/24 =", 25/24)
