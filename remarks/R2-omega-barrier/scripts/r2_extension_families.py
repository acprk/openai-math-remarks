#!/usr/bin/env python3
"""
R2 extension check: the balanced assignment
    lambda*(C(a,b)^pi) = P*(a,b)^(3/4),  lambda*(B_L(m)) = m^(3/4),  lambda*(T_m) = m^(9/4)
against the additional VALID inequality families found in the Claude <-> Codex
(gpt-6-astra) rounds (CODEX-ROUND1/2/3).  This is an independent re-implementation
(it shares no code with the round transcripts).  Every comparison is exact:
rationals, or integer 4th-root brackets for the (3/4)-powers.

Environment: Python >= 3.8, standard library only.
Run:  /usr/bin/python3 r2_extension_families.py
Expected output: expected_output/r2_extension_families.txt  (runtime ~2-4 min, 1 thread)

These finite checks SUPPLEMENT the all-size proofs listed in note.tex / VERIFICATION.md;
for families marked "finite only" there they are the only evidence.
"""
from fractions import Fraction as Fr
from itertools import product
from math import isqrt
import sys

PHI = [None, Fr(1)]
for m in range(1, 2000):
    PHI.append(PHI[m] * Fr(3 * m + 1, 3 * m))

def F(a, b):                       # P*(a,b)
    m, M = (a, b) if a <= b else (b, a)
    return PHI[m] * (M + Fr(m - 1, 2))

SCALE = 10 ** 30
def p34(x):
    """Return (lo, hi) rationals with lo <= x^(3/4) <= hi (x > 0 rational)."""
    c = x ** 3
    n = c.numerator * SCALE ** 4 // c.denominator
    r = isqrt(isqrt(n))            # floor of 4th root of floor(n)
    lo = Fr(r, SCALE)
    hi = Fr(r + 1, SCALE)
    assert lo ** 4 <= c <= hi ** 4
    return lo, hi

stats = {}
def check(name, lhs_lo, rhs_hi):
    s = stats.setdefault(name, [0, 0])
    s[0] += 1
    if not lhs_lo >= rhs_hi:
        s[1] += 1
        if s[1] <= 3:
            print("  VIOLATION", name, float(lhs_lo), float(rhs_hi))

# (A) transfer:  m P(a,b) >= sum_{i in I} P(a+m-1-2i, b-m+1+2i)      [R1 sec.2]
# (B) higher dot: m P(a,b) >= sum_{i<min(b,m)} P(a, b+m-1-2i)         [R1 sec.3]
for a, b in product(range(1, 41), repeat=2):
    for m in range(1, 25):
        r = m - 1
        rhsA = sum((F(a + r - 2 * i, b - r + 2 * i) for i in range(min(a, m))
                    if 0 <= r - i < min(b, m)), Fr(0))
        check("(A) transfer along a+b=const", m * F(a, b), rhsA)
        rhsB = sum((F(a, b + r - 2 * i) for i in range(min(b, m))), Fr(0))
        check("(B) B_X(m) determinant filtration", m * F(a, b), rhsB)

# (C) bivariate shared output, every k                                [R1 sec.4]
# (D) full direct sum, checked with the real 3/4 power                 [R1 sec.5]
for a, b, c, d in product(range(1, 12), repeat=4):
    lhs = F(a, b) * F(c, d)
    rX, rY = min(a, c), min(b, d)
    for k in range(rX + rY - 1):
        rhs = sum((F(a + c - 1 - 2 * i, b + d - 1 - 2 * (k - i)) for i in range(rX)
                   if 0 <= k - i < rY), Fr(0))
        check("(C) bivariate grading, shared output", lhs, rhs)
    R = min(a, b, c, d)
    terms = [F(a + c - 1 - 2 * i, b + d - 1 - 2 * i) for i in range(R)]
    if len(terms) == 1:
        check("(D) full direct sum [P P]^t >= sum P^t", lhs, terms[0])
    else:
        check("(D) full direct sum [P P]^t >= sum P^t", p34(lhs)[0], sum(p34(x)[1] for x in terms))

# Claude's earlier families (barrier_extra.py): unequal 3-chains, 5-chains, Kronecker
for a in range(1, 16):
    for h1, h2, h3 in product(range(1, 25), repeat=3):
        check("unequal 3-sector chain", F(a, h1 + h2 + h3 + a - 1), F(a, h1) + F(a, h2) + F(a, h3))
    for hs in product(range(1, 8), repeat=5):
        check("5-sector alternating chain", F(a, sum(hs) + 2 * (a - 1)), sum(F(a, h) for h in hs))
for a, b, a2, b2 in product(range(1, 10), repeat=4):
    Nn = a + b - 1
    check("Kronecker restriction", F(a + (a2 - 1) * Nn, b + (b2 - 1) * Nn), F(a, b) * F(a2, b2))

# (F) Fourier filtering with explicit cost: (2M-1) P(Ma,Mb) >= M P(a,b)  [R1 sec.6]
for M, a, b in product(range(1, 13), repeat=3):
    check("(F) Fourier filtering with cost", (2 * M - 1) * F(M * a, M * b), M * F(a, b))

# (O) outer product P <= ab ; (L) P >= max(a,b) ; (5) sandwich         [R1 sec.7, R2 sec.1.2]
for a, b in product(range(1, 151), repeat=2):
    check("(O) P(a,b) <= ab", a * b, F(a, b))
    check("(L) P(a,b) >= max(a,b)", F(a, b), max(a, b))
    check("(5) sandwich lower: P^9 >= (ab(a+b-1))^4", F(a, b) ** 9, (a * b * (a + b - 1)) ** 4)
    check("(5) sandwich upper: (3(a+b)-2)^4/4^4 >= P^3", Fr(3 * (a + b) - 2, 4) ** 4, F(a, b) ** 3)
    check("(4) P(a,b) >= a+b-1", F(a, b), a + b - 1)

for n in range(1, 400):
    f = F(n, n)
    # (13) wrap-around C(n,n) >= C[Z_n] ~ <n>:  P(n,n)^(3/4) >= n
    check("(13) wrap-around P^t >= n", f ** 3, Fr(n) ** 4)
    # (6) scalar peeling C(n,n) ~> 1 (+) C(n-1,n-1)
    if n >= 2:
        check("(6) scalar peeling P_n^t >= 1 + P_{n-1}^t", p34(f)[0], 1 + p34(F(n - 1, n - 1))[1])
    # (16) small-CW laser: P(n,n)^9 >= (27/4)^4 (n-2)^9
    if n >= 3:
        check("(16) small CW laser", f ** 9, Fr(27, 4) ** 4 * (n - 2) ** 9)
    # (24) mixed T_m (x) C(n,n), m=2 sufficient condition
    check("(24) T_2 (x) C(n,n) sufficient condition", 2 * f ** 3, (Fr(3, 8) * (4 * n - 1)) ** 4)

# (21) T_m ~> (+) C's : uses sum F^(3/4) <= (9/8) m^2 and m^(1/4) > 9/8 for m >= 2
check("(21) constant 2*8^4 > 9^4", 2 * 8 ** 4, 9 ** 4 + 1)

# zero/rest laser (R3 sec.5): U(q)^3 <= R_q, need R_q^4 < F(q+1,q+1)^9
for q in range(1, 16):
    if q == 1:
        Rq_hi = Fr(12)
    else:
        lo, hi = p34(Fr(q) ** 3)          # q^(9/4) = (q^3)^(3/4)
        Rq_hi = 16 + Fr(16, 3) * hi + Fr(16, 9) * q ** 3
    check("zero/rest laser head q<=15", F(q + 1, q + 1) ** 9, Rq_hi ** 4 + Fr(1, 10 ** 40))
check("zero/rest laser tail constant (22/9+1/256)^4 < (27/8)^3", Fr(27, 8) ** 3, (Fr(22, 9) + Fr(1, 256)) ** 4 + Fr(1, 10 ** 40))
for q in range(16, 400):   # spot-check of the tail argument's two inequalities
    check("zero/rest tail: F(q+1,q+1)^3 > (27/8) q^4", F(q + 1, q + 1) ** 3, Fr(27, 8) * q ** 4)

# Codex R1 slack formula for (A2) (proof of all-size (A)):
for a in range(2, 120):
    for b in range(a, 240):
        slack = 2 * F(a, b) - F(a + 1, b - 1) - F(a - 1, b + 1)
        exp = PHI[a] * (Fr(2, 3 * a - 2) if b in (a, a + 1) else Fr(2 * (2 * a + b - 1), 3 * a * (3 * a - 2)))
        s = stats.setdefault("(A2) slack closed form (R1)", [0, 0]); s[0] += 1
        if slack != exp:
            s[1] += 1

bad = 0
for name, (n, f) in stats.items():
    print(f"{name:52s} checks={n:8d} failures={f}")
    bad += f
print("ALL PASS" if bad == 0 else f"FAILURES: {bad}")
sys.exit(1 if bad else 0)
