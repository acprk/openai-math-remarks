#!/usr/bin/env python3
"""Unit test (toy scale; NOT evidence about the asymptotic constant): exact full-law total-variation
distance of the Thorp shuffle from uniform for n = 2^d cards, d = 2, 3, in exact rational arithmetic.

Environment: Python >= 3.8, standard library only.  Single thread, ~1 min.
Run:         python3 small_d_exact.py > expected_output/small_d_exact.txt

Model (paper, Section 2): positions F_2^d; one physical shuffle = R S with S an independent fair swap on
each pair {x, x+e_1}.  We use the cyclic-coordinate frame X_t = R^{-t} Y_t (switch t acts in direction
1 + (t mod d)); TV to uniform is the same in both frames because they differ by a fixed position
permutation.  Expected: t_mix(2) = 3 and t_mix(3) = 5 (distance <= 1/4), matching the support lower
bound 2^{tn/2} >= 3n!/4 at d = 2 and exceeding it by one at d = 3.
"""
from fractions import Fraction
from math import factorial, ceil, log2


def tv_curve(d, T):
    n = 1 << d
    law = {tuple(range(n)): Fraction(1)}         # g[label] = position, start at identity
    out = []
    for t in range(T):
        e = 1 << (t % d)
        for x in range(n):
            if x & e:
                continue
            y = x | e
            new = {}
            for g, p in law.items():
                h = tuple(y if q == x else x if q == y else q for q in g)
                half = p / 2
                new[g] = new.get(g, 0) + half
                new[h] = new.get(h, 0) + half
            law = new
        u = Fraction(1, factorial(n))
        tv = (sum(abs(p - u) for p in law.values()) + u * (factorial(n) - len(law))) / 2
        out.append((t + 1, len(law), tv))
    return out


for d, T in ((2, 5), (3, 7)):
    n = 1 << d
    lb = ceil((2 / n) * log2(3 * factorial(n) / 4))
    print(f"d={d}, n={n}: support lower bound t_mix >= {lb}")
    tmix = None
    for t, supp, tv in tv_curve(d, T):
        print(f"  t={t}: support={supp:6d}  TV={str(tv):>28s} = {float(tv):.6f}")
        if tmix is None and tv <= Fraction(1, 4):
            tmix = t
    print(f"  exact t_mix(1/4) = {tmix}")
    assert tmix == {2: 3, 3: 5}[d]
print("ALL CHECKS PASSED")
