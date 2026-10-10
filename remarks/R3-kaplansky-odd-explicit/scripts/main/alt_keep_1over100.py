# Alternative instance that keeps the paper's constant 1/100 unchanged (fewer edited lines):
# (N,k,s) = (86,42,14).  Usage: /usr/bin/python3 main/alt_keep_1over100.py   (needs sympy)
from math import comb
from fractions import Fraction
from sympy import isprime
N, k, s = 86, 42, 14
m = comb(N, k); t = sum(comb(N, j) for j in range(k - s + 1))
assert 3 * k - 3 * (s - 1) > N                 # triangle-free (3*42-3*13 = 87 > 86)
assert Fraction(t, m) < Fraction(1, 100)
p = m + 1
while not (p % 4 == 1 and isprime(p)): p += 1
print("(N,k,s)=(86,42,14): m =", m, " t =", t, " t/m =", float(Fraction(t, m)), "< 1/100")
print("least prime p>m with p=1 mod 4:", p, " p-m =", p - m, " digits", len(str(p)))
print("denominators 1..%d, all < p" % (k - s))
