#!/usr/bin/env python3
"""
Independent check (exact integer arithmetic only, no floating point, no mpmath) of the exponent of
eps_n = (n^{3/2}/2) sqrt((31/32)^{198d} + 32 e^{-n/256}), n = 2^d  (family 238 main paper, eq. (h-8)).
With e* = 99 log2(32/31) - 3/2 we show 3.03456 < e* < 3.03457 by comparing integers:
   e* > a/N  <=>  (32/31)^{99N} > 2^{N(a/N + 3/2)}  <=>  2^{5*99N - (a + 3N/2)} > 31^{99N}.
Run: python3 indep_eps.py   (~5 s)
"""
from fractions import Fraction
def exponent_greater(x):        # is 99*log2(32/31) - 3/2 > x  (x a Fraction)?
    y = x + Fraction(3, 2)      # need 99*log2(32/31) > y  <=> (32/31)^99 > 2^y
    p, q = y.numerator, y.denominator
    # (32^99/31^99)^q > 2^p  <=> 2^{495q - p} > 31^{99q}
    return 2 ** (495 * q - p) > 31 ** (99 * q)
lo, hi = Fraction(303456, 100000), Fraction(303457, 100000)
print(f"e* > {float(lo)}: {exponent_greater(lo)};  e* > {float(hi)}: {exponent_greater(hi)}")
assert exponent_greater(lo) and not exponent_greater(hi)
# also: 198 log2(32/31) > 8 (the paper's step) and > 9
for t in (8, 9):
    print(f"198*log2(32/31) > {t}: {2 ** (990 - t) > 31 ** 198}")
# the paper's O(n^{-5/2}) is the special case obtained by replacing 198*log2(32/31) with 8
print(f"=> eps_n = Theta(n^(-e*)), e* in ({float(lo)}, {float(hi)}); stated exponent 5/2 is not sharp")
