# Aggregates the outputs of pair_lp_exact.py: checks all 102400 pairs present, prints exact max Z_jk.
# Usage: python3 summarize_pairs.py pairs_all.txt
import sys
from fractions import Fraction as F
L = [l.split() for l in open(sys.argv[1]) if not l.startswith('#')]
keys = {(int(a[0]), int(a[1])) for a in L}
assert len(L) == 102400 and keys == {(j, k) for j in range(320) for k in range(320)}
Z = [(F(a[2]), int(a[0]), int(a[1])) for a in L]
m = max(Z)
arg = sorted((j, k) for z, j, k in Z if z == m[0])
print("pairs:", len(L), " pairs with Z_jk >= 0:", sum(1 for z in Z if z[0] >= 0))
print("exact max Z_jk =", m[0], "=", float(m[0]))
print("in units of eps^2 (primary script's U_jk):", float(m[0] / F(4, 25)))
print("attained at pairs:", arg)
print("COVERAGE CERTIFIED:", m[0] < 0)
