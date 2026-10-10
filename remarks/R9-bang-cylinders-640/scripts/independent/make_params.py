# Exports the exact rational parameters of the 640-cylinder cover (eps = 2/5, n = 320 sectors per family).
# Input: persector_0.4_320.npy (float per-sector margins e_j from primary/cyl_persector.py).
# Rounding rule (identical to primary/cyl_cert.py, so that the same cover is certified):
#   E_j = Fraction(e_j).limit_denominator(10^7) + 1/10^7,   T_j = 1/2 + eps^2 * (d_j + E_j),
#   d_j = max(d(q_j), d(q_{j+1})), d(q) = q^2 (1+q^2)/16, q_j = -1 + 2j/n.
# The same T_j is used for both families (x<->y symmetry).  Output: params_T.txt, one line "j T_j" (exact).
# Run: sage -python make_params.py ../primary/persector_0.4_320.npy > params_T.txt
import sys, numpy as np
from fractions import Fraction as F
eps = F(2, 5); n = 320
ev = np.load(sys.argv[1])[0]; assert len(ev) == n
q = [F(-1) + F(2 * j, n) for j in range(n + 1)]
d = lambda u: u * u * (1 + u * u) / 16
for j in range(n):
    E = F(float(ev[j])).limit_denominator(10 ** 7) + F(1, 10 ** 7)
    T = F(1, 2) + eps * eps * (max(d(q[j]), d(q[j + 1])) + E)
    print(j, T)
