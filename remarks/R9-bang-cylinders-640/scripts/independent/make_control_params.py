# Negative-control parameters (should NOT cover): same 640 axes, cutoffs T_j = 1/2 + eps^2 d_j (per-sector margins removed).
# Run: python3 make_control_params.py > params_control_nomargin.txt
from fractions import Fraction as F
eps = F(2, 5); n = 320
q = [F(-1) + F(2 * j, n) for j in range(n + 1)]
d = lambda u: u * u * (1 + u * u) / 16
for j in range(n): print(j, F(1, 2) + eps * eps * max(d(q[j]), d(q[j + 1])))
