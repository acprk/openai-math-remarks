# Independent exact coverage check for the 640-cylinder cover (does not reuse primary/cyl_cert.py).
# Environment: SageMath (sage -python), exact rational LP via the Parma Polyhedra Library backend (solver='PPL').
# Geometry (written directly from construction.tex:111-123, not from the F_i/G_i sequences):
#   point (x,y,Ht) in K: 0<=t<=1, |x|<=1-t, |y|<=t.
#   family 1, sector j: intercept t0 = t - eps*x*beta_j, y0 = y - eps*x*alpha_j;
#       wedge W_j: q_j t0 <= y0 <= q_{j+1} t0;  cylinder j covers the point iff in W_j and t0 <= T_j.
#   family 2, sector k: s0 = 1 - t + eps*y*beta_k, x0 = x + eps*y*alpha_k;
#       wedge V_k: q_k s0 <= x0 <= q_{k+1} s0;  cylinder k covers iff in V_k and s0 <= T_k.
# Every point of K lies in some W_j and some V_k (first-crossing lemma, coverage.tex:13-45, valid for all eps).
# For every pair (j,k) we compute EXACTLY  Z_jk = max { min(t0-T_j, s0-T_k) : p in K cap W_j cap V_k }.
# Z_jk < 0 for all pairs  ==>  every point is covered by its selected family-1 or family-2 cylinder.
# Usage: sage -python pair_lp_exact.py params_T.txt J0 J1 > out_J0_J1.txt   (rows j in [J0,J1))
import sys
from fractions import Fraction as F
from sage.all import MixedIntegerLinearProgram, QQ
eps = QQ(2) / 5; n = 320
T = {}
for line in open(sys.argv[1]):
    j, v = line.split(); T[int(j)] = QQ(v)
q = [QQ(-1) + QQ(2 * j) / n for j in range(n + 1)]
al = [(1 + q[j] * q[j + 1]) / 4 for j in range(n)]
be = [(q[j] + q[j + 1]) / 4 for j in range(n)]
J0, J1 = int(sys.argv[2]), int(sys.argv[3])
worst = None
for j in range(J0, J1):
    for k in range(n):
        P = MixedIntegerLinearProgram(solver='PPL', maximization=True)
        v = P.new_variable(real=True, nonnegative=False)
        x, y, t, z = v['x'], v['y'], v['t'], v['z']
        for c in (t >= 0, t <= 1, x <= 1 - t, -x <= 1 - t, y <= t, -y <= t):
            P.add_constraint(c)
        t0 = t - eps * be[j] * x; y0 = y - eps * al[j] * x
        s0 = 1 - t + eps * be[k] * y; x0 = x + eps * al[k] * y
        P.add_constraint(q[j] * t0 <= y0); P.add_constraint(y0 <= q[j + 1] * t0)
        P.add_constraint(q[k] * s0 <= x0); P.add_constraint(x0 <= q[k + 1] * s0)
        P.add_constraint(z <= t0 - T[j]); P.add_constraint(z <= s0 - T[k])
        P.set_objective(z)
        Z = QQ(P.solve())
        X, Y, Tt = (QQ(P.get_values(w)) for w in (x, y, t))
        print(j, k, Z, X, Y, Tt, flush=False)
        if worst is None or Z > worst[0]: worst = (Z, j, k)
sys.stdout.flush()
print("#WORST", worst[0], float(worst[0]), worst[1], worst[2], file=sys.stderr)
