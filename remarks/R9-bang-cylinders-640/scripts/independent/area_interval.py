# Independent rigorous area bound (SageMath interval arithmetic, 256-bit RealIntervalField).
# Base of cylinder j (either family) is the orthogonal projection of the intercept triangle
#   P_j = {(0, y0, H t0): 0<=t0<=T_j, q_j t0 <= y0 <= q_{j+1} t0}   (area H * Delta * T_j^2 / 2 = H T_j^2 / n)
# onto v_j^perp; since the intercept plane x=0 has unit normal e1 and e1.v_j = 1,
#   |B_j| = |P_j| / |v_j|,  |v_j|^2 = 1 + eps^2 (alpha_j^2 + H^2 beta_j^2) = 1 + eps^2(alpha_j^2 + 2 beta_j^2).
# Family 2 is identical (w_j has the same length, intercept plane y=0 with e2.w_j = 1).
# A_min(K) = H (construction.tex Lemma 2.1, checked numerically below).  Hence
#   S / A_min = sum_j 2 T_j^2 / (n |v_j|).
# Usage: sage -python area_interval.py params_T.txt
import sys
from sage.all import RealIntervalField, QQ, sqrt
R = RealIntervalField(256)
eps = QQ(2) / 5; n = 320
T = {}
for line in open(sys.argv[1]):
    j, v = line.split(); T[int(j)] = QQ(v)
q = [QQ(-1) + QQ(2 * j) / n for j in range(n + 1)]
S = R(0)
for j in range(n):
    a = (1 + q[j] * q[j + 1]) / 4; b = (q[j] + q[j + 1]) / 4
    S += R(2 * T[j] ** 2 / n) / R(1 + eps ** 2 * (a ** 2 + 2 * b ** 2)).sqrt()
print("S/A_min in", S.str(style='brackets'))
print("upper endpoint", S.upper().n(60), " < 1/2 :", bool(S.upper() < QQ(1) / 2))
print("1/2 - S >=", (QQ(1) / 2 - S).lower().n(40))
print("number of cylinders:", 2 * n)
# sanity: untilted Bang cover (eps=0, T=1/2) has S/A_min = 1/2 exactly by the same formula
print("Bang check (eps=0,T=1/2):", sum(QQ(2) * QQ(1) / 4 / n for j in range(n)))
