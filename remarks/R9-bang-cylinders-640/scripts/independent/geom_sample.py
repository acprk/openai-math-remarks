# Direct geometric coverage test (independent of every LP / dual argument and of the first-crossing lemma).
# For each sample point p=(x,y,Ht) of K (exact rationals) we look for a cylinder containing p and confirm
# membership with exact Fraction arithmetic.  Membership in cylinder P_j + R v_j, v_j=(1,eps*al_j,H*eps*be_j):
#   the unique intercept with the plane x=0 is (0, y0, H t0), y0 = y - eps*x*al_j, t0 = t - eps*x*be_j,
#   and p is covered iff 0 <= t0 <= T_j and q_j t0 <= y0 <= q_{j+1} t0.
# Family 2 (Q_k + R w_k, w_k=(-eps*al_k,1,H*eps*be_k)), intercept with y=0: x0 = x + eps*y*al_k,
#   s0 = 1 - t + eps*y*be_k, covered iff 0 <= s0 <= T_k and q_k s0 <= x0 <= q_{k+1} s0.
# Floats are only used to ORDER candidate cylinders; acceptance is always exact.  A point is reported
# UNCOVERED only after all 640 cylinders fail exactly.
# Control: R9_SCALE=0.04 sage -python geom_sample.py params_control_half.txt 20261009  (must report uncovered points)
# Usage: sage -python geom_sample.py params_T.txt SEED [hard_points_file ...]
#   params file: "j T_j" lines.  Prints counts per sample class and the number of uncovered points.
import sys, os, random, numpy as np
from fractions import Fraction as F
eps = F(2, 5); n = 320
T = [None] * n
for line in open(sys.argv[1]):
    j, v = line.split(); T[int(j)] = F(v)
q = [F(-1) + F(2 * j, n) for j in range(n + 1)]
al = [(1 + q[j] * q[j + 1]) / 4 for j in range(n)]
be = [(q[j] + q[j + 1]) / 4 for j in range(n)]
qf = np.array([float(v) for v in q]); af = np.array([float(v) for v in al]); bf = np.array([float(v) for v in be])
Tf = np.array([float(v) for v in T]); ef = float(eps)

def cov1(x, y, t, j):
    t0 = t - eps * x * be[j]; y0 = y - eps * x * al[j]
    return 0 <= t0 <= T[j] and q[j] * t0 <= y0 <= q[j + 1] * t0
def cov2(x, y, t, k):
    s0 = 1 - t + eps * y * be[k]; x0 = x + eps * y * al[k]
    return 0 <= s0 <= T[k] and q[k] * s0 <= x0 <= q[k + 1] * s0

def slack(xs, ys, ts):  # float slack (min of the 4 defining inequalities) for all 2n cylinders, shape (m, 2n)
    X = xs[:, None]; Y = ys[:, None]; Tt = ts[:, None]
    t0 = Tt - ef * X * bf; y0 = Y - ef * X * af
    s1 = np.minimum.reduce([t0, Tf - t0, y0 - qf[:-1] * t0, qf[1:] * t0 - y0])
    s0 = 1 - Tt + ef * Y * bf; x0 = X + ef * Y * af
    s2 = np.minimum.reduce([s0, Tf - s0, x0 - qf[:-1] * s0, qf[1:] * s0 - x0])
    return np.hstack([s1, s2])

def check(pts):
    xs = np.array([float(p[0]) for p in pts]); ys = np.array([float(p[1]) for p in pts]); ts = np.array([float(p[2]) for p in pts])
    S = slack(xs, ys, ts); order = np.argsort(-S, axis=1)
    unc = []; minslack = S.max(axis=1).min()
    for i, (x, y, t) in enumerate(pts):
        assert 0 <= t <= 1 and abs(x) <= 1 - t and abs(y) <= t  # exact membership in K
        ok = False
        for c in list(order[i, :4]) + list(range(2 * n)):
            if (cov1(x, y, t, c) if c < n else cov2(x, y, t, c - n)): ok = True; break
        if not ok: unc.append((x, y, t))
    return unc, minslack

rng = random.Random(int(sys.argv[2])); D = 2 ** 24
def rat(a, b):  # random rational in [a,b] (a,b Fractions), denominator D * den
    return a + (b - a) * F(rng.randrange(D + 1), D)
def uniform_K():
    while True:  # rejection from the box [-1,1]^2 x [0,1]
        t = rat(F(0), F(1)); x = rat(F(-1), F(1)); y = rat(F(-1), F(1))
        if abs(x) <= 1 - t and abs(y) <= t: return (x, y, t)
def band():  # |t - 1/2| <= 1/20: the region where both families are near their cutoffs
    t = rat(F(9, 20), F(11, 20)); return (rat(-(1 - t), 1 - t), rat(-t, t), t)
def face():
    t = rat(F(0), F(1)); x = rat(-(1 - t), 1 - t); y = rat(-t, t); r = rng.randrange(6)
    if r == 0: x = 1 - t
    elif r == 1: x = -(1 - t)
    elif r == 2: y = t
    elif r == 3: y = -t
    elif r == 4: t = F(0); x = rat(F(-1), F(1)); y = F(0)
    else: t = F(1); x = F(0); y = rat(F(-1), F(1))
    return (x, y, t)
def edge_or_corner():
    t = rat(F(0), F(1)); r = rng.randrange(5)
    if r == 4:  # within 1e-3 of a vertex, still in K
        V = [(F(1), F(0), F(0)), (F(-1), F(0), F(0)), (F(0), F(1), F(1)), (F(0), F(-1), F(1))][rng.randrange(4)]
        while True:
            p = (V[0] + rat(F(-1, 1000), F(1, 1000)), V[1] + rat(F(-1, 1000), F(1, 1000)), V[2] + rat(F(-1, 1000), F(1, 1000)))
            if 0 <= p[2] <= 1 and abs(p[0]) <= 1 - p[2] and abs(p[1]) <= p[2]: return p
    return ((1 if r & 1 else -1) * (1 - t), (1 if r & 2 else -1) * t, t)

SC = float(os.environ.get('R9_SCALE', '1'))  # R9_SCALE<1 shrinks the random classes (used for the control run)
classes = [("uniform", uniform_K, 600000), ("band |t-1/2|<=1/20", band, 250000), ("faces", face, 100000), ("edges+near-vertices", edge_or_corner, 50000)]
vert = [(F(1), F(0), F(0)), (F(-1), F(0), F(0)), (F(0), F(1), F(1)), (F(0), F(-1), F(1)), (F(0), F(0), F(1, 2))]
total = 0; total_unc = 0
for name, gen, N in classes:
    N = int(N * SC)
    unc = 0; ms = 1e9
    for b in range(0, N, 5000):
        u, m = check([gen() for _ in range(min(5000, N - b))]); unc += len(u); ms = min(ms, m)
        if u and unc <= len(u): print("  example uncovered:", [tuple(map(float, p)) for p in u[:2]])
    print(f"{name}: {N} points, uncovered {unc}, min over points of best float slack {ms:.3e}", flush=True); total += N; total_unc += unc
u, m = check(vert); print(f"vertices+centre: {len(vert)} points, uncovered {len(u)}"); total += len(vert); total_unc += len(u)
for fn in sys.argv[3:]:  # exact LP maximizers "j k Z x y t" (one per sector pair) = hardest points
    pts = []
    for line in open(fn):
        if line.startswith('#'): continue
        w = line.split(); pts.append((F(w[3]), F(w[4]), F(w[5])))
    unc = 0
    for b in range(0, len(pts), 5000):
        u, m = check(pts[b:b + 5000]); unc += len(u)
    print(f"LP-maximizer points from {fn.split('/')[-1]}: {len(pts)}, uncovered {unc}", flush=True); total += len(pts); total_unc += unc
print(f"TOTAL {total} points, UNCOVERED {total_unc}")
