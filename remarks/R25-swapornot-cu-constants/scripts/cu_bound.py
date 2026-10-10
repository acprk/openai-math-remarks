# Rigorous upper bound for C_u = sup_m Z_u(m), Z_u(m) = sum_{lambda |- m} f_lambda^{-u}.
# Run with:  sage -python cu_bound.py <u as p/q> [M0=55] [ME=200] [M1=1000] [K1=auto]
# Region 1 (m<=M0): exact enumeration (exact_tables_M0.pkl, integer f_lambda), sums in MPFI intervals (RIF).
# Region 2 (M0<m<=ME): 3-D DP over classes (k, lambda_1=r, ell=c) using Lemma (peel): for lambda=(L,mu), mu|-k,
#     f_lambda >= f_mu * N0(L,k),  N0(L,k) = f_{(L,min(k,L))}   (+ conjugate symmetry).
#     Float arithmetic, every vector op inflated by INF=1+2^-40 => all stored numbers are upper bounds.
# Region 3 (ME<m<=M1): 1-D recursion Z(m) <= 2[1 + sum_{L>=ceil(sqrt m)} N0(L,m-L)^{-u} R(m-L,L)].
# Tail (m>M1): monotone majorant T(m) (Lemma tail), checked T(M1) <= C in interval arithmetic.
import sys, pickle, math
import numpy as np
from math import comb, isqrt
from sage.all import RIF, QQ, Integer
u_q = QQ(sys.argv[1]); u = RIF(u_q)
M0 = int(sys.argv[2]) if len(sys.argv) > 2 else 55
ME = int(sys.argv[3]) if len(sys.argv) > 3 else 200
M1 = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
K1 = int(sys.argv[5]) if len(sys.argv) > 5 else None
INF = 1.0 + 2.0 ** -40
def up(x):  # RIF -> float upper endpoint (exact double, RIF has 53-bit endpoints)
    return float(x.upper())
def N0(L, k):
    j = min(k, L); return comb(L + j, j) - (comb(L + j, j - 1) if j else 0)
wcache = {}
def w(L, k):  # upper bound for N0(L,k)^{-u}
    key = (L, min(k, L))
    if key not in wcache: wcache[key] = up(RIF(Integer(N0(L, k))) ** (-u))
    return wcache[key]
# ---------- region 1: exact ----------
T = pickle.load(open(f"exact_tables_{M0}.pkl", "rb"))
Zex = {}
F = np.zeros((ME + 1, ME + 1, ME + 1))
F[0, 0, 0] = 1.0
for k in range(1, M0 + 1):
    acc = {}
    tot = RIF(0)
    for (r, c, f), mult in T[k].items():
        val = RIF(Integer(mult)) * RIF(Integer(f)) ** (-u)
        acc[(r, c)] = acc.get((r, c), RIF(0)) + val
        tot += val
    for (r, c), v in acc.items(): F[k, r, c] = up(v)
    Zex[k] = tot
mstar = max(Zex, key=lambda k: Zex[k].upper())
print(f"u={u_q}: exact region m<={M0}: sup at m*={mstar}, Z in [{math.floor(float(Zex[mstar].lower()) * 1e10) / 1e10:.10f}, {math.ceil(float(Zex[mstar].upper()) * 1e10) / 1e10:.10f}]")  # endpoints rounded outward (the double conversion is within 1e-15, far below the 1e-10 grid)
# ---------- region 2: 3-D DP ----------
Pc = np.cumsum(F, axis=1) * INF          # Pc[k,a,h] >= sum_{r<=a} F[k,r,h]
Zhat = {k: up(Zex[k]) for k in Zex}
for k in range(M0 + 1, ME + 1):
    F[k, k, 1] = F[k, 1, k] = 1.0
    for r in range(2, k):
        kp = k - r; a = min(r, kp); ww = w(r, kp)
        cmax = min(r, kp + 1)
        if cmax < 2: continue
        vals = (ww * Pc[kp, a, 1:cmax]) * INF       # c = 2..cmax (row peel, r >= c)
        F[k, r, 2:cmax + 1] = vals
        if cmax == r:   # diagonal c=r filled once; conjugates c<r
            F[k, 2:r, r] = vals[:-1]
        else:
            F[k, 2:cmax + 1, r] = vals
    Pc[k] = np.cumsum(F[k], axis=0) * INF
    Zhat[k] = math.nextafter(math.fsum(F[k].ravel().tolist()), math.inf)
nz = F[F > 0]; assert nz.min() > 1e-290, "underflow risk"
# restricted sums R(k,L) = sum_{r<=L, c<=L-1} F[k,r,c] for k<=ME
RR = np.cumsum(Pc, axis=2) * INF   # RR[k,a,h] >= sum_{r<=a,c<=h}
# ---------- region 3: 1-D recursion ----------
def nup(x): return math.nextafter(x, math.inf)
for m in range(ME + 1, M1 + 1):
    s = isqrt(m - 1) + 1  # ceil(sqrt(m))
    acc = 1.0
    for L in range(s, m):
        k = m - L
        R = float(RR[k, min(L, k), min(L - 1, k)]) if k <= ME else Zhat[k]
        acc = nup(acc + nup(w(L, k) * R))
    Zhat[m] = nup(2.0 * acc)
reg2 = max(range(M0 + 1, M1 + 1), key=lambda m: Zhat[m])
print(f"  DP/recursion region {M0}<m<={M1}: max bound {Zhat[reg2]:.6f} at m={reg2}; value at M1: {Zhat[M1]:.6f}")
Cval = max(Zhat.values())
Cdec = math.ceil(Cval * 1e4) / 1e4
C = RIF(QQ(round(Cdec * 1e4)) / 10000)
# ---------- tail: T(M1) <= C ----------
if K1 is None: K1 = int(math.ceil(2 / float(u_q)))
assert K1 <= M0 and u_q * K1 >= 1 + QQ(1) / M1 and K1 + 1 <= M1 // 2
s = isqrt(M1 - 1) + 1; assert s >= 4
m = M1
S1 = sum((RIF(Integer(N0(m - k, k))) ** (-u)) * Zex[k] for k in range(1, K1 + 1))
G1 = RIF(m) / 2 * (RIF(m + 1) / RIF(Integer(comb(m, K1 + 1)))) ** u
cat = comb(2 * s, s) // (s + 1)
G2 = RIF(Integer(cat)) ** (-u) / (1 - RIF(3) ** (-u))
Tval = 2 * (1 + S1 + C * (G1 + G2))
ok = Tval.upper() <= C.lower()
print(f"  tail check (K1={K1}, M1={M1}): T(M1) <= {Tval.upper():.6f} vs C={Cdec}: {'OK' if ok else 'FAIL'}"
      f"  [S1={S1.upper():.4f}, G1={G1.upper():.2e}, G2={G2.upper():.2e}]")
print(f"RESULT u={u_q}: C_u <= {Cdec}  (exact-region sup {Zex[mstar].upper():.6f} at m={mstar}; "
      f"{'sup attained in exact region' if Cval == up(Zex[mstar]) else 'bound set by DP region'}); tail {'proved' if ok else 'NOT proved'}")
