# Per-sector cutoffs for the family-100 angular cover: T^1_j = 1/2 + eps^2 (d_j + e1_j), T^2_k = 1/2 + eps^2 (d_k + e2_k).
# Pair (j,k) is safe iff no point of K cap W_j cap W_k has A_j > e1_j and B_k > e2_k, where
# A_j = (t0-1/2)/eps^2 - d_j, B_k = (s0-1/2)/eps^2 - d_k.  Coordinate descent:
#   e1_j <- max_k max{A_j : p in R_jk, B_k >= e2_k}   (and symmetrically for e2).
import numpy as np, sys
from scipy.optimize import linprog

def setup(eps, n):
    qs = np.linspace(-1, 1, n + 1); a = (1 + qs[:-1] * qs[1:]) / 4; b = (qs[:-1] + qs[1:]) / 4
    d = lambda q: q * q * (1 + q * q) / 16; dj = np.maximum(d(qs[:-1]), d(qs[1:]))
    return qs, a, b, dj

def region(eps, qs, a, b, j, k):
    A = [[0, 0, -1], [0, 0, 1], [1, 0, 1], [-1, 0, 1], [0, 1, -1], [0, -1, -1]]; B = [0, 1, 1, 1, 0, 0]
    for qq, sg in ((qs[j], 1), (qs[j + 1], -1)):
        A.append([sg * (-qq * eps * b[j] + eps * a[j]), -sg, sg * qq]); B.append(0)
    for qq, sg in ((qs[k], 1), (qs[k + 1], -1)):
        A.append([-sg, sg * (qq * eps * b[k] - eps * a[k]), -sg * qq]); B.append(-sg * qq)
    return A, B

def maxA_given_B(eps, qs, a, b, dj, j, k, e2k, which):
    # which=1: maximize A_j s.t. B_k >= e2k ; which=2: maximize B_k s.t. A_j >= e1j (e2k plays e1j)
    E = eps * eps
    A, B = region(eps, qs, a, b, j, k)
    # A_j*E = t - eps b_j x - 1/2 - E d_j ; B_k*E = 1 - t + eps b_k y - 1/2 - E d_k
    if which == 1:
        A.append([0, -eps * b[k], 1]); B.append(0.5 - E * dj[k] - E * e2k)   # B_k >= e2k
        c = [eps * b[j], 0, -1]  # minimize -(t - eps b_j x)
        r = linprog(c, A_ub=A, b_ub=B, bounds=[(None, None)] * 3, method='highs')
        return -np.inf if r.status == 2 else ((-r.fun) - 0.5 - E * dj[j]) / E
    else:
        A.append([eps * b[j], 0, -1]); B.append(-0.5 - E * dj[j] - E * e2k)   # A_j >= e1j
        c = [0, -eps * b[k], 1]  # minimize -(1 - t + eps b_k y) + 1
        r = linprog(c, A_ub=A, b_ub=B, bounds=[(None, None)] * 3, method='highs')
        return -np.inf if r.status == 2 else ((1 - r.fun) - 0.5 - E * dj[k]) / E

def area(eps, n, e1, e2, qs, a, b, dj):
    sp = np.sqrt(1 + eps ** 2 * (a ** 2 + 2 * b ** 2))
    T1 = 0.5 + eps ** 2 * (dj + e1); T2 = 0.5 + eps ** 2 * (dj + e2)
    return 0.5 * np.sum((2 / n) * (T1 ** 2 + T2 ** 2) / sp)

if __name__ == "__main__":
    eps = float(sys.argv[1]); n = int(sys.argv[2]); passes = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    qs, a, b, dj = setup(eps, n)
    e = np.full(n, float(sys.argv[4]) if len(sys.argv) > 4 else 0.05)
    for p in range(passes):
        f = np.array([max(max(maxA_given_B(eps, qs, a, b, dj, j, k, e[k], 1) for k in range(n)), -0.2) for j in range(n)]) + 1e-7
        e = 0.5 * (e + f)  # damped symmetric fixed-point iteration (e1=e2=e by the x<->y, t<->1-t symmetry)
        # feasibility check of current e: worst pair slack
        worst = max(maxA_given_B(eps, qs, a, b, dj, j, k, e[k], 1) - e[j] for j in range(n) for k in range(n))
        S = area(eps, n, e, e, qs, a, b, dj)
        print(f"pass {p}: area/Amin={S:.9f} (1/2-S={0.5-S:.3e}) e range [{e.min():.5f},{e.max():.5f}] median {np.median(e):.5f} worst-slack {worst:.2e} {'FEASIBLE' if worst<0 else 'infeasible'}", flush=True)
    e1 = e2 = e
    np.save(f"persector_{eps}_{n}.npy", np.vstack([e1, e2]))
    print("e1 profile (every n/16):", np.round(e1[::max(1, n // 16)], 5))
