# Numerical sanity check of A_min(K) = sqrt(2) (construction.tex Lemma 2.1; proved there, this is only a check).
# Shadow area in direction u = area of the convex hull of the 4 projected vertices; minimized over a
# Fibonacci grid of 200k directions followed by Nelder-Mead refinement.  Also reports the shadow areas
# in the 640 cylinder directions (each >= A_min, so the directionwise ratio is also < 1/2).
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
H = np.sqrt(2); V = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, H], [0, -1, H]])
def shadow(u):
    u = np.asarray(u, float); u = u / np.linalg.norm(u)
    a = np.cross(u, [1, 0, 0] if abs(u[0]) < 0.9 else [0, 1, 0]); a /= np.linalg.norm(a); b = np.cross(u, a)
    P = np.c_[V @ a, V @ b]
    return ConvexHull(P).volume
N = 200000; i = np.arange(N) + 0.5; phi = np.arccos(1 - 2 * i / N); th = np.pi * (1 + 5 ** 0.5) * i
U = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
vals = np.array([shadow(u) for u in U]); k = vals.argmin()
r = minimize(shadow, U[k], method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14})
print(f"grid min {vals.min():.10f}  refined min {r.fun:.12f}  sqrt(2) = {H:.12f}")
eps = 0.4; n = 320; q = np.linspace(-1, 1, n + 1); al = (1 + q[:-1] * q[1:]) / 4; be = (q[:-1] + q[1:]) / 4
sh = [shadow([1, eps * a, H * eps * b]) for a, b in zip(al, be)] + [shadow([-eps * a, 1, H * eps * b]) for a, b in zip(al, be)]
print(f"shadow areas in the 640 cylinder directions: min {min(sh):.6f}, max {max(sh):.6f} (all >= sqrt 2: {min(sh) >= H - 1e-12})")
T = np.array([float(eval(l.split()[1])) for l in open('params_T.txt')])
B = H * T ** 2 / n / np.sqrt(1 + eps ** 2 * (al ** 2 + 2 * be ** 2))   # |B_j|, same for both families
Rdir = np.sum(B / np.array(sh[:n])) + np.sum(B / np.array(sh[n:]))
print(f"directionwise relative cost R = sum |B_i| / |pi_(u_i^perp) K| = {Rdir:.6f} (float; <= S/A_min)")
