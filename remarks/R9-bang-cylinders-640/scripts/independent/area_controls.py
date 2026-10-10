# Float area of comparison covers with the same 640 axes (eps=2/5, n=320), formula (area) of the note:
#   per-sector cutoffs (certified cover), uniform margin max_j E_j (the preprint's style of cutoff), no margin (not certified).
# Usage: sage -python area_controls.py params_T.txt
import sys, numpy as np
from fractions import Fraction as F
eps = 0.4; n = 320; q = np.linspace(-1, 1, n + 1); al = (1 + q[:-1] * q[1:]) / 4; be = (q[:-1] + q[1:]) / 4
d = lambda u: u * u * (1 + u * u) / 16; dj = np.maximum(d(q[:-1]), d(q[1:]))
T = np.array([float(F(l.split()[1])) for l in open(sys.argv[1])]); E = (T - 0.5) / eps ** 2 - dj
sp = np.sqrt(1 + eps ** 2 * (al ** 2 + 2 * be ** 2)); A = lambda T: np.sum(2 * T ** 2 / n / sp)
print(f"margins E_j in [{E.min():.4e}, {E.max():.4e}]; argmax at sectors {np.argsort(-E)[:2]}, argmin {np.argsort(E)[:2]}")
print(f"per-sector cutoffs (certified): {A(T):.10f}")
print(f"uniform margin max E_j:          {A(0.5 + eps**2 * (dj + E.max())):.10f}")
print(f"no margin (coverage unknown):    {A(0.5 + eps**2 * dj):.10f}")
