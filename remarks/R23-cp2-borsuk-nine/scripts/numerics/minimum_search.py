#!/usr/bin/env python3
"""ILLUSTRATION ONLY (not used in any proof).

Estimates t* = min_{x in CP^2} max_i |<x,p_i>|^2/(|x|^2|p_i|^2) for the nine points in P9.txt
by random sampling (seeded) followed by Nelder-Mead polishing, and prints it to 6 decimals.
The certificate (cert_interval.py, independent/cert_exact.py) proves t* > 263/500 = 0.526;
this script shows the certified threshold is within ~6e-5 of the numerical minimum.
Also prints the largest pairwise overlap |<p_i,p_j>|^2 (the nine points are far from a SIC: 1/4).
Needs numpy + scipy.  Single thread.
"""
import sys
import numpy as np
from scipy.optimize import minimize

rows = [[float.fromhex(t) for t in l.split()] for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
P = np.array([[r[0] + 1j * r[1], r[2] + 1j * r[3], r[4] + 1j * r[5]] for r in rows])
P = P / np.linalg.norm(P, axis=1, keepdims=True)
g = lambda X: (np.abs(X @ P.conj().T) ** 2).max(-1)
rng = np.random.default_rng(20261009)
X = rng.normal(size=(1000000, 3)) + 1j * rng.normal(size=(1000000, 3))
X /= np.linalg.norm(X, axis=1, keepdims=True)
gv = g(X)
print("sampled min over 1e6 random points  : %.4f" % gv.min())
def f(v):
    x = v[:3] + 1j * v[3:]
    return g((x / np.linalg.norm(x))[None])[0]
vals = [minimize(f, np.r_[X[i].real, X[i].imag], method='Nelder-Mead',
                 options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 8000}).fun for i in np.argsort(gv)[:100]]
print("best local minimum (Nelder-Mead)    : %.6f" % min(vals))
ov = max(abs(np.vdot(P[i], P[j])) ** 2 for i in range(9) for j in range(i))
print("max pairwise overlap |<p_i,p_j>|^2  : %.4f" % ov)
