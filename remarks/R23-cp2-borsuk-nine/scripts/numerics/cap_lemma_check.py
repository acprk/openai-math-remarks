#!/usr/bin/env python3
"""ILLUSTRATION ONLY: random test of the cap lemma used in the note.

Lemma: unit x, y, p in C^3 with |<x,p>|^2 >= T and |<y,p>|^2 >= T (T > 1/2) satisfy
|<x,y>| >= 2T - 1, hence ||xx* - yy*||_F^2 = 2 - 2|<x,y>|^2 <= 2 - 2(2T-1)^2.
We sample pairs near the boundary of the cap and report the smallest |<x,y>| - (2T-1) (should be >= 0
up to rounding) (the slack is small because sampled points sit near the cap boundary).
"""
import numpy as np
rng = np.random.default_rng(7)
T = 21 / 40
p = np.array([1, 0, 0], dtype=complex)
worst = np.inf
for _ in range(20):
    n = 200000
    def cap_pts():
        a2 = rng.uniform(T, T + 1e-3, n)  # |<x,p>|^2 just above T (near the cap boundary)
        w = rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))
        w /= np.linalg.norm(w, axis=1, keepdims=True)
        ph = np.exp(2j * np.pi * rng.uniform(size=n))
        return np.c_[np.sqrt(a2) * ph, np.sqrt(1 - a2)[:, None] * w]
    x, y = cap_pts(), cap_pts()
    s = np.abs((x * y.conj()).sum(1)) - (2 * T - 1)
    worst = min(worst, s.min())
print("T = 21/40, 2T-1 = %.4f; min over 4e6 sampled pairs of |<x,y>| - (2T-1) = %.2e (>= 0 expected)"
      % (2 * T - 1, worst))
print("diameter bound sqrt(2 - 2(2T-1)^2) = sqrt(399/200) = %.6f  (vs sqrt 2 = %.6f)"
      % (np.sqrt(399 / 200), np.sqrt(2)))
