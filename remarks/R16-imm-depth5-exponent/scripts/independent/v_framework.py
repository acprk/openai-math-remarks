#!/usr/bin/env python3
"""Continuum framework value (units sqrt(n) log n).  Saving of an unexpanded degree-y sqrt(n)
factor = dist(sigma*y, Z)/4 per factor (adversary picks parity), expanded factor: gain y/4, cost 1 S.
Single-degree adversary: value(sigma) = min_y max(dist(sigma y)/(4y), y/4)  (x->0 or x->all optimal,
the per-S objective is linear-fractional in x).  Two-degree mixtures: adversary mixes fraction p of
degree y1 and 1-p of y2; prover chooses expansion subset among {none, y1, y2, both}."""
import numpy as np
def dist(z): return np.abs(z - np.round(z))
def single(sig, Y):
    return np.min(np.maximum(dist(sig*Y)/(4*Y), Y/4))
Y = np.linspace(1e-3, 6, 600001)
best = max((single(s, Y), s) for s in np.linspace(0.3, 1.5, 1201))
print("single-degree adversary: max_sigma value = %.6f at sigma=%.4f; 1/(4sqrt2)=%.6f" % (*best, 1/(4*2**.5)))
for s in [0.5, 2**-.5, 1.0]:
    print("  sigma=%.4f single=%.6f  rule-lower min(s/4,1/(8s))=%.6f" % (s, single(s, Y), min(s/4, 1/(8*s))))
# mixtures at sigma=1/sqrt2: mass fractions p (of total degree) on y1, 1-p on y2
s = 2**-.5; Yg = np.linspace(0.02, 4, 400); worst = 1
for y1 in Yg:
    for y2 in Yg[::4]:
        for p in (0.25, 0.5, 0.75):
            # per unit total degree: counts p/y1, (1-p)/y2 (times sqrt n); value in units sqrt n log n
            sav = lambda e1, e2: ((0 if e1 else p/y1*dist(s*y1)/4) + (0 if e2 else (1-p)/y2*dist(s*y2)/4)
                                   + (p/4 if e1 else 0) + ((1-p)/4 if e2 else 0))
            cnt = lambda e1, e2: (p/y1 if e1 else 0) + ((1-p)/y2 if e2 else 0)
            # sqrt n -> inf: denominator 1+cnt*sqrt n ~ cnt*sqrt n unless cnt=0
            v = max(sav(0, 0), *(sav(a, b)/cnt(a, b) * 1 for a, b in ((1,0),(0,1),(1,1))))
            # (expanded options: value = sav/(cnt) only if unexpanded part counted in sqrt n units;
            #  sav is in units n, cnt in units sqrt n -> per-S value sav*sqrt n/(cnt sqrt n)... )
            worst = min(worst, v)
print("two-degree mixture adversary at sigma=1/sqrt2 (coarse): %.6f" % worst)
