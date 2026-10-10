"""Independent recomputation of the cubic constant for family 108 (does not import constants_box.py).
Pipeline of the manuscript (labels of sections/*.tex):
  prop:selector     : G_D is an affine projection of per_m for every m >= 11k-2q, 0<=q<=floor(k/2)
  prop:gradient-codimension + lem:characteristic-transfer : codim Gamma_D >= Delta  (Delta>=1)
  prop:smooth-slice : smooth d-dim section if 2<=d<k^2 and d-1<Delta   (integers: d<=Delta)
  thm:smooth-border : bdc >= (r-1)(d-1)/(4e), r=k-2q>=2
Hence bdc(per_m) >= P/(4e) with P=max (r-1)(d-1);  c^{-1}=4 m^3/P  (in units of e).
Models for Delta:
  paper : C(q+2,2)-k, needs 3q<=k-1 (q<=c)         [OpenAI's Prop. adjugate-length]
  BOX   : (q+1)^2-k, needs 2q<=k-2 (r>=2)          [Lemma BOX]
"""
import sys
from fractions import Fraction
from math import comb

def P(model, k, q):
    r = k - 2*q
    if r < 2: return 0
    if model == 'paper':
        if 3*q > k - 1: return 0
        D = comb(q+2, 2) - k
    else:
        D = (q+1)**2 - k
    d = min(D, k*k - 1)
    if d < 2: return 0
    return (r-1)*(d-1)

def best_exhaustive(m, model):
    bv, arg = 0, None
    k = 1
    while 11*k - 2*(k//2) <= m:          # smallest possible size for this k
        for q in range(0, k//2 + 1):
            if 11*k - 2*q > m: continue
            v = P(model, k, q)
            if v > bv: bv, arg = v, (k, q)
        k += 1
    return bv, arg

if __name__ == '__main__':
    ms = [int(a) for a in sys.argv[1:]] or [1408, 5000, 20000]
    print("stated: c^-1 = 5529600 e   (thm:main, m0=1408)")
    print("asymptotic optimum: paper-tuned 4*2*31^3 = %d ;  BOX 4*29700 = %d" % (4*2*31**3, 4*29700))
    for m in ms:
        out = []
        for mod in ['paper', 'BOX']:
            v, a = best_exhaustive(m, mod)
            k, q = a
            r = k-2*q
            out.append((mod, Fraction(4*m**3, v), a))
            print(f"m={m} {mod:5s}: c^-1 = {float(Fraction(4*m**3, v)):.1f} e   at k={k} q={q} r={r} Delta={(comb(q+2,2) if mod=='paper' else (q+1)**2)-k}  (r-1)(d-1)={v}")
        pa, bo = out[0][1], out[1][1]
        print(f"   stated/paper-tuned = {5529600/float(pa):.3f}   paper-tuned/BOX = {float(pa/bo):.4f}   stated/BOX = {5529600/float(bo):.3f}")
