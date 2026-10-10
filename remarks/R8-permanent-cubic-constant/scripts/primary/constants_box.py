"""Constants for 108 with the BOX lemma (length Z_X >= (q+1)^2, char-free, any 1<=q, r=k-2q>=2).
Chain (paper): bdc(per_m) >= (r-1)(d-1)/(4e), m >= 11k-2q, 2<=d, d-1<Delta, d<k^2, Delta>=1.
Models (Delta = lower bound on codim of gradient locus):
  paper : C(q+2,2)-k  (q<=c needed)           L2  : T'(q,r)-k (Hilbert-function max; q<=c not needed)
  BOX   : (q+1)^2-k   [PROVED here]            BOXB: 2(q+1)^2-k [Lemma B + BOX, conditional]
  A     : (2q+2)^2    [Lemma A, conditional; = rank-family ceiling]
Searches all 1<=q<=(k-2)/2.  Prints c^{-1}/e = 4 m^3/((r-1)(d-1)) minimized over (k,q,d)."""
from math import comb, e
import sys
def Tprime(q, r):
    ns = [min(s, 2*q - s) + 1 for s in range(2*q + 1)]
    pre = [0]
    for v in ns: pre.append(pre[-1] + v)
    best = 0
    for E in range(0, 2*q + r):
        lo = max(0, E - r + 1); hi = min(2*q, E)
        if lo <= hi: best = max(best, pre[hi+1] - pre[lo])
    return best
def Delta(model, k, q):
    r = k - 2*q
    if model == 'paper' and q > r - 1: return -1
    return {'paper': comb(q+2, 2) - k, 'L2': Tprime(q, r) - k, 'BOX': (q+1)**2 - k,
            'BOXB': 2*(q+1)**2 - k, 'A': (2*q+2)**2}[model]
def best(m, model):
    bv, bp = 0, None
    for k in range(4, m//9 + 2):
        for q in range(1, (k-2)//2 + 1):
            if 11*k - 2*q > m: continue
            r = k - 2*q
            D = Delta(model, k, q)
            if D < 2: continue
            d = min(D, k*k - 1)       # d-1 < Delta  <=>  d <= Delta
            v = (r - 1)*(d - 1)
            if v > bv: bv, bp = v, (k, q, r, d)
    return 4*m**3/bv, bp
if __name__ == '__main__':
    ms = [int(a) for a in sys.argv[1:]] or [1408, 5000, 20000]
    print("paper stated: c^-1 = 5529600 e ;  asymptotic closed forms: BOX 4*29791=119164, BOXB 59582, A 29791")
    for m in ms:
        res = {mod: best(m, mod) for mod in ['paper', 'L2', 'BOX', 'BOXB', 'A']}
        print(f"m={m}: " + " | ".join(f"{mod}: {v:.0f}e {p}" for mod, (v, p) in res.items()), flush=True)
        print(f"   ratios: paper_opt/BOX={res['paper'][0]/res['BOX'][0]:.3f}  L2/BOX={res['L2'][0]/res['BOX'][0]:.3f}  stated/BOX={5529600/res['BOX'][0]:.1f}  BOX/BOXB={res['BOX'][0]/res['BOXB'][0]:.3f}  BOX/A={res['BOX'][0]/res['A'][0]:.3f}")
