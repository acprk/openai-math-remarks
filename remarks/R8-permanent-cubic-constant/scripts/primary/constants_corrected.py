"""108 constants after the 10/09 Lemma-A/B audit.
bdc(per_m) >= (r-1)(d-1)/(4e), m>=11k-2q, 2<=d<=Delta, d<k^2;  printed number X means c^{-1} = X*e (=4m^3/((r-1)(d-1))).
 BOX      : Delta=(q+1)^2-k                                         PROVED (box lemma, any char, any 1<=q<=(k-2)/2)
 L2       : Delta=T'(q,r)-k                                         (previous rigorous best)
 A_old    : Delta=(2q+2)^2                                          Lemma A -- FALSE (block families), kept for comparison
 A_ceil   : Delta=min((2q+2)^2, famcodim)                           CEILING for every route bounding codim of the char-2 locus
 PER_ceil : Delta=min((2q+2)k, DAGcodim)                            ceiling for the complex per-locus route
famcodim(k,q): smallest codim of a coordinate block-upper-triangular family (one ranked block, rank r-1-m) inside Gamma_D."""
import sys
from math import comb
def Tprime(q, r):
    ns = [min(s, 2*q - s) + 1 for s in range(2*q + 1)]
    pre = [0]
    for v in ns: pre.append(pre[-1] + v)
    return max(pre[min(2*q,E)+1]-pre[max(0,E-r+1)] for E in range(0, 2*q + r) if max(0,E-r+1) <= min(2*q,E))
def famdim(k, q):
    r = k - 2*q; best = k*k - (k - r + 2)**2
    for m in range(2, r):
        rho = r - 1 - m; t = m - 1
        c = (k + 2*rho*t)/m
        for n1 in {int(c)-1, int(c), int(c)+1, int(c)+2}:
            if n1 < max(rho,1) or k - n1 < t: continue
            a, b = divmod(k - n1, t)
            sumsq = n1*n1 + b*(a+1)**2 + (t-b)*a*a
            d = (k*k - sumsq)//2 + rho*(2*n1 - rho)
            best = max(best, d)
    return best
def dagcodim(k, q):
    r = k - 2*q; m = r - 1
    a, b = divmod(k, m); sumsq = b*(a+1)**2 + (m-b)*a*a
    return k*k - (k*k - sumsq)//2
def Delta(model, k, q):
    r = k - 2*q
    if model == 'BOX': return (q+1)**2 - k
    if model == 'L2': return Tprime(q, r) - k
    if model == 'A_old': return (2*q+2)**2
    if model == 'A_ceil': return min((2*q+2)**2, k*k - famdim(k, q))
    if model == 'PER_ceil': return min((2*q+2)*k, dagcodim(k, q))
def best(m, model, qlo=0.2, qhi=0.42):
    bv, bp = 0, None
    for k in range(int(m/11.2), m//10 + 2):
        for q in range(max(1,int(qlo*k)), min((k-2)//2, int(qhi*k)) + 1):
            if 11*k - 2*q > m: continue
            r = k - 2*q
            D = Delta(model, k, q)
            if D < 2: continue
            d = min(D, k*k - 1)
            v = (r - 1)*(d - 1)
            if v > bv: bv, bp = v, (k, q, r, d)
    return 4*m**3/bv, bp
if __name__ == '__main__':
    ms = [int(a) for a in sys.argv[1:]] or [1408, 5000]
    for m in ms:
        res = {mod: best(m, mod) for mod in ['L2', 'BOX', 'A_old', 'A_ceil', 'PER_ceil']}
        print(f"m={m}: " + " | ".join(f"{mod}: {v:.0f}e {p}" for mod, (v, p) in res.items()), flush=True)
        print(f"   stated/BOX={5529600/res['BOX'][0]:.1f}x  L2/BOX={res['L2'][0]/res['BOX'][0]:.3f}x  BOX/A_ceil={res['BOX'][0]/res['A_ceil'][0]:.2f}x  A_ceil/A_old={res['A_ceil'][0]/res['A_old'][0]:.3f}  stated/A_ceil={5529600/res['A_ceil'][0]:.0f}x  stated/PER_ceil={5529600/res['PER_ceil'][0]:.0f}x", flush=True)
