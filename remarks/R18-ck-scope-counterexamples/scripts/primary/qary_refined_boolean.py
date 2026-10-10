#!/usr/bin/env python3
"""Does the output-bias refinement I(f(X);Y) <= psi(lam * psi^{-1}(H(m)))  (OpenAI 119 Cor.,
true on the uniform binary cube) survive on [q]^n with q-ary symmetric noise?  Test single-coordinate
Boolean f with k ones out of q (n=1 suffices: any f of one coordinate) and the exhaustive optima."""
import mpmath as mp
mp.mp.dps = 30
def psi(r):
    p=(1+r)/2; q=1-p
    return mp.log(2) + (p*mp.log(p) if p>0 else 0) + (q*mp.log(q) if q>0 else 0)
def psi_inv(I):
    if I >= mp.log(2) - mp.mpf("1e-28"): return mp.mpf(1)
    return mp.findroot(lambda r: psi(r)-I, (mp.mpf('1e-25'), 1-mp.mpf('1e-25')), solver='bisect')
worst = -1
for q in range(3, 10):
    for k in range(1, q):
        w = mp.mpf(k)/q; m = 2*w-1
        Hm = psi(m)*0  # placeholder
        H = -( (w*mp.log(w)) + ((1-w)*mp.log(1-w)) )   # H(f) nats
        for lam in [mp.mpf(x)/10 for x in range(1,10)]:
            I1 = w*psi(lam*1+(1-lam)*m) + (1-w)*psi(-lam+(1-lam)*m) - psi(m)
            rhs = psi(lam*psi_inv(H))
            g = I1 - rhs
            worst = max(worst, g)
            if g > 0:
                print(f"VIOLATION q={q} k={k} lam={float(lam)}: I={float(I1):.6f} refined bound={float(rhs):.6f} gap={float(g):+.3e}")
print("max gap (refined Boolean, single-coordinate):", float(worst))
