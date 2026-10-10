"""Search the smallest k for which the block family violates Lemma B:
famdim(k,q) > k^2 + k - 2*C(2q+2,2), using the rigorous bound len Z_X >= C(2q+2,2) for every X in the family."""
from math import comb
from family_exact import famdim
for k in range(8, 140):
    best=None
    for q in range(1, (k-2)//2+1):
        d, info = famdim(k, q)
        marg = d - (k*k + k - 2*comb(2*q+2, 2))
        if best is None or marg > best[0]: best = (marg, q, info, d)
    if best[0] > 0:
        print("first violation: k=%d q=%d family=%s dim=%d  bound=%d  margin=%d" % (k, best[1], best[2], best[3], k*k+k-2*comb(2*best[1]+2,2), best[0])); break
    if k % 10 == 0: print(k, best)
