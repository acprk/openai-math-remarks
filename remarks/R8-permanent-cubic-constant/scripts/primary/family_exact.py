"""Exact (integer) best block family inside Gamma for given (k,q): coordinate partition into m blocks,
block-upper-triangular X, one 'big' block of size n1 with rank <= r-1-m, the other m-1 blocks rank 0
(sizes as equal as possible).  All almost-principal (r-1)-minors vanish => family subset Gamma_D
(every theta, every characteristic).  Returns max dim.  (Optimal shape by the Lagrange computation:
ranked blocks equal size, best with one ranked block; unranked blocks as small/equal as possible.)
Also: rigorous facts for every X in the family: adj(M) in m^{2q+1}, F in m^{2q+1+m}  =>  len Z_X >= C(2q+2,2).
Lemma B (dim{len Z>=l} <= k^2+k-2l) is violated as soon as famdim > k^2+k-2*C(2q+2,2)."""
import sys
from math import comb
def famdim(k, q):
    r = k - 2*q; best = (k*k - (k - r + 2)**2, (1, k, r-2))   # m=1: rank family
    for m in range(2, r):
        rho = r - 1 - m
        if rho < 0: break
        for n1 in range(max(rho,1), k - (m-1) + 1):
            rest = k - n1; t = m - 1
            a, b = divmod(rest, t)                     # b blocks of size a+1, t-b of size a
            if a == 0: continue
            sumsq = n1*n1 + b*(a+1)**2 + (t-b)*a*a
            d = (k*k - sumsq)//2 + rho*(2*n1 - rho)
            if d > best[0]: best = (d, (m, n1, rho))
    return best
if __name__ == '__main__':
    print("k q r | rank-fam codim (2q+2)^2 | best family codim [m,n1,rho] | ratio | LemmaB margin famdim-(k^2+k-2C(2q+2,2))")
    for k in [int(a) for a in sys.argv[1:]] or [22, 31, 61, 100, 151, 211, 301, 485]:
        q = (k-1)//3; r = k-2*q
        d, info = famdim(k, q)
        cod = k*k - d
        lb = d - (k*k + k - 2*comb(2*q+2, 2))
        print(f"{k} {q} {r} | {(2*q+2)**2} | {cod} {info} | {cod/(2*q+2)**2:.4f} | {lb}")
