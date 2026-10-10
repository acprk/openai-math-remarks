# R3 parameter search (heuristic: largest admissible s for each (N,k)). Run: /usr/bin/python3 main/tk_opt.py (<1 s)
# Expected output: expected_output/tk_opt.txt.  Exhaustive minimality is independent/minN.py.
# T-OPT for 197 odd-char: smallest threshold-Kneser instance (N,k,s):
# vertices k-subsets of [N], edge iff |I∩J|<s; triangle-free iff 3k-3(s-1)>N;
# fit rank t = sum_{j<=k-s} C(N,j) (multilinearized degree-(k-s) polynomial); need t/m < c.
from math import comb, log10
for c,name in [(1/100,'paper 1/100'),(1/96,'1/96'),(1/48,'1/48')]:
    best=None
    for N in range(3,1300):
        for k in range(1,N+1):
            m=comb(N,k)
            # triangle-free: 3k-3(s-1)>N  -> s < (3k-N)/3+1 ; largest s maximizes nothing? smaller degree k-s wants s large
            smax=max(s for s in range(0,k+1) if 3*k-3*(s-1)>N) if 3*k+3>N else None
            if smax is None or smax<1: continue
            s=smax; t=sum(comb(N,j) for j in range(0,k-s+1))
            if t*1.0/m < c and t<m:
                Delta=sum(comb(k,j)*comb(N-k,k-j) for j in range(0,s))
                best=(N,k,s,m,t,Delta); break
        if best: break
    N,k,s,m,t,D=best
    print(f"{name}: N={N} k={k} s={s}  m=C(N,k)≈10^{log10(m):.2f}  t/m={t/m:.5f}  maxdeg≈10^{log10(D):.2f}")
print("paper: N=1200 k=600 s=200 m≈10^%.1f"%log10(comb(1200,600)))
