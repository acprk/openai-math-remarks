# R3 exhaustive minimality over all threshold-Kneser (N,k,s), N<=120, with 3k-3(s-1)>N. Run: /usr/bin/python3 independent/minN.py (~75 s)
# (N>120 cannot help: k-s >= floor(N/3) so t >= C(N,floor(N/3)) and m > 96 t >= 96*C(121,40) > 10^34.)  Expected output: expected_output/minN.txt
from math import comb
from fractions import Fraction as Fr
for bound,name in ((Fr(1,96),'1/96'),(Fr(1,48),'1/48')):
    best=None
    for N in range(2,121):
        for k in range(1,N):
            for s in range(1,k+1):
                if 3*k-3*(s-1)<=N: continue
                t=sum(comb(N,j) for j in range(k-s+1)); m=comb(N,k)
                if Fr(t,m)<bound:
                    if best is None or comb(N,k)<best[0]: best=(comb(N,k),N,k,s,float(Fr(t,m)))
    print(name,best)
