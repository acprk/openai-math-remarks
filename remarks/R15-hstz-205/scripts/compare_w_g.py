#!/usr/bin/env python3
"""compare_w_g.py -- compare the 205 lower bound w_m (from w205.sage, file w_m{m}.txt)
with the exact Kronecker coefficients g (from crt_exact.py inputs, i.e. the kron_stair
.full files for several primes).  Reports:
  * w <= g for every lambda (consistency of both computations),
  * whether the 205 bound alone certifies HSTZ at this m:
        for every lambda with a < 0:  max(w(lambda), w(lambda')) >= |a| + 1,
        and max(w(rho)) >= 2  (then g(rho,rho,rho) = 1 mod 4 forces >= 5),
  * the lambdas where the 205 bound is one-sided (w(lambda) < |a|+1 <= w(lambda')).
Usage: python3 compare_w_g.py m w_m{m}.txt ks_m{m}_P*.full
"""
import sys
from math import prod
m = int(sys.argv[1]); wf = sys.argv[2]; files = sys.argv[3:]
def conj(l):
    return tuple(sum(1 for x in l if x > j) for j in range(l[0])) if l else ()
W = {}
for line in open(wf):
    p, v = line.split(); W[tuple(map(int, p.split(',')))] = int(v)
primes = [int(f.split('_P')[1].split('.')[0]) for f in files]
M = prod(primes); res = {}
for p, f in zip(primes, files):
    for line in open(f):
        q, g, a = line.split(); lam = tuple(map(int, q.split(',')))
        res.setdefault(lam, [0, int(a)])
        Mi = M // p; res[lam][0] += int(g) * Mi * pow(Mi, -1, p)
G = {l: (v[0] % M, v[1]) for l, v in res.items()}
assert set(G) == set(W)
viol = [l for l in G if W[l] > G[l][0]]
print(f"m={m}: #lambda={len(G)}  #(w > g) = {len(viol)}")
fails = []; onesided = []
for l, (g, a) in G.items():
    if a < 0:
        need = -a + 1; wl, wc = W[l], W[conj(l)]
        if max(wl, wc) < need: fails.append((l, wl, wc, g, a))
        elif wl < need: onesided.append(l)
rho = tuple(range(m, 0, -1))
print("205 bound certifies all lambda with a<0 (two-sided max):", not fails, " failures:", fails[:5])
print("w(rho) =", W[rho], " (>=2 needed)")
print("#lambda with a<0 where w(lambda) alone is too small but w(lambda') suffices:", len(onesided))
