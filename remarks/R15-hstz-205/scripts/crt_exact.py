#!/usr/bin/env python3
"""crt_exact.py -- combine kron_stair runs for several primes into EXACT
Kronecker coefficients g(rho_m, rho_m, lambda) (valid because
0 <= g <= dim S^rho_m < product of the primes, checked below), then report

  * the HSTZ criterion  g > -eps*a  for every lambda (eps = +1 for m != 2 mod 4),
  * g(rho,rho,rho) (must be >= 5, and = 1 mod 4 by Prop. 1),
  * the minimal margin  g - |a|  over lambda with eps*a < 0, and where it is attained,
  * exact g on all hooks.

Usage: python3 crt_exact.py m file_p1.full file_p2.full ...
(the .full files are written by kron_stair for N <= 45).
"""
import sys
from math import factorial, prod

def dim(lam):
    lt = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    h = 1
    for i, r in enumerate(lam):
        for j in range(r):
            h *= (r - j) + (lt[j] - i) - 1
    return factorial(sum(lam)) // h

def read(fn):
    out = {}
    for line in open(fn):
        part, g, a = line.split()
        out[tuple(map(int, part.split(',')))] = (int(g), int(a))
    return out

m = int(sys.argv[1]); N = m * (m + 1) // 2
files = sys.argv[2:]
primes = [int(f.split('_P')[1].split('.')[0]) for f in files]
data = [read(f) for f in files]
M = prod(primes)
rho = tuple(range(m, 0, -1))
drho = dim(rho)
assert drho < M, "product of primes too small"
d = (m + 1) // 2
eps = 1 if ((N - d) // 2) % 2 == 0 else -1
G = {}
for lam in data[0]:
    x = 0
    for p, D in zip(primes, data):
        r, a = D[lam]
        Mi = M // p
        x += r * Mi * pow(Mi, -1, p)
    G[lam] = (x % M, data[0][lam][1])
# exact sanity: sum g dim = dim(rho)^2 over the integers
s = sum(g * dim(l) for l, (g, a) in G.items())
print(f"m={m} N={N} eps={eps} primes={primes}")
print("exact identity sum_lambda g*dim(lambda) == dim(rho)^2 :", s == drho * drho)
print("g(rho,rho,rho) =", G[rho][0], " (mod 4 =", G[rho][0] % 4, ")")
bad = [(l, g, a) for l, (g, a) in G.items() if g + eps * a <= 0]
print("lambdas violating g + eps*a > 0 :", len(bad), bad[:5])
neg = [(g - abs(a), l, g, a) for l, (g, a) in G.items() if eps * a < 0]
neg.sort()
print("#lambda with eps*a<0 :", len(neg))
for t in neg[:5]:
    print("  smallest margins g-|a|:", t[0], "lambda", t[1], "g", t[2], "a", t[3])
zero = [(g, l) for l, (g, a) in G.items() if a == 0]
zero.sort()
print("min g over a=0 lambdas:", zero[0])
print("hooks (k = leg length): g(rho,rho,(N-k,1^k)), a")
for k in range(N):
    lam = (N - k,) + (1,) * k
    g, a = G[lam]
    print(f"  k={k} g={g} a={a}")
