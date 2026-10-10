#!/usr/bin/env python3
# R13 primary script. Env: /usr/bin/python3 + mpmath, single thread.
# Expected: expected_output/primary_main.txt (args 6.15 1.4 0.4 6), primary_control_paper_params.txt (args 50 10 1 16 2/3).
"""Rigorous (mpmath interval + exact rational binomial) check of the 176 constant.
Usage: kk176_verify.py L al be g [target]   (th = 1/2 fixed)
Every quantity below is an UPPER bound for the corresponding term in the generalised proof:
  lam  = (g+1)^(g+1)/g^g             [2C((g+1)r,r) < lam^r for all r>=1, by C(n,k)<=sqrt(n/(2pi k(n-k))) n^n/(k^k (n-k)^(n-k))]
  K    = prod_{v>=0} (1-e^{-be 2^v})^{-1}
  B    >= P(root bad)   (sup_m psi0 + sum_{j>=2} sup_{m>=g^{j-1}} psi)
  N(s) = min N with P(Bin(N,p_lo)>=s) >= target  (exact rationals, s<=S0); N(s)<=ceil(s/p_lo) beyond (median>=floor(Np))
  C1   = e^L K lam sup_T N(T+1)/(T+1);  C2 = e^L K lam sup_n N(s(n))/log2 n,  s(n)=floor(log2 C(n,2))+1
"""
import sys
from mpmath import iv, mp
from fractions import Fraction as Fr
from math import comb, floor, ceil
iv.dps = 40

def run(L, al, be, g, target=Fr(1, 2), M0=400, S0=150, quiet=False):
    L, al, be = iv.mpf(L), iv.mpf(al), iv.mpf(be)
    x = iv.exp(-L)
    c = al + be + iv.log(1 + x)
    # --- local lemma + pruning bound B
    psi0, psi = [None], [None]
    for m in range(1, M0+1):
        s = iv.mpf(0)
        for r in range(m//2 + 1, m+1):
            s += comb(m, r) * x**r
        d = iv.exp(-al*m) + iv.exp(al*m) * s
        v0 = iv.exp(be*m) * d
        psi0.append(v0)
        psi.append(v0 * iv.exp(c * (m//g) * iv.mpf(g)/(g-1)))
    # analytic tail for m>M0: psi(m) <= e^{-r1 m} + e^{-r2 m}, with
    # r1 = al-be-c/(g-1), r2 = L/2 - (al+be+log2+c/(g-1))  [x^{floor(m/2)+1} <= x^{m/2}, sum C(m,r)<=2^m]
    r1 = al - be - c/(g-1)
    r2 = L/2 - (al + be + iv.log(2) + c/(g-1))
    assert r1.a > 0 and r2.a > 0, ('tail rates', r1, r2)
    tailsup = lambda m: iv.exp(-r1*m) + iv.exp(-r2*m)
    def sup_from(m0):
        best = iv.mpf(0)
        for m in range(m0, M0+1):
            best = iv.mpf([max(best.a, psi[m].a), max(best.b, psi[m].b)])
        t = tailsup(M0+1)
        return iv.mpf([max(best.a, t.a), max(best.b, t.b)])
    B = iv.mpf(0)
    for m in range(1, M0+1):
        pass
    sup0 = max(p_.b for p_ in psi0[1:])
    t0 = tailsup(M0+1).b
    B = iv.mpf(max(sup0, t0))
    j = 2
    while g**(j-1) <= M0:
        B += sup_from(g**(j-1)).b
        j += 1
    # remaining j: sum_{j} (e^{-r1 g^{j-1}} + e^{-r2 g^{j-1}}) <= 2*(that at first j) (super-geometric)
    B += 2*tailsup(g**(j-1)).b
    Bu = B.b
    p_lo = 1 - Bu
    assert p_lo > 0
    # --- spread inflation K
    lK = iv.mpf(0); v = 0
    while True:
        t = iv.exp(-be * 2**v)
        if t.b < iv.mpf(10)**-30: break
        lK += -iv.log(1 - t); v += 1
    lK += 2*iv.exp(-be*2**v)  # tail: sum_{v'>=v} -log(1-e^{-be 2^v'}) <= 2 e^{-be 2^v} (tiny)
    K = iv.exp(lK)
    lam = iv.mpf(g+1)**(g+1) / iv.mpf(g)**g
    # --- attempts: exact rationals with rational p_lo' <= p_lo
    pq = Fr(int(mp.floor(mp.mpf(p_lo.a)*10**12)), 10**12)
    def tail_ge(N, s):
        return sum(Fr(comb(N, k)) * pq**k * (1-pq)**(N-k) for k in range(s, N+1))
    Ns = {}
    N = 1
    for s in range(1, S0+1):
        N = max(N, s)
        while tail_ge(N, s) < target:
            N += 1
        Ns[s] = N
    def Nb(s):
        if s in Ns: return Fr(Ns[s])
        assert target <= Fr(1, 2)
        return Fr(ceil(Fr(s)/pq))   # median of Bin(N,p) >= floor(Np) (Kaas-Buhrman)
    # C1: worst h=2^T, s=T+1, 1+log2 h = T+1; for s>S0 ratio <= ceil(s/pq)/s <= 1/pq + 1/s
    ar1 = max(max(Fr(Ns[s], s) for s in Ns), 1/pq + Fr(1, S0+1)) if target <= Fr(1, 2) else max(Fr(Ns[s], s) for s in Ns)
    base = iv.exp(L) * K * lam
    C1 = base * float(ar1)
    # C2: sup over n>=2 of N(s(n))/log2 n, s(n)=floor(log2 C(n,2))+1
    ar2 = 0.0; argn = None
    for n in range(2, 1 << 16):
        h = n*(n-1)//2; s = h.bit_length()  # floor(log2 h)+1
        nb = Nb(s); val = (iv.mpf(nb.numerator)/nb.denominator) / (iv.log(n)/iv.log(2))
        if val.b > ar2: ar2, argn = val.b, n
    # n >= 2^16: s <= 2 log2 n, ratio <= (2log2n/pq+1)/log2n <= 2/pq + 1/16
    ar2 = max(ar2, float(2/pq + Fr(1, 16)))
    C2 = base * ar2
    out = dict(B=float(Bu), p_lo=float(p_lo.a), K=float(K.b), lam=float(lam.b), ar1=float(ar1), ar2=float(ar2), argn=argn,
               C1=float(C1.b), C2=float(C2.b), N1to6=[Ns[s] for s in range(1, 7)])
    if not quiet: print(out)
    return out

if __name__ == '__main__':
    a = sys.argv[1:]
    tgt = Fr(a[4]) if len(a) > 4 else Fr(1, 2)
    run(a[0], a[1], a[2], int(a[3]), tgt)
