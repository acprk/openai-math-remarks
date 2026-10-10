#!/usr/bin/env python3
"""Independent integer re-check of the finite colour budget of Section 8
(08-support.tex, lines 326-510) of OpenAI's preprint
'A high-arity counterexample to Pixton completeness in Chow'.

Only the paper's own inequalities are used; every constant (3g, 10K+2, 2K C_K,
11g, 99g^2, 30g^2, 140g^2, t_b=25(b-1)^2-1, e in [m0/2-2, m0/2], 3K+5 / 4K+20
blocks) is copied from the text.  Written without reference to pixton_budget.py.

Conditions for a parameter point (g, K, m0, nblocks), with C_K = binom(2K,3),
M_K = max(10K+3, C_K+1), S0 = g - nblocks*m0, e = 2*floor(m0/4):
 (P1) K >= 5, g >= 3                     [b>=5 for Lemma large-vertex-test, triples exist]
 (P2) g > 2K                             [unique large vertex: b > g-K >= g/2]
 (P3) 4(g-K-2) > 3g-3                    [at most 3 nonterminal expensive advancements]
 (P4) 3g-2 < (g-K)^2 - 1                 [dim S < (b-1)^2-1 with b >= g-K+1]
 (P5) loss + t_max <= d_e/2,  loss = 3g*M_K + g(10K+2) + 2K*C_K,
      t_max = 25(g-1)^2 - 1,  d_e = binom(e,3)   [eq. small-loss & large-test-threshold]
 (P6) C_K + 1 < 11g                      [ordinary costs <= 11g]
 (P7) S0 > 0 and S0*(g(K-4)+2) > 280 g^2 [contradiction in large-deficit regime, l>=K]
 (P8) nblocks >= 3K+5                    [paper: 'at most 3K+5 < 4K+20 blocks are needed']
"""
from math import comb
import sys

def ok(g, K, m0, nb):
    if K < 5 or g < 3: return False
    CK = comb(2*K, 3); MK = max(10*K + 3, CK + 1)
    if not g > 2*K: return False
    if not 4*(g - K - 2) > 3*g - 3: return False
    if not 3*g - 2 < (g - K)**2 - 1: return False
    e = 2 * (m0 // 4)
    if not (e >= 4 and 2*e >= m0 - 4 and 2*e <= m0): return False
    loss = 3*g*MK + g*(10*K + 2) + 2*K*CK
    tmax = 25*(g - 1)**2 - 1
    if not 2*(loss + tmax) <= comb(e, 3): return False
    if not CK + 1 < 11*g: return False
    S0 = g - nb*m0
    if not (S0 > 0 and S0*(g*(K - 4) + 2) > 280*g*g): return False
    if nb < 3*K + 5: return False
    return True

def icbrt_ceil(x):
    r = int(round(x ** (1/3))) if x < 10**300 else 1 << ((x.bit_length() + 2)//3)
    while r**3 < x: r += 1
    while r > 0 and (r-1)**3 >= x: r -= 1
    return r

def min_m0(g, K):
    """least m0 (necessarily a multiple of 4) with (P5)."""
    CK = comb(2*K, 3); MK = max(10*K + 3, CK + 1)
    need = 2*(3*g*MK + g*(10*K + 2) + 2*K*CK + 25*(g - 1)**2 - 1)
    # smallest even e with binom(e,3) >= need
    e = max(4, icbrt_ceil(6*need)); e += e % 2
    while e > 4 and comb(e - 2, 3) >= need: e -= 2
    while comb(e, 3) < need: e += 2
    return 2*e   # m0 = 2e gives e = 2*floor(m0/4) exactly

def h_ok(g, K, nb, j):
    """(P7) with m0 = 4j:  (g - 4 nb j)(g(K-4)+2) > 280 g^2."""
    S0 = g - 4*nb*j
    return S0 > 0 and S0*(g*(K - 4) + 2) > 280*g*g

def gmin(K, nb, j):
    """least g satisfying (P7) with m0=4j and the m0-free conditions (P2-P4,P6);
    all of these are monotone increasing in g (P7 because K-4 > 280)."""
    CK = comb(2*K, 3)
    def free(g):
        return (g > 2*K and 4*(g-K-2) > 3*g-3 and 3*g-2 < (g-K)**2-1 and CK+1 < 11*g)
    hi = 8
    while not (free(hi) and h_ok(hi, K, nb, j)): hi *= 2
    lo = 0
    while lo + 1 < hi:
        mid = (lo + hi)//2
        if free(mid) and h_ok(mid, K, nb, j): hi = mid
        else: lo = mid
    return hi

def Gmax(K, j):
    """largest g satisfying (P5) with e = 2j (need(g) is increasing in g)."""
    CK = comb(2*K, 3); MK = max(10*K + 3, CK + 1); de = comb(2*j, 3)
    def p5(g): return 2*(3*g*MK + g*(10*K+2) + 2*K*CK + 25*(g-1)**2 - 1) <= de
    if not p5(1): return 0
    lo, hi = 1, 2
    while p5(hi): hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi)//2
        if p5(mid): lo = mid
        else: hi = mid
    return lo

def least_g(K, nb):
    """Feasible set = union over j of [gmin(j), Gmax(j)] (m0 = 4j; a minimal m0
    is always a multiple of 4).  gmin is ~linear and Gmax ~ j^{3/2} in j, so
    f(j) = Gmax(j) - gmin(j) changes sign once; least feasible g = gmin(j*)
    for the least j* with f(j*) >= 0 (both monotone in j)."""
    f = lambda j: Gmax(K, j) - gmin(K, nb, j)
    hi = 2
    while f(hi) < 0: hi *= 2
    lo = 1
    while lo + 1 < hi:
        mid = (lo + hi)//2
        if f(mid) >= 0: hi = mid
        else: lo = mid
    for j in range(max(1, hi - 200), hi): assert f(j) < 0
    return gmin(K, nb, hi), 4*hi

if __name__ == "__main__":
    g0 = 10**60; K0 = 10**5; m00 = 2000 * 10**40
    print("positive control (paper point g=1e60, K=1e5, m0=2000 g^(2/3), 4K+20 blocks):",
          ok(g0, K0, m00, 4*K0 + 20))
    print("  paper's displayed ratio (4K+20)m0/g =", (4*K0 + 20)*m00 / g0)
    # paper-point sub-claims
    CK = comb(2*K0, 3)
    print("  C_K < 1.334e15:", CK < 1334*10**12, " binom(m0,3) > 1e9 g^2:", comb(m00, 3) > 10**9 * g0**2)
    for rule, nbf in (("4K+20", lambda K: 4*K + 20), ("3K+5", lambda K: 3*K + 5)):
        best = None
        Kmin = 5
        for K in range(Kmin, 4000):
            nb = nbf(K)
            # quick skip: S0 condition needs K-4 > 280 g/S0 >= 280
            if K - 4 <= 280: continue
            gK, _ = least_g(K, nb)
            if best is None or gK < best[0]: best = (gK, K)
            # once C_K/11 alone exceeds best, larger K cannot help
            if comb(2*K, 3) // 11 > best[0]: break
        gK, K = best
        _, m0 = least_g(K, nbf(K))
        assert ok(gK, K, m0, nbf(K)) and m0 == min_m0(gK, K)
        print(f"[{rule}] least g = {gK} (~{gK:.4e}) at K = {K}, m0 = {m0}, blocks = {nbf(K)}, "
              f"S0 = {gK - nbf(K)*m0}; ok(g-1)={ok(gK-1, K, min_m0(gK-1, K), nbf(K))}; n = 3*C(g,3) ~ {3*comb(gK,3):.3e}")
