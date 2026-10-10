#!/usr/bin/env python3
# R14 primary script (copied from reviews/R7-skill/rescan/topt-lean/). Env: sage -python or python3+numpy.
# Expected: expected_output/primary_spectral.txt
"""196 Kaplansky zero-divisor: squared-word decay needs rho(M)<1, M_{t,u}=w(t,u)^2.
Builds PG(2,q) for q=2^a, the alphabet T = points + 3 extras, involution pairing,
and computes rho(M) (numeric) for paper's bound and for several pairings.
Positive control: q=128 must give the paper's row bounds 500731261911/504183783481 and 820350863/1363395845."""
import numpy as np, itertools, random, sys
from fractions import Fraction as Fr

def gf_mul(a, b, a_bits, poly):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> a_bits: a ^= poly
    return r

POLYS = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011}

def pg2(abits):
    q = 1 << abits; poly = POLYS[abits]
    mul = lambda a, b: gf_mul(a, b, abits, poly)
    pts = []
    for x in range(q):
        for y in range(q): pts.append((x, y, 1))
    for x in range(q): pts.append((x, 1, 0))
    pts.append((1, 0, 0))
    idx = {p: i for i, p in enumerate(pts)}
    lines = []
    for (a, b, c) in pts:  # line a x + b y + c z = 0
        lines.append(frozenset(i for i, (x, y, z) in enumerate(pts) if (mul(a, x) ^ mul(b, y) ^ mul(c, z)) == 0))
    return q, pts, lines

def build_W(q, lines):
    v = q*q+q+1; p = Fr(q+1, v)
    n = v+3  # letters 0..v-1 ordinary, v,v+1,v+2 extra
    W = [[Fr(0)]*n for _ in range(n)]
    for a in range(n):
        for u in range(n):
            if u == a: continue
            ao, uo = a < v, u < v
            if ao and uo: w = Fr(1, q+1)
            elif ao != uo: w = p
            else: w = Fr(1)
            W[a][u] = w*w
    return W, v, p

def rho_for_pairing(W, bar):
    n = len(W)
    M = np.array([[float(W[bar[t]][u]) if u != bar[t] else 0.0 for u in range(n)] for t in range(n)])
    return max(abs(np.linalg.eigvals(M))), M

def pairing(v, extras_to, rest_order):
    bar = {}
    for k, o in enumerate(extras_to):
        bar[v+k] = o; bar[o] = v+k
    rest = [i for i in rest_order if i not in extras_to]
    for i in range(0, len(rest), 2):
        bar[rest[i]] = rest[i+1]; bar[rest[i+1]] = rest[i]
    return [bar[i] for i in range(v+3)]

if __name__ == '__main__':
    # positive control on paper's row bounds (q=128) -- symbolic only
    q = 128; v = q*q+q+1; p = Fr(q+1, v)
    r1 = Fr(q, q+1) + 3*p*p + Fr(12, (q+1)**2)
    r2 = (2 + v*p*p + 12*p*p)/5
    print('control q=128 row bounds', r1, r2, r1 == Fr(500731261911, 504183783481), r2 == Fr(820350863, 1363395845))
    for abits in [1, 2, 3, 4, 5]:
        q, pts, lines = pg2(abits)
        W, v, p = build_W(q, lines)
        best = None
        rng = random.Random(1)
        # pairings: extras paired to collinear / non-collinear triples, random rest
        cands = []
        L0 = sorted(lines[0])
        noncol = None
        for tri in itertools.combinations(range(v), 3):
            if not any(set(tri) <= l for l in lines): noncol = tri; break
        for ext in [tuple(L0[:3]), noncol]:
            for trial in range(20):
                order = list(range(v)); rng.shuffle(order)
                bar = pairing(v, ext, order)
                r, M = rho_for_pairing(W, bar)
                cands.append((r, ext, bar))
        cands.sort(key=lambda t: t[0])
        print(f'q={q} v={v} |T|={v+3} rose edges={(v+3)//2}: rho(M) min={cands[0][0]:.6f} max={cands[-1][0]:.6f}', flush=True)
