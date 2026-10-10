#!/usr/bin/env python3
"""T-OPT search for family 161 (Sidorenko counterexample).
Find the smallest triangulated 2-sphere satisfying Proposition prop:complex of
the source (sections/complex.tex):
 (1) every pair in exactly two faces, face-neighbour graph connected  [automatic]
 (2) a labelling beta: faces -> {0,1}^3 s.t. for each bit position and value the
     class of faces is shared-pair connected, covers all V points, and has an
     ordering adding exactly one new point per face  (<=> |class| = V-2,
     dual-connected, covers all points)
 (3) the two faces on every pair have different labels; some pair differs in >=2 bits
 (4) colour refinement on the incidence graph (start: one colour per part) is discrete.
Input: plantri -a output on stdin.  Usage: plantri -a V | python3 sid_complex_search.py
Also: python3 sid_complex_search.py --paper   (sanity check on the source's complex)
"""
import sys, itertools
PAPER = [(0,1,3),(0,1,9),(0,2,3),(0,2,9),(1,3,10),(1,7,10),(1,7,12),(1,9,12),(2,3,4),(2,4,9),
         (3,4,11),(3,10,11),(4,8,9),(4,8,11),(5,6,7),(5,6,11),(5,7,10),(5,10,11),(6,7,12),(6,8,11),
         (6,8,12),(8,9,12)]
PAPER_BITS = ['000','010','001','110','001','000','100','110','101','100','111','101','101','011',
              '000','110','010','111','001','010','011','111']

def faces_from_plantri(line):
    n, rest = line.split(None, 1)
    n = int(n)
    adj = [[ord(c) - ord('a') for c in w] for w in rest.strip().split(',')]
    F = set()
    for v in range(n):
        d = adj[v]
        for i in range(len(d)):
            a, b = d[i], d[(i + 1) % len(d)]
            F.add(tuple(sorted((v, a, b))))
    return n, sorted(F)

def wl_discrete(n, F):
    # colour refinement on incidence graph; points 0..n-1, faces n..n+|F|-1
    N = n + len(F)
    nb = [[] for _ in range(N)]
    for j, f in enumerate(F):
        for i in f:
            nb[i].append(n + j); nb[n + j].append(i)
    col = [0] * n + [1] * len(F)
    while True:
        sig = [(col[v], tuple(sorted(col[u] for u in nb[v]))) for v in range(N)]
        m = {s: k for k, s in enumerate(sorted(set(sig)))}
        new = [m[s] for s in sig]
        if len(set(new)) == len(set(col)):
            return len(set(new)) == N
        col = new

def dual_edges(F):
    where = {}
    for j, f in enumerate(F):
        for e in itertools.combinations(f, 2):
            where.setdefault(e, []).append(j)
    assert all(len(v) == 2 for v in where.values())
    return [tuple(v) for v in where.values()]

def good_class(mask, F, n, DE):
    idx = [j for j in range(len(F)) if mask >> j & 1]
    if len(idx) != n - 2: return False
    if len(set(itertools.chain(*[F[j] for j in idx]))) != n: return False
    S = set(idx); seen = {idx[0]}; st = [idx[0]]
    adj = {}
    for a, b in DE:
        if a in S and b in S:
            adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    while st:
        x = st.pop()
        for y in adj.get(x, []):
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == len(S)

def search(n, F, want_all=False):
    DE = dual_edges(F); nF = len(F); full = (1 << nF) - 1
    splits = []
    for comb in itertools.combinations(range(1, nF), n - 3):   # fix face 0 in class 0
        mask = 1 | sum(1 << c for c in comb)
        if good_class(mask, F, n, DE) and good_class(full ^ mask, F, n, DE):
            cut = sum(1 << k for k, (a, b) in enumerate(DE) if (mask >> a & 1) != (mask >> b & 1))
            splits.append((mask, cut))
    allE = (1 << len(DE)) - 1
    cuts = sorted(set(c for _, c in splits))
    sols = []
    for i in range(len(cuts)):
        for j in range(i, len(cuts)):
            u = cuts[i] | cuts[j]
            for k in range(j, len(cuts)):
                if u | cuts[k] == allE:
                    two = (cuts[i] & cuts[j]) | (cuts[i] & cuts[k]) | (cuts[j] & cuts[k])
                    if two:
                        sols.append((cuts[i], cuts[j], cuts[k]))
                        if not want_all: return len(splits), sols
    return len(splits), sols

if __name__ == '__main__':
    if '--paper' in sys.argv:
        n, F = 13, [tuple(sorted(f)) for f in PAPER]
        DE = dual_edges(F)
        print('paper complex: WL discrete', wl_discrete(n, F))
        ok = True
        for k in range(3):
            m = sum(1 << j for j in range(len(F)) if PAPER_BITS[j][k] == '0')
            ok &= good_class(m, F, n, DE) and good_class(((1 << len(F)) - 1) ^ m, F, n, DE)
        diff = [sum(PAPER_BITS[a][k] != PAPER_BITS[b][k] for k in range(3)) for a, b in DE]
        print('paper bits: classes good', ok, 'min diff', min(diff), 'max diff', max(diff))
        sys.exit()
    tot = disc = hit = 0
    for line in sys.stdin:
        if not line.strip(): continue
        n, F = faces_from_plantri(line); tot += 1
        if not wl_discrete(n, F): continue
        disc += 1
        ns, sols = search(n, F)
        if sols:
            hit += 1
            print('FOUND V=%d faces=%s splits=%d' % (n, F, ns), flush=True)
    print('V=%d triangulations=%d WL-discrete=%d satisfying all=%d' % (n, tot, disc, hit))
