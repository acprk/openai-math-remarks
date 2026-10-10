"""Common complex utilities for the 161 (Sidorenko) relax-P4 study.
A complex = list of faces (sorted 3-tuples of points 0..V-1).
H = incidence graph (points I, faces J, corners C = (i,j) with i in face j).
"""
import itertools, sys
from collections import defaultdict, Counter


def parse_line(line):
    part = line.split('|')[0].split()
    faces = [tuple(sorted(int(x) for x in t.split(','))) for t in part]
    V = 1 + max(max(f) for f in faces)
    return V, faces


class Complex:
    def __init__(self, faces):
        self.faces = [tuple(sorted(f)) for f in faces]
        self.V = 1 + max(max(f) for f in self.faces)
        self.F = len(self.faces)
        self.corners = [(i, j) for j, f in enumerate(self.faces) for i in f]
        self.cidx = {c: k for k, c in enumerate(self.corners)}
        pf = defaultdict(list)
        for j, f in enumerate(self.faces):
            for a, b in itertools.combinations(f, 2):
                pf[(a, b)].append(j)
        assert all(len(v) == 2 for v in pf.values()), "pair not in exactly 2 faces"
        self.pairs = sorted(pf)
        self.pidx = {e: k for k, e in enumerate(self.pairs)}
        self.pfaces = {e: tuple(pf[e]) for e in self.pairs}
        # U_e : the four corners of pair e (as corner indices)
        self.U = []
        for e in self.pairs:
            v, w = e
            j, l = self.pfaces[e]
            self.U.append(frozenset(self.cidx[(x, y)] for x in (v, w) for y in (j, l)))
        # corner -> its two pairs (indices)
        self.cpairs = []
        for (i, j) in self.corners:
            others = [x for x in self.faces[j] if x != i]
            self.cpairs.append(tuple(self.pidx[tuple(sorted((i, o)))] for o in others))
        # M_kappa and Gamma_kappa
        self.M, self.G = [], []
        for k, (e, f) in enumerate(self.cpairs):
            assert self.U[e] & self.U[f] == {k}
            self.M.append(self.U[e] | self.U[f])
            self.G.append(self.U[e] ^ self.U[f])

    def H_edges(self):
        return [(i, self.V + j) for (i, j) in self.corners]

    def nvert(self):
        return self.V + self.F

    def degrees(self):
        return [sum(1 for f in self.faces if i in f) for i in range(self.V)]


def wl_classes(cx):
    """Color refinement on H starting from the two part colors.
    Returns stable color list over V+F vertices."""
    n = cx.nvert()
    adj = defaultdict(list)
    for u, v in cx.H_edges():
        adj[u].append(v); adj[v].append(u)
    col = [0] * cx.V + [1] * cx.F
    while True:
        sig = [(col[u], tuple(sorted(Counter(col[v] for v in adj[u]).items()))) for u in range(n)]
        m = {s: k for k, s in enumerate(sorted(set(sig)))}
        new = [m[s] for s in sig]
        if len(set(new)) == len(set(col)):
            return new
        col = new


def automorphisms(cx):
    """All point permutations preserving the face set (= Aut(H) part-preserving,
    since a face is determined by its point set; H connected bipartite with
    unequal parts, so Aut(H) = these)."""
    fs = set(cx.faces)
    nbr = defaultdict(set)
    for a, b in cx.pairs:
        nbr[a].add(b); nbr[b].add(a)
    deg = cx.degrees()
    out = []
    V = cx.V
    def rec(p, used):
        k = len(p)
        if k == V:
            if all(tuple(sorted(p[x] for x in f)) in fs for f in cx.faces):
                out.append(tuple(p))
            return
        for t in range(V):
            if t in used or deg[t] != deg[k]:
                continue
            ok = True
            for x in range(k):
                if (x in nbr[k]) != (p[x] in nbr[t]):
                    ok = False; break
            if ok:
                p.append(t); used.add(t)
                rec(p, used)
                p.pop(); used.discard(t)
    rec([], set())
    return out


def load(path):
    out = []
    for line in open(path):
        if line.strip():
            V, faces = parse_line(line)
            out.append(Complex(faces))
    return out
