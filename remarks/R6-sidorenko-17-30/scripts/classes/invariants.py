#!/usr/bin/env python3
"""Graph invariants of H = point-face incidence graph of the V=7 complex, as used in the
"known Sidorenko classes" table of README/note (Section 5).  Pure computation; the literature
criteria themselves are stated in the note, this script only supplies the numbers.
Environment: Python >= 3.8, networkx (tested 3.4.2).  Runtime < 2 s.
Run: python3 invariants.py > ../expected_output/invariants.txt
"""
import itertools
import networkx as nx

FACES = [(0,1,2),(0,1,6),(0,2,3),(0,3,4),(0,4,5),(0,5,6),(1,2,5),(1,5,6),(2,3,4),(2,4,5)]
V, F = 7, len(FACES)
G = nx.Graph()
P = [('p', i) for i in range(V)]
Q = [('f', j) for j in range(F)]
G.add_nodes_from(P + Q)
for j, f in enumerate(FACES):
    for i in f:
        G.add_edge(('p', i), ('f', j))
n, m = G.number_of_nodes(), G.number_of_edges()
print(f"|V(H)| = {n} (prime: {all(n % k for k in range(2, n))}), |E(H)| = {m}, parts {len(P)} + {len(Q)}")
print("edge list (point i -- face j), faces j=0..9 =", ["".join(map(str, f)) for f in FACES])
for j, f in enumerate(FACES):
    print(f"   f{j}: " + ", ".join(f"p{i}" for i in f))
dp = sorted(G.degree(v) for v in P); dq = sorted(G.degree(v) for v in Q)
print("point-side degrees", dp, " face-side degrees", dq)
print("biregular (both sides regular):", len(set(dp)) == 1 and len(set(dq)) == 1)
print("connected:", nx.is_connected(G), " tree:", nx.is_tree(G), " min degree:", min(dp + dq),
      " (min degree >= 3 => not a subdivision with degree-2 vertices)")
print("girth:", nx.girth(G), " (H is not a cycle: degrees >= 3)")
dom = [v for v in P if G.degree(v) == len(Q)] + [v for v in Q if G.degree(v) == len(P)]
print("vertices complete to the other side (Conlon-Fox-Sudakov):", dom)
print("max degree on face side vs |P|-1:", max(dq), "vs", len(P) - 1,
      "; max degree on point side vs |Q|-1:", max(dp), "vs", len(Q) - 1)
twins = [(u, v) for u, v in itertools.combinations(G.nodes, 2) if set(G[u]) == set(G[v])]
print("twin pairs (identical neighbourhoods; a nontrivial blow-up has them):", twins)
GM = nx.algorithms.isomorphism.GraphMatcher(
    G, G, node_match=lambda a, b: True)
auts = [a for a in GM.isomorphisms_iter() if all(a[v][0] == v[0] for v in G)]
print("|Aut(H)| preserving parts:", len(auts))
allauts = list(nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
print("|Aut(H)| (all):", len(allauts))
for a in allauts:
    if any(a[v] != v for v in G):
        fix = [v for v in G if a[v] == v]
        rest = G.subgraph([v for v in G if v not in fix])
        print("nontrivial automorphism fixes", fix, "; H - Fix connected:", nx.is_connected(rest))
print("4-cycles exist (two faces sharing a pair):", any(len(set(f) & set(g)) == 2 for f, g in itertools.combinations(FACES, 2)))
