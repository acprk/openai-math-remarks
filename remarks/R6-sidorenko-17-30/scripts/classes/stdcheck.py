"""Necessary condition for strong tree-decomposability (CKLL 1510.06533, Def. before Thm 1.2):
a leaf bag X with parent Y gives Z=X\\Y nonempty with N(Z) subset of X, H[X] an induced tree,
S=X∩Y=N(Z)\\Z independent, and H-Z must contain an induced copy of the minimal subtree T_X(S)
of H[X] spanning S, via an isomorphism fixing S pointwise.  (Degenerate leaf bags X⊆Y are
single vertices and can be deleted.)  If no Z works, H is not strongly tree-decomposable, hence
neither a reflection tree (Li-Szegedy) nor tree-arrangeable (Kim-Lee-Lee).
Applied iteratively this is only the first step; failure at step one already suffices."""
import itertools, networkx as nx
FACES=[(0,1,2),(0,1,6),(0,2,3),(0,3,4),(0,4,5),(0,5,6),(1,2,5),(1,5,6),(2,3,4),(2,4,5)]
import sys
if len(sys.argv)>1 and sys.argv[1]=='orig':
    FACES=[(0,1,3),(0,1,9),(0,2,3),(0,2,9),(1,3,10),(1,7,10),(1,7,12),(1,9,12),(2,3,4),(2,4,9),(3,4,11),(3,10,11),(4,8,9),(4,8,11),(5,6,7),(5,6,11),(5,7,10),(5,10,11),(6,7,12),(6,8,11),(6,8,12),(8,9,12)]
V=1+max(max(f) for f in FACES)
G=nx.Graph(); G.add_nodes_from(range(V+len(FACES)))
for j,f in enumerate(FACES):
    for i in f: G.add_edge(i,V+j)
n=G.number_of_nodes(); nodes=list(G.nodes)
def steiner(T,S):
    # minimal subtree of tree T containing S: repeatedly strip leaves not in S
    T=T.copy()
    while True:
        lv=[v for v in T if T.degree(v)<=1 and v not in S]
        if not lv or T.number_of_nodes()==len(S): break
        T.remove_nodes_from(lv)
    return T
found=[]; tested=0
for mask in range(1,1<<n):
    Z={v for v in nodes if mask>>v&1}
    if len(Z)>n-2: continue
    X=set(Z)
    for z in Z: X|=set(G[z])
    if len(X)>=n: continue
    HX=G.subgraph(X)
    if not nx.is_tree(HX): continue
    S=X-Z
    if not S: continue
    tested+=1
    T=steiner(HX,S); inner=[v for v in T if v not in S]
    Hm=G.subgraph(set(nodes)-Z)
    ok=False
    cand=[v for v in Hm if v not in S]
    for img in itertools.permutations(cand,len(inner)):
        mp=dict(zip(inner,img)); mp.update({s:s for s in S})
        if all(Hm.has_edge(mp[a],mp[b]) for a,b in T.edges()):
            W=set(mp.values())
            if Hm.subgraph(W).number_of_edges()==T.number_of_edges(): ok=True;break
    if ok: found.append((sorted(Z),sorted(S)))
print('candidate leaf sets Z with H[Z∪N(Z)] a tree:',tested,' passing the copy test:',len(found))
for x in found[:20]: print(x)
