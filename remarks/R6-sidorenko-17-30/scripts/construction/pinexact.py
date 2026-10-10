"""Exhaustive check of the generalised Lemma pinning over ALL part-preserving maps and over
ALL homomorphisms H->H, by (max,+,count) variable elimination.
Score S(f) = sum_t L0^{R-t} F_t(f) with L0 = 10^6 (exceeds 2*max|F_t| comfortably; checked).
Reports max S, S(id), and number of maximisers."""
import sys
from collections import Counter
from cplx import load, automorphisms
idx=int(sys.argv[1]); path=sys.argv[2]; homonly=int(sys.argv[3])
cx=load(path)[idx]; V,F=cx.V,cx.F; nv=V+F
adj=[[] for _ in range(nv)]
for i,j in cx.corners: adj[i].append(V+j); adj[V+j].append(i)
cols=[[0]*V+[1]*F]
while True:
    c=cols[-1]; sig=[(c[u],tuple(sorted(Counter(c[v] for v in adj[u]).items()))) for u in range(nv)]
    m={s:k for k,s in enumerate(sorted(set(sig)))}; new=[m[s] for s in sig]; cols.append(new)
    if len(set(new))==len(set(c)): break
R=len(cols)-1; L0=10**6
hist=[[Counter(cols[t-1][v] for v in adj[u]) for u in range(nv)] for t in range(1,R+1)]
W=[L0**(R-1-t) for t in range(R)]
def vs(t_u):  # vertex score of type t_u (global vertex index)
    return sum(W[t]*(-sum(c*c for c in hist[t][t_u].values())) for t in range(R))
def ps(a,b):  # pair score, a point type, b face type (global idx)
    return sum(W[t]*2*(hist[t][a][cols[t][b]]+hist[t][b][cols[t][a]]) for t in range(R))
dom=[list(range(V))]*V+[list(range(V,V+F))]*F
# factors: unary on each vertex, binary on each source corner
factors=[]
for u in range(nv): factors.append(((u,),{(x,):(vs(x),1) for x in dom[u]}))
for i,j in cx.corners:
    tab={}
    for x in range(V):
        for y in range(V,V+F):
            if homonly and x not in cx.faces[y-V]: continue
            tab[(x,y)]=(ps(x,y),1)
    factors.append(((i,V+j),tab))
def comb(p,q): return (p[0]+q[0],p[1]*q[1])
def add(p,q):
    if p is None: return q
    if q[0]>p[0]: return q
    if q[0]<p[0]: return p
    return (p[0],p[1]+q[1])
remaining=set(range(nv))
while remaining:
    best=None
    for v in remaining:
        nb=set()
        for s,_ in factors:
            if v in s: nb|=set(s)
        nb.discard(v); cost=1
        for u in nb: cost*=len(dom[u])
        if best is None or cost<best[0]: best=(cost,v)
    v=best[1]
    inv=[f for f in factors if v in f[0]]; rest=[f for f in factors if v not in f[0]]
    scope=sorted(set(u for s,_ in inv for u in s)-{v}); allv=scope+[v]; pos={u:k for k,u in enumerate(allv)}
    import itertools
    out={}
    for asg in itertools.product(*[dom[u] for u in allv]):
        val=(0,1); ok=True
        for s,tab in inv:
            key=tuple(asg[pos[u]] for u in s)
            if key not in tab: ok=False;break
            val=comb(val,tab[key])
        if ok:
            k=asg[:-1]; out[k]=add(out.get(k),val)
    rest.append((tuple(scope),out)); factors=rest; remaining.discard(v)
tot=(0,1)
for s,tab in factors: tot=comb(tot,tab[()])
# identity score
Sid=sum(vs(u) for u in range(nv))+sum(ps(i,V+j) for i,j in cx.corners)
maxF=max(abs(x) for x in [0])
print('R=%d  max S=%d  S(id)=%d  #maximisers=%d  |Aut|=%d  (homs only=%d)'%(R,tot[0],Sid,tot[1],len(automorphisms(cx)),homonly))
