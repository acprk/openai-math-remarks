"""Independent re-check (adversarial verifier) of the V=7 candidate for OAI 161.
Written from the tex (kernel.tex, complex.tex), not from the other agent's code.
(1)-(3) of prop:complex, stable colour refinement, automorphisms,
generalised pinning: maximisers of the Lemma-pinning score among ALL part-preserving
homomorphisms H->H (non-homs have C=0 exactly), and exact C(f) as polynomial in theta."""
import itertools, sys
from collections import Counter
from fractions import Fraction

FACES = [(0,1,2),(0,1,6),(0,2,3),(0,3,4),(0,4,5),(0,5,6),(1,2,5),(1,5,6),(2,3,4),(2,4,5)]
BITS  = ['000','001','001','011','111','011','100','110','100','110']   # from verify3 output, re-verified below
V=7; F=len(FACES); fs=[frozenset(f) for f in FACES]
pairs=sorted({frozenset(p) for f in FACES for p in itertools.combinations(f,2)}, key=sorted)
pf={e:[j for j in range(F) if e<=fs[j]] for e in pairs}
# ---- (1)
assert all(len(pf[e])==2 for e in pairs); print('(1) every pair in exactly 2 faces; |E|=',len(pairs),' Euler',V-len(pairs)+F)
def conn(S):
    S=list(S); seen={S[0]}; st=[S[0]]
    while st:
        j=st.pop()
        for k in S:
            if k not in seen and len(fs[j]&fs[k])==2: seen.add(k); st.append(k)
    return len(seen)==len(S)
assert conn(range(F)); print('    face-neighbour graph connected')
# ---- (2): exposure order with exactly one new point per later face, sharing a pair with an earlier face
def exposure(S):
    S=list(S)
    for root in S:
        order=[root]; pts=set(fs[root]); rest=set(S)-{root}; ok=True
        while rest:
            nxt=[k for k in rest if any(len(fs[k]&fs[o])==2 for o in order) and len(fs[k]-pts)==1]
            if not nxt: ok=False; break
            k=min(nxt); order.append(k); pts|=fs[k]; rest.discard(k)
        if ok and len(pts)==V: return order
    return None
for pos in range(3):
    for val in '01':
        S=[j for j in range(F) if BITS[j][pos]==val]
        o=exposure(S); assert len(S)==V-2 and o is not None, (pos,val)
        print('(2) bit %d=%s class'%(pos+1,val),S,'order',o)
# ---- (3)
ae={e:sum(BITS[pf[e][0]][p]!=BITS[pf[e][1]][p] for p in range(3)) for e in pairs}
assert min(ae.values())>=1 and max(ae.values())>1
print('(3) 3a_e values', sorted(Counter(ae.values()).items()))
# point-graph degrees (transverse needs max<=6 with D>=12)
deg=Counter(x for e in pairs for x in e); print('point-graph degrees',[deg[i] for i in range(V)])
# ---- H, colour refinement
nv=V+F; adj=[[] for _ in range(nv)]
corners=[(i,j) for j in range(F) for i in FACES[j]]
for i,j in corners: adj[i].append(V+j); adj[V+j].append(i)
cols=[[0]*V+[1]*F]; hists=[]
while True:
    c=cols[-1]; h=[Counter(c[w] for w in adj[u]) for u in range(nv)]
    sig=[(c[u],tuple(sorted(h[u].items()))) for u in range(nv)]
    m={s:k for k,s in enumerate(sorted(set(sig)))}; new=[m[s] for s in sig]
    hists.append(h); cols.append(new)
    if len(set(new))==len(set(c)): break
R=len(hists); stable=cols[-1]
print('rounds',R,'stable classes points',[[u for u in range(V) if stable[u]==k] for k in sorted(set(stable[:V]))],
      'faces',[[FACES[u-V] for u in range(V,nv) if stable[u]==k] for k in sorted(set(stable[V:]))])
# automorphisms
auts=[]
for p in itertools.permutations(range(V)):
    img=[frozenset(p[x] for x in f) for f in FACES]
    if set(img)==set(fs): auts.append((p,tuple(fs.index(s) for s in img)))
print('|Aut|=',len(auts), auts)
# ---- pinning score (Lemma pinning), lexicographic over rounds via L0 huge
def Ft(t,a,b):   # a: points->points, b: faces->faces ; t index 0..R-1 uses hists[t], cols[t]
    h=hists[t]; cp=cols[t]
    f=lambda u: a[u] if u<V else V+b[u-V]
    s=0
    for u in range(nv):
        s-=sum(x*x for x in h[f(u)].values())
        for w in adj[u]: s+=2*h[f(u)][cp[f(w)]]
    return s
ida=tuple(range(V)); idb=tuple(range(F))
L0=10**6
def Q(a,b): return sum(L0**(R-1-t)*Ft(t,a,b) for t in range(R))
Qid=Q(ida,idb)
# decomposed: Q = sum_u vs(f(u)) + sum_corners ps(a(i),b(j))
def vs(x): return sum(L0**(R-1-t)*(-sum(c*c for c in hists[t][x].values())) for t in range(R))
def ps(x,y): return sum(L0**(R-1-t)*2*(hists[t][x][cols[t][y]]+hists[t][y][cols[t][x]]) for t in range(R))
VS=[vs(x) for x in range(nv)]; PS=[[ps(x,V+j) for j in range(F)] for x in range(V)]
assert Qid==sum(VS)+sum(PS[i][j] for i,j in corners)
best=None; maxers=[]; nhom=0
for a in itertools.product(range(V),repeat=V):
    tot=sum(VS[a[i]] for i in range(V)); choices=[]; ok=True
    for j in range(F):
        img=frozenset(a[x] for x in FACES[j]); cand=[y for y in range(F) if img<=fs[y]]
        if not cand: ok=False; break
        sc=[(VS[V+y]+sum(PS[a[x]][y] for x in FACES[j]),y) for y in cand]
        mx=max(s for s,_ in sc); tot+=mx; choices.append([y for s,y in sc if s==mx])
    if not ok: continue
    nhom+=1
    if best is None or tot>best: best=tot; maxers=[(a,choices)]
    elif tot==best: maxers.append((a,choices))
print('point maps admitting a hom:',nhom,' max hom score - Q(id) =',best-Qid)
mh=[(a,b) for a,ch in maxers for b in itertools.product(*ch)]
print('maximising homs:',len(mh)); [print('   ',x) for x in mh]
# ---- exact C(f) as polynomial in theta (dict deg->Fraction), sigma: one negative corner (0,face0)
K=corners; kidx={k:n for n,k in enumerate(K)}
def U(e): j,l=pf[e]; return frozenset(kidx[(x,y)] for y in (j,l) for x in e)
def two_pairs(i,j): return [frozenset({i,x}) for x in FACES[j] if x!=i]
Mset={};Gset={}
for (i,j) in K:
    e,f=two_pairs(i,j); Mset[kidx[(i,j)]]=U(e)|U(f); Gset[kidx[(i,j)]]=U(e)^U(f)
    assert len(Mset[kidx[(i,j)]])==7 and len(Gset[kidx[(i,j)]])==6
def A(tau,sig,al,be):  # returns dict poly
    Qs=al|be; Rs=al^be; out=Counter()
    for half,S,sg in ((Fraction(1,2),Mset[tau],1),(Fraction(1,2),Gset[tau],-1)):
        if Qs<=S:
            if not Rs: out[0]+=half
            elif Rs==Gset[tau]: out[1]+=half*sg*sig
    return out
def pmul(p,q):
    r=Counter()
    for a,x in p.items():
        for b,y in q.items(): r[a+b]+=x*y
    return Counter({k:v for k,v in r.items() if v})
def Cpoly(a,b,sigma):
    tgt={kidx[(i,j)]:kidx[(a[i],b[j])] for (i,j) in K}
    cand={}
    for e in pairs:
        inter=frozenset(range(len(K)))
        for kk in U(e): inter&=Mset[tgt[kk]]
        L=sorted(inter); cand[e]=[frozenset(s) for r in range(1,len(L)+1) for s in itertools.combinations(L,r)]
    order=sorted(pairs,key=lambda e:len(cand[e]))
    cornerpairs={kidx[(i,j)]:two_pairs(i,j) for (i,j) in K}
    total=Counter()
    def rec(n,lab,poly):
        nonlocal total
        if n==len(order):
            for k,v in poly.items(): total[k]+=v
            return
        e=order[n]
        for s in cand[e]:
            lab[e]=s; p=poly; ok=True
            for kk in U(e):  # corners of e : check those whose both pairs labelled and include e
                e1,e2=cornerpairs[kk]
                if e1 in lab and e2 in lab and (e in (e1,e2)):
                    other=e2 if e==e1 else e1
                    if order.index(other)<=n:
                        fac=A(tgt[kk],sigma[tgt[kk]],lab[e1],lab[e2])
                        if not fac: ok=False;break
                        p=pmul(p,fac)
            if ok: rec(n+1,lab,p)
            del lab[e]
    rec(0,{},Counter({0:Fraction(1)}))
    return {k:v for k,v in total.items() if v}
sigma={n:1 for n in range(len(K))}; sigma[kidx[(0,0)]]=-1
tot=Counter()
for a,b in mh:
    c=Cpoly(a,b,sigma); print('C',a,b,'=',{k:str(v) for k,v in c.items()}, ' (theta/2)^30 coeff check', c.get(30,0)*2**30)
    for k,v in c.items(): tot[k]+=v
print('sum over maximisers:',{k:str(v*2**30) for k,v in tot.items()},'x (theta)^k/2^30')
