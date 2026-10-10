# Independent edge-aware pinning check. Homs f=(g,b): H->H, g: points->points, b: faces->faces with g(j) ⊆ b(j).
# Score Q(f)=sum_i s_P[g(i)] + sum_j s_F[b(j)] + sum_{i in j} s_I[(g(i),b(j))].
# Want integer scores with Q(id)-Q(f) >= 1 for every hom f that is not an automorphism.
# Column generation (HiGHS LP) to find scores; final check in exact integer arithmetic over all V^V point maps.
# Usage: python3 pin_margin.py <candidate index in construction/cands.txt>   (3 = K_7, 6 = V=8 #6)
import sys, itertools, numpy as np
from scipy.optimize import linprog
from fractions import Fraction as Fr
def load(i, path='../construction/cands.txt'):
    line=open(path).read().splitlines()[i]
    return [tuple(map(int,f.split(','))) for f in line.split('|')[0].split()]
F=load(int(sys.argv[1])); V=1+max(map(max,F)); nF=len(F)
inc=[(i,j) for j,f in enumerate(F) for i in f]; incidx={c:k for k,c in enumerate(inc)}
nvar=V+nF+len(inc)
G=np.array(list(itertools.product(range(V),repeat=V)),dtype=np.int8)   # all point maps
fs=[set(f) for f in F]
def automorphisms():
    out=[]
    for p in itertools.permutations(range(V)):
        if all(set(p[i] for i in f) in fs for f in F): out.append(p)
    return out
AUT=automorphisms(); autset=set(AUT)
print('V',V,'|Aut|',len(AUT))
def vec(g,b):
    x=np.zeros(nvar)
    for i in range(V): x[g[i]]+=1
    for j in range(nF):
        x[V+b[j]]+=1
        for i in F[j]: x[V+nF+incidx[(g[i],b[j])]]+=1
    return x
xid=vec(range(V),range(nF))
def oracle(s, need=200):
    """For every point map g, best face images; return list of (gap, g, b) for non-automorphisms, sorted."""
    sP,sF,sI=s[:V],s[V:V+nF],s[V+nF:]
    # table T[j][tuple g(j)] -> (best value, best face) ; precompute per face over V^3 image triples
    base=sP[G].sum(1).astype(object) if s.dtype==object else sP[G].sum(1)
    tot=base.copy(); choice=[]
    ok=np.ones(len(G),bool)
    for j,f in enumerate(F):
        a,b_,c=f
        best=np.full((V,V,V),-np.inf if s.dtype!=object else None,dtype=s.dtype)
        arg=np.full((V,V,V),-1)
        for x,y,z in itertools.product(range(V),repeat=3):
            cand=[k for k in range(nF) if {x,y,z}<=fs[k]]
            if not cand: continue
            vals=[sF[k]+sI[incidx[(x,k)]]+sI[incidx[(y,k)]]+sI[incidx[(z,k)]] for k in cand]
            m=max(range(len(cand)),key=lambda t:vals[t]); best[x,y,z]=vals[m]; arg[x,y,z]=cand[m]
        idx=(G[:,a],G[:,b_],G[:,c])
        ok&=arg[idx]>=0
        tot=tot+np.where(arg[idx]>=0,best[idx],0)
        choice.append(arg[idx])
    choice=np.array(choice).T
    Qid=xid@s
    rows=[]
    order=np.argsort(-np.where(ok,tot,-1e18).astype(float))
    for r in order:
        if not ok[r]: break
        g=tuple(int(t) for t in G[r])
        if g in autset: continue
        rows.append((Qid-tot[r],g,tuple(int(t) for t in choice[r])))
        if len(rows)>=need: break
    return rows
# LP: max delta s.t. (xid-x)·s >= delta for x in rows, |s|<=1
rows=[]; s=np.zeros(nvar)
for it in range(200):
    viol=oracle(s)
    new=[vec(g,b) for gap,g,b in viol if gap< (delta if it else 1e9)-1e-9] if it else [vec(g,b) for _,g,b in viol]
    if it and not new: break
    rows+=new
    A=np.array([x-xid for x in rows]); n=len(rows)
    # vars s (nvar), delta ; minimize -delta ; A s + delta <= 0
    res=linprog(np.r_[np.zeros(nvar),-1],A_ub=np.c_[A,np.ones(n)],b_ub=np.zeros(n),bounds=[(-1,1)]*nvar+[(None,None)],method='highs')
    s=res.x[:nvar]; delta=-res.fun
    
    if delta<1e-9: print('NO PINNING (delta*=0)'); sys.exit()
# exact verification with integer scores
den=10**4; si=np.array([round(v*den) for v in s],dtype=np.int64)
viol=oracle(si.astype(np.float64),need=5)
print('LP margin delta* = %.4f (scores in [-1,1])' % delta)
print('integer scores (x10^4): every non-automorphism homomorphism scores strictly below id:', viol[0][0] > 0)
