"""Independent checks for candidate #3 (V=7):
 (a) brute force over ALL maps preserving the stable WL colouring: which are homomorphisms;
 (b) prop:complex (1)-(3) with explicit bit labels;
 (c) WL classes listed."""
import itertools, sys
from cplx import load, wl_classes, automorphisms
idx = int(sys.argv[1]) if len(sys.argv)>1 else 3
cx = load('cands.txt')[idx]
col = wl_classes(cx)
print('faces', cx.faces)
pcls = {}; fcls = {}
for u in range(cx.V): pcls.setdefault(col[u], []).append(u)
for j in range(cx.F): fcls.setdefault(col[cx.V+j], []).append(j)
print('point classes', list(pcls.values())); print('face classes', [[cx.faces[j] for j in c] for c in fcls.values()])
pa = [pcls[col[u]] for u in range(cx.V)]; fb = [fcls[col[cx.V+j]] for j in range(cx.F)]
nmaps = 0; nh = []
fs=[set(f) for f in cx.faces]
for a in itertools.product(*pa):
    for b in itertools.product(*fb):
        nmaps += 1
        if all(a[i] in fs[b[j]] for (i,j) in cx.corners): nh.append((a,b))
print('colour-preserving maps', nmaps, 'homomorphisms among them', len(nh))
for a,b in nh: print('  ', a, b)
print('automorphisms', automorphisms(cx))
# (b) find bit labels: need 3 bipartitions into two good halves with a_e>=1/3 all, >1/3 some
def good(S):
    pts=set(x for j in S for x in cx.faces[j])
    if len(pts)!=cx.V or len(S)!=cx.V-2: return False
    S=list(S); seen={S[0]}; st=[S[0]]
    while st:
        j=st.pop()
        for k in S:
            if k not in seen and len(fs[j]&fs[k])==2: seen.add(k); st.append(k)
    return len(seen)==len(S)
cuts=[]
for S in itertools.combinations(range(cx.F), cx.V-2):
    if 0 in S and good(S) and good([j for j in range(cx.F) if j not in S]): cuts.append(set(S))
print('number of good splits', len(cuts))
sol=None
for c3 in itertools.combinations(cuts,3):
    ok=True; strict=False
    for e in cx.pairs:
        j,l=cx.pfaces[e]
        d=sum((j in c)!=(l in c) for c in c3)
        if d==0: ok=False;break
        if d>1: strict=True
    if ok and strict: sol=c3; break
print('bits found:', sol is not None)
if sol:
    for j in range(cx.F): print('   face',j,cx.faces[j],'bits',''.join('0' if j in c else '1' for c in sol))
    for e in cx.pairs:
        j,l=cx.pfaces[e]; print('   pair',e,'faces',cx.faces[j],cx.faces[l],'a_e=%d/3'%sum((j in c)!=(l in c) for c in sol))
