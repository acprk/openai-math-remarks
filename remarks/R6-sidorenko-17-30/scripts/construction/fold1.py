"""#1 (V=6) obstruction certificate:
 the 8 non-automorphism WL-face homs with lowest term theta^22 split into 4 pairs
 (f,f') with vec(f)+vec(f') = 2 vec(id)  (vec = point, face, corner image counts),
 and their C is +theta^22/2^24 for random sigma with prod sigma = -1 (C(aut)= -theta^24/2^24).
 AM-GM: w(f)+w(f') >= 2 w(id) for any priors/thinning (w = monomial x^vec)."""
import random, itertools
from cplx import load, wl_classes, automorphisms
from homs import homs
from ccoef import C_poly
from vec import vec
cx=load('cands.txt')[1]; col=wl_classes(cx); nv=cx.nvert()
allowed=[set(w if w<cx.V else w-cx.V for w in range(nv) if col[w]==col[u]) for u in range(nv)]
H=list(homs(cx,cx,allowed)); auts=set(automorphisms(cx))
idv=vec(cx,tuple(range(cx.V)),tuple(range(cx.F)))
folds=[h for h in H if C_poly(cx,cx,*h)==[0]*22+[1]]
print('fold maps', len(folds))
used=set(); pairs=[]
for f,g in itertools.combinations(folds,2):
    s=vec(cx,*f); s.update(vec(cx,*g))
    if all(s[k]==2*idv[k] for k in set(s)|set(idv)): pairs.append((f,g))
print('midpoint pairs', len(pairs)); 
for p in pairs: print('  ',p)
random.seed(7); bad=0
for t in range(20):
    sig={k:random.choice((1,-1)) for k in range(len(cx.corners))}
    if eval('*'.join(str(v) for v in sig.values()))!=-1: sig[0]*=-1
    for f in folds:
        if C_poly(cx,cx,*f,sigma=sig)!=[0]*22+[1]: bad+=1
    for a in auts:
        b=tuple(cx.faces.index(tuple(sorted(a[x] for x in f))) for f in cx.faces)
        if C_poly(cx,cx,a,b,sigma=sig)!=[0]*24+[-1]: bad+=1
print('sigma-dependence violations', bad)
