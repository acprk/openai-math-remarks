"""Fractional (Shearer-weight) version of Prop. complex (2)(3): is there a probability
distribution w on good splits with  P_w(split cuts e) >= 1/3 for every pair e ?
(The tail.tex Prop global-cost argument works verbatim with b bit positions / rational
weights; b=3 is the integral case.)  For infeasible cases emit an exact dual certificate:
y >= 0 on pairs, sum y = 1, with y(cut) < 1/3 for EVERY good split (checked in Fractions).
stdin: output of DUMP=1 ./sidsearch_dump"""
import sys
from fractions import Fraction as Fr
from scipy.optimize import linprog
import numpy as np
best=0; nfeas=0; ncert=0; tot=0
for line in sys.stdin:
    if not line.startswith('T '): continue
    p=line.split(); idx,V,E=int(p[1]),int(p[2]),int(p[3]); cuts=[int(x,16) for x in p[4:]]
    tot+=1
    if not cuts: ncert+=1; continue
    S=len(cuts)
    A=np.array([[(c>>e)&1 for c in cuts] for e in range(E)],dtype=float)
    # dual: min s  s.t. sum_e y_e A[e,s] <= s for all splits, sum y=1, y>=0
    c=np.zeros(E+1); c[-1]=1
    Aub=np.hstack([A.T, -np.ones((S,1))]); bub=np.zeros(S)
    Aeq=np.hstack([np.ones((1,E)),np.zeros((1,1))]); beq=[1]
    r=linprog(c,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*E+[(None,None)],method='highs')
    tstar=r.fun; best=max(best,tstar)
    if tstar>=1/3-1e-9:
        nfeas+=1; print('FRACTIONAL-FEASIBLE V=%d idx=%d t*=%.6f'%(V,idx,tstar),flush=True); continue
    y=[Fr(v).limit_denominator(10**6) if v>1e-12 else Fr(0) for v in r.x[:E]]
    s=sum(y); y=[v/s for v in y]
    m=max(sum(y[e] for e in range(E) if c_>>e&1) for c_ in cuts)
    assert m<Fr(1,3), (idx,m)
    ncert+=1
print('SUMMARY V=%d triangulations=%d fractional-feasible=%d exact-dual-certified-infeasible=%d max t*=%.6f'%(V,tot,nfeas,ncert,best))
