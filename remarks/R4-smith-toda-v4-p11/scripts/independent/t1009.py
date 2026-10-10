# R4 independent positive control at p=1009. Run from scripts/independent/: python t1009.py (z3-solver).
# Expected output: expected_output/t1009.txt
from mayenum import *
import time
p=1009; q=2*p-2
for lab,(s,w) in {'k1 (q+1,u4+1)':(q+1,u(p,4)+1),'(q+1,u5+1)':(q+1,u(p,5)+1),'k=2':(2*q+1,u(p,4)+2),'k=p-1':(q*(p-1)+1,u(p,4)+p-1),'k=p':(q*p+1,u(p,4)+p),'k=p+1':(q*(p+1)+1,u(p,4)+p+1),'k=3':(3*q+1,u(p,4)+3)}.items():
    t=time.time(); S=enum_z3(p,s,w); print(lab,len(S),[name(p,m) for m in S], '%.1fs'%(time.time()-t))
