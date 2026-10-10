# R4 independent d1 / quotient-table / page-bound checks. Run from scripts/independent/: python d1check.py 11 (z3-solver).
# Expected output: expected_output/d1check_p{7,11,13,17}.txt
# Independent d1 (splitting) computations + Frobenius quotient-group table + inequality page bound, at prime p.
import sys
from mayenum import enum_z3, u, name
p=int(sys.argv[1]); q=2*p-2
def iv(g): return (g[2], g[2]+g[1])          # h_{i,j} -> [j, j+i)
def hg(a,b): return ('h', b-a, a)
def d1_mono(mono):
    ext, pol = mono
    ext = sorted(ext, key=iv)
    out = {}
    for k,g in enumerate(ext):
        a,b = iv(g)
        for l in range(a+1,b):
            new = ext[:k]+[hg(a,l),hg(l,b)]+ext[k+1:]
            if len(set(new))<len(new): continue
            sign = (-1)**k
            # sort with permutation sign
            arr=list(new); sg=1
            for i in range(len(arr)):
                for j in range(len(arr)-1-i):
                    if iv(arr[j])>iv(arr[j+1]): arr[j],arr[j+1]=arr[j+1],arr[j]; sg=-sg
            key=(tuple(arr), tuple(sorted(pol)))
            out[key]=(out.get(key,0)+sign*sg)%p
    return {k:v for k,v in out.items() if v}
def rank(rows):
    rows=[dict(r) for r in rows]; r=0; piv=[]
    keys=sorted({k for row in rows for k in row}, key=str)
    M=[[row.get(k,0)%p for k in keys] for row in rows]
    rk=0; cols=len(keys)
    for c in range(cols):
        pr=next((i for i in range(rk,len(M)) if M[i][c]),None)
        if pr is None: continue
        M[rk],M[pr]=M[pr],M[rk]; inv=pow(M[rk][c],p-2,p)
        M[rk]=[x*inv%p for x in M[rk]]
        for i in range(len(M)):
            if i!=rk and M[i][c]:
                f=M[i][c]; M[i]=[(x-f*y)%p for x,y in zip(M[i],M[rk])]
        rk+=1
    return rk
def canon(m): return (tuple(sorted(m[0],key=iv)), tuple(sorted(m[1])))
# (C) n5 table: d1 injective on span of m=3,5 entries
S5=[canon(m) for m in enum_z3(p,q+1,u(p,5)+1)]
low=[m for m in S5 if len(m[0]) in (3,5)]
print('(C) #entries',len(S5),' m-profile',sorted(len(m[0]) for m in S5),' rank d1 on m<=5 span =',rank([d1_mono(m) for m in low]),'of',len(low))
# also k=1 target: acb not cycle, acybz cycle
for m in enum_z3(p,q+1,u(p,4)+1): print('    k=1',name(p,m),' d1 nonzero:',bool(d1_mono(canon(m))))
# (E) (6,W+1) d1-cohomology
W=2*p*p+4*p+1
src=[canon(m) for m in enum_z3(p,5,W+1)]; mid=[canon(m) for m in enum_z3(p,6,W+1)]
rin=rank([d1_mono(m) for m in src]); rout=rank([d1_mono(m) for m in mid])
print('(E) (6,W+1): dim',len(mid),' rank in',rin,' rank out',rout,' => H =',len(mid)-rout-rin)
# (D) quotient groups at rescaled weights (P_r, r<=3)
tab=[(p*(p-1),2*p-4,0),(p*(p-1),2*p-3,1),(p*(p-2),2*p-4,1),(p*(p-3),2*p-5,0),(p+2,1,0),(p+2,3,0),(p+2,2,1),(2*p+2,2,0),(2*p+2,3,1),(3*p+2,3,0),(2,1,0),(2,3,0)]
okD=True
for w,s,c in tab:
    L=[m for m in enum_z3(p,s,w) if all(g[1]<=3 for g in m[0]) and all(g[0][1]<=3 for g in m[1])]
    if len(L)!=c: okD=False; print('   D mismatch',w,s,[name(p,m) for m in L])
print('(D) quotient table matches:',okD)
N=p*(p-1)-2; K=(p-1)**2; S=5+2*K
print('(F) S=%d < qp+1=%d : %s ; N-K=%d'%(S,q*p+1,S<q*p+1,N-K))
