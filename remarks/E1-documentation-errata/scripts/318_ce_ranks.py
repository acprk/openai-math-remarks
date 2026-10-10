"""Recomputation of the trivial-coefficient Chevalley-Eilenberg ranks quoted in
family 318, "A rational obstruction to strong chromatic splitting at height
three" (SUPPORT-COVERAGE.md: support/r-ce/ce_certificate.py gives ranks
(0,3,0,0) for sl2 and (0,8,20,35) for sl3; that script is not in the repo).
Basis E_ij (i!=j), H_i = E_ii - E_{i+1,i+1}.  rank(d_k: C^k -> C^{k+1}) equals
the rank of the chain boundary Lambda^{k+1} g -> Lambda^k g.  Exact Fractions."""
from fractions import Fraction as F
from itertools import combinations

def sl(n):
    def E(i,j):
        M=[[0]*n for _ in range(n)]; M[i][j]=1; return M
    B=[E(i,j) for i in range(n) for j in range(n) if i!=j]
    for i in range(n-1):
        M=[[0]*n for _ in range(n)]; M[i][i]=1; M[i+1][i+1]=-1; B.append(M)
    return B
def mul(A,B):
    n=len(A); return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def coords(M,B):
    # solve M = sum c_b B_b over Q (B is a basis of sl_n)
    n=len(M); vecs=[[b[i][j] for i in range(n) for j in range(n)] for b in B]
    target=[M[i][j] for i in range(n) for j in range(n)]
    # least-squares-free exact solve via elimination on augmented system
    m=len(B); rows=[[F(vecs[c][r]) for c in range(m)]+[F(target[r])] for r in range(n*n)]
    piv=[];r=0
    for c in range(m):
        p=next((i for i in range(r,len(rows)) if rows[i][c]!=0),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        rows[r]=[x/rows[r][c] for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]!=0:
                f=rows[i][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[r])]
        piv.append(c); r+=1
    sol=[F(0)]*m
    for i,c in enumerate(piv): sol[c]=rows[i][m]
    assert all(all(x==0 for x in row[:m]) and row[m]==0 for row in rows[r:])
    return sol
def rank(rows):
    rows=[r[:] for r in rows]; rk=0; ncol=len(rows[0]) if rows else 0
    for c in range(ncol):
        p=next((i for i in range(rk,len(rows)) if rows[i][c]!=0),None)
        if p is None: continue
        rows[rk],rows[p]=rows[p],rows[rk]
        for i in range(rk+1,len(rows)):
            if rows[i][c]!=0:
                f=rows[i][c]/rows[rk][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[rk])]
        rk+=1
    return rk
def ce_ranks(n):
    B=sl(n); m=len(B)
    br={}
    for a in range(m):
        for b in range(m):
            br[a,b]=coords([[x-y for x,y in zip(r1,r2)] for r1,r2 in zip(mul(B[a],B[b]),mul(B[b],B[a]))],B)
    ranks=[0]  # d_0 = 0 for trivial coefficients
    for k in range(1,4):
        src=list(combinations(range(m),k+1)); tgt={t:i for i,t in enumerate(combinations(range(m),k))}
        rows=[]
        for s in src:
            v=[F(0)]*len(tgt)
            for i in range(len(s)):
                for j in range(i+1,len(s)):
                    rest=[x for t,x in enumerate(s) if t not in (i,j)]
                    for c,coef in enumerate(br[s[i],s[j]]):
                        if coef==0 or c in rest: continue
                        lst=[c]+rest; sign=(-1)**(i+j)
                        # sort with sign
                        perm=sorted(range(len(lst)),key=lambda t:lst[t]); inv=sum(1 for x in range(len(perm)) for y in range(x+1,len(perm)) if perm[x]>perm[y])
                        v[tgt[tuple(sorted(lst))]]+=sign*(-1)**inv*coef
            rows.append(v)
        ranks.append(rank(rows) if rows else 0)
    return ranks
print("sl2 ranks d0..d3:", ce_ranks(2), "(paper/support: (0,3,0,0))")
print("sl3 ranks d0..d3:", ce_ranks(3), "(paper/support: (0,8,20,35))")
