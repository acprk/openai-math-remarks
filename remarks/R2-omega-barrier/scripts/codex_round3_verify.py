# Verification program extracted VERBATIM from proofs/CODEX-ROUND3.md
# Author: OpenAI Codex session (model gpt-6-astra, Codex CLI v0.157.0), 2026-10-09,
# in a Claude<->Codex collaboration round directed by [AUTHOR NAME].
# Environment: /usr/bin/python3 with sympy (round 1), stdlib (round 2), numpy+scipy+mpmath (round 3).
# Expected output: expected_output/codex_round3_verify.txt
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import log, comb
import sys
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.iv.dps=40
iv=mp.iv
DEN=10**9
stats=dict(states=0,partitions=0,optimizations=0,nonunique=0)
records={}
@lru_cache(None)
def phi(a):return Q(1) if a==1 else phi(a-1)*Q(3*a-2,3*a-3)
def F(a,b):
 a,b=sorted((a,b));return phi(a)*(b+Q(a-1,2))
def I(q):
 q=Q(q);return iv.mpf(q.numerator)/q.denominator
def dyadic(t):
 sign,man,e,bc=t
 assert bc>=0
 return (-1 if sign else 1)*Q(man)*(Q(2)**e)
def enclose(v):return tuple(dyadic(t) for t in v._mpi_)
@lru_cache(None)
def canonical(a,b,c,s):
 ds=(a,b,c);lo=[max(0,s-(sum(ds)-ds[i]-2)) for i in range(3)]
 hi=[min(d-1,s) for d in ds]
 if any(l>h for l,h in zip(lo,hi)):return None
 ds=tuple(sorted(h-l+1 for l,h in zip(lo,hi)));s-=sum(lo)
 return (*ds,min(s,sum(ds)-3-s))
def blocks(T,cuts):
 a,b,c,s=T;D=(a,b,c)
 out={}
 for bits in product(range(2),repeat=3):
  lows=[cuts[j] if bits[j] else 0 for j in range(3)]
  lens=[D[j]-cuts[j] if bits[j] else cuts[j] for j in range(3)]
  B=canonical(*lens,s-sum(lows))
  if B:out[bits]=B
 return out
@lru_cache(None)
def faces(S):
 ans=set()
 for w in product((-2,-1,1,2),repeat=3):
  scores=[sum(x*y for x,y in zip(w,s)) for s in S];m=min(scores)
  A=tuple(s for s,v in zip(S,scores) if v==m)
  ans.add(A)
 return tuple(sorted(ans))
def exact_rank(A):
 A=[[Q(x) for x in row] for row in A];r=0
 for j in range(len(A[0])):
  piv=next((i for i in range(r,len(A)) if A[i][j]),None)
  if piv is None:continue
  A[r],A[piv]=A[piv],A[r];d=A[r][j];A[r]=[x/d for x in A[r]]
  for i in range(r+1,len(A)):
   d=A[i][j];A[i]=[x-d*y for x,y in zip(A[i],A[r])]
  r+=1
 return r
def entropy_solver(S,values):
 stats['optimizations']+=1
 A=np.array([[int(s[j]==k) for s in S] for j in range(3) for k in sorted(set(s[j] for s in S))],dtype=float)
 # Check injectivity of the marginal map by exact rational elimination.
 if exact_rank(A.astype(int).tolist())!=len(S):
  stats['nonunique']+=1;raise ValueError(S)
 l=np.array([log(float((a+b)/2)) for a,b in values])
 def fun(p):
  q=A@p;return sum(x*log(x) for x in q if x>0)/3-l@p
 def jac(p):return A.T@(np.log(np.maximum(A@p,1e-150))+1)/3-l
 k=len(S)
 if k==1:p=np.ones(1)
 else:
  r=minimize(fun,np.full(k,1/k),jac=jac,bounds=[(1e-12,1)]*k,constraints={'type':'eq','fun':lambda p:p.sum()-1,'jac':lambda p:np.ones(k)},method='SLSQP',options={'ftol':1e-12,'maxiter':150})
  p=r.x
 counts=[max(1,round(float(v)*DEN)) for v in p];total=sum(counts);P=[Q(v,total) for v in counts]
 M=[[sum(P[i] for i,s in enumerate(S) if s[j]==k) for k in sorted(set(s[j] for s in S))] for j in range(3)]
 H=-sum(I(q)*iv.ln(I(q)) for row in M for q in row)/3
 LB=iv.exp(H+sum(I(p)*iv.ln(I(v[0])) for p,v in zip(P,values)))
 ub=[]
 for i,s in enumerate(S):
  prodq=Q(1)
  for j in range(3):prodq*=M[j][sorted(set(v[j] for v in S)).index(s[j])]
  ub.append(I(values[i][1])/iv.exp(iv.ln(I(prodq))/3))
 lower=enclose(LB)[0];upper=max(enclose(v)[1] for v in ub)
 return lower,upper
@lru_cache(None)
def solve(T,mode='character'):
 stats['states']+=1
 a,b,c,s=T
 if a==1:return enclose(iv.exp(iv.ln(I(b**3))/4))
 seeds=[(Q(1),Q(1))]
 if c==a+b-1:seeds.append(enclose(iv.exp(iv.ln(I(a*b*c))/3)))
 if a==b==c and s==a-1 and a>=3:
  CW=((2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0))
  dot=enclose(iv.exp(iv.ln(I((a-2)**3))/4))
  seeds.append(entropy_solver(CW,[(Q(1),Q(1))]*3+[dot]*3))
 if c==a+b-1 and mode=='character':seeds.append(enclose(iv.exp(iv.ln(I(F(a,b)**3))/4)))
 bestlo=max(v[0] for v in seeds);besthi=max(v[1] for v in seeds)
 # Consecutive interval halves; all floor/ceil choices on all three legs.
 cutchoices=[sorted(set((d//2,(d+1)//2))) for d in (a,b,c)]
 seen=set();newlo=Q(1);newhi=Q(1)
 for cuts in product(*cutchoices):
  bs=blocks(T,cuts)
  for S in faces(tuple(bs)):
   kids=tuple(bs[x] for x in S)
   key=(S,kids)
   if key in seen:continue
   seen.add(key);stats['partitions']+=1
   if len(S)==1:L,U=solve(kids[0],mode)
   else:L,U=entropy_solver(S,[solve(k,mode) for k in kids])
   bestlo=max(bestlo,L);besthi=max(besthi,U);newlo=max(newlo,L);newhi=max(newhi,U)
 records[(T,mode)]=(newlo,newhi)
 return bestlo,besthi
def floor4(q):return f'{(q.numerator*10000//q.denominator)/10000:.4f}'
def ceil4(q):return floor4(Q(-((-q.numerator*10000)//q.denominator),10000))
def all_cuts():
 for n in range(2,13):
  T=(n,n,2*n-1,2*n-2);seen=set();L=U=Q(1);ncuts=0
  for cuts in product(range(1,n),range(1,n),range(1,2*n-1)):
   ncuts+=1;bs=blocks(T,cuts)
   for S in faces(tuple(bs)):
    kids=tuple(bs[x] for x in S);key=(S,kids)
    if key in seen:continue
    seen.add(key)
    if len(S)==1:l,u=solve(kids[0])
    else:l,u=entropy_solver(S,[solve(k) for k in kids])
    L=max(L,l);U=max(U,u)
  f=enclose(iv.exp(iv.ln(I(F(n,n)**3))/4))[0]
  assert U<f
  print(n,ncuts,len(seen),floor4(L),ceil4(U),floor4(f),flush=True)
 print(stats)

if __name__=='__main__':
 if '--all-cuts' in sys.argv:
  all_cuts();sys.exit(0)
 Z=((0,0,0),(0,1,1),(1,0,1),(1,1,2))
 for n in range(1,31):
  T=(n,n,2*n-1,2*n-2);C=solve(T);O=solve(T,'operational')
  fl,fu=enclose(iv.exp(iv.ln(I(F(n,n)**3))/4))
  new=records.get((T,'character'),(Q(1),Q(1)))
  if n==1:zlo=zhi=zupper=Q(1)
  else:
   q=n-1;D=enclose(iv.exp(iv.ln(I(q**3))/4))
   tr=solve((q,q,q,q-1),'operational')
   zlo,zhi=entropy_solver(Z,[(Q(1),Q(1)),D,D,tr])
   _,zupper=entropy_solver(Z,[(Q(1),Q(1)),D,D,(Q(q),Q(q))])
   assert new[1]<fl and zupper<fl
  print(n,floor4(max(O[0],zlo)),ceil4(new[1]),floor4(zlo),ceil4(zupper),floor4(fl),2*n-1,flush=True)
 print(stats)
 # All-size zero/rest proof: finite head and rational tail.
 for q in range(1,16):
  if q==1:R=Q(12)
  else:R=enclose(I(16)+I(Q(16,3))*iv.exp(iv.ln(I(q**9))/4)+I(Q(16,9)*q**3))[1]
  assert R**4<F(q+1,q+1)**9
 assert (Q(22,9)+Q(1,256))**4 < Q(27,8)**3
 print('zero/rest finite-head and all-size tail: PASS')

 for a,b in product(range(1,31),repeat=2):
  p={(i,j):Q(comb(a-1,i)*comb(b-1,j),
               (a+b-1)*comb(a+b-2,i+j))
     for i in range(a) for j in range(b)}
  assert all(sum(p[i,j] for j in range(b))==Q(1,a) for i in range(a))
  assert all(sum(p[i,j] for i in range(a))==Q(1,b) for j in range(b))
  assert all(sum(p[i,j] for i in range(a) for j in range(b) if i+j==k)
             ==Q(1,a+b-1) for k in range(a+b-1))
 print('Uniform marginals: 900 exact checks PASS')
 assert 2**8*3**5 > F(2,2)**9
 print('Quantum/group combined bound does not exclude F at n=2: PASS')
 for q in range(2,30):
  assert solve((q,q,q,q-1),'character')[1]<q
 print('Truncated character recipes q=2..29: all upper endpoints below q PASS')
