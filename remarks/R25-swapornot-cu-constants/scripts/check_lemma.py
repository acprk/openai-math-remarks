# falsification test: f_lambda >= f_mu * N0(L,k) for lambda=(L,mu), all lambda |- m <= 22 (exact integers)
from math import comb, factorial
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest
def dim(lam):
    if not lam: return 1
    conj=[sum(1 for r in lam if r>j) for j in range(lam[0])]
    p=1
    for i,r in enumerate(lam):
        for j in range(r): p*=r-j-1+conj[j]-i
    return factorial(sum(lam))//p
def N0(L,k):
    j=min(k,L); return comb(L+j,j)-(comb(L+j,j-1) if j else 0)
bad=0; tight=0; tot=0
for m in range(1,23):
    for lam in partitions(m):
        tot+=1
        L=lam[0]; mu=lam[1:]; k=m-L
        lhs=dim(lam); rhs=dim(mu)*N0(L,k)
        if lhs<rhs: bad+=1; print("VIOLATION",lam,lhs,rhs)
        if lhs==rhs: tight+=1
print("checked",tot,"violations",bad,"equality cases",tight)
