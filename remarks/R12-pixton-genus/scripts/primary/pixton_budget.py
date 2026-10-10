# Exact integer re-optimisation of the colour budget in 053 (Pixton counterexample), Section 8.
# Inequalities exactly as in 08-support.tex (lines 326-510) and Lemma 6.x (large-vertex-test), parameters g,K,m0,e,nblocks.
from math import comb
def need_m0(g,K):
    CK=comb(2*K,3); MK=max(10*K+3,CK+1)
    loss=3*g*MK + g*(10*K+2) + 2*K*CK      # paper's bound on lost internal entries (eq. small-loss)
    tb=25*(g-1)**2-1                       # t_b <= 25(b-1)^2-1 with b<=g
    # need loss + t_b <= d_e/2,  d_e=binom(e,3), e = largest even <= m0/2 (also >= m0/2-2 automatically)
    lo,hi=8,1
    def ok(m0):
        e=(m0//2)//2*2
        return e>=4 and 2*(loss+tb)<=comb(e,3)
    while not ok(hi): hi*=2
    lo=hi//2
    while lo+1<hi:
        mid=(lo+hi)//2
        if ok(mid): hi=mid
        else: lo=mid
    return hi
def feasible(g,K,nblocks_rule):
    if K<5: return False
    CK=comb(2*K,3)
    if not (CK+1<11*g): return False
    if not (4*(g-K-2)>3*g-3): return False
    if not (3*g-2<(g-K)**2-1): return False     # dim S < (b-1)^2-1 with b>=g-K+1
    m0=need_m0(g,K)
    nb = 4*K+20 if nblocks_rule=='paper' else 3*K+5
    S0=g-nb*m0
    if S0<=0: return False
    if not (S0*(g*(K-4)+2) > 2*140*g*g): return False   # large-deficit contradiction
    return (m0,nb,S0)
def mingK(K,rule):
    lo,hi=1,2
    while not feasible(hi,K,rule): hi*=2
    lo=hi//2
    while lo+1<hi:
        mid=(lo+hi)//2
        if feasible(mid,K,rule): hi=mid
        else: lo=mid
    return hi
# sanity: paper's parameters
g=10**60; K=10**5
print("paper point feasible:",bool(feasible(g,K,'paper')))
for rule in ('paper','3K+5'):
    best=None
    for K in list(range(290,2000,10)):
        gK=mingK(K,rule)
        if best is None or gK<best[0]: best=(gK,K)
    gK,K=best
    # refine K locally
    for K2 in range(K-10,K+11):
        g2=mingK(K2,rule)
        if g2<best[0]: best=(g2,K2)
    gK,K=best
    m0,nb,S0=feasible(gK,K,rule)
    print(rule,"min g = %d (~%.3e) at K=%d m0=%d blocks=%d |S0|=%d ; g-1 feasible? %s"%(gK,gK,K,m0,nb,S0,bool(feasible(gK-1,K,rule))))
