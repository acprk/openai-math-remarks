"""Constants obtainable if the dimension input is only what the openai/math Lean development proves
(Main.lean `critical_dimension_effective`: krullDim <= size^2 - floor(C(q+2,2)/18), needs 3q <= size-1),
fed into the Lean conditional pipeline (ConditionalBound.lean `selector_lower_of_critical_dimension`:
bdc >= (r-1)(d-1)/(4e) when B+(d-1) < k^2, m >= 11k-2q).
  lean    : D = floor(C(q+2,2)/18)            (Lean-certified input, parameters re-chosen)
  leanBOX : D = floor((q+1)^2/18)             (hypothetical: box lemma inserted into Lean's T/18 route; not formalized)
  lean-fixed : Lean's own choice k=m//11, q=(k-1)//3, d=floor(C(q+2,2)/18).
Usage: python3 lean_dimension_model.py   (~2 s, single thread)
"""
from math import comb
from fractions import Fraction
def P(model,k,q):
    r=k-2*q
    if r<2: return 0
    if model=='lean':
        if 3*q>k-1: return 0
        D=comb(q+2,2)//18
    elif model=='leanBOX':
        D=(q+1)**2//18
    d=min(D,k*k-1)
    if d<2: return 0
    return (r-1)*(d-1)
def best(m,model):
    bv,arg=0,None
    k=1
    while 11*k-2*(k//2)<=m:
        for q in range(k//2+1):
            if 11*k-2*q>m: continue
            v=P(model,k,q)
            if v>bv: bv,arg=v,(k,q)
        k+=1
    return bv,arg
for m in [1408,5000,20000]:
    for mod in ['lean','leanBOX']:
        v,a=best(m,mod); c=4*m**3/v
        print(m,mod,round(c),a,'stated/this=%.3f'%(5529600/c))
# Lean's own fixed choice: k=m//11, q=(k-1)//3, d=C(q+2,2)//18
for m in [1408,5000,20000]:
    k=m//11;q=(k-1)//3;d=comb(q+2,2)//18;r=k-2*q
    print('lean-fixed',m,round(4*m**3/((r-1)*(d-1))))
# asymptotic
# asymptotic: q=a*k, m=(11-2a)k, (r-1)(d-1) ~ (1-2a)*g(a)*k^3, g=a^2/36 (lean) or a^2/18 (leanBOX)
def asym(g,amax):
    best=min(4*(11-2*a)**3/((1-2*a)*g(a)) for a in [i/10**5*amax for i in range(1,10**5+1)] if a<0.5)
    return best
print('asym lean %.0f' % asym(lambda a:a*a/36,1/3), 'asym leanBOX %.0f' % asym(lambda a:a*a/18,0.4999))
