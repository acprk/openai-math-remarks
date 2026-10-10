# Usage (from scripts/): /usr/bin/python3 main/hw319_detail.py  -- per-M breakdown and the C1v variant; ~2 s.
# Expected output: expected_output/hw319_detail.txt
import sys; sys.argv=['x']
import os; src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'hw319_pmin.py')).read().split("for Rmin,tag in [(4,'literal")[0]
exec(src)
def c1(p,M): return F(p+1,2*M)-F(p*(M+1),p-1) > 1
# variant: lemma op:coefficient-control hypothesis b>=X/(2M) replaced by b>=X(1/M-2/p) (first band lower end) and base depth d=floor(p/M)
def c1v(p,M): return (p+1)*(F(1,M)-F(2,p)) - F(p*(M+1),p-1) > 1 and F((p+1)*(p//M),p) - F(p*(M+1),p-1) > 1
for M in (7,9,11,13,15):
    t=next(p for p in primerange(11,10**6) if c1(p,M)); tv=next(p for p in primerange(11,10**6) if c1v(p,M))
    print('M=%2d Rmax=%d B0max=%2d B0^3=%4d  C1 first prime=%d  C1v first prime=%d'%(M,(M-3)//2,2*((M-3)//2),(2*((M-3)//2))**3,t,tv))
for p,M,R0 in [(293,11,4),(199,9,3)]:
    print(p,M,R0,'CS(+3):',(F(2*R0*p,M)+3)*(p+1)<=p*p,' dev<=2/p(+3):',F(1,p)+(F(2*R0*p,M)+3)/p**2<=F(2,p),
          ' C1 margin:',F(p+1,2*M)-F(p*(M+1),p-1)-1, ' B0^3-p:',(2*R0)**3-p)
# paper-original params at p=293 fail where?
print('paper at p=293: need M>=81 and p>16M^2 =',16*81*81)
# full search with C1 replaced by C1v
_c1_orig = ok
def okv(p,M,R0,Rmin):
    # same as ok but C1 -> C1v
    if M%2==0 or M>=p or R0<Rmin: return False
    s0=(M-1)//2
    if not (s0-R0+1>=1 and s0+R0<=M-2): return False
    if not c1v(p,M): return False
    if not ((2*R0)**3 > p): return False
    if not ((F(2*R0*p,M)+3)*(p+1) <= p*p): return False
    if not F(1,p)+(F(2*R0*p,M)+3)/p**2<=F(2,p): return False
    if not (F(1,M)+F(2,p) < F(2,M)-F(2,p)): return False
    if not (F(1,M)-F(2,p) >= F(1,2*M)): return False
    d=p//M
    if not (2*M*d>=p and (p+1)*d<p*p and d<=(p-1)//2): return False
    cent=[(2*i+1)*p//(2*M) for i in range(M)]
    return len(set(cent))==M and cent[s0]==(p-1)//2
for Rmin,tag in [(4,'C1v literal R0>=4'),(2,'C1v sharpened R0>=2')]:
    for p in primerange(11,10**4):
        sols=[(M,R0) for M in range(3,p,2) for R0 in range(Rmin,(M-3)//2+1) if okv(p,M,R0,Rmin)]
        if sols: print(tag,'p_min =',p,sols[:4]); break
