# Usage (from scripts/): /usr/bin/python3 main/hw319_pmin.py   -- exact Fraction arithmetic, ~10 s, no third-party packages.
# Expected output: expected_output/hw319_pmin.txt
# 319 Hahn-Wilson height 2: minimal prime for operator.tex conditions, exact arithmetic.
# Conditions (operator.tex line refs):
#  C1 (l.322-326 op:correction-bound, kappa<=p/(p-1) l.91): (p+1)/(2M) - p(M+1)/(p-1) > 1
#  CB (l.503 beta=log_p B0 > 1/3, needed obstruction.tex l.1082/1123): B0^3 > p, B0=2R0
#  CS (l.590-597 same-cluster diff < X/q): (2R0 p/M + 2)(p+1) <= p^2   [paper: p/4+3 with R0=M/8]
#  CE (l.507 I meets neither endpoint): 1 <= s0-R0+1 and s0+R0 <= M-2  => R0 <= (M-3)/2
#  CBand (l.603-606 first band ±2/p, separation from k1=2, b>=X/(2M)):
#     dev = 1/p + (2R0 p/M+2)/p^2 <= 2/p ; 1/M+2/p < 2/M-2/p ; 1/M-2/p >= 1/(2M)
#  CJ (l.692-694 'within three places of s0 ... inside I'): R0 >= 4 (literal), R0>=2 (sharpened |i-s0|<=1)
#  CD base stage l.546-562: d=floor(p/M) >= p/(2M), (p+1)d < p^2, d <= (p-1)/2 ; M<p distinct centers
#  global: p>8 (comparison.tex l.19,41), p>5 (geometry.tex l.14)
from fractions import Fraction as F
def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def primerange(a,b):
    for n in range(a,b):
        if isprime(n): yield n
def nextprime(n):
    n+=1
    while not isprime(n): n+=1
    return n
def ok(p,M,R0,Rmin):
    if M%2==0 or M>=p or R0<Rmin: return False
    s0=(M-1)//2
    if not (s0-R0+1>=1 and s0+R0<=M-2): return False
    if not (F(p+1,2*M)-F(p*(M+1),p-1) > 1): return False
    if not ((2*R0)**3 > p): return False
    if not ((F(2*R0*p,M)+2)*(p+1) <= p*p): return False
    dev=F(1,p)+(F(2*R0*p,M)+2)/p**2
    if not dev<=F(2,p): return False
    if not (F(1,M)+F(2,p) < F(2,M)-F(2,p)): return False
    if not (F(1,M)-F(2,p) >= F(1,2*M)): return False
    d=p//M
    if not (2*M*d>=p and (p+1)*d<p*p and d<=(p-1)//2): return False
    cent=[ (2*i+1)*p//(2*M) for i in range(M)]
    if len(set(cent))<M or cent[s0]!=(p-1)//2: return False
    return True
def paper(p,M):  # original parameters (op:parameters l.265-266)
    return M%2==1 and M>=81 and p>16*M*M and (2*(M//8))**3>p
for Rmin,tag in [(4,'literal R0>=4'),(2,'sharpened R0>=2')]:
    for p in primerange(11,10**5):
        sols=[(M,R0) for M in range(3,p,2) for R0 in range(Rmin,(M-3)//2+1) if ok(p,M,R0,Rmin)]
        if sols:
            print(tag,'p_min =',p,'witnesses',sols[:6]); break
# original mechanism minimum prime:
a=129; M=8*a+1
print('paper: need a=floor(M/8) with 8a^3 > 1024a^2 => a>=129; M=',M,'16M^2=',16*M*M,'(2a)^3=',(2*a)**3)
ps=[p for p in primerange(16*M*M+1,(2*a)**3) if paper(p,M)]
print('paper p_min =',ps[0] if ps else None,' count in window',len(ps))
# check no smaller prime works for paper conditions with any M
best=None
for M in range(81,1200,2):
    lo=16*M*M+1; hi=(2*(M//8))**3
    if hi>lo:
        q=nextprime(lo-1)
        if q<hi and (best is None or q<best): best=q
print('paper p_min over all M:',best, ' 16^6=',16**6)

