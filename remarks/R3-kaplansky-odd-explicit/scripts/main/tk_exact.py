# R3 main computation. Run from scripts/: /usr/bin/python3 main/tk_exact.py  (needs sympy; <1 s)
# Expected output: expected_output/tk_exact.txt
from math import comb
from fractions import Fraction
from sympy import isprime, nextprime
N,k,s=86,41,13
m=comb(N,k); t=sum(comb(N,j) for j in range(0,k-s+1))
assert 3*k-3*(s-1)>N            # triangle-free
assert Fraction(t,m)<Fraction(1,96)
print("m =",m," t =",t," t/m =",float(Fraction(t,m))," 1/96 =",1/96)
p=m+1
while not (p%4==1 and isprime(p)): p+=1
print("least prime p>m with p≡1 mod 4:",p," digits",len(str(p)))
# exact check that the fit polynomial values vanish for every nonadjacent intersection size and =1 on diagonal
for a in range(0,k+1):
    val=Fraction(1)
    for kk in range(s,k): val*=Fraction(a-kk,k-kk)
    if a==k: assert val==1
    elif a>=s: assert val==0
print("fit identities ok; max degree Delta =",sum(comb(k,j)*comb(N-k,k-j) for j in range(s)))
# paper comparison
M0=comb(1200,600); print("paper m digits:",len(str(M0)))
