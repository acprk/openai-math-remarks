# R4: every p-size inequality written in the p=1009 manuscript as a predicate in p. Run: /usr/bin/python3 main/ineq308.py (needs sympy)
# Expected output: expected_output/ineq308.txt
#!/usr/bin/env python3
"""Every p-size inequality written in the V(4)@1009 tex, as an exact predicate in p.
Prints, for each, the least odd prime from which it holds for all larger primes up to 10^4, and its value at p=11."""
from fractions import Fraction as F
from sympy import primerange
u = lambda p, i: (p**i-1)//(p-1)
I = [
 ("01-framework:250", "2u4<p^4<u5", lambda p: 2*u(p,4) < p**4 < u(p,5)),
 ("02a-cobar:61",  "p-1>1 (carry lowers filtration)", lambda p: p-1 > 1),
 ("02-lower:72",   "max w (product n=3) < u4  [additive comparison]", lambda p: (2*p*p+4*p+6) + (2*p*p+4*p+13)//(2*p-3) < u(p,4)),
 ("02-lower:83-88","max w (lower tables) < p^3", lambda p: (2*p*p+4*p+6) + (2*p*p+4*p+13)//(2*p-3) < p**3),
 ("02-lower:109-111","3(p^2-p-1) > 2p^2+21p/2+6  (=> k<=2)", lambda p: p*p - F(27,2)*p - 9 > 0),
 ("02-lower:136",  "F0+k<=8<p", lambda p: 8 < p),
 ("02-lower:139-142","2p^2+4p+8<3p^2", lambda p: p*p-4*p-8 > 0),
 ("02-lower:154-155","[-14,3] contains no nonzero multiple of p", lambda p: p > 14),
 ("02-lower:256",  "rho<=8<q", lambda p: 8 < 2*p-2),
 ("03-may:33",     "q/(q-1) u4 < 4p^3 < p^4", lambda p: F(2*p-2, 2*p-3)*u(p,4) < 4*p**3 < p**4),
 ("03-may:54",     "p^2-9p-4>0 (=> k<=p+2)", lambda p: p*p-9*p-4 > 0),
 ("03-may:64-67",  "p^2-9p-3>0 (no position 3 for k>=2)", lambda p: p*p-9*p-3 > 0),
 ("03-may:129",    "H0<=n=5<p", lambda p: 5 < p),
 ("03-may:146",    "p+8<2p (carries<=1)", lambda p: p+8 < 2*p),
 ("03-may:152-153","4 < p-6 (some carry)", lambda p: 4 < p-6),
 ("03-may:161-163","2p-19 > p-2 (single high block)", lambda p: 2*p-19 > p-2),
 ("03-may:273",    "minor 3 != 0 mod p", lambda p: 3 % p != 0),
 ("03-may:325",    "F7-F* = -p^2-11p+9 < 0", lambda p: -p*p-11*p+9 < 0),
 ("03-may:409-410","2p-2 > L=2p-3 bound", lambda p: 2*p-2 > 2*p-3),
 ("03-may:435",    "4+f<=6<p-1", lambda p: 6 < p-1),
 ("03-may:440",    "[-3,9] strictly inside (-p,p)", lambda p: 9 < p),
 ("04-frob:262",   "p^2+p > p^2-p+1 (no E)", lambda p: p*p+p > p*p-p+1),
 ("04-frob:337",   "3p+2 < p(p-1)", lambda p: 3*p+2 < p*(p-1)),
 ("04-frob:464",   "N0=(p-1)^2>=3", lambda p: (p-1)**2 >= 3),
 ("04-frob:501-502","6 != 0 mod p", lambda p: 6 % p != 0),
 ("04-frob:554",   "2 != 0 mod p", lambda p: 2 % p != 0),
 ("05-top:33",     "W=2p^2+4p+1<u4, u4+p<u5", lambda p: 2*p*p+4*p+1 < u(p,4) and u(p,4)+p < u(p,5)),
 ("05-top:142-149","w<=W+1+jmax < p^3 < u4", lambda p: 2*p*p+4*p+2+(2*p*p+4*p-4)//(2*p-3) < p**3 < u(p,4)),
 ("05-top:166-169","3(p(p-1)-1) > W+1 (=> j<=2)", lambda p: 3*(p*(p-1)-1) > 2*p*p+4*p+2),
 ("05-top:181",    "w < 3p^2", lambda p: 2*p*p+4*p+2+2 < 3*p*p),
 ("05-top:187-188","[-3,8] contains no nonzero multiple of p", lambda p: p > 8),
 ("05-top:256",    "rho... (j=0 source m=1: 1+2<4)", lambda p: True),
 ("05-top:839",    "qp+1-S=2p-6>0, N-K=p-3>0", lambda p: 2*p-6 > 0 and p-3 > 0),
]
P = list(primerange(3, 10000))
for tag, desc, f in I:
    thr = next(p for p in P if all(f(r) for r in P if r >= p))
    print(f"{tag:18s} {desc:52s} holds for all p>={thr:<3d}  at p=11: {bool(f(11))}")
print("literal threshold (all as written):", max(next(p for p in P if all(f(r) for r in P if r >= p)) for _, _, f in I))
