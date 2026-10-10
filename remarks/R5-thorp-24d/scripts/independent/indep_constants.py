#!/usr/bin/env python3
"""Independent re-derivation (R7 adversarial review of R5) of the Thorp constants.
Does NOT reuse R5's constants_certificate.py.

kappa(b) = 2 / log2(2b/(2b-1)),  e(b) = 1 - 2/b  (integer b >= 3; b<=2 has e<=0 and positive cost).
Route H: minimise sum kappa(b_i) s.t. sum e(b_i) > 1 ; Route O: s.t. sum e(b_i) > 2.
Admissibility exact (Fraction); costs in mpmath.iv (outward rounded).
Completeness: L(b)=log2(1+1/(2b-1)) <= 1/((2b-1) ln2)  =>  kappa(b) >= 2(2b-1) ln2,
and kappa is increasing in b, kappa(b) >= kappa(3) for every part.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import mpmath
from mpmath import iv
iv.dps = 40; mpmath.mp.dps = 40
lo = lambda x: mpmath.mpf(x.a)
hi = lambda x: mpmath.mpf(x.b)

def kap(b):  # interval
    return iv.mpf(2) / (iv.log(iv.mpf(2*b) / iv.mpf(2*b-1)) / iv.log(2))
def e(b): return 1 - F(2, b)

def search(tau, Ubound):
    # Ubound: rigorous upper bound for the optimum (cost of some admissible multiset)
    k3lo = lo(kap(3))
    rmax = int(mpmath.floor(Ubound / k3lo))
    # part b allowed only if 2(2b-1)ln2 <= Ubound (lower bound for kappa(b))
    bmax = 3
    while 2*(2*(bmax+1)-1)*mpmath.log(2) <= Ubound: bmax += 1
    res = []
    for r in range(1, rmax+1):
        for ms in combinations_with_replacement(range(3, bmax+1), r):
            if sum(e(b) for b in ms) > tau:
                c = sum((kap(b) for b in ms), iv.mpf(0))
                res.append((c, ms))
    res.sort(key=lambda t: lo(t[0]))
    return res, rmax, bmax

for name, tau, ub in [("H", 1, 23.6), ("O", 2, 44.31)]:
    res, rmax, bmax = search(tau, ub)
    best, second = res[0], res[1]
    assert hi(best[0]) < lo(second[0]), "optimum not separated"
    # every multiset not enumerated has cost > ub >= best upper end
    assert hi(best[0]) <= ub
    print(f"route {name}: parts<= {bmax}, r<= {rmax}, #admissible enumerated={len(res)}")
    print(f"  optimum {best[1]}: [{mpmath.nstr(lo(best[0]),15)}, {mpmath.nstr(hi(best[0]),15)}]")
    print(f"  runner-up {second[1]}: >= {mpmath.nstr(lo(second[0]),10)}")
    for c, ms in res[:5]: print("   ", ms, mpmath.nstr(lo(c), 10), "sum e =", sum(e(b) for b in ms))
    # sigma window: need sum(1-(2+s)/b) > tau  <=>  s < (r - tau)/sum(1/b) - 2
    ms = best[1]; S = sum(F(1, b) for b in ms)
    print("  sigma window: sigma <", F(len(ms) - tau) / S - 2)

# fractional lower bound min_b kappa(b)/e(b) over real b>2
f = lambda b: (2/ mpmath.log(2*b/(2*b-1), 2)) / (1 - 2/b)
bstar = mpmath.findroot(lambda b: mpmath.diff(f, b), 3.9)
print("real minimiser b* =", mpmath.nstr(bstar, 8), " min kappa/e =", mpmath.nstr(f(bstar), 10))
# rigorous: check on integer b (the relevant set) min kappa/e
print("integer min kappa/e:", min(lo(kap(b)/iv.mpf(e(b).numerator)*e(b).denominator) for b in range(3, 60)))

# ledger (route H)
def L2(x): return iv.log(x)/iv.log(2)
steps = []
L = L2(iv.mpf(32)/31); steps.append(("paper 8x200d", iv.mpf(1600)))
steps.append(("paper inequalities minimal T, 8 factors", 8*(2 + 3/L)))
steps.append(("two HS factors, b=8", 2*(2 + 3/L)))
L = L2(iv.mpf(16)/15); steps.append(("theta->0, beta=1/8", 2*(2 + 3/L)))
steps.append(("KL chain", 2*(2 + 2/L)))
steps.append(("look-back", 2*(2/L)))
steps.append(("blocks {4,5}", kap(4) + kap(5)))
for nme, v in steps: print(f"  {nme:45s} {mpmath.nstr(lo(v),8)}")
# route O ledger key points
L = L2(iv.mpf(32)/31)
print("route O: paper ineq, 5 factors (Lemma5.1 sigma=1, b=8: need 5*(1-3/8)>2 -> r>=4? check)")
for r in range(2, 9):
    # op exponent with b=8, sigma=1: D^2 prod D^{-1+3/8} -> 2 - r*5/8 <0
    print("   r=",r," exponent 2 - r*5/8 =", F(2) - r*F(5,8))
print("  5 factors x (2+3/L) =", mpmath.nstr(lo((5*(2+3/L))), 8))
L = L2(iv.mpf(16)/15)
print("  theta->0:", mpmath.nstr(lo((5*(2+3/L))),8), " KL:", mpmath.nstr(lo((5*(2+2/L))),8), " look-back:", mpmath.nstr(lo((5*(2/L))),8), " 3 factors b=8 sigma->0:", mpmath.nstr(lo((3*(2/L))),8))
