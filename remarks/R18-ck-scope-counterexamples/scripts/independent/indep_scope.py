#!/usr/bin/env python3
"""Independent re-implementation (does not import the primary scripts).

Setting: a finite probability space (Omega, mu), n = 1, a set A with mu(A) = w, and the Boolean
function f = 2*1_A - 1 with mean m = 2w - 1.  Resampling channel R_lam: Y = X with probability lam,
otherwise Y is a fresh independent sample of mu.  Then E[f(X) | Y] = lam*f(Y) + (1-lam)*m.
For mu uniform on {-1,1} (w = 1/2) this is the binary symmetric channel with correlation rho = lam.
Examples with w != 1/2: Omega = [q] uniform, |A| = k (w = k/q); or the p-biased bit (w = p).

  I(f(X);Y) = H(m) - E H(lam f + (1-lam) m)           (nats; H(r) = entropy of a sign of mean r)
  Hellinger loss  L = sqrt(1-m^2) - E sqrt(1-(lam f + (1-lam) m)^2)

Claims checked (they concern ONLY these non-uniform-binary settings; family 119 does not claim them):
 (A) refined Boolean bound  I <= psi(lam * r_m),  r_m = psi^{-1}(H(m)),  psi = log 2 - H,
     fails for small lam whenever w != 1/2 (this is also Theorem 1.1 with u = f, since I(f) = H(m));
 (B) Hellinger bound  L <= 1 - sqrt(1 - lam^2)  fails for small lam whenever w != 1/2;
 (C) the plain Courtade-Kumar curve  I <= psi(lam)  is NOT violated by these single-coordinate
     examples (consistent with it being open, not refuted, in the general setting).

Part 1: rigorous interval-arithmetic (mpmath.iv) certificates of three witnesses.
Part 2: Taylor coefficients in lam (exact, sympy).
Part 3: high-precision scan of the violation window lam in (0, lam*(w)).
Part 4: control - f = AND on the uniform 2-cube under product bit flips satisfies both bounds.
Environment: python3, mpmath >= 1.3, sympy >= 1.10. Single thread, < 1 min.
"""
import mpmath as mp
from mpmath import iv
import sympy as sp

iv.dps = 50
mp.mp.dps = 50

# ---------------- Part 1: interval certificates ----------------
def H_iv(r):
    a = (1 + r) / 2; b = (1 - r) / 2
    return -(a * iv.log(a) + b * iv.log(b))

def psi_iv(r):
    return iv.log(2) - H_iv(r)

def psi_inv_upper(target_iv):
    """rigorous upper bound r_hi in (0,1) with psi(r_hi) >= target (psi is increasing on [0,1))."""
    lo, hi = mp.mpf(0), mp.mpf(1) - mp.mpf(10) ** -40
    for _ in range(160):
        mid = (lo + hi) / 2
        if psi_iv(iv.mpf(mid)).a > target_iv.b:     # psi(mid) certainly > target
            hi = mid
        else:
            lo = mid
    assert psi_iv(iv.mpf(hi)).a > target_iv.b
    return iv.mpf(hi)

def certify(label, w, lam):
    w = iv.mpf(w[0]) / w[1]; lam = iv.mpf(lam[0]) / lam[1]      # exact rationals, enclosed
    m = 2 * w - 1
    up = lam * 1 + (1 - lam) * m          # posterior mean when f(Y) = +1
    dn = lam * (-1) + (1 - lam) * m
    Hm = H_iv(m)
    I = Hm - (w * H_iv(up) + (1 - w) * H_iv(dn))
    r_hi = psi_inv_upper(Hm)
    refined = psi_iv(lam * r_hi)          # >= psi(lam r_m) since psi increasing
    ck = psi_iv(lam)
    hel = iv.sqrt(1 - m * m) - (w * iv.sqrt(1 - up * up) + (1 - w) * iv.sqrt(1 - dn * dn))
    hel_rhs = 1 - iv.sqrt(1 - lam * lam)
    lo = lambda x: mp.nstr(mp.mpf(x.a.a), 12)
    hi = lambda x: mp.nstr(mp.mpf(x.b.b), 12)
    print(f"[{label}] w={hi(w)} lam={hi(lam)}  m={hi(m)}")
    print(f"   I(f;Y)            in [{lo(I)}, {hi(I)}]")
    print(f"   psi(lam*r_m)      <= {hi(refined)}   -> refined bound (A) violated: {I.a > refined.b}")
    print(f"   psi(lam) (CK)     >= {lo(ck)}   -> plain CK curve respected: {I.b < ck.a}")
    print(f"   Hellinger loss    in [{lo(hel)}, {hi(hel)}]")
    print(f"   1-sqrt(1-lam^2)   in [{lo(hel_rhs)}, {hi(hel_rhs)}]  -> (B) violated: {hel.a > hel_rhs.b}")

print("Part 1: interval certificates")
certify("q=4 uniform, |A|=1  (= 1/4-biased bit)", (1, 4), (1, 10))
certify("q=3 uniform, |A|=1", (1, 3), (1, 10))
certify("p-biased bit, p=0.3", (3, 10), (1, 5))

# ---------------- Part 2: Taylor coefficients ----------------
print("Part 2: Taylor expansions in lam (exact)")
lam, m = sp.symbols('lam m', real=True)
Hs = lambda r: -((1 + r) / 2 * sp.log((1 + r) / 2) + (1 - r) / 2 * sp.log((1 - r) / 2))
w = (1 + m) / 2
I = Hs(m) - (w * Hs(lam + (1 - lam) * m) + (1 - w) * Hs(-lam + (1 - lam) * m))
cI = sp.simplify(sp.series(I, lam, 0, 3).removeO().coeff(lam, 2))
L = sp.sqrt(1 - m**2) - (w * sp.sqrt(1 - (lam + (1 - lam) * m)**2) + (1 - w) * sp.sqrt(1 - (-lam + (1 - lam) * m)**2))
cL = sp.simplify(sp.series(L, lam, 0, 3).removeO().coeff(lam, 2))
r = sp.symbols('r', positive=True)
psi = sp.log(2) - Hs(lam * r)
cpsi = sp.simplify(sp.series(psi, lam, 0, 3).removeO().coeff(lam, 2))
crhs = sp.simplify(sp.series(1 - sp.sqrt(1 - lam**2), lam, 0, 3).removeO().coeff(lam, 2))
print("   I(f;Y)            = (", cI, ") lam^2 + O(lam^3)")
print("   psi(lam r)        = (", cpsi, ") lam^2 + O(lam^4)")
print("   Hellinger loss    = (", cL, ") lam^2 + O(lam^3)")
print("   1-sqrt(1-lam^2)   = (", crhs, ") lam^2 + O(lam^4)")
print("   => (A) fails for small lam iff r_m < 1 iff H(m) < log 2 iff m != 0;")
print("      (B) fails for small lam iff 1/sqrt(1-m^2) > 1 iff m != 0.")

# ---------------- Part 3: violation windows ----------------
print("Part 3: largest lam with a violation (high precision root finding), for w = mu(A)")
def Hm_(x): return -((1 + x) / 2 * mp.log((1 + x) / 2) + (1 - x) / 2 * mp.log((1 - x) / 2))
def psi_(x): return mp.log(2) - Hm_(x)
def psi_inv_(t): return mp.findroot(lambda x: psi_(x) - t, (mp.mpf('1e-30'), 1 - mp.mpf('1e-30')), solver='bisect')
for wv in [mp.mpf(1) / 4, mp.mpf(1) / 3, mp.mpf(2) / 5, mp.mpf(1) / 10]:
    mm = 2 * wv - 1
    rm = psi_inv_(Hm_(mm))
    def gA(l):
        return Hm_(mm) - (wv * Hm_(l + (1 - l) * mm) + (1 - wv) * Hm_(-l + (1 - l) * mm)) - psi_(l * rm)
    def gB(l):
        return (mp.sqrt(1 - mm**2) - (wv * mp.sqrt(1 - (l + (1 - l) * mm)**2) + (1 - wv) * mp.sqrt(1 - (-l + (1 - l) * mm)**2))
                - (1 - mp.sqrt(1 - l**2)))
    def last_pos(g):
        grid = [mp.mpf(k) / 1000 for k in range(1, 1000)]
        last = None
        for a, b in zip(grid, grid[1:]):
            if g(a) > 0 and g(b) <= 0:
                last = mp.findroot(g, (a, b), solver='bisect')
        return last
    lA, lB = last_pos(gA), last_pos(gB)
    print(f"   w={mp.nstr(wv,4)}: (A) violated for 0<lam<{mp.nstr(lA,8)} ;  (B) violated for 0<lam<{mp.nstr(lB,8)}")

# ---------------- Part 4: control (product bit-flip noise) ----------------
print("Part 4: control - f = AND on the uniform {-1,1}^2 cube (same law as q=4, |A|=1) under product bit flips")
import itertools
rho = mp.mpf(1) / 10
pts = list(itertools.product([-1, 1], repeat=2))
f = {x: (1 if x == (1, 1) else -1) for x in pts}
m = mp.mpf(sum(f.values())) / 4
def Tf(y):
    tot = 0
    for x in pts:
        w = 1
        for xi, yi in zip(x, y):
            w *= (1 + rho * xi * yi) / 2
        tot += w * f[x]
    return tot
Ivals = [Tf(y) for y in pts]
I_and = Hm_(m) - sum(Hm_(v) for v in Ivals) / 4
rm = psi_inv_(Hm_(m))
Lh = mp.sqrt(1 - m**2) - sum(mp.sqrt(1 - v**2) for v in Ivals) / 4
print(f"   rho=0.1: I = {mp.nstr(I_and,6)} <= psi(rho r_m) = {mp.nstr(psi_(rho*rm),6)} : {I_and <= psi_(rho*rm)}")
print(f"            Hellinger loss = {mp.nstr(Lh,6)} <= 1-sqrt(1-rho^2) = {mp.nstr(1-mp.sqrt(1-rho**2),6)} : {Lh <= 1-mp.sqrt(1-rho**2)}")
