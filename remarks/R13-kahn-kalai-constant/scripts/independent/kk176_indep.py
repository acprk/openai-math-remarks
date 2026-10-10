#!/usr/bin/env sage -python
"""Independent re-implementation (adversarial review, 2026-10-09) of the constant
bookkeeping for OpenAI preprint 'The second Kahn--Kalai conjecture' (family 176).
Written from the paper, NOT from kk176_verify.py.  Arithmetic: Sage RealBallField
(Arb ball arithmetic, 256 bits) for transcendental quantities; exact integers /
rationals for binomial tails.

Generalised proof (paper constant in brackets):
  x = a/rho = e^{-L}                         [L = 50]
  local failure: Z < e^{-al m} or |T| > m/2  [al = 10]
  pruning threshold eps = e^{-be m}          [be = 1]
  capacity separation m_{j+1} >= g m_j       [g = 16]; hierarchy |S| < |I|/(g+1)
  sigma = lam q, lam = (g+1)^(g+1)/g^g       [paper: 128]
  K = prod_{v>=0} (1 - e^{-be 2^v})^{-1}     [paper: prod over all m >= 1, < 4]
  attempts N(s): P(Bin(N, p) >= s) >= 1/2    [paper: N = 4s via Markov, prob 2/3]
Usage: sage -python kk176_indep.py L al be g
"""
import sys
from sage.all import RealBallField, QQ, ZZ, binomial, floor, ceil, log as slog
R = RealBallField(256)

def up(b):            # rigorous upper endpoint as a ball with zero radius
    return R(b.upper())

def rootbad(L, al, be, g, M0=600, variant='floor'):
    """Upper bound B >= P(root bad) for every disjoint tree with capacities
    m_1>=1, m_{j+1}>=g m_j and every level a-spread with a/rho=e^{-L}.
    E D_root <= sum_j exp(c * sum_{h<j} m_h) * phi0(m_j),
    phi0(m) = e^{be m}[e^{-al m} + e^{al m} sum_{r>m/2} C(m,r) x^r],
    c = al + be + log(1+x).  Prefix sum bounded by S(m_j):
      variant 'floor':  S(m) = sum_{i>=1} floor(m/g^i)      (m_{j-i} <= floor(m_j/g^i))
      variant 'paper' : S(m) = floor(m/g) * g/(g-1)           (relaxation used by the other side)"""
    L, al, be = R(L), R(al), R(be)
    x = (-L).exp()
    c = al + be + (1 + x).log()
    def S(m):
        if variant == 'floor':
            s, d = 0, g
            while d <= m:
                s += m // d; d *= g
            return R(s)
        return R(m // g) * R(g) / R(g - 1)
    phi0, phi = {}, {}
    for m in range(1, M0 + 1):
        tail = R(0)
        for r in range(m // 2 + 1, m + 1):
            tail += R(binomial(m, r)) * x**r
        v0 = (be * m).exp() * ((-al * m).exp() + (al * m).exp() * tail)
        phi0[m] = up(v0)
        phi[m] = up(v0 * (c * S(m)).exp())
    # m > M0:  S(m) <= m/(g-1);  sum_{r>m/2} C(m,r) x^r <= 2^m x^{m/2}
    r1 = al - be - c / (g - 1)
    r2 = L / 2 - (al + be + R(2).log() + c / (g - 1))
    assert r1 > 0 and r2 > 0, (r1, r2)
    def tail_sum(m0):   # sum_{m >= m0} (e^{-r1 m} + e^{-r2 m})
        return up((-r1 * m0).exp() / (1 - (-r1).exp()) + (-r2 * m0).exp() / (1 - (-r2).exp()))
    def tail_sup(m0):   # sup_{m >= m0} of the same envelope
        return up((-r1 * m0).exp() + (-r2 * m0).exp())
    def sup_ge(table, m0):
        best = tail_sup(M0 + 1)
        for m in range(max(m0, 1), M0 + 1):
            if table[m].upper() > best.upper(): best = table[m]
        return best
    B = sup_ge(phi0, 1)                       # j = 1
    j = 2
    while g ** (j - 1) <= M0:                 # 2 <= j, g^{j-1} <= M0
        B = up(B + sup_ge(phi, g ** (j - 1)))
        j += 1
    B = up(B + tail_sum(M0 + 1))              # all j with g^{j-1} > M0: distinct m_j > M0
    return B

def Kfac(be):
    be = R(be); lk = R(0); v = 0
    while True:
        t = (-be * 2**v).exp()
        if t.upper() < 1e-40: break
        lk += -(1 - t).log(); v += 1
    # remaining v' >= v: -log(1-t) <= 2t for t<=1/2, and sum_{v'>=v} e^{-be 2^v'} <= 2 e^{-be 2^v}
    lk += 4 * (-be * 2**v).exp()
    return up(lk.exp())

def lam_check(g, rmax=3000):
    lam = QQ(g + 1)**(g + 1) / QQ(g)**g
    for r in range(1, rmax + 1):
        assert 2 * binomial((g + 1) * r, r) < lam**r, r
    # analytic for all r: C(n,k) <= sqrt(n/(2 pi k (n-k))) n^n/(k^k (n-k)^(n-k))  (Robbins)
    #   => 2C((g+1)r,r) <= 2 sqrt((g+1)/(2 pi g r)) lam^r, need 2 sqrt((g+1)/(2 pi g)) < 1
    assert (2 * (R(g + 1) / (2 * R.pi() * g)).sqrt()) < 1
    return lam

def lower_tail_ok(N, s, a, b):
    """exact: P(Bin(N, a/b) <= s-1) <= 1/2  <=>  2 * sum_{k<s} C(N,k) a^k (b-a)^(N-k) <= b^N"""
    tot = ZZ(0)
    for k in range(0, min(s, N + 1)):
        tot += binomial(N, k) * a**k * (b - a)**(N - k)
    return 2 * tot <= b**N

def kl_lower_ok(N, s, pq, target_fail):
    """rigorous Chernoff: P(Bin(N,p) <= s-1) <= exp(-N KL((s-1)/N || p)) when (s-1)/N < p"""
    k = s - 1
    if N * pq <= k: return False
    p = R(pq)
    if k == 0:
        val = (R(N) * (1 - p).log()).exp()
    else:
        t = R(k) / N
        kl = t * (t / p).log() + (1 - t) * ((1 - t) / (1 - p)).log()
        val = (-R(N) * kl).exp()
    return val.upper() <= target_fail

def run(L, al, be, g, S0=200, smax=400, target=QQ(1)/2, S1=5000):
    B = rootbad(L, al, be, g)
    Bp = rootbad(L, al, be, g, variant='paper')
    plo = 1 - B
    den = 10**15
    a = ZZ(floor(QQ(plo.lower().exact_rational()) * den)); b = ZZ(den)
    pq = QQ(a) / b
    assert pq > 0 and R(pq) <= plo
    K = Kfac(be)
    lam = lam_check(g)
    fail = 1 - target
    def exact_ok(N, s):
        tot = ZZ(0)
        for k in range(0, min(s, N + 1)):
            tot += binomial(N, k) * a**k * (b - a)**(N - k)
        return tot <= fail * b**N
    Ns = {}; N = 1
    for s in range(1, S0 + 1):
        N = max(N, s)
        while not exact_ok(N, s): N += 1
        Ns[s] = N
        if target <= QQ(1) / 2:
            assert N <= ceil(QQ(s) / pq), (s, N)
    if target <= QQ(1) / 2:
        def Nb(s):
            return QQ(Ns[s]) if s in Ns else QQ(ceil(QQ(s) / pq))   # Kaas--Buhrman
        ar1 = max(max(QQ(Ns[s]) / QQ(s) for s in Ns), 1 / pq + QQ(1) / (S0 + 1))
        smax_eff = smax
        def tail2():   # s > smax: ratio <= 2/pq + 1/log2 n, log2 n >= s/2
            return (R(2 / pq) + R(2) / smax_eff).upper()
    else:
        # S0 < s <= S1 : Chernoff-KL (rigorous, ball arithmetic); s > S1 : Hoeffding with fixed eps
        N = Ns[S0]
        for s in range(S0 + 1, S1 + 1):
            N = max(N, s)
            while not kl_lower_ok(N, s, pq, fail): N += 1
            Ns[s] = N
        # Hoeffding: P(X <= s-1) <= exp(-2 (Np-s+1)^2 / N) <= fail if Np - s >= sqrt(N log(1/fail)/2).
        # Take N = ceil((1+e) s / p): Np - s >= e s, N <= (1+e)s/p + 1. Sufficient for all s > S1:
        #   e^2 S1 >= ((1+e)/p + 1/S1) * log(1/fail) / 2.
        e = QQ(1) / 1000
        while not (R(e)**2 * S1 >= (R((1 + e) / pq) + R(1) / S1) * (1 / R(fail)).log() / 2):
            e += QQ(1) / 1000
        def Nb(s):
            return QQ(Ns[s]) if s in Ns else QQ(ceil((1 + e) * QQ(s) / pq))
        ar1 = max(max(QQ(Ns[s]) / QQ(s) for s in Ns), (1 + e) / pq + QQ(1) / (S1 + 1))
        smax_eff = S1
        def tail2():
            return (R(2 * (1 + e) / pq) + R(2) / smax_eff).upper()
        print('hoeffding eps =', e)
    def nmin(s):     # least n>=2 with C(n,2) >= 2^(s-1)
        lo, hi = 2, 4
        while hi * (hi - 1) // 2 < 2**(s - 1): hi *= 2
        while lo < hi:
            mid = (lo + hi) // 2
            if mid * (mid - 1) // 2 >= 2**(s - 1): hi = mid
            else: lo = mid + 1
        return lo
    ar2 = R(0); arg = None
    for s in range(1, smax_eff + 1):
        n = nmin(s)
        val = R(Nb(s)) / (R(n).log() / R(2).log())
        if val.upper() > ar2.upper(): ar2, arg = up(val), (s, n)
    ar2 = max(ar2.upper(), tail2())
    base = up(R(L).exp() * K * R(lam))
    C1 = up(base * R(ar1)); C2 = up(base * R(ar2))
    out = dict(target=str(target), B=B.upper(), B_paper_relax=Bp.upper(), p_lo=float(pq), K=K.upper(), lam=float(lam),
               ar1=float(ar1), ar2=float(ar2), ar2_arg=(arg[0], 'n~2^%.2f' % float(R(arg[1]).log() / R(2).log())),
               N_1_8=[Ns[s] for s in range(1, 9)],
               base=base.upper(), C1=C1.upper(), C2=C2.upper(),
               ratio_C1=(R(2048) * R(50).exp() / C1).lower(), ratio_C2=(R(6144) * R(50).exp() / C2).lower())
    for k, v in out.items(): print(f'{k:>14} = {v}')
    return out

if __name__ == '__main__':
    from fractions import Fraction as F
    Q = lambda t: QQ(F(t).numerator) / F(t).denominator
    L, al, be = sys.argv[1:4]; g = int(sys.argv[4])
    tgt = Q(sys.argv[5]) if len(sys.argv) > 5 else QQ(1) / 2
    run(Q(L), Q(al), Q(be), g, target=tgt)
