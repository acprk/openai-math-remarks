#!/usr/bin/env python3
"""Family 135: explicit lower bound on log S / (sqrt(n) log n) for the modified proof.

Modified proof (see note.pdf, Section 2):
  s  = largest integer <= sqrt(n/2) with parity of n   (paper: <= sqrt(n))
  t1 = n/(2s)                                           (paper: t0 = n/(4s), single expansion)
  every factor occurrence of degree >= t1 is expanded simultaneously (S^{|X|} terms),
  unexpanded factors keep their saving q^{s e/(2n)}.
  =>  log S >= min( Lq*s/2 - K ,  t1*c' ),
      c' = Lq (n-2s)/(2n) - log 2,  K = C_I + Cost_even + Cost_odd + log 2.

Every quantity below is replaced by an explicit majorant/minorant that is monotone
in n (worst case of s in (sqrt(n/2)-2, sqrt(n/2)]).  B1(n), B2(n) are increasing
lower bounds for the two branches divided by sqrt(n) log n, so B(n0) >= c implies
B(n) >= c for all n >= n0.

Positive control: --paper reproduces the paper's chain (low 1/20, high 3/160,
Prop 1/100, Thm 1/200 -> 1/400) and the single-expansion optimum 1/8.
"""
import sys
from mpmath import mp, mpf, sqrt, log, exp, findroot

mp.dps = 60
LOG2 = log(2)


def pieces(ell):
    """ell = log n.  Return (B1, B2, K/sqrt(n)) with all majorants explicit."""
    L = mpf(ell)
    n = exp(L)
    rn = sqrt(n)
    s_hi = sqrt(n / 2)
    s_lo = s_hi - 2
    k_lo = (n - s_hi) / 2
    # alpha in [n^{-1/2}, n^{-1/2} e^{eps_a}],  eps_a = 1/a0 = (sqrt n - 1)/(k n^2)
    eps_a = (rn - 1) / (k_lo * n**2)
    inv_b0 = 1 / (((n + s_lo) / 2) * n**2 * n**(-mpf(1) / 2))   # b0 >= u*alpha^rho >= m n^2 n^{-1/2}
    # Lq = log(1/q) >= (n/(n+s)) (L/2 - eps_a) >= (1 - 1/sqrt(2n)) (L/2 - eps_a)
    Lq_lo = (1 - 1 / sqrt(2 * n)) * (L / 2 - eps_a)
    theta = exp(L / (2 * sqrt(2 * n)) + eps_a)          # q <= theta / sqrt(n)
    q_hi = theta / rn
    ns_n = exp(L / sqrt(2 * n))                          # n^{s/n}
    # ---- circuit-side costs
    cost_even = (rn / 2 * theta * ns_n + sqrt(2) * sqrt(theta) * n**(mpf(1) / 4)
                 + rn / 2 * theta) / (1 - q_hi)
    cost_odd = 4 * s_hi * (sqrt(theta) * n**(-mpf(1) / 4) + q_hi / (1 - q_hi))
    # ---- IMM side  C_I = 2 loss1 + log W + m*beta
    loss1 = 2 / rn
    logW = (2 * k_lo / n**2 + 2 / n**2) + (L / 2 + 1) + rn / 2 + 2 * n**(-mpf(3) / 2) \
        + (rn + 2 * L + log(7))
    mbeta = rn / 2 * (1 + 1 / sqrt(2 * n)) * exp(L / sqrt(2 * n) + eps_a + inv_b0)
    C_I = 2 * loss1 + logW + mbeta
    K = C_I + cost_even + cost_odd + LOG2
    B1 = Lq_lo * s_lo / 2 / (rn * L) - K / (rn * L)
    t1_lo = sqrt(n / 2)
    cprime = Lq_lo * (1 - sqrt(2 / n)) / 2 - LOG2
    B2 = t1_lo * cprime / (rn * L)
    # sanity: q <= 0.1 needed for odd low factors x<=1/8
    assert q_hi <= mpf('0.1')
    return B1, B2, K / rn, dict(Kc=K / rn, CI=C_I / rn, ce=cost_even / rn, co=cost_odd / rn)


def B(ell):
    b1, b2, _, _ = pieces(ell)
    return min(b1, b2)


def threshold(c):
    lo, hi = mpf(14), mpf(20000)          # n >= e^14 ~ 1.2e6 for side conditions
    if B(lo) >= c:
        return lo
    assert B(hi) >= c
    for _ in range(200):
        mid = (lo + hi) / 2
        if B(mid) >= c:
            hi = mid
        else:
            lo = mid
    return hi


def paper_control():
    print("== positive control: paper's own chain (03-circuits.tex l.250-299, 06-completion.tex l.9-29)")
    # low case: q^{s/4} with q<=n^{-2/5}, s>=sqrt n/2 -> n^{-sqrt n/20}
    low = mpf(2) / 5 * mpf(1) / 2 / 4
    # high case: q^{lambda t0/2}, q<=n^{-2/5}, lambda>=3/8, t0>=sqrt n/4
    high = mpf(2) / 5 * mpf(3) / 8 * mpf(1) / 4 / 2
    print("  low exponent  =", low, "(paper 1/20 =", mpf(1) / 20, ")")
    print("  high exponent =", high, "(paper 3/160 =", mpf(3) / 160, ")")
    print("  min           =", min(low, high), ">= 1/100 used in Prop 3.1")
    print("  Theorem: S^2 >= n^{sqrt n/100} e^{-O(sqrt n)} -> S >= n^{(1/200-o(1)) sqrt n}; stated 1/400")
    print("  the proof as written actually yields any c < 3/320 =", mpf(3) / 320)
    # single-expansion variant: q~n^{-1/2}, lambda~1/2, threshold n/(2s), s=sigma sqrt n
    #   low: sigma/4 ; high (one factor, S^2): 1/(16 sigma)  -> sigma=1/2, c=1/8
    f = lambda sg: min(sg / 4, 1 / (16 * sg))
    print("  single-expansion optimum: sigma=1/2 ->", f(mpf(1) / 2), "(=1/8)")
    # new: multi-expansion: min(sigma/4, 1/(8 sigma)) -> sigma=1/sqrt2 -> 1/(4 sqrt2)
    g = lambda sg: min(sg / 4, 1 / (8 * sg))
    print("  new multi-expansion optimum: sigma=1/sqrt2 ->", g(1 / sqrt(2)), " 1/(4sqrt2) =", 1 / (4 * sqrt(2)))


if __name__ == "__main__":
    paper_control()
    print("\n== explicit monotone lower bound B(n) = min(B1,B2) for the new proof")
    for ell in [14, 20, 50, 100, 200, 300, 500, 1000, 5000, 20000]:
        b1, b2, kc, d = pieces(ell)
        print(f"  log n={ell:6d}: B1={float(b1):.5f} B2={float(b2):.5f}  K/sqrt n={float(kc):.4f}"
              f" (C_I {float(d['CI']):.3f}, even {float(d['ce']):.3f}, odd {float(d['co']):.4f})")
    # monotonicity spot check on a dense grid
    prev = None
    ok = True
    for i in range(0, 4000):
        ell = 14 + i * 5
        v = B(ell)
        if prev is not None and v < prev:
            ok = False
            print("  non-monotone at", ell)
        prev = v
    print("  B monotone on grid log n in [14,20009] step 5:", ok)
    for c in [mpf(1) / 400, mpf(3) / 320, mpf(1) / 8, mpf(3) / 20, mpf(1) / 6, mpf('0.17')]:
        t = threshold(c)
        print(f"  c={float(c):.6f}: holds for all log n >= {float(t):.2f}  (n >= 10^{float(t / log(10)):.1f})")
    print("  limit as n->inf:", 1 / (4 * sqrt(2)))
