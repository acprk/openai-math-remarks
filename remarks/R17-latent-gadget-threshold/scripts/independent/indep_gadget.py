#!/usr/bin/env python3
"""Independent re-implementation for the latent-gadget remark (does not import the primary script).

Part A (lower bound, Prop. 1 of the note).  Random small K-equivariant maps C: F_2^n -> K = F_2^l,
C(k, w) = k xor f(w), and random rational noise laws mu.  Everything is computed in exact rational
arithmetic (fractions.Fraction) by direct summation (no Hadamard transform):
   p        = Pr[C(x+a) != C(x)],          p' = p 2^l / (2^l - 1)
   nu_s(h)  = (E_x (-1)^{s(C(x)) + h(x)})^2,  supported on {h : h|_K = s}, total mass 1
   E_r      = E over s_1..s_r iid uniform nonzero, h_i ~ nu_{s_i} of Pr_a[exists i: h_i(a) = 1]
 and checks  E_1 == p'/2  exactly  and  E_r <= r p'/2  for r = 2, 3.
 (These maps are not gadgets; this is a unit test of the inequality chain, which holds for any C, mu.)

Part B (upper bound, Prop. 2 of the note).  For the worst-case rank chain allowed by the paper's
generic-state estimate (loss of exactly one with probability min(1, 3 r theta)), compute the exact
extinction probability after s steps in exact rationals, and check
   (i)  W_j = e^{-3 theta j}/(r_j + 1) has non-increasing expectation (supermartingale), and
   (ii) Pr[r_s = 0] <= e^{3 theta s}/(r_0 + 1)  (the bound used in the note),
 also for a second chain whose loss probability is a fixed fraction of the allowed maximum,
 and compare with the paper's harmonic/Markov bound 3 s theta / H_{r_0}.

Part C (parameters).  For p in {1e-1, 1e-2, 1e-4, 1e-8}, the threshold certified by
 (A) the paper's main-text choice  exp(1 - H_{r_0}/8) <= p,  eps_0 = 1/(H_{r_0} 2^{r_0}), r_* = 3d+2;
 (B) the paper's effective recipe  2^{-m} <= p, H_{r_0} >= 8(m+1), r_* = 3d+2;
 (C) the paper's harmonic analysis without the enlargement step (r_* = r_0 + 2);
 (D) this note, r_* = ceil(1.5 p^{-3}) + 1;  lower bound 1/(4p);  route limit 2 ln(8/7) p^{-3}.
 For (A),(B), d is the least integer with  [binom(2^r-1, 3r+1)(6r+3) 2^{-d} + binom(2^r-1,2) 4^{-d} + 2^r 8^{-d}]
 / (1 - 2^r 8^{-d}) <= eps_0  at r = r_0 (the paper's Lemma 2.4 union bound; the term is increasing
 in r, so r = r_0 dominates all ranks 1..r_0), using log2 binom(N,k) <= k log2 N - log2 k!.
Part D. Exact (float) worst-case rank chain at the note's parameters, vs the supermartingale bound
 and the independent-particle approximation (1 - e^{-3 theta s})^{r_0}.
Environment: python3, mpmath >= 1.3. Single thread; ~2 min.
"""
import random, math
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 100

# ---------------- Part A ----------------
def parity(x):
    return bin(x).count("1") & 1

def part_A(n, l, seed, rmax=3):
    rnd = random.Random(seed)
    N, L = 1 << n, 1 << l
    f = [rnd.randrange(L) for _ in range(N >> l)]
    C = [(x & (L - 1)) ^ f[x >> l] for x in range(N)]
    supp = rnd.sample(range(1, N), 5)
    wts = [rnd.randint(1, 9) for _ in supp]
    mu = {a: Fr(w, sum(wts)) for a, w in zip(supp, wts)}
    p = sum(w * Fr(sum(1 for x in range(N) if C[x ^ a] != C[x]), N) for a, w in mu.items())
    pp = p * L / (L - 1)
    nu = {}
    for s in range(1, L):
        dist = {}
        for h in range(N):
            c = Fr(sum(1 - 2 * (parity(s & C[x]) ^ parity(h & x)) for x in range(N)), N)
            if c:
                assert (h & (L - 1)) == s          # equivariance forces h|_K = s
                dist[h] = c * c
        assert sum(dist.values()) == 1
        nu[s] = dist
    # h_i ~ nu_{s_i} with s_i uniform nonzero is the same as h_i ~ rho, rho(h) = nu_{h|_K}(h)/(L-1),
    # independently; hence E_r = sum_a mu(a) (1 - (1 - pi_a)^r),  pi_a = Pr_{h~rho}[h(a) = 1].
    rho = {h: w / (L - 1) for s_ in nu for h, w in nu[s_].items()}
    assert sum(rho.values()) == 1
    pi = {a: sum(w for h, w in rho.items() if parity(h & a)) for a in mu}
    res = [sum(mu[a] * (1 - (1 - pi[a]) ** r) for a in mu) for r in range(1, rmax + 1)]
    return p, pp, res

print("Part A: exact rational checks of E_r <= r p'/2 (r=1 must be an identity)")
worst = 0.0
for seed in range(12):
    n, l = (6, 2) if seed % 2 == 0 else (5, 3)
    p, pp, E = part_A(n, l, seed)
    if p == 0:
        continue
    assert E[0] == pp / 2, "r=1 identity fails"
    for r, Er in enumerate(E, 1):
        assert Er <= r * pp / 2, "inequality fails"
    worst = max(worst, max(float(Er / (r * pp / 2)) for r, Er in enumerate(E, 1) if r > 1))
    print(f"  seed={seed} n={n} l={l}: p={p} ({float(p):.4f})  E1 = p'/2 exactly;  "
          f"E2/(2p'/2)={float(E[1]/pp):.4f}  E3/(3p'/2)={float(E[2]/(3*pp/2)):.4f}")
print(f"  max ratio E_r/(r p'/2) for r=2,3 over all trials: {worst:.4f} (<= 1 required)")

# ---------------- Part B ----------------
def chain(r0, theta, s, frac):
    dist = {r0: Fr(1)}
    EW = []
    e3 = mp.e ** (3 * mp.mpf(theta.numerator) / theta.denominator)
    for j in range(s + 1):
        ew = sum(w / (r + 1) for r, w in dist.items())
        EW.append(mp.mpf(ew.numerator) / ew.denominator / e3 ** j)
        if j == s:
            break
        new = {}
        for r, w in dist.items():
            if r == 0:
                new[0] = new.get(0, 0) + w; continue
            pl = min(Fr(1), frac * 3 * r * theta)
            new[r - 1] = new.get(r - 1, 0) + w * pl
            new[r] = new.get(r, 0) + w * (1 - pl)
        dist = new
    ext = dist.get(0, Fr(0))
    mono = all(EW[j + 1] <= EW[j] * (1 + mp.mpf(10) ** -50) for j in range(s))
    return ext, mono

print("Part B: exact rank chains (theta = 1/200), extinction vs bounds")
theta = Fr(1, 200)
for p_target in (Fr(1, 2), Fr(3, 10)):
    s = math.ceil(math.log(1 / float(p_target)) / float(theta))
    for r0 in (4, 16, 40):
        for frac in (Fr(1), Fr(1, 2)):
            ext, mono = chain(r0, theta, s, frac)
            bound = mp.e ** (3 * mp.mpf(1) / 200 * s) / (r0 + 1)
            H = sum(mp.mpf(1) / j for j in range(1, r0 + 1))
            markov = 3 * s * mp.mpf(1) / 200 / H
            ok = mp.mpf(ext.numerator) / ext.denominator <= bound
            print(f"  p={p_target} s={s} r0={r0} loss=x{frac}: Pr[extinct]={float(ext):.4e}  "
                  f"bound e^(3 theta s)/(r0+1)={float(bound):.4e} holds={ok}  E[W_j] nonincreasing={mono}  "
                  f"paper Markov 3 s theta/H={float(markov):.3f}")

# ---------------- Part C ----------------
EG = mp.euler
def H_approx(r):     # harmonic number, exact for small r, asymptotic (error < 1/(120 r^4)) for large
    r = mp.mpf(r)
    if r < 10 ** 6:
        return mp.harmonic(r)
    return mp.log(r) + EG + 1 / (2 * r) - 1 / (12 * r ** 2)

def smallest_r_with_H_at_least(T):
    x = mp.exp(T - EG)                      # H_r = ln r + gamma + 1/(2r) - 1/(12 r^2) + O(r^-4)
    if x > 10 ** 7:                         # unit steps are below working precision: solve directly
        return mp.ceil(mp.findroot(lambda r: H_approx(r) - T, x))
    r = mp.floor(x)
    while H_approx(r) >= T: r -= 1
    while H_approx(r) < T: r += 1
    return r

def log2_binom_upper(log2N, k):
    return k * log2N - mp.log(mp.factorial(k), 2) if k < 10 ** 6 else k * log2N - (mp.loggamma(k + 1) / mp.log(2))

def d_needed(r, log2_eps0):
    # binom(2^r - 1, 3r+1) <= (2^r)^(3r+1)/(3r+1)!
    lb = log2_binom_upper(r, 3 * r + 1)
    d = mp.ceil(lb + mp.log(6 * r + 3, 2) + 1 - log2_eps0)     # main term <= eps0/2
    # remaining terms are <= eps0/4 for this d (checked):
    assert 2 * r - 2 * d < log2_eps0 - 3 and r - 3 * d < log2_eps0 - 3
    return d

print("Part C: thresholds certified by each argument (log2 values; p^-e gives the exponent e)")
print("  A = paper, main text  (exp(1-H_r0/8)<=p, r_* = 3d+2, d from the paper's union bound)")
print("  B = paper, effective recipe (2^-m<=p, H_r0>=8(m+1), r_* = 3d+2)")
print("  C = paper's harmonic analysis with K = B directly / padding (r_* = r0 + 2, r0 as in A)")
print("  D = this note: r_* = ceil(1.5 p^-3) + 1      L = lower bound 1/(4p)   M = 2 ln(8/7) p^-3 (route limit)")
for p in (mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('1e-4'), mp.mpf('1e-8')):
    lp = mp.log(1 / p, 2)
    r0A = smallest_r_with_H_at_least(8 * (1 + mp.log(1 / p)))
    dA = d_needed(r0A, -(mp.log(H_approx(r0A), 2) + r0A)); A = 3 * dA + 2
    m = mp.ceil(mp.log(1 / p, 2))
    r0B = smallest_r_with_H_at_least(8 * (m + 1))
    dB = d_needed(r0B, -(mp.log(H_approx(r0B), 2) + r0B)); B = 3 * dB + 2
    C = r0A + 2
    D = mp.ceil(mp.mpf(3) / 2 * p ** -3) + 1
    L = 1 / (4 * p)
    M = 2 * mp.log(mp.mpf(8) / 7) * p ** -3
    cells = []
    for name, v in (("A", A), ("B", B), ("C", C), ("D", D), ("L", L), ("M", M)):
        cells.append(f"{name}=2^{float(mp.log(v,2)):6.1f}(e={float(mp.log(v,2)/lp):5.2f})")
    print(f"  p={mp.nstr(p,2):6s} " + " ".join(cells))
print("  asymptotic exponents: A -> 16, B -> 16/ln 2 = %.3f, C -> 8, D -> 3, L -> 1" % (16 / math.log(2)))

# ---------------- Part D ----------------
print("Part D: worst-case chain at the note's parameters (theta = 1/4000, floats)")
print("  r0 = ceil(1.5 p^-3) - 1, s = ceil(ln(1/p)/theta); proof needs bound <= 3/4 - 1/60 (theta <= 0.01 in the note)")
def ext_float(r0, th, s):
    dist = [0.0] * (r0 + 1); dist[r0] = 1.0
    for _ in range(s):
        new = [0.0] * (r0 + 1); new[0] = dist[0]
        for r in range(1, r0 + 1):
            pl = min(1.0, 3 * r * th)
            new[r - 1] += dist[r] * pl; new[r] += dist[r] * (1 - pl)
        dist = new
    return dist[0]
th = 1 / 4000
for p in (0.5, 0.4, 0.3, 0.2):
    r0 = math.ceil(1.5 * p ** -3) - 1
    s = math.ceil(math.log(1 / p) / th)
    ext = ext_float(r0, th, s)
    bnd = math.exp(3 * th * s) / (r0 + 1)
    indep = (1 - math.exp(-3 * th * s)) ** r0
    print(f"  p={p}: r0={r0} s={s}  exact={ext:.4f}  supermartingale bound={bnd:.4f} (<= {3/4-1/60:.4f}: {bnd <= 3/4-1/60})  "
          f"(1-e^(-3 theta s))^r0={indep:.4f}")
print("  route limit (method only): r_s ~ Bin(r0, p^3), detection ~ 1 - E 2^(-r_s) ~ 1 - exp(-r0 p^3/2) >= 1/8 needs")
print("  r0 >= 2 ln(8/7) p^-3 (1+o(1)); 2 ln(8/7) = %.4f vs the note's 1.5 (ratio %.2f)" % (2*math.log(8/7), 1.5/(2*math.log(8/7))))
