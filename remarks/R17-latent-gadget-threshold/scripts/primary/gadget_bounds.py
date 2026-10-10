#!/usr/bin/env python3
"""Evidence for the new latent-gadget bounds (REPORT.md section 4, Theorems N1/N2).

Part 1 (Theorem N1, lower bound r* >= 2*kappa*Pr[indep]/p' -> 1/(4p) for kappa=1/8):
  Random K-equivariant maps C : V=F_2^n -> K=F_2^l and random noise laws mu
  on V.  We compute exactly p = Pr[C(x+a)!=C(x)], and run the Fourier-sampling
  construction of the proof: sample sigma_i uniform nonzero, h_i ~ |F_sigma^(h)|^2,
  and compute the exact detection probability of G=span(h_i).  We check the key
  inequality  E[detect(G)] <= r p'/2,  p' = p 2^l/(2^l-1)  (union-bound step),
  i.e. the averaged detection never exceeds the bound.

Part 2 (Theorem N2, death-chain analysis):
  Exact DP of the pure death chain  r -> r-1 w.p. min(1,3 r theta)  for s steps,
  vs the Chernoff bound  Pr[extinct by s] <= e^{3 theta s}/(r0+1)  used in the proof,
  and vs the paper's Markov/harmonic bound  4 s theta / H_{r0}.

Part 3: parameter comparison  r*(p)  paper (harmonic, Lemma 2.1) vs new analysis,
  using the paper's own genericity union bound to choose d.
"""
import math, random, itertools
import numpy as np

# ---------------- Part 1 ----------------
def part1(n=7, l=3, trials=30, seed=3):
    """Exact (no sampling) expectation over sigma_i uniform nonzero, h_i ~ |F^_sigma_i|^2,
    for r=1 (must be EQUAL to p'/2) and r=2 (must be <= p').  Returns max ratios."""
    rng = np.random.default_rng(seed)
    N = 1 << n; L = 1 << l; W = 1 << (n - l)
    popc = np.array([bin(i).count("1") & 1 for i in range(N)], dtype=np.int8)
    H = np.array([[1]])
    for _ in range(n): H = np.block([[H, H], [H, -H]])
    xs = np.arange(N)
    r1_dev, r2_ratio = 0.0, 0.0
    for _ in range(trials):
        f = rng.integers(0, L, size=W)
        C = (xs & (L - 1)) ^ f[xs >> l]            # C(k,w)=k xor f(w): K-equivariant
        supp = rng.integers(0, N, size=6); wts = rng.random(6); wts /= wts.sum()
        mu = np.zeros(N)
        for v, w_ in zip(supp, wts): mu[v] += w_
        p = sum(mu[v] * np.mean(C[xs ^ v] != C) for v in range(N) if mu[v])
        if p == 0: continue
        pp = p * L / (L - 1)
        nu = {}
        for sigma in range(1, L):
            Fh = (H @ (1 - 2 * popc[sigma & C])) / N
            pr = Fh ** 2
            assert abs(pr.sum() - 1) < 1e-9
            # support only on lifts of sigma (equivariance check)
            for h in np.nonzero(pr > 1e-12)[0]:
                assert (h & (L - 1)) == sigma
            nu[sigma] = [(int(h), pr[h]) for h in np.nonzero(pr > 1e-12)[0]]
        def det(hs):
            return sum(mu[v] for v in range(N) if mu[v] and any(popc[h & v] for h in hs))
        E1 = np.mean([sum(w * det([h]) for h, w in nu[s1]) for s1 in nu])
        r1_dev = max(r1_dev, abs(E1 - pp / 2))
        E2 = np.mean([sum(w1 * w2 * det([h1, h2]) for h1, w1 in nu[s1] for h2, w2 in nu[s2])
                      for s1 in nu for s2 in nu])
        r2_ratio = max(r2_ratio, E2 / pp)
    return r1_dev, r2_ratio

# ---------------- Part 2 ----------------
def death_chain_extinct(r0, theta, s):
    dist = np.zeros(r0 + 1); dist[r0] = 1.0
    for _ in range(s):
        new = np.zeros_like(dist); new[0] = dist[0]
        for r in range(1, r0 + 1):
            q = min(1.0, 3 * r * theta)
            new[r] += dist[r] * (1 - q); new[r - 1] += dist[r] * q
        dist = new
    return dist[0]

def part2():
    rows = []
    for p in (0.3, 0.1, 0.05):
        theta = 1e-3
        s = math.ceil(math.log(1 / p) / theta)
        for r0 in (int(2 * p ** -3), int(4 * p ** -3)):
            exact = death_chain_extinct(r0, theta, s)
            chern = math.exp(3 * theta * s) / (r0 + 1)
            Hr = sum(1 / j for j in range(1, r0 + 1))
            markov = 4 * s * theta / Hr
            rows.append((p, r0, s, exact, chern, markov))
    return rows

# ---------------- Part 3 ----------------
def log2_binom(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2)

def d_needed(r0, eps0):
    # paper's Lemma 2.4 union bound at rank r0 (dominant term), require <= eps0/2
    lb = log2_binom(2.0 ** r0, 3 * r0 + 1) if 2 ** r0 - 1 > 3 * r0 else 0.0
    return math.ceil(lb + math.log2(6 * r0 + 3) + math.log2(2 / eps0)) + 2

def paper_rstar_log2(p):
    m = math.ceil(math.log2(1 / p))
    # smallest r0 with H_{r0} >= 8(m+1):  H_r ~ ln r + gamma + 1/(2r)
    lnr0 = 8 * (m + 1) - 0.5772156649
    r0 = math.exp(lnr0)
    Hr0 = 8 * (m + 1)
    d = 3 * r0 * r0                # leading term of d_needed (r0 astronomically large)
    return math.log2(3 * d + 2), r0

def new_rstar(p):
    r0 = math.ceil(4 * math.e ** 0.01 * p ** -3)   # theta tiny
    s_theta = math.log(1 / p) + 1e-3
    eps0 = 1 / (8 * s_theta * 2 ** r0) if r0 < 1000 else None
    if eps0 is not None:
        d = d_needed(r0, eps0)
    else:
        d = 3 * r0 * r0
    return math.log2(3 * d + 2), r0

if __name__ == "__main__":
    d1, q2 = part1()
    print("Part 1 (exact): max |E[detect], r=1  - p'/2| = %.2e (identity); max E[detect]/(2p'/2), r=2 = %.4f (must be <=1)" % (d1, q2))
    print("Part 2: p, r0, s, exact extinction, Chernoff bound (proof), Markov/harmonic bound (paper):")
    for row in part2():
        print("   p=%.2f r0=%d s=%d exact=%.3e chernoff=%.3e harmonic=%.3e" % row)
    print("Part 3: log2 r*(p): paper vs new vs lower bound log2(1/(4p))")
    for p in (1e-1, 1e-2, 1e-4, 1e-8):
        lp, r0p = paper_rstar_log2(p)
        ln_, r0n = new_rstar(p)
        print("   p=%g  paper log2 r*~%.1f (r0~2^%.1f)   new log2 r*~%.1f (r0=%.3g)   lower log2 r*>=%.1f"
              % (p, lp, math.log2(r0p), ln_, r0n, math.log2(1 / (4 * p))))
        print("      exponents: paper r* ~ p^-%.1f, new r* ~ p^-%.2f" % (lp / math.log2(1 / p), ln_ / math.log2(1 / p)))
