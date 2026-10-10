#!/usr/bin/env python3
"""Monte Carlo check of the stage-1 energy machinery (paper Sec. 3-4) in the cyclic-coordinate
frame, for a uniformly random ordered list of length k=(1-1/b)n; tagged card = last list card
(m = n-k+1 free positions, the bottleneck).  For each run we track the conditional flow w_s of
the tagged card given the exposed paths of the first k-1 cards (rules of Lemma 3.1).
Reports, per time t:
  E_t      = E||w_t||^2
  b_t      = E sum_x w_t(x)^2 1{partner of x in direction i_t occupied}
  iota_t   = E||P_t w_0||^2   (t<d; initial term of the extended energy identity)
  RHS_t    = iota_t + sum_{j=1}^{min(t,d)} 2^-j b_{t-j}   (should equal E_t up to MC error)
  b_t/E_t  (paper bound: <= 1-beta, beta -> m/n)
and the empirical per-step decay rate of E_t vs gamma_b = 1-1/(2b) (finite-d root also shown).
Usage: python3 energy_mc.py d b runs [T_over_d]
Environment: Python >= 3.8 + NumPy (tested CPython 3.10.12, NumPy 2.2.6); fixed seed, single thread.
Expected: python3 energy_mc.py 8 4 4000 > expected_output/energy_mc_d8_b4.txt  (~1 min);
  RHS_t matches E_t within MC error (the identity of Lemma 2.2 is exact), b_t/E_t ~ 0.73-0.76 < 1-beta,
  fitted per-step decay 0.8676 vs finite-d root 0.8686 (see expected_output).  Sanity check only; not part of the proof.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, numpy as np
d, b, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Tod = float(sys.argv[4]) if len(sys.argv) > 4 else 12
n = 1 << d
k = n - n // b                       # list length; m = n-k+1
m = n - k + 1
T = int(Tod * d)
rng = np.random.default_rng(12345 + 1000 * d + b)
B = 2000                             # batch
acc_E = np.zeros(T + 1); acc_b = np.zeros(T + 1); acc_i = np.zeros(T + 1)
acc_E2 = np.zeros(T + 1)
done = 0
while done < R:
    r = min(B, R - done)
    # random ordered list: permutation; exposed = first k-1 labels, tagged = k-th
    perm = np.argsort(rng.random((r, n)), axis=1)
    occ = np.zeros((r, n), bool)
    rows = np.arange(r)[:, None]
    occ[rows, perm[:, :k - 1]] = True          # initial position = label
    tag = perm[:, k - 1]
    p = np.zeros((r, n)); p[np.arange(r), tag] = 1.0
    w0 = p - (~occ) / m
    for t in range(T + 1):
        w = p - (~occ) / m
        E = (w * w).sum(1)
        acc_E[t] += E.sum(); acc_E2[t] += (E * E).sum()
        if t < d:   # iota_t = ||P_t w0||^2, P_t averages coordinates 0..t-1 (bits)
            blk = 1 << t
            # positions grouped by high bits (bits >= t): x >> t
            ww = w0.reshape(r, n // blk, blk).mean(2)
            acc_i[t] += (ww * ww * blk).sum()
        if t == T:
            break
        bit = 1 << (t % d)
        x0 = np.array([x for x in range(n) if not x & bit]); x1 = x0 | bit
        part_occ = occ[:, np.r_[0:n]] if False else None
        # b_t: partner occupied
        partner = np.arange(n) ^ bit
        acc_b[t] += ((w * w) * occ[:, partner]).sum()
        c = rng.random((r, n // 2)) < 0.5               # swap coins per pair
        o0, o1 = occ[:, x0], occ[:, x1]
        p0, p1 = p[:, x0], p[:, x1]
        no0 = np.where(c, o1, o0); no1 = np.where(c, o0, o1)
        both = (~o0) & (~o1)
        avg = 0.5 * (p0 + p1)
        # one free: mass of free endpoint goes to the free output
        mass = np.where(~o0, p0, 0.0) + np.where(~o1, p1, 0.0)
        one = o0 ^ o1
        np0 = np.where(both, avg, np.where(one & ~no0, mass, 0.0))
        np1 = np.where(both, avg, np.where(one & ~no1, mass, 0.0))
        p[:, x0], p[:, x1] = np0, np1
        occ[:, x0], occ[:, x1] = no0, no1
    done += r
E = acc_E / R; bb = acc_b / R; io = acc_i / R
seE = np.sqrt(np.maximum(acc_E2 / R - E ** 2, 0) / R)
beta = m / n; gam = 1 - beta / 2
# finite-d root of (1-beta) sum_{j=1}^d (2g)^-j = 1
lo, hi = 0.5, 1.0
for _ in range(200):
    mid = (lo + hi) / 2
    f = (1 - beta) * sum((2 * mid) ** -j for j in range(1, d + 1)) - 1
    lo, hi = (mid, hi) if f > 0 else (lo, mid)
print(f"d={d} n={n} b={b} k={k} m={m} runs={R}  beta=m/n={beta:.4f} gamma_inf={gam:.5f} gamma_d={lo:.5f}")
print(" t     E_t        se       RHS_t      b_t/E_t   E_t/E_{t-1}")
for t in range(T + 1):
    rhs = (io[t] if t < d else 0) + sum(2.0 ** -j * bb[t - j] for j in range(1, min(t, d) + 1))
    ratio = E[t] / E[t - 1] if t else float('nan')
    br = bb[t] / E[t] if t < T else float('nan')
    print(f"{t:3d} {E[t]:.4e} {seE[t]:.1e} {rhs:.4e} {br:8.4f} {ratio:8.4f}")
# fit decay on the tail where E > 50*se
ts = [t for t in range(d, T + 1) if E[t] > 20 * seE[t]]
if len(ts) > 5:
    sl = np.polyfit(ts[len(ts) // 3:], np.log(E[ts[len(ts) // 3:]]), 1)[0]
    print(f"fitted decay per step exp(slope)={np.exp(sl):.5f} vs gamma_d={lo:.5f}; log2 rate {-sl / np.log(2):.4f} vs L={-np.log2(lo):.4f}")
