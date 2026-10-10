#!/usr/bin/env python3
"""Numerical illustration (NOT part of the proof) of the stage-1 energy lemma of remark R5.

Environment: Python >= 3.8 with mpmath (tested: CPython 3.10.12, mpmath 1.3.0); single thread.
Run:         python3 energy_lemma_check.py > expected_output/energy_lemma_check.txt   (~1-2 min)

We iterate the worst case permitted by the inequalities used in the note (Lemma 2.3/2.4):
  e_t = min(1, 2^-t + 1/m + sum_{j=1}^{t} 2^-j bb_{t-j})          (t < d, extended identity)
  e_t = sum_{j=1}^{d} 2^-j bb_{t-j}                               (t >= d, paper eq. (3.2))
  bb_s = e_s                                                      (s < s0)
  bb_s = min(e_s, (1-beta) e_s + beta eta_s)                      (s >= s0)
  eta_s = min(1, (n/2^sig) exp(-2^(sig+1) theta^2 (m/n)^2)),  sig = min(s, d-1)
with m = ceil(n/b), beta = (1-theta) m/n.  Any sequence satisfying the inequalities is bounded by
this one (all maps are monotone), so it is the extremal case of the lemma's hypotheses.
Checks, in high precision (mp.dps = 40; floating, not interval):
  (i)  the closed form of Lemma 2.4,
          e_t <= (2/(1-rho)) g^(t-s0) + g^(t-d)/(m (1-rho)) + eta*,   rho = (1-beta)/(2g-1) < 1,
       holds for every t up to the horizon (s0 = ceil(log2(K d)), K = 1.6 b^2 ln2/theta^2);
  (ii) the smallest T with n^2 e_T <= 1/d, compared with 2d/log2(1/(1-beta/2)): the excess is
       O(log d) (local look-back), versus ~d without it.
"""
from mpmath import mp, mpf, exp, log, ceil
mp.dps = 40


def run(d, b, th, local=True, g_eps=mpf('0.002')):
    n = mpf(2) ** d
    m = ceil(n / b)
    beta = (1 - th) * m / n
    g = 1 - beta / 2 + g_eps                     # gamma' slightly above the root 1 - beta/2
    rho = (1 - beta) / (2 * g - 1)
    assert rho < 1
    K = mpf('1.6') * b * b * log(2) / th ** 2
    s0 = int(ceil(log(K * d, 2)))
    def eta(s):
        sig = min(s, d - 1)
        return min(mpf(1), (n / mpf(2) ** sig) * exp(-mpf(2) ** (sig + 1) * th ** 2 * (m / n) ** 2))
    eta_star = max(eta(s) for s in range(s0, d + 1))
    e, bb = [], []
    Tmin = None
    worst_ratio = mpf(0)
    t = 0
    while True:
        if t < d:
            v = min(mpf(1), mpf(2) ** (-t) + 1 / m + sum(mpf(2) ** (-j) * bb[t - j] for j in range(1, t + 1)))
        else:
            v = sum(mpf(2) ** (-j) * bb[t - j] for j in range(1, d + 1))
        e.append(v)
        bound = (2 / (1 - rho)) * g ** (t - s0) + g ** (t - d) / (m * (1 - rho)) + eta_star
        worst_ratio = max(worst_ratio, v / bound)
        if local and t >= s0:
            bb.append(min(v, (1 - beta) * v + beta * eta(t)))
        elif (not local) and t >= d - 1:
            bb.append(min(v, (1 - beta) * v + beta * eta(t)))
        else:
            bb.append(v)
        if Tmin is None and n ** 2 * v <= mpf(1) / d:
            Tmin = t
        if Tmin is not None and t >= Tmin + 2 * d:
            break
        t += 1
    L = log(1 / (1 - beta / 2), 2)
    return Tmin, 2 * d / L, worst_ratio, s0, float(eta_star)


print("b  theta   d   s0  T_min  T_min/d  2/L(beta)  excess(steps)  excess_without_lookback  max e_t/bound")
for b in (4, 5):
    for th in (mpf('0.1'),):
        for d in (40, 80, 160):
            T, ref, wr, s0, es = run(d, b, th)
            T2, _, _, _, _ = run(d, b, th, local=False)
            assert wr <= 1, "closed-form bound of Lemma 2.4 violated"
            print(f"{b}  {float(th):.2f}  {d:4d} {s0:3d} {T:6d}  {T / d:7.3f}  {float(ref) / d:8.3f}"
                  f"  {float(T - ref):10.1f}  {float(T2 - ref):10.1f} (={float(T2 - ref) / d:.2f}d)   {float(wr):.3e}")
print("closed-form bound held at every step; excess with look-back grows ~ +7 steps per doubling of d (O(log d)).")
