"""Family 140, 'Subsphere methods for memory-sample lower bounds in noiseless Gaussian
regression', build/sections/07-beta-mixture.tex, Lemma s:mixture (l.21-48).
NUMERICAL ILLUSTRATION ONLY (the remark relies on the paper's own proof, l.49-106).

The lemma reduces to a comparison of the laws of Q = |P_{V^perp} y|^2:
   actual:  Beta(alpha, beta0)      nu0: Beta(alpha_*, beta0)     nu1: Beta(alpha, beta1)
with alpha=(n-j)/2, alpha_*=(D'-j)/2, beta0=j/2, beta1=(D'-n+j-1)/2, b=D'+1-n >= 2.
The paper states Law(v) <= 2^{D'} (nu0+nu1); its proof bounds the two density ratios by
2^{(b-1)/2} (q >= 1/2) and 2^{(b-2)/2} (q < 1/2).  For every 1 <= j < n < D' <= 40 we
compute on a grid of ~6000 points of (0,1) (uniform + log-spaced towards both ends,
binary64, log-densities via math.lgamma)
   L = log2 sup_q f_{alpha,beta0}(q) / (f_{alpha_*,beta0}(q) + f_{alpha,beta1}(q))
and check L <= (b-1)/2 (+1e-9 float tolerance).  Since (b-1)/2 <= (D'-1)/2 < D', this is
the sharper constant.  Environment: python3, numpy. A few seconds, single thread.
"""
import math
import numpy as np

u = np.linspace(0, 1, 4001)[1:-1]
e = 10.0 ** -np.arange(3, 15, 0.01)
q = np.unique(np.concatenate([u, e, 1 - e]))
lq, l1q = np.log(q), np.log1p(-q)

def logpdf(A, B):
    return (A - 1) * lq + (B - 1) * l1q - (math.lgamma(A) + math.lgamma(B) - math.lgamma(A + B))

def L(Dp, n, j):
    a, a_s, b0, b1 = (n - j) / 2, (Dp - j) / 2, j / 2, (Dp - n + j - 1) / 2
    f, g0, g1 = logpdf(a, b0), logpdf(a_s, b0), logpdf(a, b1)
    m = np.maximum(g0, g1)
    r = f - (m + np.log(np.exp(g0 - m) + np.exp(g1 - m)))
    return float(r.max() / math.log(2))

count = viol = 0; worst = None
for Dp in range(3, 41):
    for n in range(2, Dp):
        for j in range(1, n):
            b = Dp + 1 - n
            if j == 1:
                assert (Dp - n + j - 1) / 2 == (Dp - j) / 2 - (n - j) / 2   # beta1 = alpha_* - alpha
            l2 = L(Dp, n, j); count += 1
            slack = (b - 1) / 2 - l2
            if slack < -1e-9: viol += 1
            if worst is None or slack < worst[0]: worst = (slack, Dp, n, j, l2)
print("grid points:", q.size)
print("parameter triples (D',n,j), 1<=j<n<D'<=40:", count)
print("violations of log2 C <= (b-1)/2:", viol)
s, Dp, n, j, l2 = worst
print(f"smallest slack (b-1)/2 - log2 C = {s:.4f} at D'={Dp}, n={n}, j={j} (log2 C = {l2:.4f})")
for (Dp, n, j) in [(40, 38, 1), (40, 38, 20), (40, 30, 5), (40, 20, 10), (40, 3, 1)]:
    print(f"D'={Dp:2d} n={n:2d} j={j:2d}:  log2 C = {L(Dp, n, j):8.4f}   (b-1)/2 = {(Dp - n) / 2:5.1f}   stated: D' = {Dp}")
