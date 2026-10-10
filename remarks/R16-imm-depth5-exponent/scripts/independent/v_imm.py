#!/usr/bin/env python3
"""Independent IMM-side check (Sec.4-5) for the modified s<=sqrt(n/2).
(i) exact word sum W_exact = sum_tau n^{-w} alpha^{-dV} beta^{dU} prod_{N_V}(1+alpha^{-1} n^{-h_t})
    by a transfer matrix over (type_{t-1},type_t), compared with the paper's bound
    (1+2/n^2)^k (1+1/alpha)(1+1/(alpha n))^k * e^{m/b0} (2+6 n^2 e^{sqrt n});
(ii) explicit C_I(n) = 2*loss1 + m log(1+beta) + log(prefactor) + log(Zbound) and
     the constant C_I - 2 sqrt n - c log n."""
import math, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v_params_H import params
from mpmath import mpf, log, sqrt, exp, log1p

def schedule(P):
    n, s, k = P['n'], P['s'], P['k']
    lay = []
    for i in range(1, k + 1):
        l = 1 + (i * s) // k - ((i - 1) * s) // k
        assert l in (1, 2)
        lay += ['V'] + ['U'] * l
    assert len(lay) == n and lay.count('V') == k
    assert all(not (lay[i] == 'V' and lay[i + 1] == 'V') for i in range(n - 1))
    return lay

def lse(a, b):
    m = max(a, b); return m + math.log(math.exp(a - m) + math.exp(b - m))

def word_sum(P):
    n = P['n']; lay = schedule(P)
    ia = float(1 / P['al']); be = float(P['be']); ln = math.log(n)
    # state: (type_{t-1}, type_t) with weight of layers < t finalised; layer t's weight
    # needs type_{t+1}.  None = outside.
    st = {(None, 'N'): 0.0, (None, 'D'): 0.0}
    for t in range(n):
        nxt_opts = ('N', 'D') if t < n - 1 else (None,)
        new = {}
        for (p, c), lw in st.items():
            for x in nxt_opts:
                if c == 'D':
                    wt = ia if lay[t] == 'V' else be
                else:
                    if lay[t] == 'V':
                        h = (p == 'N') + (x == 'N')   # internal vertex with both sides N
                        wt = 1 + ia * n ** (-h)
                    else:
                        wt = 1.0   # exact geometric weight is 1+beta*1{i=j}; paper bounds by 1+beta outside
                val = lw + math.log(wt) - (ln if (x is not None and x != c) else 0.0)
                key = (c, x)
                new[key] = lse(new[key], val) if key in new else val
        st = new
    v = list(st.values()); r = v[0]
    for y in v[1:]: r = lse(r, y)
    return r

def CI_parts(P):
    n, k, m, a, b, v, u = P['n'], P['k'], P['m'], P['a'], P['b'], P['v'], P['u']
    al, be = P['al'], P['be']
    loss1 = mpf(k)**2 / a + mpf(k)**2 / (2 * v) + mpf(m)**2 / b + mpf(m)**2 / (2 * u)
    b0 = u * al ** P['rho'] / (1 - al ** P['rho'])
    logpref = k * log1p((1 + al) / n**2) + log1p(1 / al) + k * log1p(1 / (al * n))
    logZb = m / b0 + log(2 + 6 * mpf(n)**2 * exp(sqrt(mpf(n))))
    mlogb = m * log1p(be)
    CI = 2 * loss1 + mlogb + logpref + logZb
    return dict(CI=CI, loss1=loss1, mlogb=mlogb, logpref=logpref, logZb=logZb)

if __name__ == "__main__":
    for n in [2000, 5001, 20000, 60001]:
        P = params(n)
        ws = word_sum(P)
        d = CI_parts(P)
        bound = float(d['logpref'] + d['logZb'])
        print(f"n={n}: log exact word sum={ws:.2f}  paper bound={bound:.2f}  ok={ws <= bound}  "
              f"(exact/sqrt n={ws/math.sqrt(n):.3f})")
        assert ws <= bound
    print("explicit C_I:  C_I - 2sqrt n - 2.5 log n  and  - (2.5+1/(2sqrt2)) log n")
    for L in [20, 40, 57.94, 100, 200, 296.7]:
        n = int(math.exp(L)) if L < 600 else None
        P = params(n)
        d = CI_parts(P)
        rn = sqrt(mpf(n)); Ln = log(mpf(n))
        r1 = d['CI'] - 2 * rn - mpf(5) / 2 * Ln
        r2 = r1 - Ln / (2 * sqrt(2))
        print(f" log n={L}: C_I/sqrt n={float(d['CI']/rn):.6f}  C_I-2sqrt n-2.5logn={float(r1):.3f}"
              f"  minus (log n)/(2sqrt2) -> {float(r2):.3f}   mlog(1+b)-sqrt n/2={float(d['mlogb']-rn/2):.3f}")
