#!/usr/bin/env python3
"""Family 135, IMM side (Sec. 4-5): exact transfer-matrix evaluation of the two-letter
word sum (eq. moment-reduction) with *refined* weights (more favourable than the paper):
   switch n^{-1};  D_V alpha^{-1};  D_U beta;
   N_V 1+alpha^{-1} n^{-h};  N_U 1+beta n^{-h}   (paper bounds N_U by 1+beta).
Schedule (4.1): V U^{l_1} ... V U^{l_k}, s = sigma sqrt n, rho = k/m, beta = alpha^rho.
log(word sum) = log( T_2 D / T^2 ) up to e^{O(1)}  (i!=j relaxation and occupation
comparison are e^{O(1)}), so it is the exact loss of the trace method.

Test: alpha = n^{-gamma}.  gamma=1/2: loss O(sqrt n) (paper).  gamma>1/2 would scale the
circuit-side saving by 2 gamma, gain (2gamma-1)/(4 sqrt2) * sqrt n log n; the loss must
stay below it.  Prediction from run factorisation d(I)=1+rho(phi_end-phi_start):
log Z ~ 2 sqrt(x), x ~ alpha^{-2}/(8 log^2(1/alpha))  => ~ n^gamma/(sqrt2 gamma log n).
"""
import math, sys

def schedule(n, s):
    k = (n - s) // 2
    lay = []
    for i in range(1, k + 1):
        l = 1 + (i * s) // k - ((i - 1) * s) // k
        lay.append('V'); lay += ['U'] * l
    assert len(lay) == n, (len(lay), n)
    return lay, k

def logZ(n, gamma, sigma=1 / math.sqrt(2)):
    s = int(math.floor(sigma * math.sqrt(n))); s -= (s - n) % 2
    lay, k = schedule(n, s)
    m = n - k; rho = k / m
    al = n ** (-gamma); be = al ** rho
    ia = 1 / al
    def w(t, prev, cur, nxt):
        # weight of layer t of type cur, neighbours prev/nxt (None = external vertex)
        V = lay[t] == 'V'
        if cur == 'D':
            return ia if V else be
        h = (prev == 'N') + (nxt == 'N')   # internal end vertices between two N layers
        return 1 + (ia if V else be) * n ** (-h)
    # DP over (tau_{t-1}, tau_t), log-scaled
    T = ('N', 'D')
    st = {(None, a): 0.0 for a in T}   # log weight before finalising layer 0
    for t in range(0, n):
        new = {}
        nxts = T if t < n - 1 else (None,)
        for (p, c), lw in st.items():
            for x in nxts:
                lw2 = lw + math.log(w(t, p, c, x))
                if x is not None and x != c:
                    lw2 -= math.log(n)
                key = (c, x)
                if key in new:
                    a, b = new[key], lw2
                    mx = max(a, b); new[key] = mx + math.log(math.exp(a - mx) + math.exp(b - mx))
                else:
                    new[key] = lw2
        st = new
    vals = list(st.values()); mx = max(vals)
    return mx + math.log(sum(math.exp(v - mx) for v in vals))

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [10**4, 10**5, 10**6]
    for n in ns:
        n += n % 2
        L = math.log(n)
        for g in [0.5, 0.55, 0.6, 0.7, 0.75]:
            lz = logZ(n, g)
            gain = (2 * g - 1) / (4 * math.sqrt(2))
            pred = n ** g / (math.sqrt(2) * g * L)
            print(f"n={n:8d} gamma={g:.2f}: logZ={lz:10.2f}  logZ/sqrt n={lz/math.sqrt(n):7.3f}"
                  f"  logZ/(sqrt n log n)={lz/(math.sqrt(n)*L):7.4f}  needed<= gain {gain:.4f}"
                  f"  heuristic n^g/(sqrt2 g log n)={pred:9.1f}", flush=True)
