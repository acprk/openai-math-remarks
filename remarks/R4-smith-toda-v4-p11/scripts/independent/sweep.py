# R4 independent re-implementation (z3 exact enumeration + int64 DP mod 2^61-1).
# Run from scripts/independent/: python sweep.py 11   with python = a venv having z3-solver and numpy (tested z3 4.x/5.1.0).
# Expected output: expected_output/sweep_p{7,11,13,17}.txt
# Sweep every May-E1 enumeration claim of family 308 at a given prime p.
# usage: python sweep.py p
import sys, itertools, numpy as np
from mayenum import u, letters, enum_z3, name

MOD = (1 << 61) - 1

def dp_table(p, smax, wmax):
    H, B = letters(p, wmax)
    f = np.zeros((smax + 1, wmax + 1), dtype=np.int64)
    f[0, 0] = 1
    for _, wt in H:
        g = f.copy()
        f[1:, wt:] = (f[1:, wt:] + g[:-1, :-wt]) % MOD
    for _, wt in B:
        for t in range(2, smax + 1):
            f[t, wt:] = (f[t, wt:] + f[t - 2, :-wt]) % MOD
    return f

def main(p):
    q = 2 * p - 2
    U = [u(p, i) for i in range(7)]
    rep = {}
    # ---------- (1) Section 2, Lemma low:algebraic-vanishing ----------
    jobs = []
    for n in range(0, 4):
        for e in itertools.product(range(3), repeat=n + 1):
            rho = sum(e)
            if rho < 2: continue
            Uv = sum(ei * U[i] for i, ei in enumerate(e))
            kc = (Uv + rho - 1) // (q - 1)
            for k in range(1, kc + 1):
                s, w = q * k + 1 - rho, Uv + k
                jobs.append(('prod n=%d e=%s' % (n, e), s, w, w < U[n + 1]))
    for n in range(1, 4):
        kc = (U[n] - 1) // (q - 1)
        for k in range(1, kc + 1):
            jobs.append(('v%d target k=%d' % (n, k), q * k + 1, U[n] + k, None))
    # ---------- (2) v4 targets, all k up to normalized cutoff ----------
    kc4 = (U[4] - 1) // (q - 1)
    for k in range(1, kc4 + 1):
        jobs.append(('v4 target k=%d' % k, q * k + 1, U[4] + k, None))
    # ---------- (3) Section 5 small stems ----------
    W = 2 * p * p + 4 * p + 1
    jobs.append(('src (5,W)', 5, W, None))
    for delta in (0, 1):
        j = 1 - delta if delta == 0 else 0
        while (q - 1) * j <= W + delta - 6:
            jobs.append(('stem d=%d j=%d' % (delta, j), 6 + q * j, W + delta + j, None))
            j += 1
    smax = max(s for _, s, _, _ in jobs); wmax = max(w for _, _, w, _ in jobs)
    f = dp_table(p, smax, wmax)
    bad = []
    out = []
    for lab, s, w, flag in jobs:
        c = int(f[s, w]) if s >= 0 else 0
        if c:
            S = enum_z3(p, s, w, maxsol=50)
            out.append((lab, s, w, c, [name(p, m) for m in S], flag))
    # k=1 incoming-source bidegree (q+1,u5+1): small s, separate table
    f5 = dp_table(p, q + 1, U[5] + 1)
    c5 = int(f5[q + 1, U[5] + 1])
    S5 = enum_z3(p, q + 1, U[5] + 1, maxsol=200)
    print('p=%d q=%d u4=%d kcut4=%d W=%d jobs=%d' % (p, q, U[4], kc4, W, len(jobs)))
    for o in out:
        print('NONEMPTY', o)
    print('(q+1,u5+1): dpcount=%d z3=%d' % (c5, len(S5)))
    for m in S5: print('   ', name(p, m))
    # A^N divisibility of k=p target
    N = p * (p - 1) - 2
    St = enum_z3(p, q * p + 1, U[4] + p)
    print('k=p target A^N-divisible:', all(dict(m[1]).get(('b', 1, 0), 0) >= N for m in St))
    # product-case additive comparison flag w<u_{n+1}
    print('product w<u_{n+1} all:', all(fl for _, _, _, fl in jobs if fl is not None))

if __name__ == '__main__':
    main(int(sys.argv[1]))
