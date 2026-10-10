# Shared module for sweep.py / d1check.py / t1009.py.
# Independent May E1 enumerator for OAI family 308 (V(4) at p=1009) re-check.
# E1 = Lambda(h_{i,j}) (x) F_p[b_{i,j}],  h_{i,j}: (1, p^j u_i),  b_{i,j}: (2, p^{j+1} u_i).
# Method 1 (any p): z3 exact integer enumeration with blocking clauses.
# Method 2 (small p): DP count over (s,w) -- independent cross-check of completeness.
import z3, itertools, functools

def u(p, i): return (p**i - 1)//(p - 1)

def letters(p, w):
    H, B = [], []
    i = 1
    while u(p, i) <= w:
        j = 0
        while p**j*u(p, i) <= w:
            H.append((('h', i, j), p**j*u(p, i)))
            j += 1
        j = 0
        while p**(j+1)*u(p, i) <= w:
            B.append((('b', i, j), p**(j+1)*u(p, i)))
            j += 1
        i += 1
    return H, B

def enum_z3(p, s, w, maxsol=200):
    if s < 0 or w < 0: return []
    H, B = letters(p, w)
    e = [z3.Int('e%d' % k) for k in range(len(H))]
    n = [z3.Int('n%d' % k) for k in range(len(B))]
    S = z3.Solver()
    for v in e: S.add(v >= 0, v <= 1)
    for v in n: S.add(v >= 0)
    S.add(z3.Sum([z3.IntVal(0)] + e) + 2*z3.Sum([z3.IntVal(0)] + n) == s)
    S.add(z3.Sum([z3.IntVal(0)] + [wt*v for (_, wt), v in zip(H, e)])
          + z3.Sum([z3.IntVal(0)] + [wt*v for (_, wt), v in zip(B, n)]) == w)
    sols = []
    while S.check() == z3.sat and len(sols) < maxsol:
        M = S.model()
        ev = [M.eval(v).as_long() for v in e]
        nv = [M.eval(v).as_long() for v in n]
        mono = (tuple(H[k][0] for k in range(len(H)) if ev[k]),
                tuple((B[k][0], nv[k]) for k in range(len(B)) if nv[k]))
        sols.append(mono)
        S.add(z3.Or([v != val for v, val in zip(e + n, ev + nv)]))
    return sols

def dp_count(p, s, w):
    """exact number of E1 monomials in (s,w), by DP (small w only)."""
    H, B = letters(p, w)
    # f[t][x] = number of ways; t <= s, x <= w
    f = [[0]*(w+1) for _ in range(s+1)]
    f[0][0] = 1
    for _, wt in H:          # exterior: 0/1
        for t in range(s, 0, -1):
            ft, fp = f[t], f[t-1]
            for x in range(w, wt-1, -1):
                if fp[x-wt]: ft[x] += fp[x-wt]
    for _, wt in B:          # polynomial: unbounded, s-step 2
        for t in range(2, s+1):
            ft, fp = f[t], f[t-2]
            for x in range(wt, w+1):
                if fp[x-wt]: ft[x] += fp[x-wt]
    return f[s][w]

NAMES = {}
def name(p, mono):
    ext, pol = mono
    al = {('h',1,0):'x',('h',2,0):'a',('h',3,0):'c',('h',1,1):'y',('h',2,1):'b',('h',1,2):'z',
          ('b',1,0):'A',('b',2,0):'B',('b',1,1):'C'}
    out = []
    for g in ext:
        out.append(al.get(g, 'h[%d,%d)' % (g[2], g[2]+g[1])))
    for g, k in pol:
        out.append(al.get(g, 'b%d%d' % (g[1], g[2])) + ('^%d' % k if k > 1 else ''))
    return ' '.join(out)
