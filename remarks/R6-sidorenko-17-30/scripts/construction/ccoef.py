"""Exact C(a,b) of kernel.tex eq. (kernel-coefficient) for the base laws
(kernel.tex eq. kernel-base-moment), K = corners of the TARGET complex.

C(a,b) = sum_{alpha_e != {} (e in E)} prod_{source corners kappa=(i,j)}
           A^base_{(a(i),b(j))}(alpha_{e(kappa)}, alpha_{f(kappa)}).

A^base_tau(al,be) = 1/2 [Q<=M_tau]([R={}] + s th [R=G_tau])
                  + 1/2 [Q<=G_tau]([R={}] - s th [R=G_tau]),   Q=al|be, R=al^be,
and A=0 at a nonincident type pair (S = {} deterministically) unless al=be={}.

We return 2^m * C as an integer polynomial in theta (list of coefficients),
m = number of corners.  sigma: dict target-corner -> +-1 (default: one -1 at corner 0).
Exact integer arithmetic throughout.
"""
from itertools import combinations


def padd(p, q):
    if len(p) < len(q):
        p, q = q, p
    r = list(p)
    for k, c in enumerate(q):
        r[k] += c
    return r


def pmul(p, q):
    if not p or not q:
        return []
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    r[i + j] += a * b
    return r


def ptrim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def subsets_nonempty(s):
    s = sorted(s)
    out = []
    for r in range(1, len(s) + 1):
        for c in combinations(s, r):
            out.append(frozenset(c))
    return out


def twoA(cx, tau, al, be, sigma):
    """2*A^base_tau(al,be) as integer poly in theta; tau = target corner index."""
    M, G = cx.M[tau], cx.G[tau]
    Q, R = al | be, al ^ be
    s = sigma.get(tau, 1)
    c0 = 0
    c1 = 0
    if Q <= M:
        if not R:
            c0 += 1
        elif R == G:
            c1 += s
    if Q <= G:
        if not R:
            c0 += 1
        elif R == G:
            c1 -= s
    return ptrim([c0, c1])


def C_poly(src, tgt, a, b, sigma=None, order=None):
    """src, tgt: Complex; a: points(src)->points(tgt); b: faces(src)->faces(tgt).
    Returns integer poly P with C = 2^{-m} P(theta)."""
    if sigma is None:
        sigma = {0: -1}
    # type corner for each source corner
    tc = []
    for (i, j) in src.corners:
        ai, bj = a[i], b[j]
        if ai not in tgt.faces[bj]:
            return []          # nonincident type pair: zero factor for nonempty labels
        tc.append(tgt.cidx[(ai, bj)])
    nE = len(src.pairs)
    # domains
    dom = []
    for e in range(nE):
        X = None
        for k in src.U[e]:
            X = set(tgt.M[tc[k]]) if X is None else X & tgt.M[tc[k]]
        if not X:
            return []
        dom.append(subsets_nonempty(X))
    # factors: (vars tuple, table dict{assignment tuple: poly})
    factors = []
    for k, (e, f) in enumerate(src.cpairs):
        tab = {}
        for x, al in enumerate(dom[e]):
            for y, be in enumerate(dom[f]):
                p = twoA(tgt, tc[k], al, be, sigma)
                if p:
                    tab[(x, y)] = p
        if not tab:
            return []
        factors.append(((e, f), tab))
    return eliminate(factors, nE, [len(d) for d in dom])


def eliminate(factors, nvars, dsz):
    remaining = set(range(nvars))
    while remaining:
        # min-fill-ish: choose var minimizing product of neighbor domain sizes
        best, bcost = None, None
        for v in remaining:
            nb = set()
            for vs, _ in factors:
                if v in vs:
                    nb |= set(vs)
            nb.discard(v)
            cost = 1
            for u in nb:
                cost *= dsz[u]
            if bcost is None or cost < bcost:
                best, bcost = v, cost
        v = best
        inv = [fc for fc in factors if v in fc[0]]
        rest = [fc for fc in factors if v not in fc[0]]
        scope = sorted(set(u for vs, _ in inv for u in vs) - {v})
        # sparse join: iterate over assignments built factor by factor
        partial = {(): []}  # assignment over 'cur' vars (dict as tuple of (var,val))
        cur = []
        # represent partial as dict: frozen tuple of values for vars in cur + v
        allv = scope + [v]
        pos = {u: k for k, u in enumerate(allv)}
        part = {tuple([None] * len(allv)): [1]}
        for vs, tab in inv:
            newp = {}
            for asg, pv in part.items():
                for key, fv in tab.items():
                    ok = True
                    na = list(asg)
                    for u, val in zip(vs, key):
                        cv = na[pos[u]]
                        if cv is None:
                            na[pos[u]] = val
                        elif cv != val:
                            ok = False; break
                    if ok:
                        t = tuple(na)
                        pr = pmul(pv, fv)
                        if t in newp:
                            newp[t] = padd(newp[t], pr)
                        else:
                            newp[t] = pr
            part = newp
        # sum out v (v must be assigned: it's in every inv factor)
        out = {}
        for asg, pv in part.items():
            key = asg[:-1]
            if None in key:
                raise RuntimeError('unassigned scope var')
            out[key] = padd(out.get(key, []), pv)
        out = {k: ptrim(p) for k, p in out.items() if ptrim(p)}
        remaining.discard(v)
        if scope:
            rest.append((tuple(scope), out))
        else:
            rest.append(((), {(): out.get((), [])}))
        factors = rest
    res = [1]
    for vs, tab in factors:
        assert vs == ()
        res = pmul(res, tab.get((), []))
    return ptrim(res)
