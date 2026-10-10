#!/usr/bin/python3
"""Polynomial (F_q[x]) analogue of Jacobsthal's function -- exact computation by SAT.
For squarefree M = P_1...P_k in F_q[x] (q prime here), J(M) = max m such that some coset
f + V_m (V_m = polys of degree < m) consists entirely of polynomials sharing a factor with M,
equivalently: exists residues r_i mod P_i with  for all g in V_m  exists i: g = r_i (mod P_i).
Usage: python3 fqx_jacobsthal.py   (runs the experiment table)"""
import itertools, sys, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc

def polymod(a, p, q):           # a,p coefficient lists low->high, p monic
    a = list(a)
    dp = len(p) - 1
    for i in range(len(a) - 1, dp - 1, -1):
        c = a[i] % q
        if c:
            for j in range(dp + 1):
                a[i - dp + j] = (a[i - dp + j] - c * p[j]) % q
    return [x % q for x in a[:dp]] + [0] * max(0, dp - len(a))

def irreducibles(q, d):
    """all monic irreducible polys of degree d over F_q (q prime), by trial division"""
    out = []
    small = []
    for e in range(1, d // 2 + 1):
        small += irreducibles(q, e) if e < d else []
    for coeffs in itertools.product(range(q), repeat=d):
        p = list(coeffs) + [1]
        if d == 1: out.append(p); continue
        if p[0] == 0: continue
        if all(any(polymod(p, s, q)) for s in small): out.append(p)
    return out

def residue_table(P, q, m):
    """index of g mod P for every g in V_m (g enumerated as base-q integer, low coeff first)"""
    d = len(P) - 1
    basis = [polymod([0] * j + [1], P, q) if j >= d else ([0] * j + [1] + [0] * (d - j - 1)) for j in range(m)]
    tab = []
    for coeffs in itertools.product(range(q), repeat=m):
        r = [0] * d
        for j, c in enumerate(reversed(coeffs)):   # coeffs[::-1] low first
            if c:
                b = basis[j]
                r = [(r[t] + c * b[t]) % q for t in range(d)]
        idx = 0
        for t in reversed(range(d)): idx = idx * q + r[t]
        tab.append(idx)
    return tab

def coverable(primes, q, m, timeout=None):
    """True iff some choice of one residue class per prime covers V_m."""
    var = {}
    nv = 0
    for i, P in enumerate(primes):
        for c in range(q ** (len(P) - 1)):
            nv += 1; var[(i, c)] = nv
    s = Cadical153()
    for i, P in enumerate(primes):           # at most one class per prime (at-most-one suffices)
        lits = [var[(i, c)] for c in range(q ** (len(P) - 1))]
        enc = CardEnc.atmost(lits=lits, bound=1, top_id=nv)
        nv = max(nv, enc.nv)
        for cl in enc.clauses: s.add_clause(cl)
    tabs = [residue_table(P, q, m) for P in primes]
    for gi in range(q ** m):
        s.add_clause([var[(i, tabs[i][gi])] for i in range(len(primes))])
    ok = s.solve()
    model = s.get_model() if ok else None
    s.delete()
    if ok:
        mset = set(x for x in model if x > 0)
        choice = [next(c for c in range(q ** (len(P) - 1)) if var[(i, c)] in mset) if any(var[(i, c)] in mset for c in range(q ** (len(P) - 1))) else None for i, P in enumerate(primes)]
        return True, choice
    return False, None

def J(primes, q, mmax=40):
    m = 0
    while m + 1 <= mmax and coverable(primes, q, m + 1)[0]:
        m += 1
    return m

if __name__ == "__main__":
    print("E1: M = product of all monic irreducibles of degree <= D  (deg M, k, J)")
    for q, Dmax in [(2, 5), (3, 3), (5, 2), (7, 1)]:
        primes = []
        for D in range(1, Dmax + 1):
            primes += irreducibles(q, D)
            t = time.time()
            j = J(primes, q, mmax=sum(len(p) - 1 for p in primes))
            print(f"  q={q} D={D}: k={len(primes)} degM={sum(len(p)-1 for p in primes)}  J={j}  ({time.time()-t:.1f}s)", flush=True)
