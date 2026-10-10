# R4 main computation. Run from scripts/: /usr/bin/python3 main/may308.py 11 13 17   (stdlib only; ~2 s per prime)
# Negative control: main/may308.py 7 (must print FAILS ['A']).  Expected output: expected_output/may308_p*.txt
#!/usr/bin/env python3
"""308 / V(4) at p: exact re-verification of every finite May-page claim of
'A finite Smith-Toda complex V(4) at the prime 1009' for a general odd prime p.
May E1 page (Lemma may:page): exterior h_{i,j} (s=1, w=p^j u_i), polynomial b_{i,j} (s=2, w=p^{j+1} u_i).
All arithmetic is exact integer arithmetic."""
import sys, itertools
from functools import lru_cache

def setup(p):
    u = lambda i: (p**i - 1)//(p - 1)
    return u, 2*p - 2

def gens(p, wmax, rmax=None):
    u, q = setup(p)
    H, B = [], []
    for i in range(1, 12):
        if rmax is not None and i > rmax: break
        for j in range(0, 12):
            if p**j*u(i) <= wmax: H.append(((i, j), p**j*u(i)))
            if p**(j+1)*u(i) <= wmax: B.append(((i, j), p**(j+1)*u(i)))
    return H, B

def monomials(p, s, w, rmax=None, limit=10**6):
    """all May E1 monomials of bidegree (s,w): returns list of (frozenset ext, tuple poly-counts dict)"""
    H, B = gens(p, w, rmax)
    B = sorted(B, key=lambda t: -t[1])
    out = []
    nH = len(H)
    for mask in range(1 << nH):
        ext = [H[t] for t in range(nH) if mask >> t & 1]
        m = len(ext)
        if m > s or (s - m) % 2: continue
        L = (s - m)//2
        W = w - sum(e[1] for e in ext)
        if W < 0: continue
        # distribute L letters among B with total weight W
        def rec(idx, Lr, Wr, acc):
            if len(out) > limit: return
            if idx == len(B):
                if Lr == 0 and Wr == 0: out.append((tuple(e[0] for e in ext), tuple(acc)))
                return
            g, wt = B[idx]
            if idx == len(B) - 1:
                if Lr*wt == Wr: rec(idx+1, 0, 0, acc + [(g, Lr)] if Lr else acc)
                return
            minw = B[-1][1]
            # need Wr - n*wt >= (Lr-n)*minw and <= (Lr-n)*wt_next
            wn = B[idx+1][1]
            for n in range(min(Lr, Wr//wt), -1, -1):
                Lr2, Wr2 = Lr - n, Wr - n*wt
                if Wr2 < Lr2*minw: continue
                if Wr2 > Lr2*wn: break
                rec(idx+1, Lr2, Wr2, acc + [(g, n)] if n else acc)
        if L == 0:
            if W == 0: out.append((tuple(e[0] for e in ext), tuple()))
        elif B:
            rec(0, L, W, [])
    return out

NAMES = {(1,0):'x',(2,0):'a',(3,0):'c',(1,1):'y',(2,1):'b',(1,2):'z'}
PN = {(1,0):'A',(2,0):'B',(1,1):'C'}
def show(mono):
    ext, poly = mono
    e = ''.join(NAMES.get(g, 'h%d%d' % g) for g in sorted(ext, key=lambda g:(g[1], g[0])))
    pp = ' '.join('%s^%d' % (PN.get(g, 'b%d%d' % g), n) for g, n in poly)
    return (e or '1') + ' ' + pp

# ---------- splitting differential d1 on exterior parts (Lemma may:page) ----------
def interval(g): i, j = g; return (j, j+i)
def tog(a, b): return (b-a, a)
def d1_ext(ext, p):
    """ext: tuple of generators in a fixed order; returns dict {tuple(sorted canonical): coeff mod p}"""
    order = lambda g: (interval(g)[0], interval(g)[1])
    res = {}
    for pos, g in enumerate(ext):
        a, b = interval(g)
        for l in range(a+1, b):
            g1, g2 = tog(a, l), tog(l, b)
            new = list(ext[:pos]) + [g1, g2] + list(ext[pos+1:])
            if len(set(new)) < len(new): continue
            sign = (-1)**pos
            # sort with sign
            arr = new[:]; sg = 1
            for i in range(len(arr)):
                for j in range(len(arr)-1-i):
                    if order(arr[j]) > order(arr[j+1]): arr[j], arr[j+1] = arr[j+1], arr[j]; sg = -sg
            key = tuple(arr)
            res[key] = (res.get(key, 0) + sign*sg) % p
    return {k: v for k, v in res.items() if v}

def canon(ext):
    order = lambda g: (interval(g)[0], interval(g)[1])
    return tuple(sorted(ext, key=order))

def rank_mod_p(rows, p):
    rows = [r[:] for r in rows]; rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p-2, p)
        rows[rk] = [x*inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(x - f*y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def d1_cohomology_dim(p, s, w, rmax=None):
    """dim of H(d1) at (s,w) on the May E1 page (d1 acts on exterior part only; poly parts are d1-cycles)"""
    def basis(ss): return [(canon(e), pol) for e, pol in monomials(p, ss, w, rmax)]
    Bm, B0, Bp = basis(s-1), basis(s), basis(s+1)
    idx0 = {b: i for i, b in enumerate(B0)}; idxp = {b: i for i, b in enumerate(Bp)}
    def mat(src, tgt_idx, ntgt):
        rows = []
        for e, pol in src:
            r = [0]*ntgt
            for k, v in d1_ext(e, p).items():
                r[tgt_idx[(k, pol)]] = (r[tgt_idx[(k, pol)]] + v) % p
            rows.append(r)
        return rows
    r_in = rank_mod_p(mat(Bm, idx0, len(B0)), p) if Bm and B0 else 0
    r_out = rank_mod_p(mat(B0, idxp, len(Bp)), p) if B0 and Bp else 0
    return len(B0) - r_out - r_in, len(B0), r_in, r_out

def check(p, verbose=True):
    u, q = setup(p)
    ok = True; fails = []
    log = lambda *a: verbose and print(*a)
    log(f"===== p={p}  q={q}  u=(0,1,{u(2)},{u(3)},{u(4)})")
    # ---- (A) Lemma low:algebraic-vanishing ----
    bad = []
    for n in range(1, 4):  # v_n targets, coefficients R/I_n
        k = 1
        while (q-1)*k <= u(n) - 1:
            ms = monomials(p, q*k+1, u(n)+k)
            if ms: bad.append(('v%d-target' % n, k, [show(m) for m in ms][:4], len(ms)))
            k += 1
    for n in range(0, 4):  # product obstructions, coefficients R/I_{n+1}
        for e in itertools.product(range(3), repeat=n+1):
            rho = sum(e)
            if rho < 2: continue
            U = sum(ei*u(i) for i, ei in enumerate(e))
            k = 1
            while (q-1)*k <= U + rho - 1:
                s, w = q*k+1-rho, U+k
                if w >= u(n+1): bad.append(('product-comparison-fails', n, e, k))
                ms = monomials(p, s, w)
                if ms: bad.append(('product n=%d' % n, e, k, (s, w), [show(m) for m in ms][:4], len(ms)))
                k += 1
    log("(A) lower obstruction groups (low:algebraic-vanishing) empty May page:", not bad)
    for b in bad[:8]: log("    FAIL", b)
    if bad: ok = False; fails.append('A')
    # ---- (B) v_4 targets (may:target-bound, may:first-list k=1, may:remaining-list) ----
    N0, N = (p-1)**2, p*(p-1)-2
    exp = {
        1: {('acb', (('B', p-2),)), ('acybz', (('B', p-3),))},
        2: {('xac', (('B', 1), ('A', p-2), ('C', p-2))), ('xacyz', (('A', p-2), ('C', p-2)))},
        p-1: {('y', (('B', 3), ('A', N0-3))), ('b', (('B', 2), ('A', N0-2))), ('ybz', (('B', 1), ('A', N0-2)))},
        p: {('cyb', (('A', N+1),)), ('ayb', (('B', 1), ('A', N)))},
    }
    def norm(m):
        ext, pol = m
        e = ''.join(sorted((NAMES.get(g, '?%d%d' % g) for g in ext), key='xacybz'.index)) if all(g in NAMES for g in ext) else 'OTHER:' + str(ext)
        pl = tuple(sorted(((PN.get(g, 'b%d%d' % g), n) for g, n in pol), key=lambda t: 'BAC'.index(t[0]) if t[0] in 'BAC' else 9))
        return (e, pl)
    expn = {k: {(e, tuple(sorted(pl, key=lambda t: 'BAC'.index(t[0])))) for e, pl in v} for k, v in exp.items()}
    found = {}
    k = 1
    import os
    krange = range(1, (u(4)-1)//(q-1)+1) if not os.environ.get('SKIPB') else [1, 2, 3, p-2, p-1, p, p+1, p+2]
    for k in krange:
        ms = monomials(p, q*k+1, u(4)+k)
        if ms: found[k] = {norm(m) for m in ms}
    good = set(found) == set(expn) and all(found[k] == expn[k] for k in found)
    log("(B) v_4 targets: nonempty k =", sorted(found), " matches paper tables:", good)
    if not good:
        ok = False; fails.append('B')
        for k in sorted(set(found) | set(expn)):
            if found.get(k) != expn.get(k):
                log("    k=%d found %s" % (k, sorted(found.get(k, []))[:6]), "\n         expected", sorted(expn.get(k, [])))
    # k=1 cycle structure: acb not d1-cycle, acybz d1-cycle
    # ---- (C) Table may:n5-table at (q+1, u5+1) and injectivity of splitting ----
    ms5 = monomials(p, q+1, u(5)+1)
    iv = lambda ext: ','.join('%d%d' % interval(g) for g in canon(ext))
    lst5 = sorted((iv(e), tuple(pol)) for e, pol in ms5)
    exp5 = sorted([('02,04,14', (((3,0), p-2),)), ('03,04,12,14,34', (((3,0), p-3),)), ('02,04,13,14,34', (((3,0), p-3),)),
                   ('02,04,12,14,24', (((3,0), p-3),)), ('02,04,12,13,14,24,34', (((3,0), p-4),))])
    goodC = lst5 == exp5
    # injectivity of d1 on span of the length-3 entry and on span of the three length-5 entries
    inj = True
    for grp in ([m for m in ms5 if len(m[0]) == 3], [m for m in ms5 if len(m[0]) == 5]):
        imgs = [d1_ext(canon(e), p) for e, pol in grp]
        keys = sorted({k for d in imgs for k in d})
        rows = [[d.get(k, 0) for k in keys] for d in imgs]
        if rank_mod_p(rows, p) != len(grp): inj = False
    log("(C) (q+1,u5+1) list = Table may:n5-table:", goodC, "; d1 injective on m=3 and m=5 spans:", inj)
    if not (goodC and inj): ok = False; fails.append('C'); log("    found", lst5)
    # filtration comparison F7-F* with M=3p (Prop may:injectivity)
    M = 3*p
    Fs = 3*p*p - 9*p + 14 + M*(2 + p*(p-3)); F7 = 5*p*p - 20*p + 23 + M*(2 + p*(p-4))
    log("    F7-F* =", F7 - Fs, "(<0 needed)")
    if not F7 < Fs: ok = False; fails.append('C-filt')
    # ---- (D) Lemma frob:quotient-groups (P_3 May page at rescaled weights) ----
    tab = [((p*(p-1), 2*p-4), []), ((p*(p-1), 2*p-3), ['h11 b10^%d' % (p-2)]), ((p*(p-2), 2*p-4), ['b10^%d' % (p-2)]),
           ((p*(p-3), 2*p-5), []), ((p+2, 1), []), ((p+2, 3), []), ((p+2, 2), ['h10 h20']), ((2*p+2, 2), []),
           ((2*p+2, 3), ['h10 h11 h20']), ((3*p+2, 3), []), ((2, 1), []), ((2, 3), [])]
    goodD = True
    for (w, s), e in tab:
        ms = monomials(p, s, w, rmax=3)
        got = sorted(' '.join(sorted('h%d%d' % g for g in ext)) + ('' if not pol else (' ' if ext else '') + ' '.join('b%d%d^%d' % (g[0], g[1], n) for g, n in pol)) for ext, pol in ms)
        if got != sorted(e): goodD = False; log("    D mismatch at", (w, s), got, e)
        if not w < p*p: goodD = False; log("    D weight not < p^2", (w, s))
    log("(D) frob:quotient-groups table:", goodD, "; 3p+2<p(p-1):", 3*p+2 < p*(p-1), "; units 2,3,6 mod p:", all(c % p for c in (2, 3, 6)))
    if not goodD: ok = False; fails.append('D')
    # Replacing A by E: p^2+p exceeds residual weights p(p-1), p^2-p+1 (always)
    # ---- (E) Lemma top:small-stems ----
    W = 2*p*p + 4*p + 1
    goodE = W < u(4)
    src = sorted(show(m) for m in monomials(p, 5, W))
    badE = []
    for delta in (0, 1):
        j = 0 if delta else 1
        while (q-1)*j <= W + delta - 6:
            s, w = 6 + q*j, W + delta + j
            if w >= u(4): badE.append(('coef', j, delta))
            ms = monomials(p, s, w)
            if ms and not (j == 0 and delta == 1): badE.append((j, delta, [show(m) for m in ms][:4]))
            j += 1
    four = sorted(show(m) for m in monomials(p, 6, W + 1))
    h, n0, rin, rout = d1_cohomology_dim(p, 6, W + 1)
    log("(E) (5,W) source list:", src, "; other stems empty:", not badE, "; (6,W+1) list:", four, "; d1-cohomology dim:", h)
    if badE or h != 0 or not goodE: ok = False; fails.append('E'); log("    ", badE[:5])
    # ---- (F) page bound top:early-disappearance ----
    K = (p-1)**2; S = 5 + 2*K
    goodF = (q*p + 1 - S > 0) and (N - K > 0)
    log("(F) qp+1-S=%d, N-K=%d :" % (q*p+1-S, N-K), goodF)
    if not goodF: ok = False; fails.append('F')
    log("RESULT p=%d: %s %s" % (p, 'ALL MAY-PAGE / FINITE CLAIMS HOLD' if ok else 'FAILS', fails))
    return ok, fails

if __name__ == '__main__':
    for p in map(int, sys.argv[1:]):
        check(p)
