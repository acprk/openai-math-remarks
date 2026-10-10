#!/usr/bin/python3
"""Independent re-verification for remark R6 (OpenAI family 161, Sidorenko counterexample).

Checks, for every simple triangulation of the 2-sphere on V = 4..12 points (and V = 13 for
the classification / positive control), the four conditions of Proposition prop:complex
(complex.tex L73-89 of the source) in the following form:
  (1) every pair in exactly two faces, face-neighbour graph connected   [automatic here]
  (2) three "bit splits": partitions of the faces into two classes, each class having an
      exposure order (every face after the first shares a pair with an earlier face and
      adds exactly one new point)
  (3) every pair separated by at least one bit, some pair by at least two
  (4) colour refinement of the point-face incidence graph is discrete.
Fractional version of (2)+(3): a probability vector w on admissible splits with
P_w(split separates e) >= 1/3 for all pairs e.  Fractional infeasibility implies integral
infeasibility (take w uniform on the three bit splits).

Implementation choices deliberately differ from the primary computation (primary/):
  * plantri 5.8 in binary planar_code (primary: plantri 5.5, ASCII -a output);
  * faces obtained by dart face-tracing of the rotation system (primary: consecutive
    neighbours around each vertex);
  * own colour-refinement code;
  * admissible splits enumerated as Hamiltonian cycles of the point graph (Lemma A.1 of
    note.tex), cross-checked against brute-force face-subset enumeration for V <= 10;
  * LP solved in primal form; the infeasibility certificate is read off the dual values and
    then verified in exact rational arithmetic (fractions.Fraction) against EVERY cycle.
    Floating point is used only to *find* certificates, never to accept them.

Environment: Python 3.10, numpy, scipy>=1.9 (HiGHS).  Needs plantri 5.8 built by
../get_plantri.sh 58.  Single-threaded.  Runtime: about 3-6 minutes on one core.
Usage:  python3 r6_verify.py [--vmax 12] [--no13]
Expected output (with --vmax 10 --no13): ../expected_output/r6_verify_vmax10.txt
"""
import sys, os, subprocess, itertools, argparse
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
PLANTRI = os.path.join(HERE, '..', '_build', 'plantri58', 'plantri')
OEIS_A000109 = {4: 1, 5: 1, 6: 2, 7: 5, 8: 14, 9: 50, 10: 233, 11: 1249, 12: 7595, 13: 49566}

# The 13-point complex of the source, transcribed from Table tab:complex (complex.tex L29-50)
PAPER_FACES = [(0,1,3),(0,1,9),(0,2,3),(0,2,9),(1,3,10),(1,7,10),(1,7,12),(1,9,12),(2,3,4),
               (2,4,9),(3,4,11),(3,10,11),(4,8,9),(4,8,11),(5,6,7),(5,6,11),(5,7,10),(5,10,11),
               (6,7,12),(6,8,11),(6,8,12),(8,9,12)]
PAPER_BITS = ['000','010','001','110','001','000','100','110','101','100','111','101','101',
              '011','000','110','010','111','001','010','011','111']


def planar_code(V):
    """Yield rotation systems (list of neighbour lists, 0-based) from plantri's binary output."""
    data = subprocess.run([PLANTRI, str(V)], capture_output=True, check=True).stdout
    hdr = b'>>planar_code<<'
    assert data.startswith(hdr)
    p = len(hdr)
    while p < len(data):
        n = data[p]; p += 1
        assert n == V
        rot = []
        for _ in range(n):
            nb = []
            while data[p] != 0:
                nb.append(data[p] - 1); p += 1
            p += 1
            rot.append(nb)
        yield rot


def faces_by_tracing(rot):
    """Trace faces of the rotation system: dart (u,v) -> (v, w), w = successor of u at v."""
    pos = [{x: i for i, x in enumerate(nb)} for nb in rot]
    seen = set(); faces = []
    for u in range(len(rot)):
        for v in rot[u]:
            if (u, v) in seen: continue
            cyc = []; a, b = u, v
            while (a, b) not in seen:
                seen.add((a, b)); cyc.append(a)
                nb = rot[b]; w = nb[(pos[b][a] + 1) % len(nb)]
                a, b = b, w
            assert len(cyc) == 3, 'not a triangulation'
            faces.append(tuple(sorted(cyc)))
    return faces


def check_closed(V, faces):
    """Condition (1): distinct faces, each pair in exactly two faces, face graph connected."""
    assert len(set(faces)) == len(faces)
    inc = {}
    for j, f in enumerate(faces):
        for e in itertools.combinations(f, 2):
            inc.setdefault(e, []).append(j)
    if any(len(x) != 2 for x in inc.values()): return False
    adj = {j: [] for j in range(len(faces))}
    for a, b in inc.values(): adj[a].append(b); adj[b].append(a)
    st, seen = [0], {0}
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == len(faces) and len(set(itertools.chain(*faces))) == V


def refine(V, faces):
    """Colour refinement on the incidence graph, parts coloured differently at start.
    Colours are canonical (ranks of sorted signatures).  Returns final colour list."""
    N = V + len(faces)
    nb = [[] for _ in range(N)]
    for j, f in enumerate(faces):
        for i in f: nb[i].append(V + j); nb[V + j].append(i)
    col = [0] * V + [1] * len(faces)
    ncol = 2
    while True:
        sig = [(col[x],) + tuple(sorted(col[y] for y in nb[x])) for x in range(N)]
        rank = {s: r for r, s in enumerate(sorted(set(sig)))}
        new = [rank[s] for s in sig]
        if len(rank) == ncol:
            return new
        col, ncol = new, len(rank)


def discrete(V, faces):
    c = refine(V, faces)
    return len(set(c)) == len(c)


def canonical(V, faces):
    """Canonical form for a complex with discrete refinement (colours are invariant)."""
    c = refine(V, faces)
    lab = {i: c[i] for i in range(V)}
    return tuple(sorted(tuple(sorted(lab[i] for i in f)) for f in faces))


def edge_index(faces):
    E = sorted({e for f in faces for e in itertools.combinations(f, 2)})
    return E, {e: k for k, e in enumerate(E)}


def ham_cycles(V, E, eidx):
    """All Hamiltonian cycles of the point graph, each once, as edge bitmasks."""
    adj = [[] for _ in range(V)]
    for a, b in E: adj[a].append(b); adj[b].append(a)
    out = []
    path = [0]
    def dfs(u, visited, mask):
        if len(path) == V:
            if 0 in adj[u] and path[1] < path[-1]:          # orientation dedupe
                out.append(mask | (1 << eidx[(min(u, 0), max(u, 0))]))
            return
        for w in adj[u]:
            if not visited >> w & 1:
                path.append(w)
                dfs(w, visited | (1 << w), mask | (1 << eidx[(min(u, w), max(u, w))]))
                path.pop()
    dfs(0, 1, 0)
    return out


def brute_splits(V, faces, E, eidx):
    """Admissible splits by brute force over face subsets (independent of Lemma A.1)."""
    F = len(faces); inc = {}
    for j, f in enumerate(faces):
        for e in itertools.combinations(f, 2): inc.setdefault(e, []).append(j)
    def has_order(S):
        S = list(S)
        for root in S[:1]:
            have = set(faces[root]); used = {root}
            progress = True
            while progress:
                progress = False
                for j in S:
                    if j in used: continue
                    f = faces[j]
                    if any(set(e) <= have and any(k in used for k in inc[e]) for e in itertools.combinations(f, 2)) \
                            and len(set(f) - have) == 1:
                        used.add(j); have |= set(f); progress = True
            return len(used) == len(S) and len(have) == V
    cuts = set()
    for comb in itertools.combinations(range(1, F), V - 3):
        S = (0,) + comb
        if len(S) != V - 2: continue
        T = tuple(j for j in range(F) if j not in S)
        if has_order(S) and has_order(T):
            Sset = set(S)
            cuts.add(sum(1 << eidx[e] for e, (a, b) in inc.items() if (a in Sset) != (b in Sset)))
    return cuts


def integral_cover(cuts, nE):
    """Exists C1,C2,C3 (repetition allowed) covering all pairs (=> some pair doubly cut,
    since 3V > |E| = 3V-6)."""
    full = (1 << nE) - 1
    by_edge = [[c for c in cuts if c >> e & 1] for e in range(nE)]
    for i, c1 in enumerate(cuts):
        for c2 in cuts[i:]:
            miss = full & ~(c1 | c2)
            if not miss: return True
            e = (miss & -miss).bit_length() - 1
            for c3 in by_edge[e]:
                if c3 & miss == miss: return True
    return False


def lp_certificate(cuts, nE):
    """Return ('cert', y, tmax_exact) with exact y, sum y=1, max_C y(C) < 1/3, or
    ('feasible', t_float, w) if the LP optimum is >= 1/3 (up to tolerance)."""
    H = len(cuts)
    A = np.array([[(c >> e) & 1 for c in cuts] for e in range(nE)], dtype=float)  # nE x H
    cost = np.zeros(H + 1); cost[-1] = -1.0
    A_ub = np.hstack([-A, np.ones((nE, 1))]); b_ub = np.zeros(nE)
    A_eq = np.hstack([np.ones((1, H)), np.zeros((1, 1))]); b_eq = [1.0]
    r = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                bounds=[(0, None)] * H + [(None, None)], method='highs-ds')
    assert r.status == 0
    t = -r.fun
    if t >= 1 / 3 - 1e-9:
        return ('feasible', t, r.x[:H])
    yf = -r.ineqlin.marginals
    for D in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 6, 10 ** 9):
        y = [Fr(v).limit_denominator(D) if v > 1e-13 else Fr(0) for v in yf]
        s = sum(y)
        if s <= 0: continue
        y = [v / s for v in y]
        m = max(sum((y[e] for e in range(nE) if c >> e & 1), Fr(0)) for c in cuts)
        if m < Fr(1, 3):
            return ('cert', y, m, t)
    raise RuntimeError('could not rationalise a certificate')


def verify_cert(y, cuts, nE):
    """Independent exact re-check of a dual certificate."""
    assert all(isinstance(v, Fr) and v >= 0 for v in y) and sum(y) == 1
    return max(sum((y[e] for e in range(nE) if c >> e & 1), Fr(0)) for c in cuts) < Fr(1, 3)


def analyse(V, faces, brute=False):
    E, eidx = edge_index(faces)
    nE = len(E); assert nE == 3 * V - 6
    cuts = ham_cycles(V, E, eidx)
    if brute:
        assert set(cuts) == brute_splits(V, faces, E, eidx), 'Lemma 1 cross-check failed'
    deg = [0] * V
    for a, b in E: deg[a] += 1; deg[b] += 1
    return E, eidx, cuts, max(deg)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--vmax', type=int, default=12)
    ap.add_argument('--no13', action='store_true'); args = ap.parse_args()
    src = open(os.path.join(os.path.dirname(PLANTRI), 'plantri.c')).read()
    print('plantri version:', src.split('#define VERSION', 1)[1].split('\n', 1)[0].strip())
    print('V | all (OEIS A000109) | P4 discrete | discrete & maxdeg<=6 | discrete: frac-feasible / exact certs '
          '| max t* (discrete) | non-discrete: frac-feasible / integral-feasible')
    tot_disc = tot_cert = 0
    for V in range(4, args.vmax + 1):
        n = nd = nd6 = nfeas = ncert = nf_frac = nf_int = 0; tmax = 0.0; tmax_cert = Fr(0)
        for rot in planar_code(V):
            faces = faces_by_tracing(rot); n += 1
            assert check_closed(V, faces)
            disc = discrete(V, faces)
            E, eidx, cuts, md = analyse(V, faces, brute=(V <= 10))
            nE = len(E)
            if disc:
                nd += 1; nd6 += md <= 6
                res = lp_certificate(cuts, nE) if cuts else ('cert', [Fr(1, nE)] * nE, Fr(0), 0.0)
                if res[0] == 'feasible':
                    nfeas += 1; tmax = max(tmax, res[1]); print('  FEASIBLE (discrete) V=%d faces=%s' % (V, faces))
                else:
                    assert verify_cert(res[1], cuts, nE)
                    assert not integral_cover(cuts, nE)
                    ncert += 1; tmax = max(tmax, res[3]); tmax_cert = max(tmax_cert, res[2])
            else:
                if cuts and lp_certificate(cuts, nE)[0] == 'feasible': nf_frac += 1
                if cuts and integral_cover(cuts, nE): nf_int += 1
        assert n == OEIS_A000109[V]
        tot_disc += nd; tot_cert += ncert
        print('%2d | %5d | %5d | %3d | %d / %d | %.6f (worst certified bound %s) | %d / %d'
              % (V, n, nd, nd6, nfeas, ncert, tmax, tmax_cert, nf_frac, nf_int), flush=True)
    print('TOTAL V=4..%d: discrete=%d, exact rational dual certificates verified=%d' % (args.vmax, tot_disc, tot_cert))

    # positive control: the source's 13-point complex
    pf = [tuple(sorted(f)) for f in PAPER_FACES]
    assert check_closed(13, pf) and discrete(13, pf)
    E, eidx = edge_index(pf)
    inc = {}
    for j, f in enumerate(pf):
        for e in itertools.combinations(f, 2): inc.setdefault(e, []).append(j)
    cuts = ham_cycles(13, E, eidx)
    bitcuts = []
    for k in range(3):
        c = sum(1 << eidx[e] for e, (a, b) in inc.items() if PAPER_BITS[a][k] != PAPER_BITS[b][k])
        assert c in cuts, 'paper bit split is not admissible'
        bitcuts.append(c)
    sep = [sum(c >> e & 1 for c in bitcuts) for e in range(len(E))]
    assert min(sep) >= 1 and max(sep) >= 2
    res = lp_certificate(cuts, len(E))
    print('POSITIVE CONTROL paper 13-point complex: (1) ok, (4) discrete, %d admissible splits, '
          'paper bits admissible, separation min/max = %d/%d, LP t* = %.6f (feasible: %s)'
          % (len(cuts), min(sep), max(sep), res[1] if res[0] == 'feasible' else res[3], res[0] == 'feasible'))

    if args.no13: return
    # V = 13: classification (discrete & maxdeg <= 6 & fractionally feasible), compare with paper
    pcan = canonical(13, pf); n = nd = nd6 = 0; feas = []
    for rot in planar_code(13):
        faces = faces_by_tracing(rot); n += 1
        if not discrete(13, faces): continue
        nd += 1
        E, eidx, cuts, md = analyse(13, faces)
        if md > 6: continue
        nd6 += 1
        res = lp_certificate(cuts, len(E)) if cuts else ('cert',)
        if res[0] == 'feasible':
            feas.append((faces, integral_cover(cuts, len(E)), canonical(13, faces) == pcan, res[1]))
        else:
            assert verify_cert(res[1], cuts, len(E))
    assert n == OEIS_A000109[13]
    print('V=13: all=%d discrete=%d discrete&maxdeg<=6=%d fractionally feasible=%d' % (n, nd, nd6, len(feas)))
    for faces, integ, ispaper, t in feas:
        print('  feasible class: integral (3 bits) feasible=%s, isomorphic to paper complex=%s, t*=%.6f, faces=%s'
              % (integ, ispaper, t, faces))


if __name__ == '__main__':
    main()
