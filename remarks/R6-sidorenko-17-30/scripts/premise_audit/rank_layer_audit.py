"""Independent recomputation of every complex-dependent quantity used in part (R)
(transverse.tex, orbits.tex, planted.tex, tail.tex) of openai/math 161,
for the paper's 13-point complex (calibration) and for candidate 3 (V=7).
Usage: python3 rank_layer_audit.py
No dependency on check_p123.py / cplx.py (written from scratch)."""
import itertools
from collections import Counter, deque
from math import comb, sqrt

PAPER13 = "0,1,3 0,1,9 0,2,3 0,2,9 1,3,10 1,7,10 1,7,12 1,9,12 2,3,4 2,4,9 3,4,11 3,10,11 4,8,9 4,8,11 5,6,7 5,6,11 5,7,10 5,10,11 6,7,12 6,8,11 6,8,12 8,9,12"
PAPER13_BITS = "000 010 001 110 001 000 100 110 101 100 111 101 101 011 000 110 010 111 001 010 011 111".split()
CAND3 = "0,1,2 0,1,6 0,2,3 0,3,4 0,4,5 0,5,6 1,2,5 1,5,6 2,3,4 2,4,5"


def faces_of(s):
    return [tuple(sorted(map(int, t.split(',')))) for t in s.split()]


def structure(F):
    V = 1 + max(max(f) for f in F)
    pf = {}
    for j, f in enumerate(F):
        for e in itertools.combinations(f, 2):
            pf.setdefault(e, []).append(j)
    adj = {j: set() for j in range(len(F))}
    for e, js in pf.items():
        if len(js) == 2:
            a, b = js; adj[a].add(b); adj[b].add(a)
    return V, pf, adj


def bfs(adj, s, allowed=None):
    allowed = set(adj) if allowed is None else set(allowed)
    dist = {s: 0}; q = deque([s])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y in allowed and y not in dist:
                dist[y] = dist[x] + 1; q.append(y)
    return dist


def exposure_order(F, adj, cls):
    """Greedy order: each later face shares a pair with an earlier face of the class
    and adds exactly one new point. Returns (order,newpoints) or None."""
    cls = sorted(cls)
    for start in cls:
        order = [start]; pts = set(F[start]); new = []
        rest = set(cls) - {start}
        ok = True
        while rest:
            cand = [j for j in sorted(rest) if adj[j] & set(order) and len(set(F[j]) - pts) == 1]
            if not cand:
                ok = False; break
            j = cand[0]; new.append((set(F[j]) - pts).pop()); pts |= set(F[j])
            order.append(j); rest.remove(j)
        if ok:
            return order, new, pts
    return None


def check_bits(F, bits, V, pf, adj, verbose=True):
    nb = len(bits[0])
    ok = True
    for pos in range(nb):
        for val in '01':
            cls = [j for j in range(len(F)) if bits[j][pos] == val]
            r = exposure_order(F, adj, cls)
            good = r is not None and r[2] == set(range(V)) and len(cls) == V - 2
            ok &= good
            if verbose:
                print('   bit %d=%s: %2d faces, order %s, new %s, %s' % (
                    pos + 1, val, len(cls), r[0] if r else None, r[1] if r else None, 'OK' if good else 'FAIL'))
    a = {e: sum(bits[js[0]][p] != bits[js[1]][p] for p in range(nb)) for e, js in pf.items()}
    p3 = min(a.values()) * 3 >= nb and max(a.values()) * 3 > nb
    return ok, a, p3


def find_bits(F, V, pf, adj, nb=3):
    nF = len(F)
    goods = []
    for cls in itertools.combinations(range(nF), V - 2):
        if 0 not in cls:
            continue
        comp = set(range(nF)) - set(cls)
        r1 = exposure_order(F, adj, cls); r2 = exposure_order(F, adj, comp)
        if r1 and r2 and r1[2] == set(range(V)) and r2[2] == set(range(V)):
            goods.append(frozenset(cls))
    sols = []
    for combo in itertools.combinations(goods, nb):  # distinct splits
        bits = [''.join('0' if j in s else '1' for s in combo) for j in range(nF)]
        a = [sum(bits[js[0]][p] != bits[js[1]][p] for p in range(nb)) for js in pf.values()]
        if 3 * min(a) >= nb and 3 * max(a) > nb:
            sols.append(bits)
    return goods, sols


def report(name, F, bits=None):
    V, pf, adj = structure(F)
    nF, nE = len(F), len(pf)
    print('=' * 70); print(name, ': V=%d F=%d |E|=%d |H|=%d edges=%d' % (V, nF, nE, V + nF, 3 * nF))
    print(' P1 every pair in exactly two faces:', set(len(v) for v in pf.values()) == {2})
    ecc = [max(bfs(adj, s).values()) for s in range(nF)]
    print(' face-neighbor graph connected:', len(bfs(adj, 0)) == nF, ' diameter l* =', max(ecc))
    deg = [sum(1 for e in pf if i in e) for i in range(V)]
    nfaces = [sum(1 for f in F if i in f) for i in range(V)]
    print(' point degrees in (I,E):', deg, ' (= #incident faces:', deg == nfaces, ')')
    print(' max degree s_max =', max(deg))
    # transverse joint-sign threshold: smallest even D with b(N-s)>s for all s<=smax, N=D-1
    b = lambda n: n * (n + 1) // 2
    D = 2
    while not all(b(D - 1 - s) > s for s in range(1, max(deg) + 1)) or D < 4:
        D += 2
    print(' transverse joint-sign: minimal even D with b(D-1-s)>s for s<=smax: D =', D, '(paper uses D>=12)')
    S618 = sum(comb(2 * d, 2) for d in deg)
    m = nE
    const = max(mm * mm - 4 * mm for mm in range(m + 1)) + S618
    print(' sum C(2deg,2) =', S618, '; m<=%d; m^2-4m+sum =' % m, const, '; gain<=2m=%d; nondirect remainder =' % (2 * m), const + 2 * m)
    print(' tail-budget constant: |J|*3/2 = %gM ; residue D-threshold uses D > M* + 2|J| C* = M* + %d C*' % (1.5 * nF, 2 * nF))
    print(' chart comparison exponent beta_D^{-%d}' % (V + nF))
    if bits is None:
        goods, sols = find_bits(F, V, pf, adj)
        print(' admissible splits (both halves exposure-orderable):', len(goods), '; 3-bit P2/P3 certificates:', len(sols))
        bits = sols[0]
    print(' bits:', ' '.join('%s:%s' % (''.join(map(str, F[j])), bits[j]) for j in range(nF)))
    ok, a, p3 = check_bits(F, bits, V, pf, adj)
    print(' P2 (all 6 classes exposure-orderable, V-2 faces, cover all points):', ok)
    print(' P3 (3a_e>=1, some 3a_e>1):', p3, ' 3a_e histogram:', dict(Counter(a.values())))
    print(' pairs with a_e>1/3:', sorted(e for e in a if a[e] > 1))
    print(' pairs with a_e=1/3 (only ones allowed h_e=2 in residual profiles):', sorted(e for e in a if a[e] == 1))
    # check that the residue lemma's propagation step never needs anything but connectivity:
    # every pair lies on some face; faces at distance t reachable from any j0
    return bits


if __name__ == '__main__':
    report('paper 13-point complex', faces_of(PAPER13), PAPER13_BITS)
    F3 = faces_of(CAND3)
    report('candidate 3 (V=7)', F3)
    # robustness: every 3-bit certificate for cand 3 satisfies P2/P3; list all a_e patterns
    V, pf, adj = structure(F3)
    goods, sols = find_bits(F3, V, pf, adj)
    pats = Counter(tuple(sum(s[js[0]][p] != s[js[1]][p] for p in range(3)) for js in pf.values()) for s in sols)
    print('candidate 3: number of 3-bit certificates (unordered distinct splits):', len(sols), '; distinct a-patterns:', len(pats))
    best = max(pats, key=lambda t: sum(x > 1 for x in t))
    print('  max #pairs with a_e>1/3 over certificates:', sum(x > 1 for x in best))
