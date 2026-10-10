"""Step 0: invariants of every complex V=5..9 satisfying prop:complex (1)-(3),
and a literature-based Sidorenko test via reflection (Cauchy-Schwarz gluing):
if an involution tau in Aut fixes no face and the faces split as Phi u tau(Phi)
with points(Phi) & points(tau Phi) inside Fix(tau), then
  t(H,W) = E_{x_S}[g(x_S)^2] >= t(F,W)^2,  F = incidence graph of Phi,
so H is Sidorenko whenever F is.  F is certified Sidorenko by
 (a) a part of size <= 4 [Sidorenko 1991], or
 (b) a vertex adjacent to the whole opposite part [Conlon-Fox-Sudakov 2010].
Run: /usr/bin/python3 survey.py cands.txt
"""
import sys, itertools
from cplx import load, wl_classes, automorphisms


def half_ok(cx, Phi):
    pts = sorted(set(x for j in Phi for x in cx.faces[j]))
    nf = len(Phi)
    # (a) part size <= 4
    if nf <= 4 or len(pts) <= 4:
        return 'part<=4 (faces %d, pts %d)' % (nf, len(pts))
    # (b) dominating vertex
    for p in pts:
        if all(p in cx.faces[j] for j in Phi):
            return 'point %d dominates the %d faces' % (p, nf)
    for j in Phi:
        if set(cx.faces[j]) >= set(pts):
            return 'face dominates'
    return None


def reflections(cx, auts):
    res = []
    for g in auts:
        if all(g[g[x]] == x for x in range(cx.V)) and any(g[x] != x for x in range(cx.V)):
            fix = {x for x in range(cx.V) if g[x] == x}
            fimg = [cx.faces.index(tuple(sorted(g[x] for x in f))) for f in cx.faces]
            if any(fimg[j] == j for j in range(cx.F)):
                continue
            orbits = sorted({tuple(sorted((j, fimg[j]))) for j in range(cx.F)})
            for choice in itertools.product((0, 1), repeat=len(orbits)):
                Phi = [o[c] for o, c in zip(orbits, choice)]
                Psi = [fimg[j] for j in Phi]
                P1 = set(x for j in Phi for x in cx.faces[j])
                P2 = set(x for j in Psi for x in cx.faces[j])
                if P1 & P2 <= fix:
                    res.append((g, sorted(Phi), sorted(P1 & P2)))
    return res


def main(path):
    for n, cx in enumerate(load(path)):
        auts = automorphisms(cx)
        col = wl_classes(cx)
        pc = len(set(col[:cx.V])); fc = len(set(col[cx.V:]))
        deg = cx.degrees()
        refl = reflections(cx, auts)
        verdict = None
        for g, Phi, S in refl:
            h = half_ok(cx, Phi)
            if h:
                verdict = (g, Phi, S, h); break
        print('#%d V=%d F=%d |V(H)|=%d |E(H)|=%d pointdeg=%s |Aut|=%d WLclasses pts=%d faces=%d reflections=%d'
              % (n, cx.V, cx.F, cx.nvert(), len(cx.corners), sorted(deg), len(auts), pc, fc, len(refl)))
        print('   faces:', ' '.join('%d%d%d' % f for f in cx.faces))
        if verdict:
            g, Phi, S, h = verdict
            print('   SIDORENKO (reflection) tau=%s half faces=%s glue S=%s ; half: %s'
                  % (g, [cx.faces[j] for j in Phi], S, h))
        else:
            print('   no certified reflection-halving')


if __name__ == '__main__':
    main(sys.argv[1])
