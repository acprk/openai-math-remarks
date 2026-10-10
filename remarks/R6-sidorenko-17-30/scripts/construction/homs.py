"""Enumerate homomorphisms f=(a,b): H(src)->H(tgt) (part preserving: a(i) in b(j)
for every corner), optionally restricted by allowed[vertex] sets
(vertex index: points 0..V-1, faces V..V+F-1)."""


def homs(src, tgt, allowed=None):
    V, F = src.V, src.F
    if allowed is None:
        allowed = [set(range(tgt.V))] * V + [set(range(tgt.F))] * F
    # face order: BFS by shared points
    order = [0]
    seen = {0}
    while len(order) < F:
        best = None
        for j in range(F):
            if j in seen:
                continue
            sc = len(set(src.faces[j]) & set(x for k in order for x in src.faces[k]))
            if best is None or sc > best[0]:
                best = (sc, j)
        order.append(best[1]); seen.add(best[1])
    a = [None] * V
    b = [None] * F
    tf = [set(f) for f in tgt.faces]

    def rec(k):
        if k == F:
            yield (tuple(a), tuple(b))
            return
        j = order[k]
        for bj in allowed[V + j]:
            fs = tf[bj]
            pts = src.faces[j]
            # points already assigned must land in fs
            if any(a[x] is not None and a[x] not in fs for x in pts):
                continue
            b[j] = bj
            free = [x for x in pts if a[x] is None]
            def rp(t):
                if t == len(free):
                    yield from rec(k + 1)
                    return
                x = free[t]
                for y in fs:
                    if y in allowed[x]:
                        a[x] = y
                        yield from rp(t + 1)
                        a[x] = None
            yield from rp(0)
            b[j] = None
    yield from rec(0)
