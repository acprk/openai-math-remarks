#!/usr/bin/env python3
"""Global interval-arithmetic certificate (primary computation).

Claim certified: for the nine points p_1..p_9 in C^3 listed in P9.txt and the
threshold T (default 13/25 = 0.52),

    for every nonzero x in C^3:  max_i |<x,p_i>|^2 - T*|p_i|^2*|x|^2 > 0 .

Domain: CP^2 is the union of the three charts x_j = 1, |x_k| <= 1 (k != j)
(normalize so that the coordinate of largest modulus is real positive and
equal to 1).  Each chart is contained in the real box [-1,1]^4 of
(Re x_k, Im x_k)_{k != j}; we certify the inequality on the whole box,
which covers CP^2.  Branch and bound: a box is accepted when, for some i,
a rigorous interval lower bound of h_i(x) = |<x,p_i>|^2 - T|p_i|^2|x|^2 on the
box is > 0; otherwise it is split into 16 congruent sub-boxes.

Rigour: every floating-point operation is IEEE-754 round-to-nearest and is
followed by an outward step with numpy.nextafter, so each interval provably
contains the exact real range (standard 'poor man's directed rounding').
The data p_i are exact binary64 numbers read from hex strings; T and all box
centres/radii are dyadic, hence exact.  No transcendental functions are used.

Usage: python3 cert_interval.py P9.txt [T_num T_den]
Environment: Python 3, numpy (tested numpy 1.26/2.x); single thread.
Expected final line: 'CERTIFIED ...' with the box counts recorded in
expected_output/cert_interval.txt.
"""
import sys
from fractions import Fraction
import numpy as np

np.seterr(over="raise", invalid="raise", divide="raise", under="ignore")
DN = lambda v: np.nextafter(v, -np.inf)
UP = lambda v: np.nextafter(v, np.inf)


class I:
    """Interval [lo, hi] over numpy arrays with outward rounding."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=np.float64)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=np.float64)

    def __add__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(DN(self.lo + o.lo), UP(self.hi + o.hi))

    def __sub__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(DN(self.lo - o.hi), UP(self.hi - o.lo))

    def __mul__(self, o):
        o = o if isinstance(o, I) else I(o)
        c = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return I(DN(np.minimum.reduce(c)), UP(np.maximum.reduce(c)))

    def sq(self):
        a, b = self.lo * self.lo, self.hi * self.hi
        lo = np.where(self.lo >= 0, a, np.where(self.hi <= 0, b, 0.0))
        lo = np.where(lo > 0, DN(lo), 0.0)  # squares are >= 0
        return I(lo, UP(np.maximum(a, b)))


def read_points(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        v = [float.fromhex(t) for t in line.split()]
        assert len(v) == 6
        rows.append([(v[0], v[1]), (v[2], v[3]), (v[4], v[5])])
    assert len(rows) >= 1
    return rows


def main():
    P = read_points(sys.argv[1])
    print(f"{len(P)} points read from {sys.argv[1]}")
    T = Fraction(13, 25) if len(sys.argv) < 4 else Fraction(int(sys.argv[2]), int(sys.argv[3]))
    # interval enclosure of T and of T*|p_i|^2 (computed exactly in rationals, then enclosed)
    TP = []
    for p in P:
        q = T * sum(Fraction(a) ** 2 + Fraction(b) ** 2 for a, b in p)
        f = float(q)
        TP.append(I(DN(f), UP(f)))
    total_acc = 0
    total_boxes = 0
    for chart in range(3):
        free = [k for k in range(3) if k != chart]
        # 16 initial boxes of half-width 1/2 inside [-1,1]^4
        g = np.array([-0.5, 0.5])
        C = np.array(np.meshgrid(g, g, g, g, indexing="ij")).reshape(4, -1).T.copy()
        R = 0.5
        depth = 0
        acc_chart = 0
        while len(C):
            n = len(C)
            total_boxes += n
            # coordinates as intervals: x_chart = 1, x_k = (a_k + i b_k)
            X = {chart: (I(np.ones(n)), I(np.zeros(n)))}
            for t, k in enumerate(free):
                X[k] = (I(C[:, 2 * t] - R, C[:, 2 * t] + R), I(C[:, 2 * t + 1] - R, C[:, 2 * t + 1] + R))
            # these are exact: centres are multiples of R, R a power of 2
            nx = I(np.ones(n))
            for k in free:
                nx = nx + X[k][0].sq() + X[k][1].sq()
            best = np.full(n, -np.inf)
            for i, p in enumerate(P):
                # <x,p> := sum_k x_k * conj(p_k) ; |<x,p>|^2 is what matters
                re = I(np.zeros(n))
                im = I(np.zeros(n))
                for k in range(3):
                    pr, pi = p[k]
                    xr, xi = X[k]
                    # (xr + i xi)(pr - i pi) = (xr pr + xi pi) + i (xi pr - xr pi)
                    re = re + xr * pr + xi * pi
                    im = im + xi * pr - xr * pi
                h = re.sq() + im.sq() - TP[i] * nx
                best = np.maximum(best, h.lo)
            ok = best > 0
            acc_chart += int(ok.sum())
            # failure detection at box centres (degenerate intervals)
            if not ok.all():
                Cb = C[~ok]
                m = len(Cb)
                Xc = {chart: (I(np.ones(m)), I(np.zeros(m)))}
                for t, k in enumerate(free):
                    Xc[k] = (I(Cb[:, 2 * t]), I(Cb[:, 2 * t + 1]))
                nxc = I(np.ones(m))
                for k in free:
                    nxc = nxc + Xc[k][0].sq() + Xc[k][1].sq()
                hb = np.full(m, -np.inf)
                for i, p in enumerate(P):
                    re = I(np.zeros(m)); im = I(np.zeros(m))
                    for k in range(3):
                        pr, pi = p[k]; xr, xi = Xc[k]
                        re = re + xr * pr + xi * pi
                        im = im + xi * pr - xr * pi
                    hc = re.sq() + im.sq() - TP[i] * nxc
                    hb = np.maximum(hb, hc.hi)
                if (hb <= 0).any():
                    print("FAIL: centre with max_i h_i <= 0 (certificate impossible at this T):",
                          chart, Cb[hb <= 0][0])
                    sys.exit(1)
            C = C[~ok]
            depth += 1
            print(f"chart {chart} depth {depth} R={R} boxes {n} accepted {int(ok.sum())} remaining {len(C)}",
                  flush=True)
            if len(C) == 0:
                break
            if depth > 30 or len(C) > 2e7:
                print("ABORT: too deep / too many boxes")
                sys.exit(2)
            R2 = R / 2
            sg = np.array(np.meshgrid(*[[-1.0, 1.0]] * 4, indexing="ij")).reshape(4, -1).T
            C = (C[None, :, :] + R2 * sg[:, None, :]).reshape(-1, 4)
            R = R2
        total_acc += acc_chart
        print(f"chart {chart}: accepted boxes {acc_chart}", flush=True)
    print(f"CERTIFIED: max_i |<x,p_i>|^2 > {T} |p_i|^2 |x|^2 for all x in C^3\\0 ({len(P)} points); "
          f"accepted leaf boxes {total_acc}; boxes examined {total_boxes}")


if __name__ == "__main__":
    main()
