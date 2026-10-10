#!/usr/bin/env python3
"""Independent re-implementation: exact-integer branch and bound.

Written independently of cert_interval.py (shares only the input file P9.txt).
Differences from the primary certificate:
  * arithmetic: exact integers / rationals only (no floating point at all
    after reading the hex inputs, which are converted exactly);
  * bound: exact second-order Taylor expansion of the real quadratic
        H_i(v) = den(T)*|<x,p_i>|^2 - num(T)*|p_i|^2*|x|^2 ,
    x = chart vector with x_j = 1 and free coordinates v in R^4,
    lower bound on a box  c +- r :
        H(c) - sum_k |dH/dv_k(c)| r_k + sum_k min(Q_kk,0) r_k^2
             - sum_{k != l} |Q_kl| r_k r_l          (exact, since H is quadratic);
  * subdivision: bisection of one coordinate at a time (the widest; ties ->
    lowest index), two children, depth-first with an explicit stack.

Claim: for all nonzero x in C^3, max_i H_i(x) > 0, i.e.
    max_i |<x,p_i>|^2 / (|x|^2 |p_i|^2) > T.
Domain argument: every point of CP^2 has a representative with some x_j = 1
and |x_k| <= 1, so (Re x_k, Im x_k)_{k != j} lies in [-1,1]^4.

Usage: python3 cert_exact.py P9.txt [T_num T_den]
Environment: plain Python 3 (no third-party modules).  Single thread.
"""
import sys
from fractions import Fraction


def read_points(path):
    pts = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        v = [Fraction(float.fromhex(t)) for t in line.split()]
        pts.append([(v[0], v[1]), (v[2], v[3]), (v[4], v[5])])
    assert len(pts) == 9
    return pts


def chart_forms(pts, j, tn, td):
    """Integer quadratic forms H_i(v) = c0 + g.v + v^T Q v (common positive scale)."""
    D = 1
    for p in pts:
        for a, b in p:
            D = max(D, a.denominator, b.denominator)
    free = [k for k in range(3) if k != j]
    forms = []
    for p in pts:
        P = [(int(a * D), int(b * D)) for a, b in p]          # exact: D is a power of 2
        for (a, b), (A, B) in zip(p, P):
            assert Fraction(A, D) == a and Fraction(B, D) == b
        n2 = sum(A * A + B * B for A, B in P)
        # <x,p> = sum_k x_k conj(p_k); x_j = 1, x_k = v_{2t} + i v_{2t+1}
        # Re = Re p_j + sum (pr*a + pi*b), Im = -pi_j + sum (pr*b - pi*a)
        alpha = (P[j][0], -P[j][1])
        Lre = [0] * 4
        Lim = [0] * 4
        for t, k in enumerate(free):
            pr, pi = P[k]
            Lre[2 * t], Lre[2 * t + 1] = pr, pi
            Lim[2 * t], Lim[2 * t + 1] = -pi, pr
        c0 = td * (alpha[0] ** 2 + alpha[1] ** 2) - tn * n2
        g = [td * 2 * (alpha[0] * Lre[k] + alpha[1] * Lim[k]) for k in range(4)]
        Q = [[td * (Lre[k] * Lre[l] + Lim[k] * Lim[l]) - (tn * n2 if k == l else 0)
              for l in range(4)] for k in range(4)]
        forms.append((c0, g, Q))
    return forms


def lower_bound(form, C, Rk, S):
    """2^{2S} * (lower bound of H on box centre C/2^S, half-widths Rk/2^S)."""
    c0, g, Q = form
    s = 1 << S
    QC = [sum(Q[k][l] * C[l] for l in range(4)) for k in range(4)]
    val = c0 * s * s + s * sum(g[k] * C[k] for k in range(4)) + sum(C[k] * QC[k] for k in range(4))
    grad = [s * g[k] + 2 * QC[k] for k in range(4)]            # = 2^S * dH/dv_k(c)
    lb = val - sum(abs(grad[k]) * Rk[k] for k in range(4))
    for k in range(4):
        if Q[k][k] < 0:
            lb += Q[k][k] * Rk[k] * Rk[k]
        for l in range(4):
            if l != k:
                lb -= abs(Q[k][l]) * Rk[k] * Rk[l]
    return lb


def value_at(form, C, S):
    c0, g, Q = form
    s = 1 << S
    return c0 * s * s + s * sum(g[k] * C[k] for k in range(4)) + sum(
        C[k] * Q[k][l] * C[l] for k in range(4) for l in range(4))


def main():
    pts = read_points(sys.argv[1])
    tn, td = (21, 40) if len(sys.argv) < 4 else (int(sys.argv[2]), int(sys.argv[3]))
    assert Fraction(tn, td) > Fraction(1, 2)
    total_leaves = 0
    total_nodes = 0
    max_depth = 0
    for j in range(3):
        forms = chart_forms(pts, j, tn, td)
        # box: centre C/2^S, half-widths Rk/2^S ; start [-1,1]^4 : C=0, R=1, S=0
        stack = [((0, 0, 0, 0), (1, 1, 1, 1), 0, 0)]
        leaves = 0
        while stack:
            C, Rk, S, d = stack.pop()
            total_nodes += 1
            max_depth = max(max_depth, d)
            if any(lower_bound(f, C, Rk, S) > 0 for f in forms):
                leaves += 1
                continue
            if all(value_at(f, C, S) <= 0 for f in forms):
                print("FAIL: centre point violates the claim:", j, [Fraction(c, 1 << S) for c in C])
                sys.exit(1)
            if d > 200:
                print("NOT CERTIFIED: depth limit reached (claim fails or is too tight at this T)")
                sys.exit(2)
            k = max(range(4), key=lambda q: (Rk[q], -q))
            # rescale so that the half-width along k is even, then halve it
            if Rk[k] % 2:
                C = tuple(2 * c for c in C); Rk = tuple(2 * r for r in Rk); S += 1
            h = Rk[k] // 2
            Rn = tuple(h if q == k else Rk[q] for q in range(4))
            for sgn in (-1, 1):
                Cn = tuple(C[q] + sgn * h if q == k else C[q] for q in range(4))
                stack.append((Cn, Rn, S, d + 1))
        total_leaves += leaves
        print(f"chart {j}: accepted leaf boxes {leaves}", flush=True)
    print(f"CERTIFIED (exact integers): max_i |<x,p_i>|^2 > {tn}/{td} |p_i|^2 |x|^2 on C^3\\0; "
          f"leaves {total_leaves}, nodes {total_nodes}, max depth {max_depth}")


if __name__ == "__main__":
    main()
