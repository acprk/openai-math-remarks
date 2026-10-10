#!/usr/bin/env python3
"""Independent exhaustive check of perfect binary arrays on C2 x C2 x P, |P| = 9.

Written for this remark (no code shared with the upstream artifacts).
Requires numpy; single thread.

A function f : C2 x C2 x P -> {+1,-1} is written h = sum_{a,b} X^a Y^b H_ab with
H_ab in {+-1}^P.  For alpha, beta in {+-1} put x_ab = h(X=alpha, Y=beta) / 2.
Perfectness (h h^* = 36 in Z[C2 x C2 x P]) is equivalent to x x^* = 9 in Z[P]
for each of the four x (the characters of C2 x C2 separate the four
components).  Coordinatewise, the 4-vector of x-values at g is (1/2) M s with
M the 4x4 Sylvester-Hadamard matrix and s = (H_00, H_10, H_01, H_11)(g); so every
coefficient of every x lies in {-2,...,2}.

Step 1: enumerate S = {x in {-2..2}^P : x x^* = 9}.
Step 2: count 4-tuples (x_{++}, x_{+-}, x_{-+}, x_{--}) in S^4 such that at every
        g the 4-vector equals (1/2) M s for some s in {+-1}^4.
Each such tuple corresponds to exactly one perfect array h.
P = C9 must give 0 (Theorem T1 with u = 3); P = C3 x C3 must give > 0
(Menon difference sets with parameters (36,15,6) exist there).
"""
import itertools
import sys
import numpy as np

M = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
# rows: (alpha,beta) = (+,+), (-,+), (+,-), (-,-)  acting on s = (H00,H10,H01,H11)
ALLOWED = {tuple((M @ np.array(s)) // 2) for s in itertools.product((1, -1), repeat=4)}


def group(shape):
    elems = list(itertools.product(*[range(k) for k in shape]))
    idx = {e: i for i, e in enumerate(elems)}
    n = len(elems)
    sub = np.zeros((n, n), dtype=np.int64)       # sub[g, s] = index of g - s
    for i, g in enumerate(elems):
        for j, s in enumerate(elems):
            sub[i, j] = idx[tuple((a - b) % k for a, b, k in zip(g, s, shape))]
    return n, sub


def norm_solutions(n, sub):
    vals = np.array(list(itertools.product(range(-2, 3), repeat=n)), dtype=np.int64)
    vals = vals[(vals * vals).sum(1) == n]       # coefficient at identity of x x^*
    keep = np.ones(len(vals), dtype=bool)
    for s in range(1, n):                         # (x x^*)(s) = sum_g x(g) x(g - s)
        keep &= (vals * vals[:, sub[:, s]]).sum(1) == 0
    return [tuple(v) for v in vals[keep]]


def count_arrays(S):
    n = len(S[0])
    total = 0
    first = None
    for a in S:
        for b in S:
            # per coordinate: allowed (c, d) given (a, b)
            opts = []
            for g in range(n):
                o = {(v[2], v[3]) for v in ALLOWED if v[0] == a[g] and v[1] == b[g]}
                if not o:
                    break
                opts.append(o)
            else:
                for c in S:
                    if not all(any(cc == c[g] for (cc, _) in opts[g]) for g in range(n)):
                        continue
                    for d in S:
                        if all((c[g], d[g]) in opts[g] for g in range(n)):
                            total += 1
                            if first is None:
                                first = (a, b, c, d)
    return total, first


def direct_check(first, shape):
    """Rebuild h on C2 x C2 x P from one tuple and test perfectness directly."""
    x = np.array(first)                          # rows ordered as the rows of M
    H = (M.T @ x) // 2                           # s = (1/2) M^T x, rows H00,H10,H01,H11
    assert set(np.unique(H)) <= {-1, 1}
    full = (2, 2) + shape
    f = {}
    for (k, (a1, b1)) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1))):
        for gi, g in enumerate(itertools.product(*[range(m) for m in shape])):
            f[(a1, b1) + g] = int(H[k, gi])
    keys = list(f)
    for s in keys:
        if all(v == 0 for v in s):
            continue
        ac = sum(f[g] * f[tuple((gg + ss) % m for gg, ss, m in zip(g, s, full))] for g in keys)
        if ac != 0:
            return False
    return True


def main():
    for name, shape in (("C9", (9,)), ("C3xC3", (3, 3))):
        n, sub = group(shape)
        S = norm_solutions(n, sub)
        t, first = count_arrays(S)
        print(f"P={name}: |{{x in {{-2..2}}^P : xx*=9}}| = {len(S)}; "
              f"perfect binary arrays on C2xC2x{name} = {t}")
        if first is not None:
            print(f"  direct autocorrelation check of the first array found: "
                  f"{'perfect' if direct_check(first, shape) else 'NOT perfect'}")


if __name__ == "__main__":
    sys.exit(main())
