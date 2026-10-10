#!/usr/bin/env python3
"""Adversarial re-derivation of the 269 objects FROM THE DEFINITIONS (no code shared with cert/ or engine.py).

* one-particle basis e_i = |Q/2, Q/2 - i>  (Condon-Shortley standard basis of spin Q/2)
* pair spin Q-1 states built from Clebsch-Gordan coefficients (Racah closed formula, exact rational squares)
* v_p := sigma * |Q-1, Q-1-p>  with the single global sign sigma fixed by v_0 = -e_0 ^ e_1 (paper 03-spin.tex l.10)
* Fock operators by our own Jordan-Wigner construction; B(e_x ^ e_y) = c_y c_x  (paper 02-operators.tex l.211)
"""
import itertools
from fractions import Fraction
from math import factorial, sqrt
import numpy as np
import scipy.sparse as sp


def cg_sq(j1, m1, j2, m2, J, M):
    """Clebsch-Gordan <j1 m1 j2 m2 | J M> via Racah's formula; all arguments are DOUBLED integers.
    Returns (sign, Fraction square)."""
    if m1 + m2 != M:
        return 0, Fraction(0)
    h = lambda a: a // 2  # all combinations below are even
    for a in (j1 + m1, j2 + m2, J + M, j1 + j2 - J, J + j1 - j2, J - j1 + j2):
        assert a % 2 == 0 and a >= 0
    pre = Fraction((J + 1) * factorial(h(J + j1 - j2)) * factorial(h(J - j1 + j2)) * factorial(h(j1 + j2 - J)),
                   factorial(h(j1 + j2 + J) + 1))
    pre *= factorial(h(J + M)) * factorial(h(J - M)) * factorial(h(j1 - m1)) * factorial(h(j1 + m1)) \
        * factorial(h(j2 - m2)) * factorial(h(j2 + m2))
    s = Fraction(0)
    for k in range(0, h(j1 + j2 - J) + h(j1 + j2 + J) + 2):
        args = [k, h(j1 + j2 - J) - k, h(j1 - m1) - k, h(j2 + m2) - k, h(J - j2 + m1) + k, h(J - j1 - m2) + k]
        if min(args) < 0:
            continue
        d = 1
        for a in args:
            d *= factorial(a)
        s += Fraction((-1) ** k, d)
    if s == 0:
        return 0, Fraction(0)
    return (1 if s > 0 else -1), pre * s * s


def pair_vectors_cg(Q):
    """dict p -> {(x,y): coefficient of e_x ^ e_y (x<y)} for the spin-(Q-1) summand of wedge^2 U_Q, from CG."""
    j = Q          # doubled spin of U_Q
    J = 2 * Q - 2  # doubled pair spin Q-1
    vec = {}
    for p in range(0, 2 * Q - 1):
        M = J - 2 * p
        d = {}
        for x in range(Q + 1):
            for y in range(x + 1, Q + 1):
                mx, my = j - 2 * x, j - 2 * y
                if mx + my != M:
                    continue
                sg, sq = cg_sq(j, mx, j, my, J, M)
                if sg:
                    d[(x, y)] = sg * sqrt(2 * float(sq))   # |JM> = sum_{x<y} CG(x,y)(|xy>-|yx>) = sum sqrt2 CG e_x^e_y
        vec[p] = d
    sigma = -1.0 / vec[0][(0, 1)]       # v_0 = -e_0^e_1
    assert abs(abs(sigma) - 1) < 1e-12
    return {p: {k: sigma * v for k, v in d.items()} for p, d in vec.items()}


class Fock:
    def __init__(self, Q):
        self.Q = Q
        self.nm = Q + 1
        self.bases, self.index, self._c = {}, {}, {}

    def basis(self, N):
        if N not in self.bases:
            L = list(itertools.combinations(range(self.nm), N))
            self.bases[N] = L
            self.index[N] = {s: i for i, s in enumerate(L)}
        return self.bases[N]

    def c(self, j, N):
        """annihilator c_j : wedge^N -> wedge^{N-1};  c_j e_{s1}^...^e_{sN} = (-1)^pos ... (pos = position of j)"""
        key = (j, N)
        if key not in self._c:
            src = self.basis(N); self.basis(N - 1); idx = self.index[N - 1]
            r, cc, v = [], [], []
            for col, s in enumerate(src):
                if j in s:
                    pos = s.index(j)
                    r.append(idx[s[:pos] + s[pos + 1:]]); cc.append(col); v.append((-1.0) ** pos)
            self._c[key] = sp.csr_matrix((v, (r, cc)), shape=(len(self.bases[N - 1]), len(src)))
        return self._c[key]

    def B(self, vec, N):
        """B(v) for v in wedge^2 given as {(x,y): coef}: B(e_x ^ e_y) = c_y c_x (conjugate-linear; real here)."""
        M = None
        for (x, y), a in vec.items():
            t = a * (self.c(y, N - 1) @ self.c(x, N))
            M = t if M is None else M + t
        return M

    def lowering(self, N):
        """J_- on wedge^N: e_i -> sqrt((i+1)(Q-i)) e_{i+1}"""
        Q = self.Q; idx = self.index[N]
        r, cc, v = [], [], []
        for col, s in enumerate(self.basis(N)):
            for pos, i in enumerate(s):
                if i + 1 <= Q and (i + 1) not in s:
                    t = list(s); t[pos] = i + 1
                    r.append(idx[tuple(t)]); cc.append(col); v.append(sqrt((i + 1) * (Q - i)))
        n = len(self.bases[N])
        return sp.csr_matrix((v, (r, cc)), shape=(n, n))
