#!/usr/bin/env python3
"""Independent re-implementation (does not import the primary code, uses no SAT/ILP solver).

Computes K_q(m) exactly by exhaustive branch-and-bound, with its own finite-field arithmetic
(GF(p^k) for q in {2,3,4,5,7}), its own polynomial reduction and its own irreducibility test.

Model. Fix m >= 1.  A coset f + V_m (V_m = polynomials of degree < m) is identified with V_m via
g -> f + g.  For a monic irreducible P:
  * deg P >= m: P divides at most one element of f + V_m (a difference of two elements is a
    nonzero polynomial of degree < m).  Conversely any single prescribed element can be made
    divisible by a fresh irreducible of degree >= m by the Chinese remainder theorem on f.
    Such a prime is a "wildcard" covering one point.
  * deg P < m: P divides f+g iff g == -f (mod P); as f varies, -f mod P is an arbitrary residue r,
    chosen independently for different P (CRT).  So P covers the class {g in V_m : g == r mod P}.
K_q(m) = min over choices of (number of low-degree primes used) + (number of uncovered points).

Exhaustive search over (skip | residue r) for every low-degree prime, with points stored as bit
masks; pruning uses the trivial bound cost >= current cost.  Prints one line per (q, m).
Environment: python3 >= 3.8, standard library only.  Runtime ~1 min.
"""
import itertools, sys

# ---------- GF(q) ----------
MODULI = {2: (2, [0, 1]), 3: (3, [0, 1]), 5: (5, [0, 1]), 7: (7, [0, 1]),
          4: (2, [1, 1, 1])}     # GF(4) = F_2[t]/(t^2+t+1); for prime q the "modulus" is t

class GF:
    def __init__(self, q):
        self.q = q
        p, mod = MODULI[q]
        self.p = p
        k = len(mod) - 1
        self.k = k
        elems = list(itertools.product(range(p), repeat=k))   # coefficient tuples, low first
        self.idx = {e: i for i, e in enumerate(elems)}
        self.el = elems
        def mul(a, b):
            prod = [0] * (2 * k - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    prod[i + j] = (prod[i + j] + x * y) % p
            for d in range(len(prod) - 1, k - 1, -1):      # reduce by monic modulus
                c = prod[d]
                if c:
                    for t in range(k + 1):
                        prod[d - k + t] = (prod[d - k + t] - c * mod[t]) % p
            return tuple(prod[:k])
        self.add = [[self.idx[tuple((x + y) % p for x, y in zip(self.el[a], self.el[b]))]
                     for b in range(q)] for a in range(q)]
        self.mul = [[self.idx[mul(self.el[a], self.el[b])] for b in range(q)] for a in range(q)]
        self.neg = [self.idx[tuple((-x) % p for x in self.el[a])] for a in range(q)]
        self.one = self.idx[tuple([1] + [0] * (k - 1))]
        assert all(any(self.mul[a][b] == self.one for b in range(q)) for a in range(1, q)), "not a field"

def polyrem(F, a, P):
    """remainder of a by monic P (coefficient lists, low first) over F"""
    a = list(a)
    dP = len(P) - 1
    for d in range(len(a) - 1, dP - 1, -1):
        c = a[d]
        if c:
            for t in range(dP + 1):
                a[d - dP + t] = F.add[a[d - dP + t]][F.neg[F.mul[c][P[t]]]]
    a = a[:dP] + [0] * max(0, dP - len(a))
    return tuple(a)

def monic_irreducibles(F, d):
    q = F.q
    lower = [P for e in range(1, d // 2 + 1) for P in monic_irreducibles(F, e)]
    out = []
    for low in itertools.product(range(q), repeat=d):
        P = list(low) + [F.one]
        if d > 1 and P[0] == 0:
            continue
        if all(any(polyrem(F, P, Q)) for Q in lower):
            out.append(P)
    return out

def K(q, m):
    F = GF(q)
    pts = list(itertools.product(range(q), repeat=m))        # g as coefficient tuple, low first
    full = (1 << len(pts)) - 1
    primes = [P for d in range(1, m) for P in monic_irreducibles(F, d)]
    options = []                                              # per prime: list of class masks
    for P in primes:
        classes = {}
        for i, g in enumerate(pts):
            r = polyrem(F, g, P)
            classes[r] = classes.get(r, 0) | (1 << i)
        options.append(sorted(classes.values()))
    best = [len(pts)]                                         # all wildcards
    n = len(primes)
    def dfs(i, cov, used):
        unc = len(pts) - bin(cov).count("1")
        if used + unc < best[0]:
            best[0] = used + unc
        if i == n or used + 1 >= best[0]:
            return
        dfs(i + 1, cov, used)                                 # skip prime i
        for msk in options[i]:
            if msk & ~cov:
                dfs(i + 1, cov | msk, used + 1)
    dfs(0, 0, 0)
    return best[0], len(primes)

if __name__ == "__main__":
    cases = [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3),
             (4, 1), (4, 2), (5, 1), (5, 2), (7, 1), (7, 2)]
    for q, m in cases:
        k, npr = K(q, m)
        extra = ""
        if m == 1: extra = f"   [claim K_q(1)=q={q}: {'OK' if k == q else 'FAIL'}]"
        if m == 2: extra = f"   [claim K_q(2)=2q-1={2*q-1}: {'OK' if k == 2*q-1 else 'FAIL'}]"
        print(f"q={q} m={m}: K_q(m)={k}  (low-degree primes available: {npr}){extra}", flush=True)
