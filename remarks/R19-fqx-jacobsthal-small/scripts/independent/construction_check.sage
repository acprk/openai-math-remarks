# Independent check in SageMath (separate code base from both Python scripts).
# (1) Explicit construction for K_q(2) <= 2q-1, q in {2,3,4,5,7,8,9}:
#     f = x^q - x, M = prod_{a in F_q} (x - a) * prod_{c != 0} P_c,  P_c = an irreducible factor of
#     x^q - x + c (all of degree >= 2).  Check: the P_c are pairwise distinct, omega(M) = 2q-1, and
#     every g in f + V_2 (q^2 polynomials) has gcd(g, M) != 1.
# (2) Equality structure for q in {3,4,5}: with all q linear primes used, residue vectors
#     v = (v_a)_a (prime x-a covers {c0 + c1 a = v_a}) give total cost 2q-1 exactly when the q lines
#     are concurrent, i.e. a -> v_a is affine.
# Run: sage construction_check.sage     (Sage 9.5+; single thread; ~1 min)
import itertools
for q in [2, 3, 4, 5, 7, 8, 9]:
    F = GF(q, 'w'); R = PolynomialRing(F, 'x'); x = R.gen()
    f = x^q - x
    Ps = []
    for c in F:
        if c == 0: continue
        fac = (x^q - x + c).factor()
        assert all(P.degree() >= 2 for P, e in fac)
        Ps.append(fac[0][0])
    assert len(set(Ps)) == q - 1
    M = prod(x - a for a in F) * prod(Ps)
    om = len(M.factor())
    ok = all(gcd(f + c0 + c1 * x, M) != 1 for c0 in F for c1 in F)
    print("q=%d: construction omega(M)=%d (2q-1=%d), all %d elements of f+V_2 non-coprime: %s"
          % (q, om, 2 * q - 1, q * q, ok))
for q in [3, 4, 5]:
    F = GF(q, 'w'); els = list(F)
    pts = [(c0, c1) for c0 in els for c1 in els]
    affine = set(tuple(al + be * a for a in els) for al in els for be in els)
    n_opt = 0; n_bad = 0
    for v in itertools.product(els, repeat=q):
        cov = sum(1 for (c0, c1) in pts if any(c0 + a * c1 == v[i] for i, a in enumerate(els)))
        cost = q + (q * q - cov)
        assert cost >= 2 * q - 1
        if cost == 2 * q - 1:
            n_opt += 1
            if v not in affine: n_bad += 1
    print("q=%d: residue vectors attaining 2q-1 with all q lines: %d, of which non-affine: %d (affine total %d)"
          % (q, n_opt, n_bad, len(affine)))
