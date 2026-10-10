# Verification program extracted VERBATIM from proofs/CODEX-ROUND1.md
# Author: OpenAI Codex session (model gpt-6-astra, Codex CLI v0.157.0), 2026-10-09,
# in a Claude<->Codex collaboration round directed by [AUTHOR NAME].
# Environment: /usr/bin/python3 with sympy (round 1), stdlib (round 2), numpy+scipy+mpmath (round 3).
# Expected output: expected_output/codex_round1_verify.txt
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb
from sympy import Matrix, Rational, integer_nthroot, symbols, expand

@lru_cache(None)
def phi(a):
    return Q(1) if a == 1 else phi(a-1) * Q(3*a-2, 3*a-3)

@lru_cache(None)
def F(a, b):
    a, b = sorted((a, b))
    assert a >= 1
    return phi(a) * (b + Q(a-1, 2))

counts = {}
def ge(name, lhs, rhs):
    counts[name] = counts.get(name, 0) + 1
    assert lhs >= rhs, (name, lhs, rhs)

for a, b in product(range(1, 81), repeat=2):
    assert F(a, b) == F(b, a)
    ge('rank upper (cubed)', (a+b-1)**4, F(a, b)**3)
    ge('outer product upper', a*b, F(a, b))
    ge('dot product lower', F(a, b), max(a, b))
    ge('monotonicity', F(a, b+1), F(a, b))
    if b >= 2:
        ge('separate concavity', 2*F(a, b), F(a, b-1)+F(a, b+1))
    ge('tripling', F(a, 3*b+a-1), 3*F(a, b))
    if b >= a-1:
        assert F(a, 3*b+a-1) == 3*F(a, b)
for a in range(1, 81):
    assert F(1, a) == a
    ge('concavity at zero', 2*F(a, 1), F(a, 2))

for a, b in product(range(1, 51), repeat=2):
    for m in range(1, 31):
        transfer = [F(a+m-1-2*i, b-m+1+2*i)
                    for i in range(min(a, m))
                    if 0 <= m-1-i < min(b, m)]
        higher = [F(a, b+m-1-2*i) for i in range(min(b, m))]
        ge('A transfer', m*F(a, b), sum(transfer, Q(0)))
        ge('B higher dot', m*F(a, b), sum(higher, Q(0)))
        # Exact endpoint/mean identity used in the all-size Jensen proof.
        if m < a+b:
            coords = [a+m-1-2*i for i in range(min(a, m))
                      if 0 <= m-1-i < min(b, m)]
            coords += [0]*max(m-a, 0) + [a+b]*max(m-b, 0)
            assert len(coords) == m and sum(coords) == m*a

for a, b, c, d in product(range(1, 10), repeat=4):
    for k in range(min(a, c)+min(b, d)-1):
        terms = [F(a+c-1-2*i, b+d-1-2*(k-i))
                 for i in range(min(a, c)) if 0 <= k-i < min(b, d)]
        ge('C bivariate shared output', F(a, b)*F(c, d), sum(terms, Q(0)))

SCALE = 10**24
@lru_cache(None)
def power_interval(q):
    # Bound q^(3/4) using only integer arithmetic and rational comparisons.
    cube = q**3
    n = integer_nthroot(cube.numerator*SCALE**4 // cube.denominator, 4)[0]
    lo = Q(int(n), SCALE)
    if lo**4 == cube:
        return lo, lo
    hi = Q(int(n)+1, SCALE)
    assert lo**4 < cube < hi**4
    return lo, hi

for a, b, c, d in product(range(1, 11), repeat=4):
    lhs = F(a, b)*F(c, d)
    terms = [F(a+c-1-2*i, b+d-1-2*i) for i in range(min(a, b, c, d))]
    if len(terms) == 1:
        ge('D full direct sum', lhs, terms[0])
    else:
        # Failure here would be reported, not accepted by a tolerance.
        ge('D full direct sum', power_interval(lhs)[0],
           sum((power_interval(x)[1] for x in terms), Q(0)))

for m, a, b in product(range(1, 11), repeat=3):
    ge('F Fourier with cost', (2*m-1)*F(m*a, m*b), m*F(a, b))

# Exact coefficients of adapted D-filtration bases.
# Keys (i,j) denote u^(r-i)v^i s^(s-j)w^j.
def lift(r, s, k, ell):
    i = min(ell, r-k)
    j = ell-i
    return {(i+h, j+k-h): Q((-1)**h * comb(k, h)) for h in range(k+1)}

@lru_cache(None)
def basis(r, s):
    return [(k, ell, lift(r, s, k, ell))
            for k in range(min(r, s)+1)
            for ell in range(r+s-2*k+1)]

apolar_checks = 0
for a, b, m in product(range(1, 5), repeat=3):
    r = m-1
    for i, ell, x in basis(a-1, r):
        for j, n, y in basis(b-1, r):
            if i+j < r:
                continue
            out = {}
            for (v, w), cx in x.items():
                for (vv, ww), cy in y.items():
                    if w+ww == r:
                        out[v+vv] = out.get(v+vv, Q(0)) + cx*cy*Q((-1)**w, comb(r, w))
            out = {k: v for k, v in out.items() if v}
            expected = {} if i+j > r else {ell+n: Q((-1)**i, comb(r, i))}
            assert out == expected, (a, b, m, i, j, out, expected)
            apolar_checks += 1

@lru_cache(None)
def inverse_basis(r, s):
    bs = basis(r, s)
    size = (r+1)*(s+1)
    mat = Matrix.zeros(size, size)
    for col, (_, _, p) in enumerate(bs):
        for (i, j), coeff in p.items():
            mat[i*(s+1)+j, col] = Rational(coeff.numerator, coeff.denominator)
    return mat.inv()

bivariate_checks = 0
for a, b, c, d in product(range(1, 4), repeat=4):
    rr, ss = a+b-2, c+d-2
    out_basis = basis(rr, ss)
    inv = inverse_basis(rr, ss)
    for i, ell, x in basis(a-1, c-1):
        for j, n, y in basis(b-1, d-1):
            raw = Matrix.zeros((rr+1)*(ss+1), 1)
            for (v, w), cx in x.items():
                for (vv, ww), cy in y.items():
                    coeff = cx*cy
                    raw[(v+vv)*(ss+1)+w+ww] += Rational(coeff.numerator, coeff.denominator)
            transformed = inv*raw
            for coeff, (k, deg, _) in zip(transformed, out_basis):
                if k < i+j:
                    assert coeff == 0
                elif k == i+j:
                    assert coeff == int(deg == ell+n)
            bivariate_checks += 1

# The coordinate-weight obstruction is an identity in all leg weights.
coeffs = [0]*15  # 4 X weights, 4 Y weights, 7 Z weights.
for i, j, sign in [(1,1,1), (2,2,1), (2,1,-1), (0,2,-1), (1,3,-1), (0,3,1)]:
    coeffs[i] += sign
    coeffs[4+j] += sign
    coeffs[8+i+j] += sign
assert coeffs == [0]*15

assert F(2,2) == Q(10,3) and F(2,3) == Q(14,3)
assert F(3,3) == Q(56,9) and F(4,4) == Q(770,81)
assert 2*F(2,2) - F(3,1) - F(1,3) == Q(2,3)
assert F(2,2)**2 - F(4,4) == Q(130,81)
assert (F(2,2)**2)**3 > 6**4 and F(3,3)**3 < 4**4
assert 7 < 9  # Output flattenings of the rejected carry-free degeneration.

# Symbolic certificates used in the all-size survival proof.
x = symbols('x')
aa = x+1
upper_numerator = (3*aa*(3*aa-1))**3*(2*aa+1)**4 - ((3*aa+1)*(3*aa+2))**3*(2*aa-1)**4
positive_coefficients = [9496, 81476, 297618, 603680, 743614, 570348, 266544, 69552, 7776]
assert expand(upper_numerator - sum(c*x**i for i, c in enumerate(positive_coefficients))) == 0
aa, hh = symbols('aa hh')
tripling_numerator = (3*aa+1)*(aa+2*hh)*(2*aa+hh-1) - 3*aa*(aa+2*hh-1)*(2*aa+hh+1)
assert expand(tripling_numerator - 2*(aa-hh)*(aa-hh+1)) == 0

assert counts['A transfer'] == counts['B higher dot'] == 75000
assert counts['C bivariate shared output'] == 39609
assert counts['D full direct sum'] == 10000
assert counts['F Fourier with cost'] == 1000
assert sum(counts.values()) == 239009
assert apolar_checks == 780 and bivariate_checks == 1296
for name, count in counts.items():
    print(f'{name}: {count} exact comparisons, PASS')
print(f'apolar leading/zero coefficients: {apolar_checks}, PASS')
print(f'bivariate adapted-basis products: {bivariate_checks}, PASS')
print('All displayed witnesses and the coordinate-weight identity: PASS')
print('No floating-point comparisons; no unresolved root intervals.')
