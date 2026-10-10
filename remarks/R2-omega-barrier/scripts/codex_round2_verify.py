# Verification program extracted VERBATIM from proofs/CODEX-ROUND2.md
# Author: OpenAI Codex session (model gpt-6-astra, Codex CLI v0.157.0), 2026-10-09,
# in a Claude<->Codex collaboration round directed by [AUTHOR NAME].
# Environment: /usr/bin/python3 with sympy (round 1), stdlib (round 2), numpy+scipy+mpmath (round 3).
# Expected output: expected_output/codex_round2_verify.txt
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product

@lru_cache(None)
def phi(a):
    assert a >= 1
    return Q(1) if a == 1 else phi(a-1)*Q(3*a-2, 3*a-3)

@lru_cache(None)
def F(a, b):
    a, b = sorted((a, b))
    return phi(a)*(b+Q(a-1, 2))

counts = dict(union=0, C=0, D_sufficient=0,
              phi_tangent=0, cw_uniform=0, sandwich=0, cw_support=0,
              mixed_square=0)

for a, b, c, d in product(range(1, 17), repeat=4):
    lhs = F(a, b)*F(c, d)
    for k in range(min(a, c)+min(b, d)-1):
        rhs = sum(F(a+c-1-2*i, b+d-1-2*(k-i))
                  for i in range(min(a, c))
                  if 0 <= k-i < min(b, d))
        assert lhs >= rhs, ('C', a, b, c, d, k)
        counts['C'] += 1
    r = min(a, b, c, d)
    s = sum(F(a+c-1-2*i, b+d-1-2*i) for i in range(r))
    assert lhs**3 >= r*s**3, ('D certificate', a, b, c, d)
    counts['D_sufficient'] += 1
    if a == r:
        assert lhs >= F(a, a)*F(c, b+d-a)
        counts['union'] += 1

for a in range(1, 33):
    assert phi(a)**3 >= a
    assert F(a, a) >= 2*a-1
    for z in range(1, 64):
        assert phi(z) <= phi(a)*(1+Q(z-a, 3*a))
        counts['phi_tangent'] += 1

for a, b in product(range(1, 101), repeat=2):
    r = min(a, b)
    assert r <= phi(r)**3 <= Q(3*r-1, 2)
    assert F(a, b)**3 <= Q(3*(a+b)-2, 4)**4
    assert F(a, b)**9 >= (a*b*(a+b-1))**4
    counts['sandwich'] += 1

for n in range(3, 65):
    assert F(n, n)**9 >= Q(27, 4)**4*(n-2)**9
    counts['cw_uniform'] += 1
    def w(i):
        return 0 if i == 0 else 2 if i == n-1 else 1
    original = {(i, j, i+j) for i in range(n) for j in range(n)
                if i+j < n}
    assert all(w(i)+w(j)-w(k) >= 0 for i, j, k in original)
    leading = {(i, j, k) for i, j, k in original
               if w(i)+w(j)-w(k) == 0}
    expected = ({(0, j, j) for j in range(n)}
                | {(i, 0, i) for i in range(n)}
                | {(i, n-1-i, n-1) for i in range(1, n-1)})
    assert leading == expected and len(leading) == 3*n-3
    small = {(i, j, k) for i, j, k in leading
             if i < n-1 and j < n-1 and k != 0}
    assert len(small) == 3*(n-2)
    counts['cw_support'] += 1

for n in range(1, 101):
    assert 2*F(n, n)**3 >= (Q(3, 8)*(4*n-1))**4
    counts['mixed_square'] += 1

# Polynomial identities: coefficient-vector operations, not sampling.
def mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def sub(a, b):
    return [(a[i] if i < len(a) else 0)
            - (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]

cube = mul(mul([1, 3], [1, 3]), [1, 3])
assert sub(mul([0, 1], cube), mul([0, 0, 0, 27], [1, 1])) == [0, 1, 9, 0, 0]
assert sub(mul([0, 0, 0, 27], [2, 3]), mul(cube, [-1, 3])) == [1, 6, 0, 0, 0]
assert 1000 > 729
assert 2*8**4 > 9**4
assert Q(512, 27) > 16
assert Q(27, 4) < 8
assert 1331 > 1296

assert counts == dict(union=18496, C=700416, D_sufficient=65536,
                     phi_tangent=2016, cw_uniform=62,
                     sandwich=10000, cw_support=62, mixed_square=100)
for name, count in counts.items():
    print(f'{name}: {count} exact checks, PASS')
print('Polynomial identities and rational constants: PASS')
print('No floating-point arithmetic or unresolved comparisons.')
