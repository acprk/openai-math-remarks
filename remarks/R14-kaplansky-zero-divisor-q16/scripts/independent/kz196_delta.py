#!/usr/bin/env sage -python
"""q=16: explicit delta in Lemma 2.1 from the certificate  sum_W P(W)^2 <= (sum f) lam^{h-1}.
Also illustrates growth/decay of 1^T M^{h-1} 1 via the exact 3x3 equitable quotient."""
from sage.all import QQ, RealBallField, matrix, vector
R = RealBallField(200)
lam = QQ(7138933) / 7179627; sf = 288
for delta in [QQ(1)/600, QQ(1)/400, QQ(1)/360]:
    # need sf * lam^(h-1) <= exp(-2 delta h) for all h >= h0;  LHS/RHS is decreasing in h iff log(lam) + 2 delta < 0
    assert R(lam).log() + 2 * R(delta) < 0
    h = 1
    while not (R(sf) * R(lam)**(h - 1) <= (-2 * R(delta) * h).exp()): h += 1
    print(f'delta={delta}: holds for all h >= {h}')
def quotient(q):
    v = q*q+q+1; p = QQ(q+1)/v; o = QQ(1)/(q+1)**2; e = p*p
    return matrix(QQ, [[(v-4)*o, 3*o, 3*e], [(v-3)*e, 3*e, 2], [(v-3)*o, 2*o, 3*e]]), v
for q in [8, 16]:
    Q, v = quotient(q); sizes = vector(QQ, [v-3, 3, 3])
    for h in [10, 100, 1000]:
        # 1^T M^{h-1} 1 = sizes . Q^{h-1} 1   (class-constant vectors are invariant)
        val = sizes * (Q**(h-1)) * vector(QQ, [1, 1, 1])
        print(f'q={q} h={h}: sum_W P(W)^2 = {float(val):.6e}')
