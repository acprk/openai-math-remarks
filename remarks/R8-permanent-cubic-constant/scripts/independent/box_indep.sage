# Independent check of Lemma BOX (does not reuse box_check.sing).
# For X in Gamma_D: I = saturation of (F_X, entries of adj M) in S=K[s,t,u], M=sI+tD+uX.
# Check: (i) X in Gamma_D, (ii) u is a nonzerodivisor mod I (Z avoids u=0), (iii) every element of I of degree
# <= k-1 lies in Ann(Lambda)=(s^{q+1},t^{q+1},u^{c+1}), (iv) I+(u) subset (s^{q+1},t^{q+1},u),
# (v) len Z = dim S/(I+u) >= (q+1)^2.
import sys
def adj_sym(X, theta, Kf):
    k = X.nrows()
    S = PolynomialRing(Kf, 's,t,u', order='degrevlex'); s, t, u = S.gens()
    M = s*identity_matrix(S, k) + t*diagonal_matrix(S, theta) + u*X.change_ring(S)
    A = M.adjugate()
    return S, M, A
def check(Kf, k, q, X, theta, name):
    c = k - 2*q - 1
    S, M, A = adj_sym(X, theta, Kf); s, t, u = S.gens()
    Fx = M.det()
    lam = s^q*t^q*u^c
    ingam = all(A[i, j].monomial_coefficient(lam) == 0 for i in range(k) for j in range(k))
    I0 = S.ideal([Fx] + [A[i, j] for i in range(k) for j in range(k)])
    I, _ = I0.saturation(S.ideal([s, t, u]))
    Iu = I + S.ideal([u])
    lenZ = Iu.vector_space_dimension()
    # regularity of u: (I : u) == I
    ureg = (I.quotient(S.ideal([u])) == I)
    hp = I.hilbert_polynomial()
    box = S.ideal([s^(q+1), t^(q+1), u])
    contained = all(g in box for g in Iu.gens())
    ann = S.ideal([s^(q+1), t^(q+1), u^(c+1)])
    G = I.groebner_basis()
    low = [g for g in G if g.degree() <= k-1]
    inann = all(g in ann for g in low)
    # check (I)_d subset Ann for each d<=k-1 via a K-basis of I_d: I_d is spanned by monomial multiples of G of degree d
    print(f"  {name}: k={k} q={q} c={c} | X in Gamma: {ingam} | rank X={X.rank()} | u regular mod I: {ureg} | "
          f"Hilbert poly of S/I = {hp} | len S/(I+u) = {lenZ} >= (q+1)^2={(q+1)**2}: {lenZ >= (q+1)**2} | "
          f"I+(u) in (s^(q+1),t^(q+1),u): {contained} | GB elts of deg<=k-1 in Ann: {inann} ({len(low)})")
def run(Kf, label):
    set_random_seed(11)
    print("=== field", label)
    if Kf.is_prime_field():   # distinct theta_i = 100 i + random  (GF(p).gen() would be 1)
        rnd = lambda a, b: random_matrix(ZZ, a, b, x=-30, y=31).change_ring(Kf)
        th = lambda k: [Kf(100*i + ZZ.random_element(1, 99)) for i in range(1, k+1)]
    else:
        rnd = lambda a, b: random_matrix(Kf, a, b)
        th = lambda k: [Kf.gen()^i for i in range(1, k+1)]
    # q > c regime (BOX applies, original Prop. adjugate-length does not): k=7,q=2 (c=2? no: c=k-2q-1)
    for (k, q, kind) in [(7, 2, 'block'), (7, 2, 'rank'), (8, 2, 'block3'), (8, 3, 'zero'), (9, 3, 'block'), (9, 3, 'rank')]:
        r = k - 2*q
        X = matrix(Kf, k, k)
        if kind == 'rank':
            X = rnd(k, r-2) * rnd(r-2, k) if r > 2 else X
        elif kind == 'block':      # m=2, rho = r-3 on block 1 of size k-k//2... two blocks (n1, n2)
            n1 = k//2 if r == 3 else k-2
            rho = r - 3
            if rho > 0: X[:n1, :n1] = rnd(n1, rho)*rnd(rho, n1)
            X[:n1, n1:] = rnd(n1, k-n1)
        elif kind == 'block3':     # m=3: sizes (k-4,2,2) ranks (r-4,0,0)
            n1 = k-4; rho = r-4
            if rho > 0: X[:n1, :n1] = rnd(n1, rho)*rnd(rho, n1)
            X[:n1, n1:] = rnd(n1, 4); X[n1:n1+2, n1+2:] = rnd(2, 2)
        check(Kf, k, q, X, th(k), kind)
which = sys.argv[1] if len(sys.argv) > 1 else 'both'
if which in ('QQ', 'both'):
    run(QQ, "QQ")
    run(GF(1000003), "GF(1000003)")
if which in ('GF2', 'both'): run(GF(2^10, 'a'), "GF(2^10)")
