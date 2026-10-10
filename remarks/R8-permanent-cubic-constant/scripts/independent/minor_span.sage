# Check: every partial derivative of G_D(X)=[s^q t^q]det(X+sI+tD) is a K-linear combination of the
# almost-principal (r-1)-minors of X  (row set U+{a}, column set U+{b}, |U|=r-2), r=k-2q.
import itertools, sys
def test(Kf, k, q, theta):
    r = k-2*q
    names = ['x%d_%d'%(a,b) for a in range(k) for b in range(k)]
    R = PolynomialRing(Kf, ['s','t']+names); s, t = R.gen(0), R.gen(1); xs = R.gens()[2:]
    X = matrix(R, k, k, xs)
    G = (X + s*identity_matrix(R,k) + t*diagonal_matrix(R, theta)).det()
    Gq = R(sum(cf*mon for cf, mon in zip(G.coefficients(), G.monomials()) if mon.degree(s)==q and mon.degree(t)==q)/(s^q*t^q))
    minors = set()
    for U in itertools.combinations(range(k), r-2):
        rest = [i for i in range(k) if i not in U]
        for a in rest:
            for b in rest:
                rows = sorted(U+(a,)); cols = sorted(U+(b,))
                minors.add(X.matrix_from_rows_and_columns(rows, cols).det())
    if r-1 == 1: minors = set(xs)
    minors = [mm for mm in minors if mm != 0]
    mons = sorted(set(mo for mm in minors for mo in mm.monomials()))
    idx = {mo:i for i, mo in enumerate(mons)}
    def vec(p):
        v = vector(Kf, len(mons))
        for cf, mo in zip(p.coefficients(), p.monomials()):
            if mo not in idx: return None
            v[idx[mo]] = cf
        return v
    Mspan = matrix(Kf, [vec(mm) for mm in minors])
    rk = Mspan.rank()
    ok = True
    for x in xs:
        d = Gq.derivative(x)
        v = vec(d)
        if v is None or matrix(Kf, list(Mspan.rows())+[v]).rank() != rk: ok = False; break
    print(f"  {Kf}: k={k} q={q} r={r}: {len(minors)} almost-principal {r-1}-minors; all k^2 partials in their span: {ok}")
for Kf in [QQ, GF(2^8,'a')]:
    for (k, q) in [(5,1), (6,1), (7,2), (7,1)]:
        th = [Kf(100*i+7*i*i) for i in range(1,k+1)] if Kf.is_prime_field() else [Kf.gen()^i for i in range(1,k+1)]
        test(Kf, k, q, th)
