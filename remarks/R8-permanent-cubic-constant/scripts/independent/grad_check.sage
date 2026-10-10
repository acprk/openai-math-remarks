# Independent exact check of the gradient condition  nabla G_D(X)=0,  G_D(X)=[s^q t^q] det(X+sI+tD).
# d G_D / d X_ij = [s^q t^q] adj(X+sI+tD)_{ji}.  adj(sI-B) = sum_j s^{k-1-j} sum_{i<=j} c_i B^{j-i}
# (char poly of B = sum c_i s^{k-i}); this identity is characteristic-free.  Here B = -(X+tD), entries in F[t].
import sys, itertools
def grad_coeff(X, theta, q, F):
    k = X.nrows(); R = PolynomialRing(F, 't'); t = R.gen()
    B = -(X.change_ring(R) + t*diagonal_matrix(R, theta))
    cp = B.charpoly()                       # monic, coefficients in R
    c = [cp[k-i] for i in range(k+1)]       # c_0=1, ..., c_k
    j = k-1-q
    S = matrix(R, k, k); P = identity_matrix(R, k)
    # S = sum_{i=0}^{j} c_i B^{j-i}  (Horner)
    S = c[0]*identity_matrix(R,k)
    for i in range(1, j+1):
        S = S*B + c[i]*identity_matrix(R,k)
    return matrix(F, k, k, lambda a, b: S[a,b][q])   # [t^q] of each entry ([s^q]adj, transposed irrelevant for vanishing)

def brute_grad(X, theta, q, F):
    k = X.nrows(); R = PolynomialRing(F, ['s','t'] + ['x%d_%d'%(a,b) for a in range(k) for b in range(k)])
    s, t = R.gen(0), R.gen(1); xs = R.gens()[2:]
    Xs = matrix(R, k, k, xs)
    G = (Xs + s*identity_matrix(R,k) + t*diagonal_matrix(R, theta)).det()
    Gq = sum(cf*mon for cf, mon in zip(G.coefficients(), G.monomials()) if mon.degree(s)==q and mon.degree(t)==q)
    Gq = Gq/(s^q*t^q)
    sub = {xs[a*k+b]: X[a,b] for a in range(k) for b in range(k)}
    return matrix(F, k, k, lambda a, b: F(Gq.derivative(xs[b*k+a]).subs(sub)))  # transpose to match adj_{ji}

def family_point(F, k, blocks, ranks, rnd):
    """coordinate block-upper-triangular X: diagonal block p of size blocks[p] and rank<=ranks[p], all
    strictly-upper blocks random, strictly-lower blocks zero."""
    X = matrix(F, k, k); off = [0]
    for b in blocks: off.append(off[-1]+b)
    for p in range(len(blocks)):
        n, rho = blocks[p], ranks[p]
        if rho > 0:
            Dg = rnd(n, rho) * rnd(rho, n)
            X[off[p]:off[p+1], off[p]:off[p+1]] = Dg
        for p2 in range(p+1, len(blocks)):
            X[off[p]:off[p+1], off[p2]:off[p2+1]] = rnd(blocks[p], blocks[p2])
    return X

def run(F, label):
    set_random_seed(2026)
    if F.characteristic() == 0:
        rnd = lambda a, b: random_matrix(ZZ, a, b, x=-99, y=100).change_ring(F)
        theta_of = lambda k: [F(100*i + ZZ.random_element(1, 99)) for i in range(1, k+1)]
    else:
        g = F.gen()
        rnd = lambda a, b: random_matrix(F, a, b)
        theta_of = lambda k: [g^i for i in range(1, k+1)]       # distinct since ord(g) > k
    print("=== field:", label)
    # (0) implementation cross-check against brute-force symbolic differentiation, k=5, q=1 and k=6, q=2
    for (k, q) in [(5, 1), (6, 2)]:
        th = theta_of(k)
        for trial in range(2):
            X = rnd(k, k)
            assert grad_coeff(X, th, q, F) == brute_grad(X, th, q, F)
        print(f"  [impl check] k={k} q={q}: Cayley-Hamilton adjugate formula == brute-force derivative (2 random X)")
    cases = [
        ("random X (control)", 9, 3, None),
        ("rank r-2 family", 9, 3, ([9], [1])),
        ("k=9 q=3 family [[0,Y],[0,0]] blocks (4,5)", 9, 3, ([4, 5], [0, 0])),
        ("random X (control)", 22, 7, None),
        ("rank r-2 family", 22, 7, ([22], [6])),
        ("k=22 q=7 family m=2 blocks (16,6) ranks (5,0)", 22, 7, ([16, 6], [5, 0])),
        ("k=22 q=7 m=2 blocks (16,6) ranks (6,0) [violates Prop F, control]", 22, 7, ([16, 6], [6, 0])),
    ]
    for name, k, q, fam in cases:
        th = theta_of(k); assert len(set(th)) == k
        r = k - 2*q
        X = rnd(k, k) if fam is None else family_point(F, k, fam[0], fam[1], rnd)
        Gm = grad_coeff(X, th, q, F)
        nz = sum(1 for a in range(k) for b in range(k) if Gm[a, b] != 0)
        extra = ""
        if fam is not None:
            bl, rk = fam
            dim = (k*k - sum(b*b for b in bl))//2 + sum(rho*(2*n-rho) for n, rho in zip(bl, rk))
            extra = f" | family dim={dim} codim={k*k-dim} vs (2q+2)^2={(2*q+2)**2} | sum(rho)+m={sum(rk)+len(bl)} vs r-1={r-1}"
        print(f"  k={k} q={q} r={r} {name}: rank X={X.rank()}, nonzero gradient entries={nz}{extra}")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("QQ", "both"): run(QQ, "QQ")
    if which in ("GF2", "both"): run(GF(2^20, 'a'), "GF(2^20)")
