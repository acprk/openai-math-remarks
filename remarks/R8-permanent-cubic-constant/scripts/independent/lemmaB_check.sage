# Lemma B counterexample (independent of bipartite_order.sing).  k=16, X=[[0,Y],[0,0]] with 8x8 blocks.
# Chart u=1: M = X + x I + y D.  Z_X near the vertex p=(0,0) is cut out by F=det M and the entries of adj M.
# Claim: every adj entry has order >= k-2 at p  =>  local length of Z at p >= C(k-1,2)=105 for EVERY Y,
# so dim{X : len Z_X >= 105} >= 64 > k^2+k-2*105 = 62, contradicting Lemma B.
import sys
k = int(sys.argv[1]) if len(sys.argv) > 1 else 16
n1 = k//2
set_random_seed(7)
F = QQ
Ry = PolynomialRing(F, 'y'); y = Ry.gen()
theta = [F(100*i + ZZ.random_element(1, 99)) for i in range(1, k+1)]
X = matrix(F, k, k)
X[:n1, n1:] = random_matrix(ZZ, n1, k-n1, x=-99, y=100)
B = -(X.change_ring(Ry) + y*diagonal_matrix(Ry, theta))
cp = B.charpoly(); c = [cp[k-i] for i in range(k+1)]
S = [identity_matrix(Ry, k)]
for j in range(1, k):
    S.append(S[-1]*B + c[j]*identity_matrix(Ry, k))
R2 = PolynomialRing(F, 'x,y', order='degrevlex'); xx, yy = R2.gens()
adj = [[sum(xx^(k-1-j) * S[j][a, b](yy) for j in range(k)) for b in range(k)] for a in range(k)]
Fx = prod(xx + th*yy for th in theta)
# sanity: adj * M == det(M) * I   (exact)
M2 = X.change_ring(R2) + xx*identity_matrix(R2, k) + yy*diagonal_matrix(R2, theta)
assert matrix(R2, adj) * M2 == Fx * identity_matrix(R2, k), "adjugate identity failed"
assert M2.det() == Fx
gens = [Fx] + [adj[a][b] for a in range(k) for b in range(k) if adj[a][b] != 0]
def order(p): return min(sum(e) for e in p.exponents())
mo = min(order(g) for g in gens[1:])
print(f"k={k}: adj*M == det(M)*I verified; F = prod(x+theta_i y); min order of adj entries at p = {mo} (claim >= k-2 = {k-2})")
mm = R2.ideal([xx, yy])
prev = None
for N in range(k-2, k+4):
    I = R2.ideal(gens) + mm^N
    L = I.vector_space_dimension()
    print(f"  dim_Q Q[x,y]/(F, adj, m^{N}) = {L}")
    if prev is not None and L == prev:
        print(f"  stabilised => local length of Z_X at p = {L}  (C(k-1,2) = {binomial(k-1,2)}; delta at p = C(k,2) = {binomial(k,2)})")
        break
    prev = L
print(f"  family dim = {n1*(k-n1)} ; Lemma B bound k^2+k-2*C(k-1,2) = {k*k+k-2*binomial(k-1,2)}")
