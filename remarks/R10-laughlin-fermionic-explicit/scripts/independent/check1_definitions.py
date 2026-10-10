#!/usr/bin/env python3
"""Item 1 (G2): definitions re-derived independently.
 (a) CG pair vectors vs paper eq. (pair) (03-spin.tex l.12-16)
 (b) H_{3,Q} = sum_{i<j} P_ij^(1) on the antisymmetric part of U^{(x)3} (first quantisation, 01-introduction l.28-35)
     equals sum_p B_p^dag B_p (02-operators l.221-224) -- coefficient-one normalisation
 (c) normal-ordered square, Lemma 2.1 (02-operators l.241-246), on N=3 and N=4
 (d) three-body Gram spectrum q_z(Q), Lemma 3.2 (03-spin l.48-53)
usage: check1_definitions.py Q"""
import sys, itertools
from math import comb, perm
import numpy as np, scipy.sparse as sp
from indep_fock import pair_vectors_cg, Fock, cg_sq

Q = int(sys.argv[1])
V = pair_vectors_cg(Q)

# (a) paper formula
def ff(s, p):
    r = 1
    for i in range(p): r *= (s - i)
    return r
err = 0.0
for p, d in V.items():
    nrm = sum(a * a for a in d.values())
    err = max(err, abs(nrm - 1))
    for i in range(Q + 1):
        for j in range(i + 1, Q + 1):
            pap = 0.0
            if i + j == p + 1:
                from fractions import Fraction
                pap = (i - j) * float(Fraction(ff(Q, i) * ff(Q, j) * ff(p, p), Q * ff(2 * Q - 2, p) * ff(i, i) * ff(j, j))) ** 0.5
            err = max(err, abs(d.get((i, j), 0.0) - pap))
print(flush=True); print(f"(a) Q={Q}: max |v_p(CG) - v_p(paper eq.pair)| over all p,i<j, and |norm-1| = {err:.2e}")

F = Fock(Q)
# (b) first quantisation for N=3 on U^{(x)3} restricted to wedge^3
n1 = Q + 1
# two-particle projector onto pair spin Q-1 in U(x)U (ordered basis): |JM> = sum_{x,y} CG |x>|y>
P2 = np.zeros((n1 * n1, n1 * n1))
for p, d in V.items():
    w = np.zeros(n1 * n1)
    for (x, y), a in d.items():
        w[x * n1 + y] += a / np.sqrt(2); w[y * n1 + x] -= a / np.sqrt(2)
    P2 += np.outer(w, w)
assert np.allclose(P2 @ P2, P2, atol=1e-10)
F.basis(3)
B3 = F.bases[3]
# isometry wedge^3 -> U^{(x)3}
Iso = np.zeros((n1 ** 3, len(B3)))
for col, s in enumerate(B3):
    for perm_ in itertools.permutations(range(3)):
        sg = np.linalg.det(np.eye(3)[list(perm_)])
        t = [s[k] for k in perm_]
        Iso[(t[0] * n1 + t[1]) * n1 + t[2], col] += sg / np.sqrt(6)
I1 = sp.identity(n1, format='csr')
P2s = sp.csr_matrix(P2)
P12 = sp.kron(P2s, I1)
P23 = sp.kron(I1, P2s)
# P13 = S23 P12 S23
S23 = sp.kron(I1, sp.csr_matrix(np.eye(n1 * n1)[[ (b % n1) * n1 + b // n1 for b in range(n1 * n1)]]))
P13 = S23 @ P12 @ S23
Hfirst = Iso.T @ ((P12 + P13 + P23) @ Iso)
Bs3 = {p: F.B(d, 3) for p, d in V.items()}
Hsecond = sum((Bs3[p].T @ Bs3[p] for p in Bs3)).toarray()
print(flush=True); print(f"(b) Q={Q} N=3: max|H_first-quantised - sum_p B_p^dag B_p| = {np.abs(Hfirst - Hsecond).max():.2e}  (||H||={np.abs(Hsecond).max():.3f})")

# (c) Lemma 2.1 on N=3 and N=4
def wedge_idx(N, tup):
    """sorted tuple + sign"""
    if len(set(tup)) < len(tup): return None, 0
    srt = sorted(tup); sg = np.linalg.det(np.eye(len(tup))[[tup.index(v) for v in srt]])
    return F.index[N][tuple(srt)], sg
F.basis(4); F.basis(2); F.basis(1)
# W3 : V (x) U -> wedge^3 ; columns indexed (p, j)
cols3 = [(p, j) for p in V for j in range(n1)]
W3 = np.zeros((len(F.bases[3]), len(cols3)))
for c, (p, j) in enumerate(cols3):
    for (x, y), a in V[p].items():
        i, sg = wedge_idx(3, [x, y, j])
        if i is not None: W3[i, c] += a * sg
X3 = W3 @ (W3.T @ W3 - np.eye(len(cols3))) @ W3.T          # matrix on wedge^3
H3 = Hsecond
lhs3 = H3 @ H3
rhs3 = H3 + X3          # L3 acts as itself on 3 particles, L4 vanishes
print(flush=True); print(f"(c) N=3: max|H^2 - (H + W3(W3^dag W3 - I)W3^dag)| = {np.abs(lhs3 - rhs3).max():.2e}")
# N=4: H^2 = H + L3(X3) + L4(W4 W4^dag), tested on random vectors (exact identity => residual ~ 1e-13)
#   L3(X) psi = sum_k Phi_k X Phi_k^T psi,  Phi_k(e_A) = e_A ^ e_k  (lift = adjoint-of-left-wedge convention,
#   <w, B(v) psi> = <v ^ w, psi>, 02-operators l.209-217);  L4(Y) = Y on 4 particles.
Bs4 = {p: F.B(d, 4) for p, d in V.items()}
H4 = sum((Bs4[p].T @ Bs4[p] for p in Bs4)).tocsr()
n4 = len(F.bases[4]); n3 = len(F.bases[3])
Phis = []
for k in range(n1):
    r_, c_, v_ = [], [], []
    for col, A in enumerate(F.bases[3]):
        if k in A: continue
        i, sg = wedge_idx(4, list(A) + [k])
        r_.append(i); c_.append(col); v_.append(sg)
    Phis.append(sp.csr_matrix((v_, (r_, c_)), shape=(n4, n3)))
cols4 = [(p, q) for p in V for q in V]
W4 = np.zeros((n4, len(cols4)))
for c, (p, q) in enumerate(cols4):
    for (x, y), a in V[p].items():
        for (u, w), b in V[q].items():
            i, sg = wedge_idx(4, [x, y, u, w])
            if i is not None: W4[i, c] += a * b * sg
rng = np.random.default_rng(1)
res = 0.0
for _ in range(5):
    psi = rng.standard_normal(n4)
    lhs = H4 @ (H4 @ psi)
    rhs = H4 @ psi + sum(Ph @ (X3 @ (Ph.T @ psi)) for Ph in Phis) + W4 @ (W4.T @ psi)
    res = max(res, np.abs(lhs - rhs).max() / np.abs(lhs).max())
print(flush=True); print(f"(c) N=4: relative residual of H^2 psi - (H + L3(W3(W3^dag W3-I)W3^dag) + L4(W4W4^dag)) psi, 5 random psi = {res:.2e}", flush=True)

# (d) q_z spectrum of W3^dag W3 - I on V (x) U, resolved by total spin via Casimir
Gm = W3.T @ W3 - np.eye(len(cols3))
# J_- on V(x)U: v_p -> sqrt((p+1)(2a-p)) v_{p+1}, a=Q-1 ;  e_j -> sqrt((j+1)(Q-j)) e_{j+1}
nV = 2 * Q - 1
idx3 = {c: k for k, c in enumerate(cols3)}
Jm = np.zeros((len(cols3), len(cols3)))
for k, (p, j) in enumerate(cols3):
    if p + 1 < nV: Jm[idx3[(p + 1, j)], k] += np.sqrt((p + 1) * (2 * (Q - 1) - p))
    if j + 1 <= Q: Jm[idx3[(p, j + 1)], k] += np.sqrt((j + 1) * (Q - j))
Jz = np.diag([(Q - 1 - p) + (Q / 2 - j) for (p, j) in cols3])
Cas = Jm.T @ Jm + Jz @ Jz + Jz        # J_+ J_- = J^2 - Jz^2 + Jz  -> J^2 = J_+J_- + Jz^2 - Jz ; use J_- ^T = J_+
Cas = (Jm @ Jm.T) + Jz @ Jz + Jz      # J^2 = J_- J_+ + Jz^2 + Jz
assert np.abs(Gm @ Cas - Cas @ Gm).max() < 1e-8
worst = 0.0
for z in range(0, Q + 1):
    Jt = (3 * Q - 2 - 2 * z) / 2       # spin of R_{3,z} = (2Q-2 + Q)/2 - z
    # highest vector of block z at level T=z: vectors with Jz = Jt
    sel = [k for k, (p, j) in enumerate(cols3) if abs((Q - 1 - p) + (Q / 2 - j) - Jt) < 1e-9]
    Cs = Cas[np.ix_(sel, sel)]
    w, U = np.linalg.eigh(Cs)
    hv = U[:, np.isclose(w, Jt * (Jt + 1), atol=1e-6)]
    assert hv.shape[1] == 1
    x = np.zeros(len(cols3)); x[sel] = hv[:, 0]
    qz_num = x @ Gm @ x
    qz = (-1) ** z * ff(Q, z) / ff(2 * Q - 2, z) * (3 * z - 1 - z * (z + 1) / Q)
    worst = max(worst, abs(qz_num - qz))
print(flush=True); print(f"(d) Q={Q}: max_z |<hw_z, (W3^dag W3 - I) hw_z> - q_z(Q) (eq. q)| = {worst:.2e}")
