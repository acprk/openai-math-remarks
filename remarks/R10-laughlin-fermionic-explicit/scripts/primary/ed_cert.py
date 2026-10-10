#!/usr/bin/env python3
"""269: rigorous gap certificate for N = 2..8 at Laughlin flux Q = 3(N-1) (Q <= 21, outside the SOS range).

Claim certified for each N:   H_{N,Q} >= g_N (I - P_L)  on  wedge^N U_Q.
Proof structure:
 (1) Lz=0 sector (sum of orbital labels = QN/2) meets every SU(2) multiplet (all spins are integers since QN/2 is),
     so it suffices to bound the spectrum of H restricted to it (the Laughlin zero mode is simple, paper Lemma 6.2).
 (2) Exact matrix: H_{S',S} = sqrt(beta(S) beta(S')) * Htilde_{S',S}, Htilde rational (beta(S)=prod_{a in S} C(Q,a)),
     computed with Fractions; the float matrix is the correctly rounded value (arb), so |Hf - H| <= 2^-52 |H| entrywise.
 (3) X = H - g I + psi psi^T with psi ANY float vector (we use the computed ground vector).  If X > 0 then
     v^T H v >= g |v|^2 for v orthogonal to psi, hence lambda_2(H) >= g by Courant-Fischer; lambda_1 = 0 (Laughlin).
 (4) X > 0 is certified by float Cholesky of A = fl(X) - tau I and Higham, Thm 10.3: R^T R = A + dA,
     |dA| <= gamma_{n+1} |R^T||R|, so ||dA||_2 <= gamma_{n+1} tr(R^T R); together with the entrywise rounding
     error of fl(X) (Frobenius bound) we need  tau > 10*gamma_{n+1} tr(R^T R) + ||X - fl(X)||_F.
usage: ed_cert.py Nmax"""
import sys, itertools, time, json
import numpy as np
from fractions import Fraction
from math import comb, factorial
from flint import arb, fmpq, ctx

ctx.prec = 128
U = 2.0 ** -53


def ff(s, p):
    r = 1
    for t in range(p):
        r *= (s - t)
    return r


def sgn_cc(S, x, y):
    """sign of <S\\{x,y}| c_y c_x |S>, S sorted tuple containing x,y"""
    px = S.index(x)
    S1 = S[:px] + S[px + 1:]
    py = S1.index(y)
    return (-1) ** (px + py), S1[:py] + S1[py + 1:]


def run(N):
    Q = 3 * (N - 1)
    M = Q * N // 2
    sec = [S for S in itertools.combinations(range(Q + 1), N) if sum(S) == M]
    idx = {S: k for k, S in enumerate(sec)}
    n = len(sec)
    beta = {}
    def B(S):
        if S not in beta:
            r = 1
            for a in S:
                r *= comb(Q, a)
            beta[S] = r
        return beta[S]
    Ht = {}
    for S in sec:
        for (xp, yp) in itertools.combinations(S, 2):          # removed pair (xp<yp), p = xp+yp-1
            p = xp + yp - 1
            if p > 2 * Q - 2:
                continue
            s1, T = sgn_cc(S, xp, yp)
            base = Fraction(factorial(p), Q * ff(2 * Q - 2, p) * B(T)) * s1 * (xp - yp)
            for x in range(0, (p + 1) // 2 + 1):
                y = p + 1 - x
                if x >= y or y > Q or x in T or y in T:
                    continue
                S2 = tuple(sorted(T + (x, y)))
                s2, _ = sgn_cc(S2, x, y)
                key = (idx[S2], idx[S])
                Ht[key] = Ht.get(key, 0) + base * s2 * (x - y)
    H = np.zeros((n, n))
    for (i, j), v in Ht.items():
        if v == 0:
            continue
        a = arb(fmpq(v.numerator, v.denominator)) * (arb(B(sec[i])) * arb(B(sec[j]))).sqrt()
        H[i, j] = float(a.mid())
    assert np.allclose(H, H.T, atol=1e-12 * max(1, np.abs(H).max()))
    H = 0.5 * (H + H.T)   # exact matrix is symmetric; symmetrisation changes entries by <= 1 ulp
    w, V = np.linalg.eigh(H)
    k0 = int((w < 1e-9).sum())
    assert k0 == 1, ("kernel dim", k0)
    gap = w[1]
    g = Fraction(int(gap * 10 ** 4) - 10, 10 ** 4)          # certified constant (rational, below the float gap)
    tau = (gap - float(g)) / 4
    psi = V[:, 0].copy()
    X = H - float(g) * np.eye(n) + np.outer(psi, psi)
    # entrywise error of fl(X) vs exact X(psi): H rounding (2u|H|), symmetrisation (u|H|), g*I, outer product, two adds
    err = (4 * U) * (np.abs(H) + float(g) * np.eye(n) + np.abs(np.outer(psi, psi))) + 4 * U * np.abs(X)
    errF = float(np.sqrt((err ** 2).sum())) * 1.01
    A = X - tau * np.eye(n)
    R = np.linalg.cholesky(A)            # raises LinAlgError if not PD in floating point
    trRR = float((R ** 2).sum())
    gam = (n + 1) * U / (1 - (n + 1) * U)
    bound = 10 * gam * trRR + errF
    ok = bound < tau
    return dict(N=N, Q=Q, dim=n, gap_float=float(gap), g_cert=str(g), tau=tau, chol_err_bound=bound, ok=bool(ok),
                second=float(w[2]) if n > 2 else None)


if __name__ == '__main__':
    out = []
    for N in range(2, int(sys.argv[1]) + 1):
        t = time.time()
        r = run(N)
        r['sec'] = round(time.time() - t, 1)
        print(json.dumps(r), flush=True)
        out.append(r)
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], 'w'), indent=1)
