# Independent re-implementation (does not import or share code with ../hardcore_certificate.py,
# ../chain_check.py or ../theorem_c_exhaustive.py).
#  Part 1. Rigorous upper bounds for the two filtering moments I1, I2, in python-flint ball arithmetic (arb),
#          computed from the DEFINITIONS of D, A, E (not from the paper's L/B identity), with different
#          quadrature size and strip widths (N = 192, inner strip 0.4, outer strip 2*0.17).
#  Part 2. The alias-side / preliminary / coverage inequalities, in exact rationals (Fraction) except for
#          binary entropies, which are enclosed in arb balls.
#  Part 3. Proposition C on small groups, enumerating subgroups as KERNELS of all homomorphisms
#          (Z/q)^k -> (Z/q)^k (every subgroup arises: G/K has exponent | q and rank <= k), vectorised in numpy.
# Usage: python3 indep_check.py      Environment: CPython 3.10, python-flint 0.9.0, numpy 2.2; 1 thread; ~1-2 min.
import json, os, sys
from fractions import Fraction as Fr
from itertools import product
import numpy as np
from flint import arb, ctx

ctx.prec = 160
HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, '..', 'params.json')))
A_ = lambda s: arb(str(s))
LN2 = arb(2).log()
l2 = lambda x: x.log() / LN2
p2 = lambda x: (x * LN2).exp()
ch = lambda t: (t.exp() + (-t).exp()) / 2


def ball(x, err):            # widen ball x by +-err (err an arb upper bound)
    return x + arb(0, err.upper())


def moments(rc, sc, kf, kg, eta, tau, lam1, lam2, N=192, SG=A_('0.4'), BT=A_('0.17')):
    rc = [A_(c) for c in rc]; sc = [A_(c) for c in sc]
    kf, kg, eta, tau, lam1, lam2 = map(A_, (kf, kg, eta, tau, lam1, lam2))
    Dl = eta - tau; a = Dl / (2 * tau - eta)
    pi = arb.pi()
    xs = [2 * pi * j / N for j in range(N)]
    trig = lambda co, x: sum((c * ((j + 1) * x).cos() for j, c in enumerate(co)), arb(0))
    smax = lambda co, s: sum((abs(c) * ch((j + 1) * s) for j, c in enumerate(co)), arb(0))
    err = lambda B0, s: 2 * B0 / ((N * s).exp() - 1)
    r = [trig(rc, x) for x in xs]; s = [trig(sc, x) for x in xs]
    R, S = smax(rc, SG), smax(sc, SG)
    Zf = ball(sum((p2(t) for t in r), arb(0)) / N, err(p2(R), SG))
    Zg = ball(sum((p2(t) for t in s), arb(0)) / N, err(p2(S), SG))
    lZf, lZg = l2(Zf), l2(Zg)
    # normalised densities and their logs at the nodes
    f = [p2(t) / Zf for t in r]; g = [p2(t) / Zg for t in s]
    lf = [t - lZf for t in r]; lgg = [t - lZg for t in s]
    # log2 of int f^kf and int g^kg (needed for log f_*, log g_*); error: |f^k| <= 2^{k (R + |log Zf|)}
    Ifk = ball(sum((p2(kf * t) for t in lf), arb(0)) / N, err(p2(kf * (R + abs(lZf))), SG))
    Igk = ball(sum((p2(kg * t) for t in lgg), arb(0)) / N, err(p2(kg * (S + abs(lZg))), SG))
    lIfk, lIgk = l2(Ifk), l2(Igk)
    # inner integrands, bounds on the strip |Im w| <= SG
    Bfg = p2(R + S) / (Zf * Zg).lower()
    eQ = err(Bfg, SG)
    eJ1 = err(Bfg * (R + abs(lZf)), SG); eJ2 = err(Bfg * (S + abs(lZg)), SG)
    X1 = arb(0); X2 = arb(0)
    for i in range(N):
        Q = arb(0); J1 = arb(0); J2 = arb(0)
        for k in range(N):
            w = f[(i + k) % N] * g[k]
            Q += w; J1 += w * lf[(i + k) % N]; J2 += w * lgg[k]
        Q = ball(Q / N, eQ); J1 = ball(J1 / N, eJ1); J2 = ball(J2 / N, eJ2)
        D = (J1 + J2) / Q - l2(Q)                         # int K log2 K
        Aa = (kf * J1) / Q - lIfk                         # int K log2 f_*(x+.)
        E = (kg * J2) / Q - lIgk                          # int K log2 g_*
        u = Aa + E - D; v = 2 * D - Aa - E
        M2 = 2 + 2 * xs[i].cos()
        X1 += M2 * p2(lam1 * (Dl - u)); X2 += M2 * p2(lam2 * (a * v - u))
    T1 = X1 / N; T2 = X2 / N
    # outer strip |Im z| <= 2 BT:  ||f(.+iy)||^2 = Zf^-2 int 2^{2 r_y}, r_y coefficients r_j cosh(j y), |y| <= BT
    rB = [c * ch((j + 1) * BT) for j, c in enumerate(rc)]; sB = [c * ch((j + 1) * BT) for j, c in enumerate(sc)]
    nf = ball(sum((p2(2 * trig(rB, x)) for x in xs), arb(0)) / N, err(p2(2 * smax(rB, SG)), SG)) / Zf ** 2
    ng = ball(sum((p2(2 * trig(sB, x)) for x in xs), arb(0)) / N, err(p2(2 * smax(sB, SG)), SG)) / Zg ** 2
    c = ((nf - 1) * (ng - 1)).sqrt()
    assert c.upper() < 1
    # crude strip bounds (all moduli): |Q| >= 1-c, |ln Q| <= -ln(1-c), |J1| <= ||f|| ||g|| sup|log2 f|, ...
    one_c = 1 - c
    nrm = (nf * ng).sqrt()
    J1m = nrm * (smax(rc, BT) + abs(lZf)); J2m = nrm * (smax(sc, BT) + abs(lZg))
    lQm = -(one_c.log()) / LN2
    # expand u and a v - u in J1/Q, J2/Q, log2 Q and the two constants (exact algebra, keeps the small
    # coefficients kappa - 1):  u = (kf-1) J1/Q + (kg-1) J2/Q + log2 Q - lIfk - lIgk
    #   a v - u = (a(2-kf)-(kf-1)) J1/Q + (a(2-kg)-(kg-1)) J2/Q - (1+2a) log2 Q + (1+a)(lIfk + lIgk)
    q1 = J1m / one_c; q2 = J2m / one_c
    um = abs(kf - 1) * q1 + abs(kg - 1) * q2 + lQm + abs(lIfk) + abs(lIgk)
    wm = abs(a * (2 - kf) - (kf - 1)) * q1 + abs(a * (2 - kg) - (kg - 1)) * q2 + (1 + 2 * a) * lQm \
        + (1 + a) * (abs(lIfk) + abs(lIgk))
    e1 = lam1 * (abs(Dl) + um); e2 = lam2 * wm
    M2m = 2 + 2 * ch(2 * BT)
    eo1 = err(M2m * p2(e1), 2 * BT); eo2 = err(M2m * p2(e2), 2 * BT)
    U1 = l2(T1 + arb(0, eo1.upper())); U2 = l2(T2 + arb(0, eo2.upper()))
    return U1, U2, c, e1, e2, eo1, eo2


def fmt(x, d=12):
    return x.upper().str(d, radius=False)


ok = True
print('Part 1: filtering moments (arb, N=192, strips 0.4 / 0.34)')
for name in ('main', 'dense', 'paper'):
    if name == 'paper':
        C = P['paper']; rc, sc, kf, kg = C['rc'], C['sc'], C['kf'], C['kg']; eta, tau = C['eta'], C['tau']
    else:
        C = P[name]; rc, sc, kf, kg = P['rc'], P['sc'], P['kf'], P['kg']
        eta, tau = C['eta'], C['tau_run'] if name == 'main' else C['tau']
    U1, U2, c, e1, e2, eo1, eo2 = moments(rc, sc, kf, kg, eta, tau, C['lambda1'], C['lambda2'])
    mx = max(U1.upper(), U2.upper())
    print('  %-5s eta=%s tau=%s : log2 I1 < %s , log2 I2 < %s   (c < %s, crude strip exponents %s %s, outer err %s %s)'
          % (name, eta, tau, fmt(U1), fmt(U2), fmt(c, 5), fmt(e1, 4), fmt(e2, 4), fmt(eo1, 2), fmt(eo2, 2)))
    if name == 'paper':
        good = mx < arb('0.5168')
        print('        paper bound 0.5168 holds:', good)
    else:
        need = A_(eta) - A_(C['finite_circle_slack']) - 48 * A_(C['eps'])
        good = mx < need
        print('        need log2 max I < eta - slack - 48 eps = %s : %s ; margin > %s'
              % (fmt(need, 10), good, (need - arb(mx)).lower().str(4, radius=False)))
    ok &= bool(good)

print('\nPart 2: alias side, exact rationals + arb entropies')
C = P['main']
F = lambda s: Fr(str(s))
eta, tr, cl = F(C['eta']), F(C['tau_run']), F(C['claim'])
sg, mu, de, ze, ga, tp = (F(C[k]) for k in ('sigma', 'mu', 'delta', 'zeta', 'gamma', 'theta_prime'))
K, gd, ep = int(C['K']), F(C['guard']), F(C['eps'])
fa = lambda q: arb(q.numerator) / arb(q.denominator)
H = lambda q: -(fa(q) * l2(fa(q)) + fa(1 - q) * l2(fa(1 - q)))
rep = (1 - H(sg)) * fa(mu)
items = {
    'prelim  1/2-delta+zeta+1e-5': fa(Fr(1, 2) - de + ze + F(C['prelim_rounding'])),
    'lists   rep+(1-mu)/2+mu h(s/2)-g mu': rep + fa((1 - mu) / 2) + fa(mu) * H(sg / 2) - fa(ga * mu),
    'pairs   rep+1-mu+2mu h(s/2)-g mu-eta+guard': rep + fa(1 - mu - eta + gd - ga * mu) + 2 * fa(mu) * H(sg / 2),
    'gen     rep+mu+1/4+guard': rep + fa(mu + Fr(1, 4) + gd),
    'P\' sampling rep+gamma mu/2': rep + fa(ga * mu / 2),
    'checksum tau_run+20eps': fa(tr + 20 * ep),
    'trial division eta/2': fa(eta / 2),
}
for k_, v_ in items.items():
    g_ = v_.upper() < fa(tr + 20 * ep).upper() + arb('1e-30') if k_.startswith('checksum') else v_.upper() < fa(tr)
    print('  %-44s < %s   %s' % (k_, fmt(v_, 10), 'ok' if g_ else 'FAIL')); ok &= bool(g_)
sp = fa(mu) + H(mu) + fa(gd)
print('  %-44s < %s   %s' % ('sparse  h(mu)+mu+guard  (< eta)', fmt(sp, 10), 'ok' if sp < fa(eta) else 'FAIL')); ok &= bool(sp < fa(eta))
th = de / (K * sg * mu)
c1 = K * sg * mu < Fr(1, 2); c2 = th < tp; c3 = ga < sg * (1 - tp)
print('  coverage: K sigma mu = %s < 1/2 %s ; theta = %s < theta\' %s ; gamma < sigma(1-theta\') = %s %s ; Markov 1-theta/theta\' = %s'
      % (float(K * sg * mu), c1, float(th), c2, float(sg * (1 - tp)), c3, float(1 - th / tp)))
ok &= c1 and c2 and c3
tot = max(float(v_.upper()) for v_ in items.values())
print('  max component exponent %.8f ; claim %s ; gap %.3g' % (tot, C['claim'], float(cl) - tot)); ok &= float(cl) - tot > 0

print('\nPart 3: Proposition C via kernels of all endomorphisms (numpy)')


def kernels(p, e, k):
    q = p ** e
    X = np.array(list(product(range(q), repeat=k)), dtype=np.int64)        # |G| x k
    seen = set()
    mats = np.array(list(product(range(q), repeat=k * k)), dtype=np.int64).reshape(-1, k, k)
    for B in np.array_split(mats, max(1, len(mats) // 2000)):
        Y = np.einsum('gk,mjk->mgj', X, B) % q                              # images of all x under all A
        Z = ~(Y.any(axis=2))                                                # kernel masks
        for row in np.packbits(Z, axis=1):
            seen.add(row.tobytes())
    Gp = np.all((p * X) % q == 0, axis=1)
    res = []
    for kb in seen:
        mask = np.unpackbits(np.frombuffer(kb, dtype=np.uint8))[:len(X)].astype(bool)
        sz = int(mask.sum()); M = len(X) // sz; c = int((mask & Gp).sum())
        res.append((M, c))
    return res


for (p, e, k) in [(2, 2, 2), (2, 2, 3), (2, 3, 2), (2, 4, 2), (3, 2, 2), (3, 3, 2), (5, 2, 2)]:
    res = kernels(p, e, k)
    viol = sum(1 for M, c in res if c ** e * M < p ** (k * e))
    eq = sum(1 for M, c in res if c ** e * M == p ** (k * e))
    print('  (Z/%d)^%d : %4d subgroups (as kernels) ; violations of |K cap G[p]|^e M >= p^{ke}: %d ; equality: %d'
          % (p ** e, k, len(res), viol, eq))
    ok &= viol == 0
print('\nINDEPENDENT CHECK PASSES:', ok)
