# Rigorous upper bounds for the two continuous filtering moments of family 138
# (paper Prop. 3.1, eq. (3.6)), for the 8-term densities of this remark, with two tilting
# exponents kappa_f != kappa_g.
#
#   I1 = int_T M(x)^2 2^{l1 (Delta - u(x))} dx,   I2 = int_T M(x)^2 2^{l2 (a v(x) - u(x))} dx,
#   Delta = eta - tau, Gamma = 2 tau - eta, a = Delta / Gamma, M(x)^2 = 2 + 2 cos x.
#
# Method (same structure as the paper's Appendix A, with kappa split in two):
#   * N-point trapezoid rule for every inner integral (Z_f, Z_g, Z_{f,kf}, Z_{g,kg}, the two
#     convolutions P(x_j), Y_f(x_j), Y_g(x_j)) and for the outer integral, all in mpmath interval
#     arithmetic (outward rounding);
#   * inner quadrature error: |T_N h - int h| <= 2 B0 / (e^{N sig} - 1) for h analytic on
#     |Im z| <= sig with |h| <= B0 there (Trefethen-Weideman 2014, Thm 3.2), sig = SIG_IN;
#   * outer quadrature error: the outer integrands extend analytically to |Im z| <= 2 BETA, via
#     Q(z) = int f(x+w+iy/2) g(w-iy/2) dw and |Q - 1| <= ||f(.+iy/2)-1||_2 ||g(.-iy/2)-1||_2 =: c < 1.
#   Score identities (kappa_f, kappa_g):  L = log2 Q,  B_f = Y_f/P - log2 Z_f,  B_g = Y_g/P - log2 Z_g,
#     H = [log2 Z_{f,kf} - kf log2 Z_f] + [log2 Z_{g,kg} - kg log2 Z_g],
#     u = L + (kf-1) B_f + (kg-1) B_g - H,   v = (2-kf) B_f + (2-kg) B_g - 2L + H.
#
# The quadrature/strip method is the one of the paper's appendix; the code is adapted from the
# certificate written by the adversarial reviewer of the constant-optimisation round (reviewer script
# "certify.py"), parameterised and cleaned up here.
#
# Usage: python3 hardcore_certificate.py [main|dense|paper] [N]
# Environment: CPython 3.10, mpmath 1.3.0; single thread; ~1 s per configuration at N = 128.
import sys, json, os
from mpmath import iv, mp

iv.dps = 30
mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
PRM = json.load(open(os.path.join(HERE, 'params.json')))
which = sys.argv[1] if len(sys.argv) > 1 else 'main'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 128

if which == 'paper':
    C = PRM['paper']
    rc, sc, kf, kg = C['rc'], C['sc'], C['kf'], C['kg']
    eta, tau, l1, l2 = C['eta'], C['tau'], C['lambda1'], C['lambda2']
elif which == 'main':
    C = PRM['main']
    rc, sc, kf, kg = PRM['rc'], PRM['sc'], PRM['kf'], PRM['kg']
    eta, tau, l1, l2 = C['eta'], C['tau_run'], C['lambda1'], C['lambda2']
elif which == 'dense':
    C = PRM['dense']
    rc, sc, kf, kg = PRM['rc'], PRM['sc'], PRM['kf'], PRM['kg']
    eta, tau, l1, l2 = C['eta'], C['tau'], C['lambda1'], C['lambda2']
else:
    raise SystemExit('unknown configuration')

I = lambda s: iv.mpf(str(s))
rc = [I(c) for c in rc]; sc = [I(c) for c in sc]
kf, kg, eta, tau, l1, l2 = map(I, (kf, kg, eta, tau, l1, l2))
LN2 = iv.log(2)
cosh = lambda t: (iv.exp(t) + iv.exp(-t)) / 2
log2 = lambda t: iv.log(t) / LN2
pow2 = lambda t: iv.exp(t * LN2)
Dl = eta - tau; Gm = 2 * tau - eta; a = Dl / Gm
SIG_IN = I('0.5'); BETA = I('0.19')
up = lambda x: mp.nstr(mp.mpf(x.b), 12)
nb = lambda x, d=6: mp.nstr(mp.mpf(x.b), d)   # print an upper endpoint


def absmax(coefs, sig):   # sup of |sum c_j cos(j z)| on |Im z| <= sig
    return sum(abs(c) * cosh((j + 1) * sig) for j, c in enumerate(coefs))


def qerr(B0, sig):        # trapezoid error bound, normalised measure
    return 2 * B0 / (iv.exp(N * sig) - 1)


def widen(x, e):
    e = iv.mpf(e).b
    return iv.mpf([x.a - e, x.b + e])


print('configuration =', which, ' N =', N)
print('eta =', mp.nstr(mp.mpf(eta.a), 8), ' tau =', mp.nstr(mp.mpf(tau.a), 8), ' lambda =', mp.nstr(mp.mpf(l1.a), 6), mp.nstr(mp.mpf(l2.a), 6),
      ' kappa_f, kappa_g =', mp.nstr(mp.mpf(kf.a), 6), mp.nstr(mp.mpf(kg.a), 6))
print('Delta =', mp.nstr(mp.mpf(Dl.a), 8), ' Gamma =', mp.nstr(mp.mpf(Gm.a), 8), ' a = Delta/Gamma in', mp.nstr(mp.mpf(a.a), 10), mp.nstr(mp.mpf(a.b), 10))

xs = [2 * iv.pi * j / N for j in range(N)]
nc = max(len(rc), len(sc))
cosk = [[iv.cos((k + 1) * x) for k in range(nc)] for x in xs]
r = [sum(c * cosk[j][k] for k, c in enumerate(rc)) for j in range(N)]
s = [sum(c * cosk[j][k] for k, c in enumerate(sc)) for j in range(N)]
F = [pow2(t) for t in r]; G = [pow2(t) for t in s]
Rin, Sin = absmax(rc, SIG_IN), absmax(sc, SIG_IN)
eZf, eZg = qerr(pow2(Rin), SIG_IN), qerr(pow2(Sin), SIG_IN)
eZfk, eZgk = qerr(pow2(kf * Rin), SIG_IN), qerr(pow2(kg * Sin), SIG_IN)
eP = qerr(pow2(Rin + Sin), SIG_IN)
eYf = qerr(Rin * pow2(Rin + Sin), SIG_IN); eYg = qerr(Sin * pow2(Rin + Sin), SIG_IN)
print('inner quadrature error bounds: Z_f', nb(eZf, 3), ' P', nb(eP, 3), ' Y', nb(eYf if eYf.b > eYg.b else eYg, 3))

Zf = widen(sum(F) / N, eZf); Zg = widen(sum(G) / N, eZg)
Zfk = widen(sum(pow2(kf * t) for t in r) / N, eZfk); Zgk = widen(sum(pow2(kg * t) for t in s) / N, eZgk)
H = (log2(Zfk) - kf * log2(Zf)) + (log2(Zgk) - kg * log2(Zg))
lZf, lZg = log2(Zf), log2(Zg)

X1 = []; X2 = []; minD = None
for i in range(N):
    Pi = sum(F[(i + k) % N] * G[k] for k in range(N)) / N
    Yf = sum(F[(i + k) % N] * G[k] * r[(i + k) % N] for k in range(N)) / N
    Yg = sum(F[(i + k) % N] * G[k] * s[k] for k in range(N)) / N
    Pi, Yf, Yg = widen(Pi, eP), widen(Yf, eYf), widen(Yg, eYg)
    L = log2(Pi) - lZf - lZg
    Bf = Yf / Pi - lZf; Bg = Yg / Pi - lZg
    u = L + (kf - 1) * Bf + (kg - 1) * Bg - H
    v = (2 - kf) * Bf + (2 - kg) * Bg - 2 * L + H
    M2 = 2 + 2 * iv.cos(xs[i])
    X1.append(M2 * pow2(l1 * (Dl - u))); X2.append(M2 * pow2(l2 * (a * v - u)))
T1 = sum(X1) / N; T2 = sum(X2) / N

# outer strip: ||f(. + i y)||_2^2 = Z_f^{-2} int 2^{2 r_y} is nondecreasing in |y|; r_y has coefficients r_j cosh(j y)
rb = [c * cosh((j + 1) * BETA) for j, c in enumerate(rc)]
sb = [c * cosh((j + 1) * BETA) for j, c in enumerate(sc)]
rbv = [sum(c * cosk[j][k] for k, c in enumerate(rb)) for j in range(N)]
sbv = [sum(c * cosk[j][k] for k, c in enumerate(sb)) for j in range(N)]
e_rb = qerr(pow2(2 * absmax(rb, SIG_IN)), SIG_IN); e_sb = qerr(pow2(2 * absmax(sb, SIG_IN)), SIG_IN)
nf2 = widen(sum(pow2(2 * t) for t in rbv) / N, e_rb) / Zf ** 2
ng2 = widen(sum(pow2(2 * t) for t in sbv) / N, e_sb) / Zg ** 2
c = iv.sqrt((nf2 - 1) * (ng2 - 1))
print('||f_beta||^2 <', up(nf2), ' ||g_beta||^2 <', up(ng2), ' c = sup|Q-1| on strip <', up(c))
assert c.b < 1
Rb, Sb = absmax(rc, BETA), absmax(sc, BETA)
nrm = iv.sqrt(nf2 * ng2)
Bfm = nrm * (Rb + abs(lZf)) / (1 - c)      # |B_f| on the strip
Bgm = nrm * (Sb + abs(lZg)) / (1 - c)
mlogQ = -log2(1 - c)                        # >= -Re L on the strip
E1 = l1 * (abs(Dl) + mlogQ + abs(kf - 1) * Bfm + abs(kg - 1) * Bgm + abs(H))
cf = a * (2 - kf) - (kf - 1); cg = a * (2 - kg) - (kg - 1)
E2 = l2 * ((1 + 2 * a) * mlogQ + abs(cf) * Bfm + abs(cg) * Bgm + (1 + a) * abs(H))
M2s = 2 + 2 * cosh(2 * BETA)
eo1 = qerr(M2s * pow2(E1), 2 * BETA); eo2 = qerr(M2s * pow2(E2), 2 * BETA)
print('strip exponents E1 <', nb(E1), ' E2 <', nb(E2), ' outer quadrature errors <',
      nb(eo1, 3), nb(eo2, 3))
U1 = (T1 + eo1).b; U2 = (T2 + eo2).b
lI1 = mp.mpf(log2(iv.mpf(U1)).b); lI2 = mp.mpf(log2(iv.mpf(U2)).b)
Ccont = max(lI1, lI2)
print('RIGOROUS  log2 I1 <', mp.nstr(mp.mpf(lI1), 12), '  log2 I2 <', mp.nstr(mp.mpf(lI2), 12))
if which == 'paper':
    print('paper claims both < 0.5168 :', Ccont < mp.mpf('0.5168'))
else:
    eps = mp.mpf(C['eps']); sl = mp.mpf(C['finite_circle_slack'])
    Cfin = Ccont + sl
    print('finite-circle constant (Lemma 3.2 analogue) C_fin = C_cont + %s <' % C['finite_circle_slack'], mp.nstr(Cfin, 10))
    mse = -mp.mpf(eta.a) + Cfin + 48 * eps
    print('MSE exponent  -eta + C_fin + 48 eps <', mp.nstr(mse, 8), ' (must be < 0):', mse < 0)
    # rounding losses in the two domination steps (paper eqs. (4.13)-(4.14)): need 1 < l1 <= 10, l2 (4a+3) < 30
    print('domination: 1 < lambda1 <= 10:', bool(l1.a > 1 and l1.b <= 10), '  lambda2 (4a+3) <', nb(l2 * (4 * a + 3), 8),
          '< 30:', bool((l2 * (4 * a + 3)).b < 30), '  lambda2 > 1:', bool(l2.a > 1))
    texp = mp.mpf(tau.b) + 20 * eps; clm = mp.mpf(C['claim'])
    print('checksum run-time exponent tau + 20 eps =', mp.nstr(texp, 10), '< claimed exponent', C['claim'], ':', texp < clm)
    print('CERTIFIED:', bool(texp < clm and mse < 0 and c.b < 1 and l1.b <= 10 and (l2 * (4 * a + 3)).b < 30))
