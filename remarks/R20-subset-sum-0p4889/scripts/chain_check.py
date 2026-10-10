# Exponent bookkeeping for the whole algorithm of family 138 with the modifications of this remark:
#   (1) sigma-variant of alias extraction (|J cap M| = sigma m instead of m/2, exponentially many trials);
#   (2) new checksum hard core (certified separately by hardcore_certificate.py);
#   (3) the three repairs requested by the CONST-OPT adversarial review:
#       A  caps on alias lists / visited pairs are CUMULATIVE over all trials (one budget 2^{ceil(tau_run n)}),
#       B  the visited-pair exponent carries the bit-guard term +1e-5 (b <= 2^{1e-5 n}),
#       C  the checksum is run (and certified) at tau_run = 0.48888 < 0.4889, not at 0.4889.
# Every inequality is evaluated in mpmath interval arithmetic (outward rounding); decimals are exact.
# The script also recomputes the paper's own chain (sigma = 1/2) for comparison.
#
# Usage: python3 chain_check.py      Environment: CPython 3.10, mpmath 1.3.0; < 1 s.
import json, os
from mpmath import iv, mp

iv.dps = 40; mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
PRM = json.load(open(os.path.join(HERE, 'params.json')))
I = lambda s: iv.mpf(str(s))
LN2 = iv.log(2)
lg = lambda x: iv.log(x) / LN2
h = lambda p: -p * lg(p) - (1 - p) * lg(1 - p)
hi = lambda x: mp.nstr(mp.mpf(x.b), 10)
lo = lambda x: mp.nstr(mp.mpf(x.a), 10)

C = PRM['main']
eta, tau_run, claim = I(C['eta']), I(C['tau_run']), I(C['claim'])
sig, mu, delta, zeta = I(C['sigma']), I(C['mu']), I(C['delta']), I(C['zeta'])
gamma, thp, K, guard, eps = I(C['gamma']), I(C['theta_prime']), int(C['K']), I(C['guard']), I(C['eps'])
ok = True


def check(name, cond, txt):
    global ok
    ok = ok and bool(cond)
    print('  %-34s %-58s %s' % (name, txt, 'ok' if cond else 'FAIL'))


print('== new chain: eta = %s, tau_run = %s, claim O(2^{%s n}) ==' % (C['eta'], C['tau_run'], C['claim']))
print('   sigma = %s, mu = %s, delta = %s, gamma = %s, theta\' = %s, K = %d, guard = %s'
      % (C['sigma'], C['mu'], C['delta'], C['gamma'], C['theta_prime'], K, C['guard'], ))
rep = (1 - h(sig)) * mu                       # log2 of the number of trials per (s,q), divided by n
print('   trial-repetition exponent (1 - h(sigma)) mu in [%s, %s]' % (lo(rep), hi(rep)))

# Preliminary search (paper Prop. 2.4 with 0.01004 -> delta): log|Sigma(U)| <= (1/2 - delta + zeta) n + 1/2 on the
# Chebyshev event; the cap 2^{ceil(c n)} uses c = 1/2 - delta + zeta + 1e-5 so that the +1/2 is absorbed (as in the paper)
prelim = I('0.5') - delta + zeta + I(C['prelim_rounding'])
check('preliminary cap', prelim.b < tau_run.a, '1/2 - delta + zeta + 1e-5 = %s < tau_run' % hi(prelim))

# Coverage (paper (5.10)-(5.12) with m/2 -> sigma m and 17 -> K)
check('K disjoint sigma*m blocks fit', (K * sig * mu).b < 0.5, 'K sigma mu = %s < 1/2' % hi(K * sig * mu))
theta = delta / (K * sig * mu)
check('Markov threshold', theta.b < thp.a, 'theta = delta/(K sigma mu) = %s < theta\'' % hi(theta))
pm = 1 - theta / thp
check('useful-block probability', pm.a > 0, '1 - theta/theta\' >= %s > 0' % lo(pm))
cov = sig * (1 - thp)
check('many distinct mask weights', gamma.b < cov.a, 'gamma < sigma(1-theta\') = %s' % lo(cov))

# Sizes (expected, summed over all trials; Markov against the cumulative cap 2^{ceil(tau_run n)})
lists = rep + (1 - mu) / 2 + mu * h(sig / 2) - gamma * mu
pairs = rep + 1 - mu + 2 * mu * h(sig / 2) - gamma * mu - eta + guard
gen = rep + mu + I('0.25') + guard          # poly(b) <= 2^{1e-5 n} under the bit guard
psamp = rep + gamma * mu / 2                 # 2^{rep n} poly calls to the P' prime sampler, trial division <= 2^{gamma m/2}
check('filtered lists (cumulative)', lists.b < tau_run.a, 'rep+(1-mu)/2+mu h(sigma/2)-gamma mu = %s' % hi(lists))
check('visited pairs (cumulative, +guard)', pairs.b < tau_run.a, 'rep+1-mu+2mu h(sigma/2)-gamma mu-eta+1e-5 = %s' % hi(pairs))
check('list generation', gen.b < tau_run.a, 'rep+mu+1/4+1e-5 = %s' % hi(gen))
check('sampling the primes P\' (all trials)', psamp.b < tau_run.a, 'rep+gamma mu/2 = %s' % hi(psamp))
sparse = h(mu) + mu + guard
check('sparse distinctions (global)', sparse.b < eta.a, 'h(mu)+mu+1e-5 = %s < eta' % hi(sparse))
cs = tau_run + 20 * eps
check('checksum time', cs.b < claim.a, 'tau_run + 20 eps = %s' % hi(cs))
td = eta / 2
check('trial division for P', td.b < tau_run.a, 'eta/2 = %s' % hi(td))
tot = max(x.b for x in (prelim, lists, pairs, gen, psamp, cs, td, tau_run))
print('   all-history exponent tau_tot = max(...) = %s ; claim - tau_tot = %s > 0 absorbs poly(n+b+2), b <= n^c'
      % (mp.nstr(mp.mpf(tot), 10), mp.nstr(mp.mpf(claim.a) - mp.mpf(tot), 4)))
ok = ok and (claim.a - tot) > 0
print('   checksum certificate at (eta, tau_run, lambda) = (%s, %s, %s/%s): see hardcore_certificate.py main'
      % (C['eta'], C['tau_run'], C['lambda1'], C['lambda2']))

print('\n== paper chain (sigma = 1/2, mu = .0567, delta = .01004, gamma = .48925, eta = .51772) ==')
m0, d0, g0, e0 = I('0.0567'), I('0.01004'), I('0.48925'), I('0.51772')
print('  lists  n/2+(h(1/4)-1/2-gamma)m           = %s   (paper, via h(1/4)<.81128: .489909101)' % hi(I('0.5') + m0 * (h(I('0.25')) - I('0.5') - g0)))
print('  pairs  1+(2h(1/4)-1-gamma)mu-eta+1e-5   = %s   (paper, via h(1/4)<.81128: .489848677)' % hi(1 + m0 * (2 * h(I('0.25')) - 1 - g0) - e0 + guard))
print('  prelim cap 0.489985, checksum 0.48996 + 20 eps, all-history 0.489995 -> claim 0.49')

print('\n== comparison ==')
print('  exponent: 0.49 -> %s ; saving below 1/2: 0.0100 -> %s (x%s)' %
      (C['claim'], mp.nstr(0.5 - mp.mpf(C['claim']), 4), mp.nstr((0.5 - mp.mpf(C['claim'])) / mp.mpf('0.01'), 4)))
print('  running-time ratio old/new = 2^{%s n}' % mp.nstr(mp.mpf('0.49') - mp.mpf(C['claim']), 4))
print('  heuristic framework limit ~0.4798 (CONST-OPT): share of the headroom 1/2 - 0.4798 used: paper %s, new %s' %
      (mp.nstr(mp.mpf('0.01') / mp.mpf('0.0202'), 3), mp.nstr((0.5 - mp.mpf(C['claim'])) / mp.mpf('0.0202'), 3)))
print('\nALL CHAIN CHECKS PASS:', ok)
