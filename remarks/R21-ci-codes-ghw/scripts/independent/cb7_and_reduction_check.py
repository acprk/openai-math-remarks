#!/usr/bin/env python3
"""Checks of the two ingredients of the proof of Theorem G, on random subsets of random CIs
(re-uses only the finite-field helpers of ghw_ci_indep.py):

 (1) Cayley-Bacharach (CB7 of Eisenbud-Green-Harris 1996, after Davis-Geramita-Orecchia 1985),
     in the normalisation used in the note:  for Gamma subset X,
         h_X(d) - h_{X minus Gamma}(d) = |Gamma| - h_Gamma(s - 1 - d),   s = sum(a_i - 1);
 (2) the down-set conclusion of the Artinian reduction: the first difference of h_Gamma is the
     degree profile of a down-set D in the box prod[0, a_i - 1]; i.e. there is a down-set D with
         #{beta in D : |beta| = t} = h_Gamma(t) - h_Gamma(t-1)  for all t.
     (This is what EGH(a) delivers; for n = 2 it is a known theorem.  We test existence of D by
     exhaustive search over down-sets.)

Usage: python3 cb7_and_reduction_check.py [seed]
"""
import sys, itertools, random, importlib.util, os
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('g', os.path.join(here, 'ghw_ci_indep.py'))
g = importlib.util.module_from_spec(spec)
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 7
sys.argv = [sys.argv[0], str(SEED)]
spec.loader.exec_module(g)
rng = random.Random(SEED + 1)

def hilb_points(pts, t, nv, p):
    if not pts:
        return 0
    M = [[g.ev_mono(e, pt, p) for pt in pts] for e in g.monomials(nv, t)]
    return g.rank_mod_p(M, p)

def profiles(a):
    prof = set()
    for D in g.downsets_bruteforce(tuple(a)):
        s = sum(x - 1 for x in a)
        v = [0] * (s + 1)
        for b in D:
            v[sum(b)] += 1
        prof.add(tuple(v))
    return prof

def main():
    plan = [((2, 2), 11, 3), ((2, 3), 11, 3), ((3, 3), 11, 3), ((3, 4), 13, 2), ((2, 5), 13, 2)]
    plan3 = [((2, 2, 2), 7, 3), ((2, 2, 3), 11, 2)]
    n_cb = n_cb_fail = n_ds = n_ds_fail = 0
    for a, p, want in plan + plan3:
        nv = len(a) + 1
        PTS = g.proj_points(nv, p)
        s = sum(x - 1 for x in a)
        prof = profiles(a)
        got = 0; tries = 0
        while got < want and tries < 5000:
            tries += 1
            res = g.ci_P2(a[0], a[1], p, PTS) if nv == 3 else g.ci_split_P3(a[0], a[1], a[2], p, PTS)
            if res is None:
                continue
            got += 1
            X = res[1]
            hX = [hilb_points(X, t, nv, p) for t in range(s + 2)]
            for _ in range(40):
                k = rng.randint(1, len(X) - 1)
                Gam = rng.sample(X, k)
                Rest = [x for x in X if x not in Gam]
                hG = [hilb_points(Gam, t, nv, p) for t in range(s + 2)]
                for d in range(0, s + 1):
                    lhs = hX[d] - hilb_points(Rest, d, nv, p)
                    m = s - 1 - d
                    rhs = len(Gam) - (hG[m] if m >= 0 else 0)
                    n_cb += 1
                    if lhs != rhs:
                        n_cb_fail += 1
                        print('CB7 FAIL', a, p, d, lhs, rhs)
                dh = tuple(hG[t] - (hG[t - 1] if t > 0 else 0) for t in range(s + 1))
                n_ds += 1
                if dh not in prof or hG[s + 1] != len(Gam):
                    n_ds_fail += 1
                    print('DOWNSET FAIL', a, p, dh)
        print(f'type {a} p={p}: CIs {got}', flush=True)
    print('CB7 identity checks', n_cb, 'failures', n_cb_fail)
    print('down-set profile checks', n_ds, 'failures', n_ds_fail)

if __name__ == '__main__':
    main()
