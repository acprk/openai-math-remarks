#!/usr/bin/env python3
"""Negative controls for indep_check_059.py: perturb the new relation in two ways
(one coefficient +1; one coefficient sign-flipped) and confirm that the product
identities AND the sampled mod-2 spectral test both reject it."""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", os.path.join(here, "indep_check_059.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
d, rel = m.load(sys.argv[1])
for name, i, f in (("e_b of row 5 plus one", 5, lambda e: e + 1), ("sign of row 100 flipped", 100, lambda e: -e)):
    rel2 = list(rel); b, e = rel2[i]; rel2[i] = (b, f(e))
    okA, okR, _, _ = m.products_direct(d, rel2)
    vA, vR = m.products_valuation(d, rel2)
    info = m.seed_multisets(d, rel2, True)   # variable counts still balance (4E)
    t, err = m.spectral_test(d, rel2, info, nsamples=100)
    rejected = (not okA or not okR) and (not vA or not vR) and err > 1e-3
    print("%-26s b=%d: direct A,R=%s,%s  valuations A,R=%s,%s  spectral max|diff|=%.3f over %d chars -> %s"
          % (name, b, okA, okR, vA, vR, err, t, "REJECTED" if rejected else "NOT REJECTED"))
    assert rejected
print("NEGATIVE CONTROLS OK")
