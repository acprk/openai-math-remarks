#!/usr/bin/env python3
"""Fast (<1 s) check of the three parameter points quoted in the remark, condition by
condition, using ok()/min_m0() from indep_budget_053.py.  Prints the slack of each
inequality of 08-support.tex (Section 8) at each point."""
import importlib.util, os
from math import comb
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("b", os.path.join(here, "indep_budget_053.py"))
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)

def report(name, g, K, m0, nb):
    CK = comb(2*K, 3); MK = max(10*K + 3, CK + 1); e = 2*(m0//4); S0 = g - nb*m0
    loss = 3*g*MK + g*(10*K + 2) + 2*K*CK; tmax = 25*(g-1)**2 - 1
    print(f"== {name}: g={g} (~{g:.4e}), K={K}, m0={m0}, blocks={nb}, |S0|={S0} ({S0/g:.3f} g)")
    rows = [
      ("P2  g > 2K", g - 2*K),
      ("P3  4(g-K-2) - (3g-3)", 4*(g-K-2) - (3*g-3)),
      ("P4  (g-K)^2-1 - (3g-2)", (g-K)**2 - 1 - (3*g-2)),
      ("P5  d_e/2 - (loss + t_max)   [eq. small-loss + large-test-threshold]", comb(e,3)/2 - (loss + tmax)),
      ("P6  11g - (C_K+1)", 11*g - (CK+1)),
      ("P7  S0(g(K-4)+2) - 280g^2   [large-deficit contradiction]", S0*(g*(K-4)+2) - 280*g*g),
      ("P8  blocks - (3K+5)", nb - (3*K+5)),
    ]
    for t, v in rows: print(f"   {t:70s} slack = {v:.6e}" + ("" if v > 0 or (t.startswith('P5') or t.startswith('P8')) and v >= 0 else "   <-- FAILS"))
    res = b.ok(g, K, m0, nb); print("   all conditions:", res); return res

assert report("paper (g=10^60, K=10^5, m0=2000 g^(2/3), 4K+20)", 10**60, 10**5, 2000*10**40, 4*10**5 + 20)
g, K = 226327553857834, 568
assert report("this remark, paper's 4K+20 reservation", g, K, b.min_m0(g, K), 4*K + 20)
print("   g-1 with optimal m0 feasible?", b.ok(g-1, K, b.min_m0(g-1, K), 4*K + 20))
g2, K2 = 93824084405552, 567
assert report("this remark, the paper's own count 3K+5", g2, K2, b.min_m0(g2, K2), 3*K2 + 5)
print("   g-1 with optimal m0 feasible?", b.ok(g2-1, K2, b.min_m0(g2-1, K2), 3*K2 + 5))
print("ratio 10^60 / g =", f"{10**60/g:.3e}", "; markings n = 3*C(g,3) =", f"{3*comb(g,3):.4e}", "vs", f"{3*comb(10**60,3):.4e}")
