#!/usr/bin/env python3
"""269: rigorous certificate for all Q >= 300 on x = 1/Q in [0, 1/300] (Taylor models of order 6, cert_box.py).
Adaptive bisection: a box [xa,xb] is accepted when cert_core.evaluate on the ball succeeds (every LDL pivot
certified > 0) and gamma_lo >= target; otherwise it is split (down to width wmin, else reported as a failure).
The box containing x = 0 covers every Q >= 1/xb (all formulas are analytic in x on [0,1/24]).
usage: certX.py data/gamma.json out.jsonl xa_num xa_den xb_num xb_den [target=0.04] [nstart=8] [rho_margin=1e-9]
(endpoints given as exact rationals; the union of output boxes equals [xa,xb])."""
import sys, json, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fractions import Fraction
from flint import arb, fmpq
import aengine as ae, cert_box as cb, rationalize as rz

gfile, outf = sys.argv[1], sys.argv[2]
xa = Fraction(int(sys.argv[3]), int(sys.argv[4]))
xb = Fraction(int(sys.argv[5]), int(sys.argv[6]))
target = float(sys.argv[7]) if len(sys.argv) > 7 else 0.04
nstart = int(sys.argv[8]) if len(sys.argv) > 8 else 8
rmarg = float(sys.argv[9]) if len(sys.argv) > 9 else 1e-10
WMIN = Fraction(1, 10 ** 8)
G = ae.Gamma(rz.load(gfile))


def ball(a, b):
    return arb(fmpq(a.numerator, a.denominator)).union(arb(fmpq(b.numerator, b.denominator)))


stack = [(xa + (xb - xa) * i / nstart, xa + (xb - xa) * (i + 1) / nstart) for i in range(nstart)][::-1]
fo = open(outf, 'w')
nacc = nfail = 0
gmin = 1e9
T0 = time.time()
while stack:
    a, b = stack.pop()
    t0 = time.time()
    try:
        o = cb.certify_box(a, b, G, rho_margin=rmarg)
        ok = (not o['fail']) and o['gamma_lo'] is not None and o['gamma_lo'] >= target
    except (ValueError, ZeroDivisionError) as e:
        o = {'fail': [str(e)], 'gamma_lo': None}
        ok = False
    if ok:
        nacc += 1
        gmin = min(gmin, o['gamma_lo'])
        rec = dict(xa=str(a), xb=str(b), gamma_lo=o['gamma_lo'], eta=o['eta'], tail_up=o['tail_up'],
                   cost3_up=o['cost3_up'], cost4_up=o['cost4'], rho={D: str(v) for D, v in o['rho'].items()},
                   sec=round(time.time() - t0, 2))
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        print(f"ACC [{float(a):.6e},{float(b):.6e}] gamma_lo={o['gamma_lo']:.6f} ({time.time()-t0:.1f}s)", flush=True)
    else:
        if b - a <= WMIN:
            nfail += 1
            fo.write(json.dumps(dict(xa=str(a), xb=str(b), FAIL=str(o['fail']))) + '\n'); fo.flush()
            print("FAIL", a, b, o['fail'], flush=True)
            continue
        print(f"split [{float(a):.6e},{float(b):.6e}] fail={o['fail']} g={o.get('gamma_lo')}", flush=True)
        m = (a + b) / 2
        stack.append((m, b))
        stack.append((a, m))
fo.close()
print(f"DONE accepted={nacc} failed={nfail} min gamma_lo={gmin:.6f} time={time.time()-T0:.0f}s", flush=True)
