#!/usr/bin/env python3
"""269: rigorous certificate for every integer flux Q in [Qlo,Qhi] with one rational Gamma (exact point x=1/Q).
usage: certQ.py Qlo Qhi data/gamma_X.json out.jsonl
Each output line: Q, certified lower bound gamma_lo(Q), pieces, the rational allowances rho_D (strings)."""
import sys, json, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from flint import arb, fmpq
import aengine as ae, cert_core as cc, rationalize as rz

Qlo, Qhi, gfile, outf = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
G = ae.Gamma(rz.load(gfile))
with open(outf, 'w') as fo:
    for Q in range(Qlo, Qhi + 1):
        t0 = time.time()
        assert Q >= 24
        o = cc.evaluate(arb(fmpq(1, Q)), G, Qint=Q)
        rec = dict(Q=Q, gamma_lo=o['gamma_lo'], eta=o['eta'], tail_up=o['tail_up'], cost3_up=o['cost3_up'],
                   cost4_up=o['cost4'], fail=o['fail'], rho={D: str(v) for D, v in o['rho'].items()},
                   gamma_file=gfile, sec=round(time.time() - t0, 2))
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        print(Q, o['gamma_lo'], o['fail'], flush=True)
