#!/usr/bin/env python3
"""269: rank-free re-check of every integer-Q certificate (b-space form, cert_core.bspace_check).
usage: bspaceQ.py certQ_X.jsonl out.jsonl   (re-uses Q, Gamma file and rho_D of the r-space certificate)"""
import sys, json, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fractions import Fraction
import aengine as ae, cert_core as cc, rationalize as rz
recs = [json.loads(l) for l in open(sys.argv[1])]
Gs = {}
with open(sys.argv[2], 'w') as fo:
    for r in recs:
        if r['gamma_file'] not in Gs:
            Gs[r['gamma_file']] = ae.Gamma(rz.load(r['gamma_file']))
        G = Gs[r['gamma_file']]
        rho = {int(D): Fraction(v) for D, v in r['rho'].items()}
        t = time.time()
        ok, c4, fails = cc.bspace_check(r['Q'], G, rho)
        gam_b = 1 - float(G.eta_q) - r['tail_up'] - r['cost3_up'] - c4 - 1e-14
        rec = dict(Q=r['Q'], ok=ok, fails=fails, cost4_b=c4, gamma_lo_bspace=gam_b, gamma_lo_rspace=r['gamma_lo'], sec=round(time.time() - t, 1))
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        print(rec, flush=True)
