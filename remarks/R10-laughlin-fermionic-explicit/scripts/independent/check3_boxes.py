#!/usr/bin/env python3
"""Item 3: Q >= 300 box certificate, adversarial spot checks.
 (0) exact contiguous cover of [0,1/300] by the 168 accepted boxes, no FAIL record (own re-implementation).
 (1) for selected boxes: Taylor coefficients c_0..c_6 of scalar cost functions (tail, three-body costs c_z) and of
     the 4-body matrices G_J, P (D in a set) re-derived by POLYNOMIAL INTERPOLATION of point values at 13 Chebyshev
     nodes of [x0-1.2h, x0+1.2h] (256-bit point evaluation, mpmath 90-digit solve) and compared with the
     certificate's arb_series coefficients at the exact centre; the certificate's order-6 remainder enclosure |F_6|
     (series at the ball base) is compared with the interpolated c_6 and with the sampled point-wise c_6.
 (2) point-wise re-certification inside each selected box with the box's own rho_D:
     r-space evaluate() at the two endpoints and the centre, and the RANK-FREE b-space check (no Lemma R)
     at integer fluxes Q with 1/Q in the box; gamma(point) must be >= gamma_lo(box).
usage: check3_boxes.py idx1,idx2,...   (indices into the sorted box list)"""
import sys, json, glob
from fractions import Fraction
import mpmath as mp
import os
PRIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'primary')
sys.path.insert(0, PRIM)
from flint import arb, fmpq, ctx
import aengine as ae, cert_core as cc, cert_box as cb, rationalize as rz

mp.mp.dps = 90
OUT = PRIM + '/out/'
recs = [json.loads(l) for f in sorted(glob.glob(OUT + 'certX_*.jsonl')) for l in open(f)]
assert not any('FAIL' in r for r in recs)
boxes = sorted(((Fraction(r['xa']), Fraction(r['xb']), r) for r in recs), key=lambda b: b[0])
cover = boxes[0][0] == 0 and boxes[-1][1] == Fraction(1, 300) and all(boxes[i][1] == boxes[i + 1][0] for i in range(len(boxes) - 1)) \
    and all(b[0] < b[1] for b in boxes)
print(f"(0) {len(boxes)} boxes, exact contiguous cover of [0,1/300], no FAIL records: {cover}; "
      f"box 0 = [{boxes[0][0]},{boxes[0][1]}] covers every Q >= {1/boxes[0][1]}", flush=True)
G = ae.Gamma(rz.load(PRIM + '/data/gamma_300_1000_100000000.json'))
sel = [int(v) for v in sys.argv[1].split(',')]
DSEL = [13, 19, 21, 23]

def to_mp(a):
    return mp.mpf(a.mid().str(80, radius=False).replace(' ', '')) if hasattr(a, 'mid') else mp.mpf(a)

def flat(res):
    """scalar functions of interest from a pipeline result, coefficient lists"""
    out = {'tail': res['tail']}
    for z in (2, 5, 7, 9):
        out[f'c3_{z}'] = res['c3'][z]
    return out

for bi in sel:
    xa, xb, rec = boxes[bi]
    xm, h = (xa + xb) / 2, (xb - xa) / 2
    print(f"\n=== box #{bi}: x in [{xa}, {xb}]  (Q in [{float(1/xb) if xb else 'inf'}, {float(1/xa) if xa else 'inf'}])  gamma_lo(box)={rec['gamma_lo']:.6f}", flush=True)
    A = cb.pipeline(arb(fmpq(xm.numerator, xm.denominator)), xm, G)
    Jc = {D: v['J'] for D, v in A['four'].items()}
    X = arb(fmpq(xa.numerator, xa.denominator)).union(arb(fmpq(xb.numerator, xb.denominator)))
    B = cb.pipeline(X, xm, G, J_choice=Jc)
    # Chebyshev nodes
    nn = 13
    nodes = [mp.mpf(xm.numerator) / xm.denominator + mp.mpf(1.2) * mp.mpf(h.numerator) / h.denominator * mp.cos(mp.pi * (2 * k + 1) / (2 * nn)) for k in range(nn)]
    pts = []
    for xn in nodes:
        fr = Fraction(mp.nstr(xn, 60, strip_zeros=False).replace('e', 'E'))  # exact rational node
        pts.append((fr, cb.pipeline(arb(fmpq(fr.numerator, fr.denominator)), xm, G, J_choice=Jc)))
    t = [mp.mpf(fr.numerator) / fr.denominator - mp.mpf(xm.numerator) / xm.denominator for fr, _ in pts]
    Vd = mp.matrix([[tk ** i for i in range(nn)] for tk in t])
    def interp(vals):
        return mp.lu_solve(Vd, mp.matrix(vals))
    hm = mp.mpf(h.numerator) / h.denominator
    worst_rel, worst_rem = 0, 0
    lines = []
    # scalars
    fA, fB = flat(A), flat(B)
    for name in fA:
        vals = [to_mp(cb.co(flat(r)[name], 0)) for _, r in pts]
        c = interp(vals)
        for i in range(6):
            ci = to_mp(cb.co(fA[name], i))
            scale = max(abs(ci), abs(c[i]), mp.mpf(10) ** -40)
            worst_rel = max(worst_rel, abs(ci - c[i]) * hm ** i / max(abs(vals[0]), mp.mpf(10) ** -40))
        F6 = cb.co(fB[name], 6)
        F6abs = max(abs(float(F6.upper())), abs(float(F6.lower())))
        c6pts = max(abs(to_mp(cb.co(flat(r)[name], 6))) for _, r in pts)
        ok6 = F6abs >= float(c6pts) * (1 - 1e-9) and F6abs >= float(abs(c[6])) * (1 - 1e-6)
        worst_rem = max(worst_rem, 0 if ok6 else 1)
        lines.append(f"   {name:6s}: c0={float(c[0]):+.6e} c1={float(c[1]):+.3e} c6(interp)={float(c[6]):+.3e} "
                     f"max_nodes|c6|={float(c6pts):.3e} <= |F6|enclosure={F6abs:.3e}: {ok6}")
    # matrices
    for D in DSEL:
        if D not in A['four']: continue
        for key in ('GJ', 'P'):
            MA, MB = A['four'][D][key], B['four'][D][key]
            m = MA[0].nrows()
            for (i, j) in [(0, 0), (m - 1, m - 1), (0, m - 1)]:
                vals = [to_mp(r['four'][D][key][0][i, j]) for _, r in pts]
                c = interp(vals)
                sc = max(abs(v) for v in vals) + mp.mpf(10) ** -40
                for k in range(6):
                    worst_rel = max(worst_rel, abs(to_mp(MA[k][i, j]) - c[k]) * hm ** k / sc)
                F6 = MB[6][i, j]
                F6abs = max(abs(float(F6.upper())), abs(float(F6.lower())))
                c6pts = max(abs(to_mp(r['four'][D][key][6][i, j])) for _, r in pts)
                ok6 = F6abs >= float(c6pts) * (1 - 1e-9) and F6abs >= float(abs(c[6])) * (1 - 1e-6)
                worst_rem = max(worst_rem, 0 if ok6 else 1)
            lines.append(f"   D={D} {key} entries (0,0),(m-1,m-1),(0,m-1): order-6 enclosure contains node/interp c6: {worst_rem == 0}")
    print("\n".join(lines))
    print(f"  (1) max_i |c_i(series) - c_i(interp)| h^i / |f| over all checked functions, i=0..5: {float(worst_rel):.2e};  "
          f"remainder-enclosure failures: {worst_rem}", flush=True)
    # (2) point re-certification with the box's rho
    rho = {D: Fraction(v) for D, v in rec['rho'].items()}
    rho = {int(D): v for D, v in rho.items()}
    rho_full = {D: rho.get(D, Fraction(0)) for D in range(1, 24)}
    pts2 = [xa, xm, xb] if xa > 0 else [xm, xb]
    for xp in pts2:
        o = cc.evaluate(arb(fmpq(xp.numerator, xp.denominator)), G, rho_given=rho_full)
        print(f"  (2) r-space point x={xp} : fails={o['fail']} gamma={o['gamma_lo']:.6f} >= box {rec['gamma_lo']:.6f}: "
              f"{o['gamma_lo'] is not None and o['gamma_lo'] >= rec['gamma_lo']}", flush=True)
    Qlo = -(-xb.denominator // xb.numerator)          # ceil(1/xb)
    Qs = [Qlo, Qlo + 1] + ([int(1 / xa)] if xa > 0 else [10 ** 6, 10 ** 9, 10 ** 15])
    for Q in Qs:
        assert Fraction(1, Q) >= xa and Fraction(1, Q) <= xb
        ok, c4, fails = cc.bspace_check(Q, G, rho_full)
        o = cc.evaluate(arb(fmpq(1, Q)), G, Qint=Q if Q <= 2000 else None, rho_given=rho_full)
        print(f"  (2) integer Q={Q}: b-space (rank-free) ok={ok} fails={fails}; r-space gamma={o['gamma_lo']:.6f} >= box: "
              f"{o['gamma_lo'] is not None and o['gamma_lo'] >= rec['gamma_lo']}", flush=True)
