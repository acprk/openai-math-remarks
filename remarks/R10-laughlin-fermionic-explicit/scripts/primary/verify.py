#!/usr/bin/env python3
"""269: aggregate all certificate outputs, check exact coverage, print the final table (expected_output.txt)."""
import json, glob
from fractions import Fraction

import os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out') + '/'
ok = True

# ---- integer fluxes 24..299 (r-space certificate + rank-free b-space recheck)
Qr = {}
for f in sorted(glob.glob(D + 'certQ_[ABC]*.jsonl')):
    for l in open(f):
        r = json.loads(l); Qr[r['Q']] = r
Qb = {json.loads(l)['Q']: json.loads(l) for l in open(D + 'bspaceQ.jsonl')}
miss = [Q for Q in range(24, 300) if Q not in Qr or Qr[Q]['fail'] or Qr[Q]['gamma_lo'] is None]
missb = [Q for Q in range(24, 300) if Q not in Qb or not Qb[Q]['ok']]
print(f"[int Q 24..299] r-space certified: {276-len(miss)}/276  missing/failed: {miss}")
print(f"[int Q 24..299] b-space (rank-free) recheck ok: {276-len(missb)}/276  failed: {missb}")
ok &= not miss and not missb
segs = {}
for Q, r in Qr.items():
    segs.setdefault(r['gamma_file'], []).append((min(r['gamma_lo'], Qb[Q]['gamma_lo_bspace']), Q))
for g, L in sorted(segs.items(), key=lambda kv: min(q for _, q in kv[1])):
    Qs = [q for _, q in L]
    m = min(L)
    print(f"   segment Q in [{min(Qs)},{max(Qs)}]  Gamma={g.split('/')[-1]}  min gamma_lo={m[0]:.6f} at Q={m[1]}")
m3 = min((min(Qr[Q]['gamma_lo'], Qb[Q]['gamma_lo_bspace']), Q) for Q in range(24, 300, 3))
print(f"   multiples of 3 (Laughlin flux, N=9..100): min gamma_lo={m3[0]:.6f} at Q={m3[1]} (N={m3[1]//3+1})")

# ---- boxes covering x in [0, 1/300]
boxes = []
for f in sorted(glob.glob(D + 'certX_*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        if 'FAIL' in r:
            ok = False; print('FAILED BOX', r)
        else:
            boxes.append((Fraction(r['xa']), Fraction(r['xb']), r['gamma_lo']))
boxes.sort()
cont = boxes[0][0] == 0 and boxes[-1][1] == Fraction(1, 300) and all(boxes[i][1] == boxes[i + 1][0] for i in range(len(boxes) - 1))
gb = min(b[2] for b in boxes)
wb = min(boxes, key=lambda b: b[2])
print(f"[Q>=300] {len(boxes)} boxes, exact contiguous cover of [0,1/300]: {cont};  min gamma_lo={gb:.6f} on x in [{float(wb[0]):.3e},{float(wb[1]):.3e}]")
print(f"   box at x=0 (Q -> infinity): gamma_lo={boxes[0][2]:.6f};  box at x=1/300: gamma_lo={boxes[-1][2]:.6f}")
ok &= cont

# ---- ED N=2..8
ed = json.load(open(D + 'ed_cert.json'))
for r in ed:
    print(f"[ED] N={r['N']} Q={r['Q']} dim(Lz=0)={r['dim']} gap_float={r['gap_float']:.6f} certified g_N={float(Fraction(r['g_cert'])):.4f} "
          f"(Cholesky err bound {r['chol_err_bound']:.1e} < tau {r['tau']:.1e}: {r['ok']})")
    ok &= r['ok']
# ---- independent cross-check logs
xl = [l for f in sorted(glob.glob(D + 'crosscheck_*.log')) for l in open(f) if l.startswith('Q=')]
nv = sum('VIOLATION' in l for l in xl)
print(f"[cross-check] {len(xl)} (Q,N) sectors checked (Jordan-Wigner), violations: {nv}")
ok &= (nv == 0 and len(xl) >= 20)
gall = min([float(Fraction(r['g_cert'])) for r in ed] + [m3[0], gb])
gfock = min(min(min(Qr[Q]['gamma_lo'], Qb[Q]['gamma_lo_bspace']) for Q in range(24, 300)), gb)
print(f"\nTheorem A (Fock): for every integer Q>=24, H_Q^2 >= gamma H_Q on F_Q with gamma = {gfock:.6f} (> 1/25: {gfock > 0.04})")
print(f"Theorem B (Laughlin, all N>=2): H_(N,3(N-1)) >= {gall:.6f} (I - P_L)  (> 1/25: {gall > 0.04})")
print("ALL CHECKS PASSED" if ok and gall > 0.04 else "SOME CHECK FAILED")
