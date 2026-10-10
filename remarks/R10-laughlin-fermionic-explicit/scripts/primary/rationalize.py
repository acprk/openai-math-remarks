#!/usr/bin/env python3
"""269: turn a float SDP solution (npz with G0..G7) into an exact rational Gram factor F_t (Gamma_t = F_t F_t^T).
usage: rationalize.py in.npz out.json [bits=44]"""
import sys, json
import numpy as np
from fractions import Fraction
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import aengine as ae


def rationalize(G, bits=44):
    F = {}
    for t in ae.TSET:
        A = 0.5 * (G[t] + G[t].T)
        w, V = np.linalg.eigh(A)
        keep = w > 1e-13 * max(1.0, w.max())
        Ft = V[:, keep] * np.sqrt(w[keep])
        F[t] = [[Fraction(int(round(v * 2 ** bits)), 2 ** bits) for v in row] for row in Ft]
    return F


def save(F, fn, meta=None):
    json.dump({'meta': meta or {}, 'F': {str(t): [[str(v) for v in row] for row in F[t]] for t in F}}, open(fn, 'w'))


def load(fn):
    d = json.load(open(fn))
    return {int(t): [[Fraction(v) for v in row] for row in rows] for t, rows in d['F'].items()}


if __name__ == '__main__':
    dat = np.load(sys.argv[1])
    bits = int(sys.argv[3]) if len(sys.argv) > 3 else 44
    G = {t: dat[f'G{t}'] for t in ae.TSET}
    F = rationalize(G, bits)
    save(F, sys.argv[2], {'source': sys.argv[1], 'bits': bits})
    print('saved', sys.argv[2], {t: (len(F[t]), len(F[t][0]) if F[t] and F[t][0] else 0) for t in F})
