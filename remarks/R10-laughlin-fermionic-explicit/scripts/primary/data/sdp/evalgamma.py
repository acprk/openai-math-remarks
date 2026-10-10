#!/usr/bin/env python3
"""269: a-posteriori evaluation of a candidate Gamma at flux Q (independent of the SDP solver's claims).
Gamma_t is projected onto the PSD cone (negative eigenvalues dropped) - so K is a genuine sum of squares.
gamma(Q) = 1 - eta - tail_Q - sum_z max(0,-slack_z)(z+1)d_3z/d_a - sum_D n_D rho_D,
rho_D = max(0, -lambda_min(U^T S^T C S U))  (scalar allowance rho*I_b, cost n_D rho).
usage: python evalgamma.py file.npz Qeval [pmax M tmax]"""
import sys, numpy as np
import os; sys.path.insert(0, os.environ.get('LAUGHLIN_STEP1_DIR', '.'))  # directory holding engine.py from the step-1 search (provenance only)
from engine import *


def evaluate(G, Q, sup, verbose=False):
    eng = Engine(Q, sup)
    da = eng.da
    eta = sum(G[t][0, 0] for t in sup.tset)
    cost3 = 0.0; w3 = (1e9, None)
    for z in eng.zlist:
        if z in (0, 1, 3):
            continue
        ez = sum((eng.C3[z][t] * G[t]).sum() for t in sup.tset)
        sl = eng.qz(z) - da / eng.d3(z) * ez
        w3 = min(w3, (sl, z))
        if sl < 0:
            cost3 += -sl * (z + 1) * eng.d3(z) / da
    tail = eng.delta_tail()
    cost4 = 0.0
    for D in range(1, sup.Dmax + 1):
        R = eng.Rset[D]
        chi = np.zeros((len(R), len(R))); chi[0, 0] = (4 * Q - 1 - 2 * D) / da
        C = chi - eng.E_of(D, G)
        bl, Smat, vecs = eng.bspace(D)
        U = range_basis(gram(vecs))
        if U.shape[1] == 0:
            continue
        lam = np.linalg.eigvalsh(U.T @ Smat.T @ C @ Smat @ U).min()
        cost4 += max(0.0, -lam) * eng.nD(D)
    g = 1 - eta - tail - cost3 - cost4
    if verbose:
        print(f"  Q={Q}: eta={eta:.6f} tail={tail:.6f} min3slack={w3[0]:+.3e}@z={w3[1]} cost3={cost3:.3e} cost4={cost4:.3e} -> gamma={g:.6f}", flush=True)
    return g


def load(fn, sup):
    dat = np.load(fn)
    G = {}
    for t in sup.tset:
        A = 0.5 * (dat[f'G{t}'] + dat[f'G{t}'].T)
        w, V = np.linalg.eigh(A)
        G[t] = (V * np.clip(w, 0, None)) @ V.T
    return G


if __name__ == '__main__':
    fn = sys.argv[1]
    sup = Support(*map(int, sys.argv[3:6]), ) if len(sys.argv) > 5 else Support()
    if len(sys.argv) > 5:
        sup = Support(int(sys.argv[3]), int(sys.argv[4]), range(int(sys.argv[5]) + 1))
    G = load(fn, sup)
    for Qs in sys.argv[2].split(','):
        evaluate(G, int(float(Qs)), sup, verbose=True)
