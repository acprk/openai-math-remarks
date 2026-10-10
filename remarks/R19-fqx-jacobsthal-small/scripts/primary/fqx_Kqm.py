#!/usr/bin/python3
"""Exact K_q(m) = min number of distinct monic irreducibles P_1..P_k in F_q[x] such that, for suitable
residues r_i, every g with deg g < m satisfies g = r_i (mod P_i) for some i.
(Irreducibles of degree >= m cover exactly one arbitrary point each -> modelled as 'wildcards' z_g.)
ILP solved by OR-tools CP-SAT (3 workers).  Usage: python3 fqx_Kqm.py q m [timelimit]"""
import sys, time
from ortools.sat.python import cp_model
from fqx_jacobsthal import irreducibles, residue_table
def K(q, m, tl=600):
    primes = [P for d in range(1, m) for P in irreducibles(q, d)]
    mdl = cp_model.CpModel(); y = {}; obj = []
    for i, P in enumerate(primes):
        lits = [mdl.NewBoolVar("") for _ in range(q ** (len(P) - 1))]
        mdl.AddAtMostOne(lits); obj += lits
        for c, l in enumerate(lits): y[(i, c)] = l
    tabs = [residue_table(P, q, m) for P in primes]
    for g in range(q ** m):
        z = mdl.NewBoolVar(""); obj.append(z)
        mdl.AddBoolOr([y[(i, tabs[i][g])] for i in range(len(primes))] + [z])
    mdl.Minimize(sum(obj))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 3
    st = s.Solve(mdl)
    used = {}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        for (i, c), l in y.items():
            if s.Value(l): used[len(primes[i]) - 1] = used.get(len(primes[i]) - 1, 0) + 1
    return s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound(), used
if __name__ == "__main__":
    q, m = int(sys.argv[1]), int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    t = time.time(); st, val, bd, used = K(q, m, tl)
    wild = int(val) - sum(used.values()) if val == val else None
    print(f"q={q} m={m}: K={val:.0f} [{st}, bound {bd:.0f}]  m(q-1)+1={m*(q-1)+1}  low-degree primes used by degree={dict(sorted(used.items()))} wildcards={wild}  ({time.time()-t:.1f}s)", flush=True)
