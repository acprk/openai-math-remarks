#!/usr/bin/env python3
"""Primary computation: exact K_q(m) by CP-SAT (OR-tools), using the model of fqx_Kqm.py.

Model (see README): every monic irreducible P with deg P < m chooses at most one residue class
mod P; every g in V_m = {deg g < m} must lie in a chosen class or be paid for by one
"wildcard" (an irreducible of degree >= m, which can cover exactly one arbitrary point).
K_q(m) = minimum number of chosen classes + wildcards.

Prints only the optimal value and the solver status (deterministic).
Environment: python3 >= 3.10, ortools 9.x, python-sat (imported by fqx_jacobsthal). Single worker.
"""
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
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 1
    st = s.Solve(mdl)
    return s.StatusName(st), int(round(s.ObjectiveValue())), int(round(s.BestObjectiveBound()))

if __name__ == "__main__":
    for q, m in [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2)]:
        st, val, bd = K(q, m)
        print(f"q={q} m={m}: K_q(m)={val} status={st} bound={bd}  (q, 2q-1)=({q},{2*q-1})")
