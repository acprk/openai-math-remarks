# Exact enumeration (integers only) of all partitions lambda |- k, k<=M0.
# Output: for each k, a dict (lambda_1, ell(lambda), f_lambda) -> multiplicity, f_lambda exact via hook formula.
# Usage: python3 exact_tables.py M0  -> writes exact_tables_M0.pkl
import sys, pickle, os
from math import factorial
from multiprocessing import Pool
from collections import Counter
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest
def table(k):
    fk = factorial(k); C = Counter()
    for lam in partitions(k):
        conj = [0]*lam[0]
        for r in lam:
            for j in range(r): conj[j] += 1
        p = 1
        for i, r in enumerate(lam):
            for j in range(r): p *= (r - j) + (conj[j] - i) - 1
        f, rem = divmod(fk, p); assert rem == 0
        C[(lam[0], len(lam), f)] += 1
    return k, dict(C)
if __name__ == "__main__":
    M0 = int(sys.argv[1])
    with Pool(int(os.environ.get("NPROC", "3"))) as pool:
        res = dict(pool.map(table, range(M0, 0, -1), chunksize=1))
    from sympy.functions.combinatorial.numbers import partition as npartitions
    for k in res: assert sum(res[k].values()) == npartitions(k)
    pickle.dump(res, open(f"exact_tables_{M0}.pkl", "wb"))
    print("done", M0, sum(len(v) for v in res.values()), "distinct (r,c,f) classes")
