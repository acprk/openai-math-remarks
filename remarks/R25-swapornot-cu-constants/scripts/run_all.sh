#!/usr/bin/env bash
# Reproduce every output of this remark and diff against expected_output/.
# Needs: SageMath (tested 9.5; used only by the primary cu_bound.py, for MPFI intervals RIF) with sympy,
#        python3 with numpy (>= 1.24) and mpmath (>= 1.3) for the other scripts.
# Threads: exact_tables.py uses NPROC=3 worker processes (~30 s); everything else is single-threaded.
# Total wall time ~8 min (primary cu_bound ~5 min, independent indep_Cu ~1.5 min). Peak RAM ~1 GB.
# All scripts are deterministic (exact integers, interval arithmetic, IEEE doubles with fixed seeds).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
SAGEPY=${SAGEPY:-"sage -python"}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NPROC=${NPROC:-3}
mkdir -p out
# ---- primary computation (Theorem 1 of the note; method of lemma_Cu: hook formula, ballot peeling, MPFI)
$SAGEPY exact_tables.py 55 > out/exact_tables.txt             # writes exact_tables_55.pkl (~60 MB)
$PY check_lemma.py > out/check_lemma.txt
{ $SAGEPY cu_bound.py 1    55 120 400
  $SAGEPY cu_bound.py 1/2  55 120 400
  $SAGEPY cu_bound.py 3/10 55 120 400
  $SAGEPY cu_bound.py 1/5  55 120 400
  $SAGEPY cu_bound.py 3/20 55 200 1000; } > out/cu_bound.txt
$PY eps_n_exponent.py > out/eps_n_exponent.txt
# ---- independent re-verification (separate code, different lemmas and arithmetic; see VERIFICATION.md)
$PY -u independent/indep_Cu.py | sed -E 's/ +\[[0-9]+s\]$//' > out/indep_Cu.txt
$PY independent/indep_eps.py > out/indep_eps.txt
for f in exact_tables check_lemma cu_bound eps_n_exponent; do
  if diff -q out/$f.txt expected_output/$f.txt >/dev/null; then echo "OK   $f"; else echo "DIFF $f"; fi
done
for f in indep_Cu indep_eps; do
  if diff -q out/$f.txt independent/expected_output/$f.txt >/dev/null; then echo "OK   independent/$f"; else echo "DIFF independent/$f"; fi
done
