#!/usr/bin/env bash
# Reproduce the computations of this remark and diff against expected_output/.
# Environment (tested 2026-10-09): CPython 3.10.12 + numpy 2.2.6 (independent checks);
#   SageMath 9.5 (primary checks; set SKIP_SAGE=1 to skip them).  Single-threaded throughout.
# Wall time: primary ~1 min, independent CB7/down-set check ~10 s, independent GHW check ~3-7 min.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p out && rm -f out/*.txt
# 1. primary computation (Sage; brute force over supports), seeds as in the original run
if [ "${SKIP_SAGE:-0}" != "1" ]; then
  (cd primary && sage ghw_check.sage 1)  > out/ghw_check_seed1.txt
  (cd primary && sage ghw_check2.sage 2) > out/ghw_check2_seed2.txt
fi
# 2. independent re-implementation (pure Python; column matroid, exhaustive down-sets)
$PY independent/cb7_and_reduction_check.py 7 > out/cb7_and_reduction_check.txt
$PY independent/ghw_ci_indep.py 2026        > out/ghw_ci_indep.txt
for f in out/*.txt; do b=$(basename "$f")
  if diff -q "$f" "expected_output/$b" >/dev/null; then echo "OK   $b"; else echo "DIFF $b"; fi
done
