#!/usr/bin/env bash
# Reproduce every output of this remark and diff against the expected outputs.
# Needs: python3 (>=3.10), mpmath 1.3, sympy (>=1.10), numpy (2.2.6 for byte-identical primary output).
# Single-threaded; total < 1 minute.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p out
$PY primary/qary_refined_boolean.py > out/qary_refined_boolean.txt
$PY primary/hellinger_counter.py    > out/hellinger_counter.txt
$PY independent/indep_scope.py      > out/indep_scope.txt
chk(){ if diff -q "$1" "$2" >/dev/null; then echo "OK   $3"; else echo "DIFF $3"; fi; }
chk out/qary_refined_boolean.txt expected_output/qary_refined_boolean.txt primary/qary_refined_boolean
chk out/hellinger_counter.txt    expected_output/hellinger_counter.txt    primary/hellinger_counter
chk out/indep_scope.txt independent/expected_output/indep_scope.txt independent/indep_scope
