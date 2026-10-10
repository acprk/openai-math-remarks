#!/usr/bin/env bash
# Reproduce every output of this remark and diff against the expected outputs.
# Needs: python3 (>=3.10) with ortools 9.x and python-sat (primary only); SageMath 9.5+ (construction check).
# All runs single-threaded; total < 1 minute.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}; SAGE=${SAGE:-sage}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p out
(cd primary && $PY run_table.py) > out/primary_table.txt
$PY independent/indep_Kqm.py > out/indep_Kqm.txt
$SAGE independent/construction_check.sage > out/construction_check.txt 2>&1
rm -f independent/construction_check.sage.py
chk(){ if diff -q "$1" "$2" >/dev/null; then echo "OK   $3"; else echo "DIFF $3"; fi; }
chk out/primary_table.txt expected_output/primary_table.txt primary/run_table
chk out/indep_Kqm.txt independent/expected_output/indep_Kqm.txt independent/indep_Kqm
chk out/construction_check.txt independent/expected_output/construction_check.txt independent/construction_check
