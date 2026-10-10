#!/usr/bin/env bash
# Reproduce every output of this remark and diff against expected_output/.
# Needs: python3 with mpmath (>= 1.3), numpy (>= 1.24) and python-flint (>= 0.6, for the independent check).
# Single-threaded; total ~6 min (theorem_c ~2 min, independent ~3.5 min).
# All scripts are deterministic (interval / ball arithmetic and exact enumeration).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out
$PY chain_check.py                    > out/chain_check.txt
$PY hardcore_certificate.py main      > out/hardcore_main.txt
$PY hardcore_certificate.py dense     > out/hardcore_dense.txt
$PY hardcore_certificate.py paper     > out/hardcore_paper.txt
$PY -u theorem_c_exhaustive.py        > out/theorem_c_exhaustive.txt
$PY independent/indep_check.py        > out/indep_check.txt
for f in chain_check hardcore_main hardcore_dense hardcore_paper theorem_c_exhaustive; do
  if diff -q out/$f.txt expected_output/$f.txt >/dev/null; then echo "OK   $f"; else echo "DIFF $f"; fi
done
if diff -q out/indep_check.txt independent/expected_output/indep_check.txt >/dev/null; then
  echo "OK   independent/indep_check"; else echo "DIFF independent/indep_check"; fi
