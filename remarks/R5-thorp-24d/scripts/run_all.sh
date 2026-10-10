#!/usr/bin/env bash
# Reproduce every output of remark R5 and diff against expected_output/.
# Needs: python3 with mpmath (>=1.3) and numpy (>=1.24). Single-threaded; total ~4-5 min.
# independent/indep_constants.py is the adversarial reviewer's own re-derivation (does not import
# constants_certificate.py); it is diffed against independent/expected_output/.
# energy_mc.py is a seeded Monte Carlo: its output is bit-identical only for the same NumPy version
# (tested 2.2.6); the other three scripts are deterministic.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out
$PY constants_certificate.py > out/constants_certificate.txt
$PY energy_lemma_check.py    > out/energy_lemma_check.txt
$PY small_d_exact.py         > out/small_d_exact.txt
$PY energy_mc.py 8 4 4000    > out/energy_mc_d8_b4.txt
$PY independent/indep_constants.py > out/indep_constants.txt
for f in constants_certificate energy_lemma_check small_d_exact energy_mc_d8_b4; do
  if diff -q out/$f.txt expected_output/$f.txt >/dev/null; then echo "OK   $f"; else echo "DIFF $f"; fi
done
if diff -q out/indep_constants.txt independent/expected_output/indep_constants.txt >/dev/null; then
  echo "OK   independent/indep_constants"; else echo "DIFF independent/indep_constants"; fi
