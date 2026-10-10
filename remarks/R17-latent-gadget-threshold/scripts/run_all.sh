#!/usr/bin/env bash
# Reproduce every output of this remark and diff against the expected outputs.
# Needs: python3 (>=3.10), numpy (2.2.6 for byte-identical primary output), mpmath 1.3.
# Single-threaded. Runtime: primary ~4-5 min (dominated by the death-chain DP), independent ~1-2 min.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out
$PY primary/gadget_bounds.py         > out/gadget_bounds.txt
$PY independent/indep_gadget.py      > out/indep_gadget.txt
chk(){ if diff -q "$1" "$2" >/dev/null; then echo "OK   $3"; else echo "DIFF $3"; fi; }
chk out/gadget_bounds.txt expected_output/gadget_bounds.txt primary/gadget_bounds
chk out/indep_gadget.txt independent/expected_output/indep_gadget.txt independent/indep_gadget
