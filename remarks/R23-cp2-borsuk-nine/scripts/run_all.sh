#!/usr/bin/env bash
# Reproduce every output of remark X156 (CP^2 projector set: Borsuk number 9) and diff with expected_output/.
# Needs: python3 with numpy (cert_interval.py, tested numpy 1.26 and 2.x) and scipy (numerics/ only).
# independent/cert_exact.py uses only the Python standard library.
# Single-threaded; total about 2 minutes on one core.
# The two certificates are deterministic (IEEE-754 binary64, resp. exact integers).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out
echo "07c6a2f1e5ded3d92038884883f6ff8f271dac5702e06700a72176da5cdadcab  P9.txt" | sha256sum -c --quiet
# primary certificate (outward-rounded interval arithmetic), headline threshold 21/40 and sharp threshold 263/500
$PY cert_interval.py P9.txt 21 40   > out/cert_interval_21_40.txt
$PY cert_interval.py P9.txt 263 500 > out/cert_interval_263_500.txt
# negative control: at 53/100 (> numerical minimum 0.526055) the certificate must FAIL with a witness
set +e; $PY cert_interval.py P9.txt 53 100 > out/control_53_100.txt; rc=$?; set -e
echo "exit code $rc" >> out/control_53_100.txt
# negative control 2: drop the ninth point; eight caps can never work (lower bound), so this must FAIL
set +e; $PY cert_interval.py controls_P8_drop9.txt 21 40 > out/control_drop9.txt; rc=$?; set -e
echo "exit code $rc" >> out/control_drop9.txt
# independent re-implementation (exact integers, Taylor bound, bisection)
$PY independent/cert_exact.py P9.txt 21 40   > out/cert_exact_21_40.txt
$PY independent/cert_exact.py P9.txt 263 500 > out/cert_exact_263_500.txt
# illustrations (not used in the proof)
$PY numerics/minimum_search.py P9.txt > out/minimum_search.txt
$PY numerics/cap_lemma_check.py        > out/cap_lemma_check.txt
for f in cert_interval_21_40 cert_interval_263_500 control_53_100 control_drop9 minimum_search cap_lemma_check; do
  if diff -q out/$f.txt expected_output/$f.txt >/dev/null; then echo "OK   $f"; else echo "DIFF $f"; fi
done
for f in cert_exact_21_40 cert_exact_263_500; do
  if diff -q out/$f.txt independent/expected_output/$f.txt >/dev/null; then echo "OK   independent/$f"; else echo "DIFF independent/$f"; fi
done
