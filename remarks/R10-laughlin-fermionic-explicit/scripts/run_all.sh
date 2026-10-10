#!/bin/bash
# R10 (OpenAI family 269): explicit Fock-space inequality H_Q^2 >= 0.05947 H_Q for all Q >= 24 and the
# Laughlin gap H_{N,3(N-1)} >= 0.05947 (I - P_L) for all N >= 2.
#
#   bash run_all.sh            = bash run_all.sh verify   (seconds) SHA256 check of data + shipped certificate outputs,
#                                re-aggregation (coverage of Q=24..299 and of x=1/Q in [0,1/300], ED table), diff vs expected
#   bash run_all.sh full       (~65 min wall, 4 single-threaded processes) recompute the whole primary certificate
#                                (certQ 3.5 min, bspaceQ 8 min, certX 43 min, ed_cert 3.3 min, crosscheck ~15 min), then verify
#   bash run_all.sh independent (~45 min wall, <= 4 single-threaded processes) adversarial re-derivations (independent/)
#
# Environment: Python 3.12 with python-flint 0.9.0, numpy 2.5, scipy 1.18, mpmath 1.4 (set PY=/path/to/python).
# Tested on Linux x86_64, 2026-10-09.  Uses at most 4 processes; BLAS forced single-threaded.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
export PY=${PY:-python3} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MODE=${1:-verify}
cd "$HERE"
sha256sum -c SHA256SUMS --quiet && echo "SHA256SUMS: code, data and expected outputs match"
if [ "$MODE" = "verify" ]; then
  sha256sum -c SHA256SUMS.outputs --quiet && echo "SHA256SUMS.outputs: shipped certificate outputs match"
fi
if [ "$MODE" = "full" ]; then
  cd "$HERE/primary"
  $PY certQ.py 24 48   data/gamma_24_36_48.json        out/certQ_A.jsonl  > out/certQ_A.log &
  $PY certQ.py 49 72   data/gamma_48_60_72.json        out/certQ_B.jsonl  > out/certQ_B.log &
  $PY certQ.py 73 186  data/gamma_72_100_150_300.json  out/certQ_C1.jsonl > out/certQ_C1.log &
  $PY certQ.py 187 299 data/gamma_72_100_150_300.json  out/certQ_C2.jsonl > out/certQ_C2.log &
  wait
  cat out/certQ_A.jsonl out/certQ_B.jsonl out/certQ_C1.jsonl out/certQ_C2.jsonl > out/certQ_all.jsonl
  $PY bspaceQ.py out/certQ_all.jsonl out/bspaceQ.jsonl > out/bspaceQ.log &
  for i in 0 1 2; do $PY certX.py data/gamma_300_1000_100000000.json out/certX_$i.jsonl $i 1200 $((i+1)) 1200 0.04 21 > out/certX_$i.log & done
  wait
  $PY certX.py data/gamma_300_1000_100000000.json out/certX_3.jsonl 3 1200 4 1200 0.04 21 > out/certX_3.log &
  $PY ed_cert.py 8 out/ed_cert.json > out/ed_cert.log &
  wait
  ./crosscheck_all.sh
fi
if [ "$MODE" = "independent" ]; then
  cd "$HERE/independent"
  mkdir -p out
  ( for q in 24 27; do $PY check1_definitions.py $q; done > out/check1_definitions.txt
    $PY check2_rank.py > out/check2_rank.txt ) &
  ( $PY check1e_sector.py 24 5; $PY check1e_sector.py 24 3 flip ) > out/check1e_A.txt &
  ( $PY check1e_sector.py 24 3,4; $PY check1e_sector.py 27 3,4; $PY check1e_sector.py 29 3,4 ) > out/check1e_B.txt &
  ( $PY check1e_sector.py 49 3; $PY check1e_sector.py 66 3; $PY check1e_sector.py 73 3; $PY check1e_sector.py 100 3 ) > out/check1e_C.txt &
  wait
  ( $PY check3_boxes.py 0 > out/check3_box0.txt ) &
  ( $PY check3_boxes.py 1,60 > out/check3_box1_60.txt ) &
  ( $PY check3_boxes.py 120,167 > out/check3_box120_167.txt ) &
  ( for n in 2 3 4 5 6 7 8; do $PY check5_ed.py $n; done > out/check5_ed.txt ) &
  wait
  for f in out/*.txt; do diff -q "$f" "$HERE/expected_output/independent/$(basename "$f")" >/dev/null \
      && echo "MATCH  $f" || echo "DIFF   $f  (only round-off digits of float residuals may differ; every line must end OK/True, except the deliberate negative control VIOLATION in check1e_A)"; done
  exit 0
fi
cd "$HERE/primary"
$PY verify.py > out/verify.txt
cat out/verify.txt
diff <(grep -v '^$' out/verify.txt) <(grep -v '^$' "$HERE/expected_output/verify.txt") && echo "MATCHES expected_output/verify.txt"
