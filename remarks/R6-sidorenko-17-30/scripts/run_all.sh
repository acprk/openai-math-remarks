#!/usr/bin/env bash
# Reproduce the outputs of remark R6 and diff against expected_output/.
# Environment (tested 2026-10-09): CPython 3.10.12, networkx 3.4.2, numpy, scipy 1.15.3 (HiGHS),
#   a C compiler, network access for plantri (third-party, Apache-2.0, downloaded + sha256-checked by
#   get_plantri.sh; NOT redistributed).  Single-threaded; total ~2-3 min.
# Not rerun here: numerics/wsearch*.py (heuristic kernel search, needs PyTorch, hours);
#   their recorded logs are expected_output/wsearch*_RECORDED.log.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
mkdir -p out
# 1. main certificate for K_7 (independent of construction/ code): (1)-(3), stable colouring,
#    Aut(H_7), generalised pinning over all homomorphisms, exact C(f)
(cd certificate && $PY check_v7.py) > out/check_v7.txt
(cd certificate && $PY pin_margin.py 3) > out/pin_margin_v7.txt
(cd premise_audit && $PY rank_layer_audit.py) > out/rank_layer_audit.txt
# 2. construction-side checks (primary code)
(cd construction && $PY verify3.py 3)               > out/verify3.txt
(cd construction && $PY pinexact.py 3 cands.txt 1)  > out/pinexact_v7_homs.txt
(cd construction && $PY pinexact.py 3 cands.txt 0)  > out/pinexact_v7_allmaps.txt
for i in 1 4 5; do (cd construction && $PY pinexact.py $i cands.txt 1) > out/pinexact_cand${i}_homs.txt; done
(cd construction && $PY survey.py cands.txt)        > out/survey.txt
(cd construction && $PY fold1.py)                   > out/fold1_v6.txt
# 3. location of H_7 relative to known Sidorenko classes
(cd classes && $PY invariants.py) > out/invariants.txt
(cd classes && $PY stdcheck.py)   > out/stdcheck_v7.txt
# 4. appendix: LP certificates with condition (4) kept, V <= 10, + positive control (paper's complex)
[ -x _build/plantri58/plantri ] || bash get_plantri.sh 58
$PY appendix_lp/r6_verify.py --vmax 10 --no13 > out/r6_verify_vmax10.txt
for f in out/*.txt; do b=$(basename "$f")
  if diff -q <(grep -v '^real\|^user\|^sys' "$f") <(grep -v '^real\|^user\|^sys' "expected_output/$b") >/dev/null; then echo "OK   $b"; else echo "DIFF $b"; fi
done
