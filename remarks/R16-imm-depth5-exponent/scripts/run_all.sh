#!/usr/bin/env bash
# Re-run every script and compare with expected_output/ (python3 + mpmath + numpy + scipy).
set -u
cd "$(dirname "$0")"
PY=${PY:-python3}
status=0
for f in main/explicit_bound main/hbound_check main/framework_limit main/imm_wordsum \
         independent/v_params_H independent/v_imm independent/v_threshold independent/v_framework; do
  b=$(basename "$f")
  if "$PY" "$f.py" | diff -q - "expected_output/$b.txt" >/dev/null; then
    echo "OK   $f"
  else
    echo "DIFF $f"; status=1
  fi
done
exit $status
