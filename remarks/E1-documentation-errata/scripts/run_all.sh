#!/usr/bin/env bash
# run_all.sh -- reproduce every check in E1-documentation-errata and compare with expected_output/.
#
# Usage:
#   git clone https://github.com/openai/math openai-math
#   OPENAI_MATH=$PWD/openai-math ./run_all.sh            # full run, ~15 min single thread
#   OPENAI_MATH=$PWD/openai-math SKIP_LONG=1 ./run_all.sh # skip the two long re-runs (~10 s)
# Optional: COMMIT (default fd4aeeb2e = openai/math main on 2026-10-09), PY (default python3).
# Requirements: bash, git, python3 (tested 3.10.12) with numpy (tested 2.2.6),
# mpmath (tested 1.3.0) and sympy (tested 1.14.0). Single-threaded throughout.
# The clone's working tree is never modified: files are read with `git show` / `git archive`.
set -euo pipefail
: "${OPENAI_MATH:?set OPENAI_MATH to a clone of github.com/openai/math}"
export OPENAI_MATH
export COMMIT="${COMMIT:-fd4aeeb2e}"
PY="${PY:-python3}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
HERE="$(cd "$(dirname "$0")" && pwd)"
EXP="$HERE/expected_output"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT
fail=0
cmp_out() {  # cmp_out NAME : compare $OUT/NAME with expected_output/NAME
  if diff -u "$EXP/$1" "$OUT/$1" > "$OUT/$1.diff"; then echo "OK    $1"
  else echo "DIFF  $1"; cat "$OUT/$1.diff"; fail=1; fi
}

cd "$HERE"
# 1. Every cited location / missing file at COMMIT (main check)
./check_upstream.sh > "$OUT/check_upstream.txt";                      cmp_out check_upstream.txt
# 2. Independent path check (cat-file over the whole history)
$PY independent/indep_check_paths.py "$OPENAI_MATH" "$COMMIT" > "$OUT/indep_check_paths.txt"; cmp_out indep_check_paths.txt
# 3. Family 238: the five rational inequalities of support/network_bounds.py
$PY 238_network_bounds_recheck.py > "$OUT/238_network_bounds_recheck.txt"; cmp_out 238_network_bounds_recheck.txt
$PY independent/indep_238.py      > "$OUT/indep_238.txt";                  cmp_out indep_238.txt
# 4. Family 318: Chevalley-Eilenberg ranks of support/r-ce/ce_certificate.py
$PY 318_ce_ranks.py               > "$OUT/318_ce_ranks.txt";               cmp_out 318_ce_ranks.txt
$PY independent/indep_318.py      > "$OUT/indep_318.txt";                  cmp_out indep_318.txt
# 5. Family 140: numerical illustration of the sharper constant in Lemma s:mixture
$PY 140_mixture_constant.py       > "$OUT/140_mixture_constant.txt";       cmp_out 140_mixture_constant.txt
# 6. Family 222: extract the Appendix A listing, check the published SHA-256
WORK="$OUT/work"; mkdir -p "$WORK"; ( cd "$WORK" && $PY "$HERE/extract_222_certificate.py" "$OPENAI_MATH" "$COMMIT" ) > "$OUT/222_extract.txt"
cmp_out 222_extract.txt

if [ -z "${SKIP_LONG:-}" ]; then
  # 7. Family 222: run the extracted certificate (~8 min); stdout = the 24 interval checkpoints
  ( cd "$WORK" && $PY 222_original_certificate.py ) > "$OUT/222_certificate_run.txt"
  cmp_out 222_certificate_run.txt
  # 8. Family 266-B: run the shipped verification/run.py (~6 min) from a git archive of the manuscript
  D=preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026
  git -C "$OPENAI_MATH" archive "$COMMIT" "$D/verification" | tar -x -C "$WORK"
  ( cd "$WORK/$D" && $PY verification/run.py ) > "$OUT/266B_verification_run.txt"
  cmp_out 266B_verification_run.txt
else
  echo "SKIP  222_certificate_run.txt, 266B_verification_run.txt (SKIP_LONG set)"
fi
[ $fail -eq 0 ] && echo "ALL OK" || { echo "SOME OUTPUTS DIFFER"; exit 1; }
