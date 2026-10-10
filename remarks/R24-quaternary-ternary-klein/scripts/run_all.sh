#!/usr/bin/env bash
# Reproduce the computational checks of this remark and diff them against expected_output/.
# Needs: a C compiler (gcc/cc), python3 with numpy. Single-threaded; total about 1 minute.
# The Lean part is NOT run by default (it needs the full openai/math Lean project, several GB);
# set LEAN_PROJECT=/path/to/openai-math/lean to run it (see lean/README.md).
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CC=${CC:-cc}; PY=${PY:-python3}
mkdir -p out
$CC -O2 -o out/perfect_search independent/perfect_search.c
{
  # Theorem Q (quaternary, length 2m, m odd > 1): expect 0
  for n in 6 10 14; do out/perfect_search 4 $n; done
  # positive controls: quaternary lengths 2, 4, 8 exist
  for n in 2 4 8; do out/perfect_search 4 $n; done
  # Theorem T3 (three-phase, length 3m, 3 does not divide m > 1): expect 0
  for n in 6 12 15; do out/perfect_search 3 $n; done
  # positive controls: three-phase lengths 3, 9 exist
  for n in 3 9; do out/perfect_search 3 $n; done
  # binary arrays on C2 x C2 x C_N, N odd non-square (trivially 0) and the control C2 x C2
  out/perfect_search 2 2 2 3; out/perfect_search 2 2 2 5; out/perfect_search 2 2 2
} > out/perfect_search.txt
$PY independent/klein_mitm.py > out/klein_mitm.txt
for f in perfect_search klein_mitm; do
  if diff -q out/$f.txt expected_output/$f.txt >/dev/null; then echo "OK   $f"; else echo "DIFF $f"; fi
done
if [ -n "${LEAN_PROJECT:-}" ]; then
  cp -r lean/OAI "$LEAN_PROJECT"/
  here=$(pwd)
  if ! ( cd "$LEAN_PROJECT" && lake build \
      OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremQ \
      OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremT1 \
      OAI.LinearAlgebra.CirculantHadamard.Eisenstein.TheoremT3 \
      OAI.LinearAlgebra.CirculantHadamard.RankOne.TheoremF \
    && lake env lean "$here/lean/Checks.lean" ) > out/lean_checks.txt 2>&1; then
    echo "FAIL lean (build or check failed; see out/lean_checks.txt)"; exit 1
  fi
  $PY - out/lean_checks.txt <<'PYEOF'
import re, sys
t = open(sys.argv[1]).read()
blocks = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", t)
std = {"propext", "Classical.choice", "Quot.sound"}
good = [n for n, ax in blocks if {a.strip() for a in ax.split(",")} <= std]
ok = len(blocks) == 9 and len(good) == 9 and "error" not in t
print(("OK   " if ok else "FAIL ") + f"lean: {len(good)}/{len(blocks)} printed theorems use only "
      "[propext, Classical.choice, Quot.sound] (expected 9/9)")
PYEOF
fi
