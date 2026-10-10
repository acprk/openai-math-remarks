#!/usr/bin/env bash
# check_upstream.sh -- mechanical check of every location cited in README.md.
#
# Usage:  OPENAI_MATH=/path/to/clone/of/openai/math  COMMIT=fd4aeeb2e  ./check_upstream.sh
# Reads the repository only through `git show COMMIT:path` and `git ls-tree`,
# so the working tree of the clone is never touched.
# For each item it prints the quoted line(s) at the stated line numbers and,
# for each referenced-but-missing file, the number of paths in the whole tree
# whose basename matches (expected: 0).
# Environment: bash, git >= 2.25, GNU coreutils. Runtime: a few seconds.
set -euo pipefail
: "${OPENAI_MATH:?set OPENAI_MATH to a clone of github.com/openai/math}"
COMMIT="${COMMIT:-fd4aeeb2e}"
G="git -C $OPENAI_MATH"
FULL=$($G rev-parse "$COMMIT")
echo "openai/math commit: $FULL"
TREE=$(mktemp); trap 'rm -f "$TREE"' EXIT
$G ls-tree -r --name-only "$FULL" > "$TREE"

P=preprints
line() {  # line FILE A B : print lines A..B of FILE at COMMIT, prefixed by line number
  $G show "$FULL:$1" | awk -v a="$2" -v b="$3" 'NR>=a && NR<=b {printf "  %s:%d: %s\n", F, NR, $0}' F="$1"
}
absent() { # absent BASENAME : count tree paths whose basename equals BASENAME (and list them)
  local n; n=$(awk -F/ -v b="$1" '$NF==b' "$TREE" | wc -l)
  printf "  files named %-28s in tree: %s\n" "$1" "$n"
  awk -F/ -v b="$1" '$NF==b {print "    (elsewhere: " $0 ")"}' "$TREE"
}
lsdir() { # lsdir DIR : list top level of a manuscript directory
  printf "  contents of %s/: " "$(basename "$1")"
  $G ls-tree --name-only "$FULL:$1" | tr '\n' ' '; echo
}

echo; echo "== 266-A  N(6)=3 verification archive =="
D=$P/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026
lsdir "$D"
line $D/build/sections/introduction.tex 154 157
line $D/build/sections/execution.tex 16 20
line $D/build/sections/execution.tex 119 119
line $D/build/sections/execution.tex 381 381
line $D/build/sections/execution.tex 390 390
for f in seed.cpp covout.cpp ub.cpp patchx.cpp atlas.cpp stageB.cpp stageC.cpp \
         extraction.json proof_run.py REPRODUCE.md run.json status.json build-manifest.json; do absent "$f"; done
printf "  paths under %s/verification: %s\n" "$(basename $D)" "$(command grep -c "^$D/verification/" "$TREE" || true)"
line lean/docs/266.md 11 13

echo; echo "== 266-B  Exact-Fourier listing path =="
D=$P/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026
lsdir "$D"
line $D/build/paper.tex 235 235
line $D/build/paper.tex 273 273
line $D/build/sections/04-mixed.tex 217 218
printf "  %s/verification/verify.py present: %s\n" "$(basename $D)" "$(command grep -cx "$D/verification/verify.py" "$TREE" || true)"
printf "  any path */certificate/verify.py in tree: %s\n" "$(command grep -c '/certificate/verify\.py$' "$TREE" || true)"

echo; echo "== 095  verify_certificate.py =="
D=$P/Hyperbolicity-Cones-Without-Semidefinite-Lifts-October-5-2026
lsdir "$D"
line $D/build/source/appendix-certificate.tex 536 546
absent verify_certificate.py
line lean/docs/095.md 1 5

echo; echo "== 318  support/ tree =="
D=$P/A-rational-obstruction-to-strong-chromatic-splitting-at-height-three-September-25-2026
lsdir "$D"
line $D/SUPPORT-COVERAGE.md 7 7
line $D/SUPPORT-COVERAGE.md 15 15
line $D/build/sections/08-elementary-background.tex 58 59
for f in ce_certificate.py ce_certificate.json ce_certificate.txt BUILDING.md SUPPORT-COVERAGE.json; do absent "$f"; done

echo; echo "== 222  certificate/original_certificate.py =="
D=$P/Microscopic-jamming-in-the-negative-spherical-perceptron-September-24-2026
lsdir "$D"
line $D/build/main.tex 9466 9476
line $D/build/main.tex 9481 9482
absent original_certificate.py

echo; echo "== 238  support/network_bounds.py =="
D=$P/Routing-densities-and-representation-contraction-for-Thorp-sweeps-September-26-2026
lsdir "$D"
$G show "$FULL:$D/build/sections/04-certificate.tex" | awk 'NR==294' | command grep -o 'The accompanying script.*' | sed "s|^|  $D/build/sections/04-certificate.tex:294: ...|"
absent network_bounds.py

echo; echo "== 227  numerics-*.py =="
D=$P/Cutoff-throughout-the-high-temperature-Sherrington-Kirkpatrick-phase-September-24-2026
lsdir "$D"
line $D/build/appendices/scalar-certificate.tex 291 294
absent numerics-interval.py
absent numerics-cross-errors.py

echo; echo "== 107  lean/docs/107.md =="
line lean/docs/107.md 5 6
line lean/docs/107.md 10 10
printf "  links to Matrix-Multiplication-Nine-Fourths in lean/docs/107.md: %s\n" "$($G show "$FULL:lean/docs/107.md" | command grep -c 'Nine-Fourths' || true)"
printf "  'Lean' mentions in Nine-Fourths README.md: %s\n" "$($G show "$FULL:$P/Matrix-Multiplication-Nine-Fourths-October-2-2026/README.md" | command grep -ci 'lean' || true)"
$G show "$FULL:lean/ComparatorChallenges/MatrixMultiplication.lean" | awk '/theorem complex_omega_le_nine_quarters/ {printf "  lean/ComparatorChallenges/MatrixMultiplication.lean:%d: %s\n", NR, $0}'

echo; echo "== 026  lean/docs/026.md =="
line lean/docs/026.md 9 9
line lean/docs/026.md 15 15
$G show "$FULL:lean/ComparatorChallenges/PrimeGaps.lean" | awk '/theorem large_gaps_and_ratio_density/ {s=NR} s && NR>=s && NR<=s+5 {printf "  lean/ComparatorChallenges/PrimeGaps.lean:%d: %s\n", NR, $0}'

echo; echo "== 099  lean/docs/099.md =="
line lean/docs/099.md 15 19

echo; echo "== 103  environment.tex =="
line $P/Exact-Derandomization-of-Logarithmic-Space-L-equals-RL-equals-BPL-September-23-2026/build/sections/environment.tex 90 94

echo; echo "== 142  CONTENTS.md =="
$G show "$FULL:CONTENTS.md" | awk '/^\*\*142\. / {printf "  CONTENTS.md:%d: ...%s\n", NR, substr($0, length($0)-90)}'
D=$P/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026
line $D/build/main.tex 56 57
line $D/build/sections/10-analytic.tex 9 9
printf "  lean/docs/029.md present: %s   lean/docs/142.md present: %s\n" \
  "$(command grep -cx 'lean/docs/029.md' "$TREE" || true)" "$(command grep -cx 'lean/docs/142.md' "$TREE" || true)"

echo; echo "== 140  Lemma s:mixture =="
D=$P/Subsphere-methods-for-memory-sample-lower-bounds-in-noiseless-Gaussian-regression-September-27-2026
line $D/build/sections/07-beta-mixture.tex 40 43
line $D/build/sections/07-beta-mixture.tex 85 103
echo; echo "done"
