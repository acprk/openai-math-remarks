#!/usr/bin/env bash
# Re-execution record for the exact sign certificate of family 229
# ("The Reconstruction Threshold for the Ferromagnetic Four-State Potts Model", openai/math).
#
# Usage:
#   OPENAI_MATH=/path/to/clone/of/openai/math ./run_all.sh            # default: ~3-9 min, 1 thread
#   OPENAI_MATH=... CXX=g++ OPT=-O3 ./run_all.sh                      # choose compiler / optimisation level
#   OPENAI_MATH=... FULL_B=1 ./run_all.sh                             # also run mode b in full (20-90 min)
#   OPENAI_MATH=... B_PREFIX_SECONDS=900 ./run_all.sh                 # run mode b for N seconds, check the prefix
#
# Needs: a C++17 compiler with GCC-compatible <bits/stdc++.h>, Boost.Multiprecision headers, python3.
# Steps (all deterministic; the summary does not depend on compiler or timing):
#  1. hash-check the upstream source, source input, logs and helper scripts (expected_hashes.txt);
#  2. check that the archived traces in logs/ are byte-identical to (resp. a prefix of) the upstream traces;
#  3. compile certificate.cpp with $CXX $OPT (assertions on: -UNDEBUG), run mode a, compare the full trace;
#  4. optionally mode b (full, or a time-limited prefix that must contain at least one progress line);
#  5. run upstream small_polynomials.py and compare its JSON with the recorded upstream output.
# Steps 1-3 and 5 are written to out/summary.txt and diffed against expected_output/summary.txt; any DIFF,
# a failed optional mode-b check, or a summary mismatch makes the script exit with status 1.
set -uo pipefail
cd "$(dirname "$0")"
HERE=$PWD
: "${OPENAI_MATH:?set OPENAI_MATH to a clone of github.com/openai/math (tested at commit fd4aeeb2)}"
D="$OPENAI_MATH/preprints/The-Reconstruction-Threshold-for-the-Ferromagnetic-Four-State-Potts-Model-October-5-2026/verification"
[ -f "$D/checks/certificate.cpp" ] || { echo "upstream bundle not found under $D"; exit 1; }
CXX=${CXX:-$(command -v clang++ >/dev/null 2>&1 && echo clang++ || echo g++)}
OPT=${OPT:--O2}
OUT="$HERE/out"; rm -rf "$OUT"; mkdir -p "$OUT"
S="$OUT/summary.txt"; : > "$S"
FAIL=0
ok()  { echo "OK   $*" | tee -a "$S"; }
bad() { echo "DIFF $*" | tee -a "$S"; FAIL=1; }
say() { echo "$*" | tee -a "$S"; }
sha() { sha256sum "$1" | cut -d' ' -f1; }

say "== 1. upstream hashes =="
while read -r h f; do
  if [ -f "$D/$f" ] && [ "$(sha "$D/$f")" = "$h" ]; then ok "$f"; else bad "$f"; fi
done < expected_hashes.txt

say "== 2. archived traces vs upstream =="
for t in logs/full-run-2026-10-07-gcc9.5-O3-mode-a.trace logs/spot-2026-10-09-clang17-O2-mode-a.trace \
         logs/spot-2026-10-09-gcc12-O1-mode-a.trace; do
  if cmp -s "$t" "$D/checks/logs/mode-a.log"; then ok "$(basename "$t") == upstream mode-a.log"; else bad "$(basename "$t")"; fi
done
t=logs/full-run-2026-10-07-gcc9.5-O3-mode-b.trace
if cmp -s "$t" "$D/checks/logs/mode-b.log"; then ok "$(basename "$t") == upstream mode-b.log"; else bad "$(basename "$t")"; fi
# upstream: paper-mode-*.log = mode-*.log + the authors' `time` output; mode-a-ubsan.log = mode-a.log
for m in a b; do
  n=$(wc -c < "$D/checks/logs/mode-$m.log")
  if head -c "$n" "$D/checks/logs/paper-mode-$m.log" | cmp -s - "$D/checks/logs/mode-$m.log"; then
    ok "upstream paper-mode-$m.log = mode-$m.log + timing lines"; else bad "paper-mode-$m.log prefix"; fi
done
if cmp -s "$D/checks/logs/mode-a-ubsan.log" "$D/checks/logs/mode-a.log"; then ok "upstream mode-a-ubsan.log == mode-a.log"
else bad "mode-a-ubsan.log"; fi
prefix_ok() {   # $1 = candidate prefix trace, $2 = full trace; needs >= 1 complete line
  local n l; n=$(wc -c < "$1"); l=$(wc -l < "$1")
  [ "$l" -ge 1 ] && [ "$(tail -c 1 "$1" | od -An -c | tr -d ' ')" = '\n' ] && head -c "$n" "$2" | cmp -s - "$1"
}
t=logs/spot-2026-10-09-clang17-O2-mode-b-prefix.trace
if prefix_ok "$t" "$D/checks/logs/mode-b.log"; then ok "$(basename "$t") ($(wc -l < "$t") lines) is a prefix of upstream mode-b.log"
else bad "$(basename "$t")"; fi

say "== 3. fresh run of mode a =="
echo "compiling with $CXX $OPT ..." >&2
# compile from inside the bundle with a relative source path (no absolute path is embedded in the binary)
if ( cd "$D" && "$CXX" -std=c++17 $OPT -UNDEBUG checks/certificate.cpp -o "$OUT/certificate" ); then
  "$CXX" --version | head -1 > "$OUT/compiler.txt"
  start=$(date +%s); "$OUT/certificate" 2> "$OUT/mode-a.trace"; rc=$?; end=$(date +%s)
  echo "mode a: exit $rc, wall time $((end-start)) s" >> "$OUT/compiler.txt"
  say "$(tail -1 "$OUT/mode-a.trace")"
  if [ $rc -eq 0 ] && cmp -s "$OUT/mode-a.trace" "$D/checks/logs/mode-a.log"; then ok "fresh mode-a trace == upstream mode-a.log"
  else bad "fresh mode-a trace (exit $rc)"; fi
else
  bad "compilation failed"
fi

say "== 5. small_polynomials.py =="
if python3 "$D/checks/small_polynomials.py" > "$OUT/small_polynomials.json"; then
  python3 - "$D/checks/logs/paper-small-polynomials.json" "$OUT/small_polynomials.json" <<'EOF' | tee -a "$S"
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:3])
print("OK   small_polynomials.py output equals upstream paper-small-polynomials.json" if a == b
      else "DIFF small_polynomials.py output")
print("     minima:", [(c["name"], c["minimum"]) for c in b["small_checks"]])
EOF
else
  bad "small_polynomials.py failed"
fi
grep -q '^DIFF' "$S" && FAIL=1

# optional mode b: reported separately (not part of the fixed summary), but a failure sets the exit status
if [ "${FULL_B:-0}" = 1 ] && [ -x "$OUT/certificate" ]; then
  "$OUT/certificate" b 2> "$OUT/mode-b.trace"; rc=$?
  tail -1 "$OUT/mode-b.trace"
  if [ $rc -eq 0 ] && cmp -s "$OUT/mode-b.trace" "$D/checks/logs/mode-b.log"; then echo "OK   fresh mode-b trace == upstream mode-b.log"
  else echo "DIFF fresh mode-b trace (exit $rc)"; FAIL=1; fi
elif [ -n "${B_PREFIX_SECONDS:-}" ] && [ -x "$OUT/certificate" ]; then
  timeout "$B_PREFIX_SECONDS" "$OUT/certificate" b 2> "$OUT/mode-b-prefix.raw"; rc=$?
  if [ $rc -ne 0 ] && [ $rc -ne 124 ]; then echo "DIFF mode b crashed (exit $rc)"; FAIL=1; fi
  grep -E '^1: ' "$OUT/mode-b-prefix.raw" > "$OUT/mode-b-prefix.trace" || true
  [ $rc -eq 0 ] && cp "$OUT/mode-b-prefix.raw" "$OUT/mode-b-prefix.trace"
  if prefix_ok "$OUT/mode-b-prefix.trace" "$D/checks/logs/mode-b.log"; then
    echo "OK   fresh mode-b prefix ($(wc -l < "$OUT/mode-b-prefix.trace") lines) matches upstream"
  else echo "DIFF fresh mode-b prefix (empty or mismatching; try a larger B_PREFIX_SECONDS)"; FAIL=1; fi
fi

if ! diff -q "$S" expected_output/summary.txt >/dev/null; then
  echo "SUMMARY DIFFERS from expected_output/summary.txt"; diff "$S" expected_output/summary.txt; FAIL=1
fi
if [ $FAIL -eq 0 ]; then echo "ALL OK"; else echo "FAILED"; fi
exit $FAIL
