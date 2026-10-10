#!/usr/bin/env bash
# Check the R1 remark (r(4,t) = t^3/(log t)^{2+o(1)}) against openai/math.
# Environment used for the recorded run: Linux x86_64, elan, lean4 v4.34.1,
# 8 cores.  Usage:  bash scripts/check.sh /path/to/openai-math-clone
# Expected output: see expected_output/lean_check.txt
set -euo pipefail
MATH=${1:?path to a clone of github.com/openai/math}
PIN=fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd "$MATH"
# The SharpRamsey sources must match the pinned tree.
test "$(git rev-parse HEAD:lean/OAI/Combinatorics/SharpRamsey)" = \
     "$(git rev-parse $PIN:lean/OAI/Combinatorics/SharpRamsey)"
cd lean
lake exe cache get || true          # Mathlib oleans
lake build OAI.Combinatorics.SharpRamsey.Main
lake env lean "$HERE/lean/R4Lower.lean"
# Optional, slow: full kernel replay of the library dependency closure.
# LEAN_NUM_THREADS=1 lake env leanchecker --fresh OAI.Combinatorics.SharpRamsey.Main
