#!/usr/bin/env bash
# R12 reproduction.  Python 3.8+ stdlib only, single process.
#   fast  (default): verify_points.py, < 1 s
#   full  (./run_all.sh full): + exhaustive least-g search (independent, ~6 min)
#                              + the primary session's search script (~10-20 min)
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
E=expected_output; T=$(mktemp -d)
sha256sum -c SHA256SUMS --quiet
python3 independent/verify_points.py > $T/verify_points.txt
diff -q $T/verify_points.txt $E/verify_points.txt && echo "ok  verify_points"
if [[ "${1:-}" == full ]]; then
  python3 independent/indep_budget_053.py > $T/indep_budget.txt
  diff -q $T/indep_budget.txt $E/indep_budget.txt && echo "ok  indep_budget"
  python3 primary/pixton_budget.py > $T/primary_pixton_budget.txt
  diff -q $T/primary_pixton_budget.txt $E/primary_pixton_budget.txt && echo "ok  primary_pixton_budget"
fi
