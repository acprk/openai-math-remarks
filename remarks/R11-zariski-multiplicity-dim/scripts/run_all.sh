#!/usr/bin/env bash
# R11 reproduction.  Python 3.8+ standard library only; single process; ~15 s total.
# Each output is diffed against expected_output/.  The relation search itself
# (primary/*.sage, SageMath 9.5, ~55 min) is optional and not run here.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
E=expected_output; T=$(mktemp -d)
NEW=data/relation_d3pow11_N1802976.json
OLD=data/openai_explicit-relation.json
sha256sum -c SHA256SUMS --quiet
python3 primary/openai_explicit_relation.py $NEW $T/rep_new.json > $T/openai_verifier_new.txt
python3 primary/openai_explicit_relation.py $OLD $T/rep_old.json > $T/openai_verifier_original.txt
python3 independent/indep_check_059.py $NEW --spectral > $T/indep_new.txt
python3 independent/indep_check_059.py $OLD > $T/indep_original.txt
python3 independent/negative_control.py $NEW > $T/negative_control.txt
for f in openai_verifier_new openai_verifier_original indep_new indep_original negative_control; do
  diff -q $T/$f.txt $E/$f.txt && echo "ok  $f"
done
