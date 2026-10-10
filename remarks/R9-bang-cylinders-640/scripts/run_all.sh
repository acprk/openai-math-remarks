#!/bin/bash
# Reproduce the R9 verification.  Requires SageMath (sage -python with numpy/scipy; PPL backend).
# Uses 3 processes for the pair LPs (~6 min wall), ~5 min for the geometric sampler.
set -e
cd "$(dirname "$0")/independent"
sage -python make_params.py ../primary/persector_0.4_320.npy > params_T.txt
md5sum -c ../params_T.md5
for r in "0 107" "107 214" "214 320"; do set -- $r
  sage -python pair_lp_exact.py params_T.txt $1 $2 > pairs_$1.txt 2> pairs_$1.err & done; wait
cat pairs_0.txt pairs_107.txt pairs_214.txt > pairs_all.txt
md5sum -c ../pairs_all.md5 || echo "WARNING: pairs_all.txt md5 differs (maximizer points can differ across Sage/PPL versions; the Z values are what matter)"
python3 summarize_pairs.py pairs_all.txt            | diff - ../expected_output/pair_lp_summary.txt
sage -python area_interval.py params_T.txt          | diff - ../expected_output/area_interval.txt
sage -python amin_check.py                           | diff - ../expected_output/amin_check.txt
sage -python geom_sample.py params_T.txt 20261009 pairs_all.txt | diff - <(grep -v '^real\|^user\|^sys\|^$' ../expected_output/geom_sample.txt)
echo "ALL OUTPUTS MATCH"
# Optional (slow, ~20 min): the primary certificate by the original author of the improvement:
#   cd ../primary && sage -python cyl_cert.py 2/5 320 persector_0.4_320.npy   (compare with cert_04_320_ps.log)
# Optional (~50 min): regenerate the margins: sage -python cyl_persector.py 0.4 320 3 0.002
