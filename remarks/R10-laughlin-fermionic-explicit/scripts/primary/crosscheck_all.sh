#!/bin/bash
# primary cross-check (Jordan-Wigner sectors) on 10 random fluxes (seed 269) in [24,40], N=3,4 (+ N=2,5 at Q=24)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P=${PY:-python3}
run() { while read Q GF G; do $P crosscheck.py $GF $Q $1 $G; done < data/xc_list.txt; }
run 3 > out/crosscheck_N3.log &
run 4 > out/crosscheck_N4.log &
wait
$P crosscheck.py data/gamma_24_36_48.json 24 2,5 $(awk '$1==24{print $3}' data/xc_list.txt) > out/crosscheck_Q24.log
