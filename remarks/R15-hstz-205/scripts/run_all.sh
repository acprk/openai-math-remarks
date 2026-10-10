#!/usr/bin/env bash
# run_all.sh -- reproduce every computation of remark R15.
# Environment used: Ubuntu 22.04, g++ 9.5 (OpenMP), GAP 4 (generic character tables),
# SageMath 9.5 (bundled lrcalc), Python 3.10.  Uses <= 3 threads.
# Wall times on the reference host (shared, nice'd): m<=10 seconds; m=11 ~2 min (2 thr);
# m=12 ~24 min with 3 threads and ~11 GB RAM.  Outputs go to ./out/, compare with ../expected_output/.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p out
g++ -O3 -march=native -fopenmp -o out/kron_stair kron_stair.cpp

# (a) Proposition 1 vs direct A_n decomposition, all self-conjugate lambda, 5 <= n <= 20; also (b).
nice gap -q -o 8g an_square_check.g > out/an_square_check.out 2>&1

# (c) direct certificate: g(rho_m,rho_m,.) mod p and a = chi(hhat), every lambda |- N_m
P=2147483647
for m in 3 4 5 6 7 8 9 10; do nice out/kron_stair $m $P 7 3 out/kron_stair_m$m.txt; done
nice out/kron_stair 11 $P 9 3 out/kron_stair_m11.txt
nice out/kron_stair 12 $P 9 3 out/kron_stair_m12.txt          # the long one
# exact values for m <= 9 via four primes + CRT
for m in 3 4 5 6 7 8 9; do
  for q in 2147483647 2147483629 2147483587 2147483579; do nice out/kron_stair $m $q 7 3 out/ks_m${m}_P$q.txt; done
  python3 crt_exact.py $m out/ks_m${m}_P*.full > out/crt_exact_m$m.txt
done
# GAP cross-check of the full g-vectors for m = 3..6
(cd out && nice gap -q ../gap_stair_g.g > /dev/null 2>&1)
for m in 3 4 5 6; do
  sort out/gap_g_m$m.txt > out/a.s; sort out/ks_m${m}_P$P.txt.full 2>/dev/null > out/b.s || sort out/kron_stair_m$m.txt.full > out/b.s
  cmp -s out/a.s out/b.s && echo "m=$m identical" || echo "m=$m DIFF"
done

# (d) the 205 lower bound w_m: independent Sage implementation
(cd out && nice sage ../w205.sage hooks 16 > w205_hooks_m1-16.txt)
(cd out && nice sage ../w205.sage full 8 > w205_full.txt)
for m in 3 4 5 6 7 8; do python3 compare_w_g.py $m out/w_m$m.txt out/ks_m${m}_P*.full 2>/dev/null || true; done
