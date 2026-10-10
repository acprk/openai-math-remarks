# VERIFICATION — R22 (four-state Potts sign certificate, third-party re-execution)

Host: Intel Xeon Gold 6230R, Ubuntu 22.04.5, Linux 6.8.0, Boost 1.86 headers; one thread per run.
`scripts/run_all.sh` (default: Clang 17.0.6 `-O2`) reproduced `expected_output/summary.txt` and ended with
`ALL OK`, exit 0, on 2026-10-09.

## 1. Who checked what

| Item | Primary | Second toolchain ("independent re-implementation" replaced by another compiler / level, as instructed) | Adversarial review | Status |
|---|---|---|---|---|
| Full traversal (a) | 2026-10-07, GCC 9.5.0 `-O3`: trace identical to upstream | 2026-10-09: Clang 17.0.6 `-O2` and GCC 12.3.0 `-O1`, full, identical; `run_all.sh` again with Clang | reviewer re-ran `run_all.sh` (ALL OK, 3m32s) | reproduced 4× |
| Full traversal (b) | 2026-10-07, GCC 9.5.0 `-O3`: trace identical to upstream | 2026-10-09: Clang 17.0.6 `-O2`, first 30 min (8 of 19 progress lines, 800 of 1935 cells), prefix identical | not re-run (20–90 min); the reviewer checked file times against the logged run time | reproduced 1× in full + 1 prefix |
| 2026-10-07 binary flags | not logged at the time | rebuild `g++ -std=c++17 -O3 -UNDEBUG checks/certificate.cpp` gives the same SHA-256 `a18f5b85…` | reviewer reproduced it, and showed `-O2` gives a different hash | identified |
| `small_polynomials.py` | 2026-10-09 | — | reviewer re-ran | JSON identical to upstream |
| Upstream hashes | `expected_hashes.txt` (11 files) | — | reviewer recomputed | match |

## 2. Adversarial review (separate AI session, 2026-10-09): **CONFIRMED WITH FIXES**

Core claims confirmed. The archived traces are byte-identical to upstream (SHA-256 439cc59… and 6d4befc…).
The 2026-10-07 file times fit real runs: binary 10:10:57, mode-a trace 10:26:44 (= +15m47.5s, the logged real
time), mode-b trace 11:57:06 (= +90m22.1s). The binary rebuilds bit-for-bit. `run_all.sh` ends with ALL OK.
All statement references (Theorem 1.1, Proposition 3.1, Lemma 5.1, Lemma 6.1, Proposition 6.2, Table 2) match
the compiled preprint. Findings and fixes:

1. **must-fix — cited commit not public.** `08fcd618…` is a local, unpushed commit of our clone. Fixed: the
   record now cites `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (origin/main, 2026-10-08). The bundle is
   identical between the two (`git diff` empty).
2. **must-fix — `run_all.sh` could print ALL OK after a mode-b mismatch, accepted an empty prefix, and always
   exited 0.** The reviewer demonstrated all three with a fake compiler and a tampered upstream copy. Fixed:
   the script now exits 1 on any DIFF; optional mode-b checks set the exit status; a prefix must contain at
   least one complete line; compilation and crashes are caught.
3. **must-fix — README referred to a missing VERIFICATION.md.** Fixed (this file).
4. **must-fix — local paths could ship.** `scripts/out/` contained a binary with an absolute `__FILE__` path,
   and LaTeX by-products contained the working directory. Fixed: `run_all.sh` now compiles from inside the
   bundle with a relative path, and `out/` and all LaTeX by-products are removed from the record.
5. **minor — upstream also has `mode-a-ubsan.log`** (identical to `mode-a.log`, not referenced), apparently a
   second authors' run of (a). Now mentioned in README and note, hashed, and compared by `run_all.sh`.
6. **minor — duplicate search should name `carlostoledo1891/cert-machine`.** That project lists this preprint as
   row F-433, "NOT YET DECIDED", and plans a clean-room re-implementation. Now cited in README and note, with
   the remark that it would supersede this record.
7. **minor — archived traces cannot by themselves prove a run.** Now stated in README ("What a reader can
   re-check") and note §4. Compiler banners, commands, binary hashes, timings and file times are in
   `scripts/logs/provenance.txt`.
8. **minor — wording.** "Equality of traces means the same cells were visited" → "strong evidence (though not
   a proof)". "Rules out UB" → "evidence against UB on the toolchains tested". "Different compiler (GCC 9.5.0)"
   → "another GCC version".
9. **minor — timing inconsistency** (`run_all.sh` mode a took 503 s wall on the shared host). README now says
   3–9 min wall.
10. **cosmetic.** The time file's blocks are now labelled in `provenance.txt`. Prop. 6.2 is described as the
    finite verification reached via Prop. 3.1 and Lemma 6.1. The 24 permutations are attributed to §5, not
    to Lemma 5.1.

## 3. Commands and key outputs

```
cd scripts && OPENAI_MATH=<clone of openai/math> ./run_all.sh
  == 1. upstream hashes ==            11 x OK
  == 2. archived traces vs upstream == 9 x OK (3 full mode-a traces, 1 full mode-b trace, paper-mode-* = mode-* + time,
                                      mode-a-ubsan.log == mode-a.log, 8-line mode-b prefix)
  == 3. fresh run of mode a ==        done 4479 19 lp=2074 lt=166 ; OK fresh mode-a trace == upstream mode-a.log
  == 5. small_polynomials.py ==       OK ; minima H+1/2 = 19/50, F = 0, G = 87/1000, (3-Euler)H = 0
  ALL OK
```
2026-10-07 commands (in a copy of the upstream `verification/` directory):
`g++ -std=c++17 -O3 -UNDEBUG checks/certificate.cpp -o artifacts/certificate; time ./artifacts/certificate;
time ./artifacts/certificate b` → `done 4479 19 lp=2074 lt=166` (real 15m47s, user 5m51s) and
`done 1935 20 lp=871 lt=97` (real 90m22s, user 45m06s). `python3 checks/verify.py` → all PASS, "Full
traversals were not repeated".

## 4. Remaining gaps

- **GAP-1.** Mode (b) was traversed in full only once by us (GCC 9.5.0). The second toolchain covers only an
  8-line prefix, because a full run takes more than 30 min here.
- **GAP-2.** No re-implementation of the subdivision or the sign tests, and no audit of the reduction
  (preprint Sections 3–6, Appendix A). Trace equality shows only that the program is reproducible.
- **GAP-3.** The full-run claims for 2026-10-07 rest on local records (file times, `time` output) that a reader
  cannot inspect; mode (a) can be re-checked cheaply with `run_all.sh`.
