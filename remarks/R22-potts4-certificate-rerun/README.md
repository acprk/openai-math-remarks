# R22 — Four-state Potts reconstruction threshold: a third-party re-execution of the exact sign certificate

## Premise and confidence (read first)

- **Kind of result: a reproducibility record, not mathematics.** We ran OpenAI's published verifier
  `certificate.cpp` unchanged and compared its output with the published logs. We did not re-implement the
  certificate, and we did not check the reduction from the theorem to the coefficient arrays.
- **What was already on record.** The authors' own full runs are recorded upstream in `checks/full-runs.json`
  (GCC 13.3.0, `-O3`, 98 s and 1288 s) together with their logs. The bundle also contains
  `checks/logs/mode-a-ubsan.log`, identical to `mode-a.log` and referenced nowhere, which suggests a second
  authors' run of mode (a), apparently under UBSan. The bundled lightweight checker `checks/verify.py`
  validates the logs but, in its own words, "Full traversals were not repeated". We found no third-party
  record of a full traversal (see Duplicate search). So the published counts rested on the authors' own
  executions. This remark adds executions on a different machine with other compiler versions.
- **Confidence: high, for what is claimed.** All comparisons are byte-for-byte on deterministic output of an
  exact-arithmetic program.

**One line.** Before → After: published counts `done 4479 19 lp=2074 lt=166` (a) and `done 1935 20 lp=871 lt=97`
(b), backed only by the authors' own runs → **both full traversals re-executed** on a different machine
with another GCC version (9.5.0) and Boost 1.86, and the complete stderr traces are **byte-identical** to the
published ones. Traversal
(a) is also byte-identical under Clang 17.0.6 `-O2` and GCC 12.3.0 `-O1`, and the first 8 progress lines
(800 of 1935 cells) of (b) under Clang 17.0.6 `-O2`.

## Statement (Before / After)

| | Before (upstream bundle) | After (this record) |
|---|---|---|
| Full traversal (a) | authors' run, GCC 13.3.0 -O3; logs + hashes | + GCC 9.5.0 -O3 (full), Clang 17.0.6 -O2 (full), GCC 12.3.0 -O1 (full): traces identical |
| Full traversal (b) | authors' run, GCC 13.3.0 -O3; logs + hashes | + GCC 9.5.0 -O3 (full): trace identical; Clang 17.0.6 -O2: first 30 min (8 progress lines) identical |
| `verify.py` | validates source binding and recorded logs; does not rerun | also passes here (2026-10-07) |
| `small_polynomials.py` | recorded JSON | rerun: identical JSON |
| Machine / compiler info, source hash, commands | GCC version and flags in `full-runs.json` | all of the above for each run, below |

## Attribution

The theorem (*The Reconstruction Threshold for the Ferromagnetic Four-State Potts Model*, Theorem 1.1:
nonreconstruction for dλ² ≤ 1 on d-ary and Poisson(d) trees), its reduction (Proposition 3.1 → Lemma 6.1,
Proposition 6.2, Table 2), and the verifier are OpenAI's work: family 229 of
[github.com/openai/math](https://github.com/openai/math), commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/The-Reconstruction-Threshold-for-the-Ferromagnetic-Four-State-Potts-Model-October-5-2026/`. This
record only executes the published program and compares outputs.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction;
an independent adversarial review was run by a separate AI session and is not human peer review.

## Runs

Machine: Intel Xeon Gold 6230R @ 2.10 GHz, Linux 6.8.0, Ubuntu 22.04.5 LTS; one thread per run (the host is
CPU-quota limited, so wall time exceeds user time). Boost headers 1.86. Source `checks/certificate.cpp`,
SHA-256 `9f9765fb4dee9f2fc618390df89d3abccf2602569709fd4c2133a201b6e86bff`, which equals the upstream
`full-runs.json` value. All builds keep assertions on (`-UNDEBUG`).

| Date | Compiler | Flags | Mode | User time | Result |
|---|---|---|---|---|---|
| 2026-10-07 | GCC 9.5.0 | `-std=c++17 -O3 -UNDEBUG` | a, full | 5m51s | `done 4479 19 lp=2074 lt=166`; trace == upstream `mode-a.log` |
| 2026-10-07 | GCC 9.5.0 | `-std=c++17 -O3 -UNDEBUG` | b, full | 45m06s | `done 1935 20 lp=871 lt=97`; trace == upstream `mode-b.log` |
| 2026-10-09 | Clang 17.0.6 | `-std=c++17 -O2 -UNDEBUG` | a, full | 2m54s | trace == upstream `mode-a.log` |
| 2026-10-09 | GCC 12.3.0 | `-std=c++17 -O1 -UNDEBUG` | a, full | 4m00s | trace == upstream `mode-a.log` |
| 2026-10-09 | Clang 17.0.6 | `-std=c++17 -O2 -UNDEBUG` | b, 30-min prefix | 19m40s | 8 lines, prefix of upstream `mode-b.log` |
| 2026-10-09 | Clang 17.0.6 | `-std=c++17 -O2 -UNDEBUG` | a, full (`run_all.sh`) | 3–9 min wall (shared host) | trace == upstream `mode-a.log` |

The 2026-10-07 binary was built from the upstream `verification/` directory. Its flags were not logged at the
time, but rebuilding with `g++ -std=c++17 -O3 -UNDEBUG checks/certificate.cpp` (GCC 9.5.0) in the same directory
reproduces its SHA-256 (`a18f5b85…`) exactly. `-O2` gives a different hash, so the hash does identify the
optimisation level; `-UNDEBUG` is a no-op here. The upstream `paper-mode-{a,b}.log` files are `mode-{a,b}.log`
followed by the authors' `time` output; `run_all.sh` checks this. Compiler banners, commands, binary hashes,
timings and file times for every archived trace are in `scripts/logs/provenance.txt`.

**What a reader can re-check.** The traces contain only deterministic counts, with no host, compiler or
timestamp, so an archived trace by itself cannot prove that a run took place. The local evidence (binary and
log modification times that differ by exactly the logged run times; raw `time` output) supports real runs, but
it does not travel with this record. A reader can re-check mode (a) in 3–6 CPU minutes and mode (b) in 20–45
with `run_all.sh`. The claim that mode (b) was run in full on 2026-10-07 rests on our local record.

Commands (from the upstream `verification/` directory):
```bash
g++ -std=c++17 -O3 -UNDEBUG checks/certificate.cpp -o artifacts/certificate   # or clang++ -O2, g++-12 -O1
./artifacts/certificate   2> mode-a.trace     # ~3-6 min CPU here
./artifacts/certificate b 2> mode-b.trace     # ~20-45 min CPU here
cmp mode-a.trace checks/logs/mode-a.log && cmp mode-b.trace checks/logs/mode-b.log
python3 checks/small_polynomials.py           # compare with checks/logs/paper-small-polynomials.json
```

## Reproduce

```bash
cd scripts
OPENAI_MATH=/path/to/clone/of/openai/math ./run_all.sh               # ~3-5 min, one thread
OPENAI_MATH=... CXX=g++ OPT=-O3 ./run_all.sh                         # other compiler / level
OPENAI_MATH=... FULL_B=1 ./run_all.sh                                # also mode b in full (20-90 min)
```
`run_all.sh` (1) checks the SHA-256 of the upstream source, source input, logs and helper script against
`expected_hashes.txt`; (2) checks that the archived traces in `scripts/logs/` equal the upstream traces;
(3) compiles and runs mode (a) afresh and compares the full trace; (4) optionally runs (b); (5) reruns
`small_polynomials.py`. Its compiler-independent summary is diffed against `expected_output/summary.txt`
and ends with `ALL OK`.

## What is NOT claimed

- **Not an independent verification of the theorem.** We did not re-derive Sections 3–6, re-implement the
  midpoint subdivision or the Type P / Type T tests, or audit Appendix A of the preprint. Agreement of traces
  shows that the published program gives the published output on another toolchain. It does not show that the
  program checks the right inequalities.
- **Not the first execution.** The authors ran both traversals (mode a apparently twice); this is the first
  *third-party* record we know of.
- **Mode (b) with a second compiler is only a prefix** (8 of 19 progress lines; the full run takes more than
  30 min here). The full (b) traversal was repeated only with GCC 9.5.0.
- **No new mathematics, no Lean.**

## Duplicate search (to 2026-10-09)

To our knowledge (search of GitHub and arXiv up to 2026-10-09), no third-party re-run of this verifier has
been published. GitHub repository search ("four-state Potts", "Potts reconstruction", "kesten stigum potts")
found nothing, and code search for the output string `lp=2074` finds only openai/math and the project noted
next. **Related:** the community project `carlostoledo1891/cert-machine` ("Independent exact certification of
machine-generated mathematics", updated 2026-10-10 UTC) lists this preprint as row F-433 of
`certs/openai-math-finite.json`, with status "NOT YET DECIDED" ("no decider written yet"). It plans a
clean-room re-implementation with its own subdivision. That would be a stronger check than this record, and it
may supersede it. arXiv: we know of no paper that re-checks this certificate.

## Files

- `note.tex`, `note.pdf` — 2-page record.
- `scripts/run_all.sh`, `scripts/expected_hashes.txt`, `scripts/expected_output/summary.txt`.
- `scripts/logs/` — archived traces: full runs of 2026-10-07 (GCC 9.5.0 -O3, modes a and b), spot checks of
  2026-10-09 (Clang 17 -O2 mode a, GCC 12 -O1 mode a, Clang 17 -O2 mode b prefix), the 2026-10-07 `time`
  output, and `provenance.txt`.
- `VERIFICATION.md` — checks and the adversarial review.
