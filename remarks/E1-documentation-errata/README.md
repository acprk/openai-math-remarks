# E1 — Documentation, packaging and reproducibility errata

## Premise and confidence (read first)

**No mathematical error is claimed.** Every item below is a missing companion file, a broken path,
a Lean-docs page that describes its coverage too narrowly, or a sentence whose wording could be read
the wrong way. In each case the mathematics, as far as we checked it, is unaffected. Where a missing
script certifies a finite computation and we could recompute it (families 222, 238, 318, 266-B), our
recomputation **agrees with the paper**.

**One line (Before → After).** Before: seven manuscripts cite companion files or paths that are not in the repository, and
six Lean-docs, CONTENTS or wording spots are imprecise. After: for each, a minimal diff (add the file, fix the
path, or reword one line) with locations at the current `main`, plus our re-runs of the four finite
computations we could reproduce (222, 238, 318, 266-B), all of which agree with the papers.

All locations are relative to [github.com/openai/math](https://github.com/openai/math) at commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (the October 7, 2026 update; `main` as of 2026-10-09,
re-fetched before writing). The repository has GitHub issues disabled, so these notes are collected here.

Every item was re-checked against this commit by a separate adversarial AI session instructed to refute it
(see `VERIFICATION.md`). Three points from our internal list were withdrawn, narrowed or downgraded as a result (§4).

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s
direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Summary table

| # | Family | Manuscript (under `preprints/`) or file | Issue | Kind | Our check |
|---|---|---|---|---|---|
| P1 | 266 | The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026 | `verification/` tree cited by Thm 1.1 is absent | missing artifacts (most significant) | absence only |
| P2 | 095 | Hyperbolicity-Cones-Without-Semidefinite-Lifts-October-5-2026 | `verify_certificate.py` absent | missing script | absence only |
| P3 | 318 | A-rational-obstruction-to-strong-chromatic-splitting-at-height-three-September-25-2026 | `support/`, `BUILDING.md`, `SUPPORT-COVERAGE.json` absent | missing support files | **ranks recomputed: agree** |
| P4 | 222 | Microscopic-jamming-in-the-negative-spherical-perceptron-September-24-2026 | `certificate/original_certificate.py` and run record absent | missing script | **listing extracted, SHA-256 matches, run passes** |
| P5 | 238 | Routing-densities-and-representation-contraction-for-Thorp-sweeps-September-26-2026 | `support/network_bounds.py` absent | missing script | **5 inequalities recomputed: all hold** |
| P6 | 227 | Cutoff-throughout-the-high-temperature-Sherrington-Kirkpatrick-phase-September-24-2026 | two `numerics-*.py` files absent | dangling reference | absence only |
| P7 | 266 | Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026 | listing path `certificate/verify.py` should be `../verification/verify.py` | wrong path | **shipped verifier re-run: output matches Appendix A** |
| D1 | 107 | `lean/docs/107.md` | 9/4 manuscript not in the paper list | Lean docs | — |
| D2 | 026 | `lean/docs/026.md` | Scope omits the main formalized theorem | Lean docs (understatement) | — |
| D3 | 142 | `CONTENTS.md` | family summary line does not name the 029 dependency | transparency | — |
| W1 | 103 | Exact-Derandomization-of-Logarithmic-Space-L-equals-RL-equals-BPL-September-23-2026 | "increases … by at most 1/q_fr" reads as additive | wording | — |
| W2 | 140 | Subsphere-methods-for-memory-sample-lower-bounds-in-noiseless-Gaussian-regression-September-27-2026 | Lemma states 2^{D'}, proof gives 2^{(D'−n)/2} | optional sharpening | numerical illustration |
| W3 | 099 | `lean/docs/099.md` | "The paper's" is ambiguous in a three-paper page | optional wording | — |

## 1. Missing or misnamed companion files (P1–P7)

### P1 — 266, N(6) = 3: the cited verification archive is not in the repository

- **Where.** `…/build/sections/introduction.tex` l.13–18 (Theorem 1.1: "the fresh complete execution
  documented in Section 7 certifies …"), l.154–157 (`verification/computation/README.md`);
  `…/build/sections/execution.tex` l.16–20 (`verification/code/` with `seed.cpp`, `covout.cpp`, `ub.cpp`,
  `patchx.cpp`, `atlas.cpp`, `stageB.cpp`, `stageC.cpp`, listed in `verification/code/extraction.json`),
  l.119 (`verification/code/proof_run.py`), l.380–381 (`verification/computation/README.md`), l.390
  (`verification/REPRODUCE.md`).
- **Issue.** The manuscript directory contains only `README.md`, the PDF and `build/`. None of the cited
  files exists anywhere in the tree or its history. (Paths such as `full-run-01/run.json`, l.386–389, are
  inside the execution archive, which is itself not provided.) Theorem 1.1 is stated as a certificate *of
  that execution*. The arithmetic is binary64 round-to-nearest with an a-priori error analysis
  (`rounding.tex`), and `lean/docs/266.md` l.12 states that the Lean development "does not establish the
  paper's upper bound of three". So the headline N(6) = 3 cannot at present be re-run or audited by a reader.
  The sibling manuscript *Exact Fourier certificates…* does ship `verification/`, so this looks like an
  omission in packaging.
- **Suggested minimal diff.** Add `verification/` as described in Section 7 (the seven sources with the
  listed SHA-256 digests, `extraction.json`, `proof_run.py`, `REPRODUCE.md`, `computation/README.md`, and
  the archive or a stable link to it). If that is not possible soon, a one-line note in the manuscript
  README would help, e.g.:
  ```diff
  +Note: the computation artifacts cited in Section 7 (verification/) are not yet included in this release.
  ```
- **Why.** The Lean-checked bound in this family is N(6) ≤ 5. The step from 5 to 3 rests on this
  computation alone.

### P2 — 095, no semidefinite lift: `verify_certificate.py` is absent

- **Where.** `…/build/source/appendix-certificate.tex` l.536–546 ("The accompanying source file
  `verify_certificate.py` reconstructs the monomial blocks, … every displayed minor using exact integer
  arithmetic modulo seven").
- **Issue.** No file of that name is in the manuscript directory. (Three files with that name exist, but
  under the unrelated *Snaky-in-21-Maker-moves* manuscript.) The paper says (l.544–546) that the displayed
  minors and the characteristic-zero argument establish Lemma `lem:seed` without the program. That is fair,
  but the minors are F₇ determinants of blocks up to about 66×66 inside a rank-302 / dimension-321
  computation, so in practice a reader needs the script. Unlike the sibling nonspectrahedrality result,
  this manuscript is not linked from any `lean/docs` page (`lean/docs/095.md` links only *A nonspectrahedral
  hyperbolicity cone*).
- **Suggested minimal diff.** Add `verify_certificate.py` (and its expected output) to the manuscript
  directory. Optionally add one sentence to `lean/docs/095.md`: "The no-semidefinite-lift theorem of
  *Hyperbolicity Cones Without Semidefinite Lifts* is not part of this formalization."
- **Why.** The page title of `lean/docs/095.md` is "Hyperbolicity cones without semidefinite lifts", so a
  reader may assume the no-lift theorem is the formalized one.

### P3 — 318, chromatic splitting: `support/` tree absent (arithmetic reproduced)

- **Where.** `SUPPORT-COVERAGE.md` l.7 (`support/r-ce/ce_certificate.{py,json,txt}`), l.11 (`BUILDING.md`),
  l.13 (`support/historical/r-ce/`), l.15 (link to `SUPPORT-COVERAGE.json`);
  `build/sections/08-elementary-background.tex` l.58.
- **Issue.** The directory contains `README.md`, `SUPPORT-COVERAGE.md`, `paper.pdf` and `build/` only; the
  link at l.15 is broken.
- **Our check.** Two independent exact computations give the Chevalley–Eilenberg ranks d₀…d₃
  `sl2: (0,3,0,0)`, `sl3: (0,8,20,35)`, as stated (and H⁰…H³ = 1,0,0,1, consistent with
  H*(sl₃) = Λ(x₃, x₅)). Scripts: `scripts/318_ce_ranks.py`, `scripts/independent/indep_318.py`.
- **Suggested minimal diff.** Add the files, or:
  ```diff
  -The active guarded trio is in `support/r-ce/`. ...
  +The support files described here are not included in this release.
  ```

### P4 — 222, negative perceptron: `certificate/original_certificate.py` absent (certificate re-run)

- **Where.** `…/build/main.tex` l.9466–9476 ("The accompanying file `certificate/original_certificate.py` is
  the unaltered code printed in the appendix … The accompanying reproducibility record …"), l.9481
  ("the actual output is retained with the certificate").
- **Issue.** There is no `certificate/` directory and no run record.
- **Our check.** The listing between `\begin{lstlisting}` (l.9672) and `\end{lstlisting}` (l.10354),
  newline-terminated, has SHA-256 `7bd5dbd9…a85dc7f`, **identical** to the digest printed at l.9496–9497.
  Running it (Python 3.10.12, NumPy 2.2.6, assertions on, one thread, about 8 min) prints all 24 checkpoints,
  from `0 410000` to `283000 304000`, and exits 0. So the certificate is fully recoverable from the TeX.
- **Suggested minimal diff.** Add `certificate/original_certificate.py` plus the saved output, or:
  ```diff
  -The accompanying file \texttt{certificate/original\_certificate.py} is the
  -unaltered code printed in the appendix.  From the manuscript directory,
  -run
  +The code printed in Appendix~\ref{app:code}, saved as
  +\texttt{original\_certificate.py}, is the certificate.  Run
  ```

### P5 — 238, Thorp sweeps: `support/network_bounds.py` absent (inequalities recomputed)

- **Where.** `…/build/sections/04-certificate.tex` l.294 ("The accompanying script
  `support/network_bounds.py` checks exactly these five rational inequalities").
- **Our check.** From the definitions of h_j (l.134–136) and H(ρ) (l.138–140) we recomputed
  the five bounds in two independent ways (exact `Fraction` arithmetic following the paper's own
  tail bound; and a generic fractional-knapsack with a separate rigorous tail bound plus `mpmath.iv`
  interval logarithms):

  | Claim (l.272–282) | Exact upper bound (paper's route) | Independent upper bound |
  |---|---|---|
  | H(1/8) < .791402 | 0.791401931 | 0.7914013842 |
  | H(13/50) < 1.098069 | 1.098068598 | 1.0980680509 |
  | H(1) < 1.958069 | 1.958068598 | 1.9580680509 |
  | bracket(7/25) < .904054 | 0.904053756 | 0.9040532085 |
  | bracket(1) < .515374 | 0.515373557 | 0.5153730100 |

  All five hold. The margins are small (about 1e-7) but positive.
- **Suggested minimal diff.** Ship the script, or replace "The accompanying script
  `support/network_bounds.py` checks" by "Exact rational arithmetic checks".

### P6 — 227, SK cutoff: two `numerics-*.py` files absent

- **Where.** `…/build/appendices/scalar-certificate.tex` l.291–294.
- **Issue.** `numerics-interval.py` and `numerics-cross-errors.py` are not in the repository. The text says
  that they are "not needed to specify the certificate", so this is only a dangling reference. We did not
  re-run this certificate.
- **Suggested minimal diff.** Add the two files, or delete the sentence.

### P7 — 266, Exact Fourier certificates: listing path

- **Where.** `…/build/paper.tex` l.273 (`\lstinputlisting[…]{certificate/verify.py}`);
  `…/build/sections/04-mixed.tex` l.217–218 ("also supplied in the companion file `certificate/verify.py`").
- **Issue.** The shipped file is `verification/verify.py` (and `paper.tex` l.235 and l.259 already use that
  name; the checksum at l.259–262 matches the shipped file). There is no `certificate/` directory and no copy
  under `build/`, so `build/` does not compile as distributed.
- **Our check.** `python3 verification/run.py` passes its checksum test and prints the three lines displayed at
  `paper.tex` l.244–246 verbatim (`rank table: 8871 909 897`, …), about 3 min CPU single-threaded.
- **Suggested minimal diff.**
  ```diff
  -\lstinputlisting[style=certificate,caption={Complete exact verifier},label={lst:verifier}]{certificate/verify.py}
  +\lstinputlisting[style=certificate,caption={Complete exact verifier},label={lst:verifier}]{../verification/verify.py}
  ```
  and `certificate/verify.py` → `verification/verify.py` in `04-mixed.tex` l.218. (A plain rename to
  `verification/verify.py` would still not compile from `build/`.)

## 2. Lean docs and CONTENTS (D1–D3)

### D1 — `lean/docs/107.md` does not list the 9/4 manuscript

- **Where.** `lean/docs/107.md` l.5–6 (paper list) and l.10 (Scope: "bound the complex
  matrix-multiplication exponent by ω(ℂ) ≤ 9/4"); comparator
  `lean/ComparatorChallenges/MatrixMultiplication.lean` l.118 (`complex_omega_le_nine_quarters`).
- **Issue.** The 9/4 proof is in `preprints/Matrix-Multiplication-Nine-Fourths-October-2-2026/`
  (*An Upper Bound of 9/4 for the Matrix Multiplication Exponent*). The 2.258 manuscript does not contain it.
  CONTENTS.md attaches the Lean link at the family level, but the docs page itself does not list the paper
  whose main theorem it formalizes.
- **Suggested minimal diff.**
  ```diff
   The following describes the scope of the Lean formalization related to the following accompanying paper(s):
   
  +- [An Upper Bound of 9/4 for the Matrix Multiplication Exponent](../../preprints/Matrix-Multiplication-Nine-Fourths-October-2-2026/paper.pdf)
   - [Complex Matrix Multiplication Below 2.258 and Rectangular Bounds](...)
  ```

### D2 — `lean/docs/026.md` describes only the corollary

- **Where.** `lean/docs/026.md` l.9 (Scope) and l.15 (table row "Positive lower density of prime-ratio
  increases"); comparator `lean/ComparatorChallenges/PrimeGaps.lean` l.36–40.
- **Issue.** `large_gaps_and_ratio_density` is the conjunction of (i) for every C > 0, a positive
  proportion of n ≤ N have p_{n+1} − p_n > C log p_n (`largeGapIndices`), which is the paper's headline
  theorem, and (ii) the prime-ratio corollary. The docs page mentions only (ii). This understates what is
  verified.
- **Suggested minimal diff.**
  ```diff
  -The formalized supplement proves that the indices $n$ for which ...
  +The formalized statement proves the paper's main theorem: for every fixed $C>0$, a positive
  +proportion of the indices $n\le N$ satisfy $p_{n+1}-p_n>C\log p_n$. It also proves that the indices $n$ for which ...
  ```
  and the table row "Positive lower density of large prime gaps and of prime-ratio increases".

### D3 — 142: name the 029 dependency in the family summary line

- **Where.** `CONTENTS.md` l.3427 (family 142 summary: "No randomness, integer-factorization or
  primitive-root oracle, or GRH assumption is required.").
- **Issue.** The unconditional statement uses Theorem 1.2 of family 029 (*Primitive roots for every
  admissible integer base*: uniform Hecke zero-free strip). The 142 manuscript states this clearly (abstract,
  `build/main.tex` l.56–57; `00-introduction.tex` l.21–24; `10-analytic.tex` l.9), and so does the abstract
  row at `CONTENTS.md` l.3436. Only the one-line family summary does not. Neither 029 nor 142 has a
  `lean/docs` page. This is a transparency suggestion, not a correction.
- **Suggested minimal diff.**
  ```diff
  -... or GRH assumption is required.
  +... or GRH assumption is required; the proof uses the Hecke zero-free theorem of family 029 (Theorem 1.2).
  ```

## 3. Wording in manuscripts (W1–W3)

### W1 — 103, Lemma `lem:env-ranks` (3): "increases by at most 1/q_fr"

- **Where.** `…/build/sections/environment.tex` l.90–94.
- **Issue.** "Conditioning on frame success increases this probability by at most $1/q_{\rm fr}$" can be
  read additively. The intended and correct step is multiplicative: Pr[E | F] ≤ Pr[E]/Pr[F] = Pr[E]/q_fr
  (q_fr is the exact frame probability, l.68–70), which gives C_var = 4/q_fr in the form C_var/(qπ_a)
  required at l.56–58. The additive reading would give a vacuous bound. The proof of assertion (2)
  (l.84–86) already uses the multiplicative form. The main theorem is Lean-checked (`lean/docs/103.md`).
- **Suggested minimal diff.**
  ```diff
  -by $4/(q\pi_a)$. Conditioning on frame success increases this probability
  -by at most $1/q_{\rm fr}$.
  +by $4/(q\pi_a)$. Conditioning on frame success multiplies this probability
  +by at most $1/q_{\rm fr}$.
  ```

### W2 — 140, Lemma `s:mixture`: the proof gives a sharper constant (optional)

- **Where.** `…/build/sections/07-beta-mixture.tex` l.21–48 (statement: Law(v) ≤ 2^{D'}(ν₀+ν₁)),
  l.86–103 (proof).
- **Observation.** The proof bounds the density ratios by 2^{α*−α} = 2^{(b−1)/2} for q ≥ 1/2 and by
  2^{β₁−β₀} = 2^{(b−2)/2} for q < 1/2, with b = D'+1−n, and then weakens both to 2^{D'}. So the same proof
  gives Law(v) ≤ 2^{(D'−n)/2}(ν₀+ν₁). The stated bound is correct. Its only use (l.253–256, "Since D' ≤ d
  … ≤ 2^{3d}U₀") remains valid either way, and the Lean statement (`MixtureEstimate`,
  `lean/ComparatorChallenges/SubsphereCurrent.lean` l.404) uses the same 2^{finrank E} form. Nothing needs
  to change. We mention it only because the lemma is the place where an explicit constant is lost, so
  anyone revisiting c_* = 2^{−16} may find it useful. We do **not** claim any improved c_*.
- **Suggested (optional) diff.**
  ```diff
  - \Law(v)\le2^{D'}(\nu_0+\nu_1).
  + \Law(v)\le2^{(D'-n)/2}(\nu_0+\nu_1)\le2^{D'}(\nu_0+\nu_1).
  ```
- **Check.** `scripts/140_mixture_constant.py` evaluates sup_q f_{α,β₀}/(f_{α*,β₀}+f_{α,β₁}) on a grid for
  all 9 880 triples 1 ≤ j < n < D' ≤ 40: no violation of the 2^{(b−1)/2} bound. The adversarial reviewer
  re-derived each step and found the same with mpmath at 40 digits for D' ≤ 15.

### W3 — 099: ambiguous "The paper's" in `lean/docs/099.md` (optional)

- **Where.** `lean/docs/099.md` l.19: "The selected statements are the binary lower bounds. The paper's
  constant-distortion binary conversion and the companion upper embedding theorem are outside them."
- **Issue.** On a page that lists three papers and states at l.15 that the upper embedding *is* formalized,
  this sentence can be misread as contradicting l.15. It is not a contradiction: it reuses the abstract of
  the *Tree constructions* paper and refers only to the two binary lower-bound comparators
  (`BinaryEditLower.lean`, `TreeEdit.lean`). The upper bound is formalized via `EditDistance.lean` and
  `FiniteCircle.lean`.
- **Suggested minimal diff.**
  ```diff
  -The selected statements are the binary lower bounds. The paper's constant-distortion binary conversion and the companion upper embedding theorem are outside them.
  +The binary and tree comparator statements are the binary lower bounds. The tree paper's constant-distortion binary conversion is outside them; the upper embedding is covered by the two comparators above.
  ```

## 4. Withdrawn or not included

- **099 as a "contradiction" (internal list).** Withdrawn after adversarial review; kept only as the wording
  note W3.
- **142 (narrowed).** The CONTENTS abstract row (l.3436) already names the dependency; D3 now concerns
  only the one-line family summary.
- **107: "the Nine-Fourths README has no Lean link."** Dropped. No manuscript README in the repository
  links to `lean/docs`, so this is not specific to 107.
- **238, `main.tex` l.386–388, ε_n = O(n^{−5/2}).** The same derivation gives O(n^{−3.035})
  (198·log₂(32/31) ≈ 9.07 > 8). This is a slack in a constant, not a documentation issue. It is left to a
  separate constants remark on family 238.
- Families 205 and 269 are handled in other remarks and are deliberately not included here.
- None of the items was found to be already fixed upstream (checked at `fd4aeeb2e`, the current `main`).

## Reproduce

```bash
git clone https://github.com/openai/math openai-math
cd scripts
OPENAI_MATH=/path/to/openai-math ./run_all.sh              # ~15 min CPU, single thread
OPENAI_MATH=/path/to/openai-math SKIP_LONG=1 ./run_all.sh  # ~10 s, skips the two long re-runs
```
Requirements: bash, git, python3 (tested 3.10.12) with numpy (2.2.6), mpmath (1.3.0), sympy (1.14.0).
The clone is read only through `git show`, `git cat-file` and `git archive`; its working tree is never touched.
`run_all.sh` prints `OK` for each of the 10 outputs in `scripts/expected_output/`:

- `check_upstream.sh` — prints every cited line at the stated line numbers and counts files with each
  cited name in the whole tree (the main check for P1–P7, D1–D3, W1–W3).
- `independent/indep_check_paths.py` — separate check with `git cat-file -e` over the commit and its history.
- `238_network_bounds_recheck.py`, `independent/indep_238.py` — P5.
- `318_ce_ranks.py`, `independent/indep_318.py` — P3.
- `140_mixture_constant.py` — W2 (illustration only).
- `extract_222_certificate.py` — P4: extracts the Appendix A listing and checks its SHA-256 against the paper;
  then (long) runs it.
- (long) the shipped `verification/run.py` of the Exact-Fourier manuscript — P7.

## What is NOT claimed

- No mathematical error in any of the manuscripts. No claim that any theorem is false or unproved.
- For P1, P2 and P6 we checked only that the files are absent. We did **not** re-implement the N(6) = 3
  computation, the F₇ rank certificate of 095, or the 227 interval certificate.
- For P3, P4, P5 and P7 our runs reproduce the finite computations as stated; they do not check the
  surrounding arguments.
- W2 is an observation about a constant inside a proof; no improved final constant is claimed.
- Line numbers refer to commit `fd4aeeb2e`; they will move if the files change.
- Duplicate search (to 2026-10-09; GitHub only, arXiv not relevant for packaging issues): issues and
  discussions are disabled on openai/math, and its only pull request (#1, the October 7 update) does not
  touch these items. In mahrud/openai-math-comments (one pull request per family) the threads for all
  families listed here had no general comments. GitHub repository searches ("openai math", "openai math
  audit/errata/review", and the missing file names) found no other report of these items. Related community
  checks of other families exist (e.g. sunnyspot114514/openai-math-audit, piiiico/openai-math-counterexamples);
  they do not overlap with this list. The search was not exhaustive; GitHub code search hit a rate limit
  for some file names.
