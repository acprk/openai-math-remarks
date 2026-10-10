# VERIFICATION — R24 (family 179, norm-m form; parallel generalization)

All runs were on 2026-10-09. `scripts/run_all.sh` reproduces both expected outputs byte for byte
(`OK perfect_search`, `OK klein_mitm`; about 1 minute, single thread).

## 1. Who checked what

| Item | Primary | Independent check | Adversarial review | Status |
|---|---|---|---|---|
| Prop′, Prop″, Master Lemma (ℤ[i], ℤ[ω]), Q, T1 (both forms), T3, F | Lean 4 proofs (17 files, 4085 lines) on top of openai/math `fd4aeeb2e` | Build session: `lake build` of the four top modules plus `Checks.lean` (`#check` and `#print axioms` for 9 declarations), under `flock`, 8 cores; output in `scripts/lean/expected_output/axioms.txt` | Separate AI session read every statement against README/note. Verdict: faithful and not vacuous; no `sorry`/`axiom`/`native_decide`/`implemented_by`/`extern`/`unsafe`/`opaque`; 15 upstream imports, none modified; no name collisions | holds |
| Normalizations (Q: ÷(1+i); T3: ÷(1−ω), e = 2ω²; T1: ÷2, t = −2) | note §2, Lean docstrings | — | Adversarial reviewer re-derived each identity by hand (e.g. y₁+y₂ = 2ω² + λω²(W₁+W₂); b₁+b₂+b₃ ≡ J) | holds |
| Small cases of Q and T3; controls | — | `independent/perfect_search.c` (written for this remark): exhaustive over μ_q-valued functions with a(0) normalized. Q: 0 at lengths 6, 10, 14. T3: 0 at lengths 6, 12, 15. Controls: 2/8/32 quaternary at lengths 2/4/8, 6/54 three-phase at lengths 3/9 | Reviewer's own FFT brute force (unnormalized) gives exactly q × these counts | agree |
| T1 at u = 3; control | — | `independent/klein_mitm.py`: 90 norm-9 elements of {−2..2}^{C₉} → 0 arrays on C₂×C₂×C₉; 792 → 15552 arrays on C₂×C₂×C₃×C₃, one of them re-checked directly | Reviewer's own hash-join (different algorithm): 90/792 and 0/15552 | agree |
| Classical coverage remarks (T1 trivial for non-square N; Turyn exponent bound) | README/note | — | Reviewer's computation: Turyn's bound covers every prime power u and every odd u < 39; first uncovered 39, 55, 117, 155, 165, 195 | wording corrected accordingly |

## 2. Lean environment

- Base: `openai/math` commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (`lean/`). Toolchain
  `leanprover/lean4:v4.34.1`. Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612` (upstream manifest).
- Command (from the upstream `lean/` directory, with `scripts/lean/OAI/` copied in):
  `flock <lockfile> taskset -c 0-7 lake build <4 modules>` followed by
  `lake env lean Checks.lean`.
- Result: `Build completed successfully (8983 jobs)`. Every one of the 9 `#print axioms` lines lists
  `[propext, Classical.choice, Quot.sound]` (see `scripts/lean/expected_output/axioms.txt`).
  The only build output besides this line is linter warnings (unused section variables). Wall time was
  6084 s; nearly all of it went to I/O contention on the shared host while loading Mathlib oleans.
- After the build, two comment lines in `GeneralNorm/TheoremT1.lean` were reworded ("new case" →
  "non-trivial case"). Code is unchanged.

## 3. Adversarial review (separate AI session)

Verdict: **CONFIRMED-WITH-FIXES**. No mathematical or formal claim was refuted. Fixes applied:
1. Added VERIFICATION.md and the Lean expected output.
2. T1 classical coverage was understated; README and note now state the Turyn exponent-bound coverage.
3. "Not re-audited length by length" was replaced by "makes no claim about individual lengths", and
   Arasu–de Launey–Ma 2002 and Ma–Ng 2009 are now cited.
4. The formalization scope now says ℂ-valued (Q, T3) and ℤ-valued (T1).
5. The duplicate-search wording now says that openai/math has issues disabled.
6. The Lean branch of `run_all.sh` now fails loudly and asserts 9/9.
7. "Δ = 1 for at least two primes" wording fixed.
8. The note says that Prop″ needs the 3-adic localization layer redone.
9. The note says that F is vacuous unless p ≡ q ≡ 1 (mod 4).

The placeholders `R24` and `[AUTHOR NAME]` are left on purpose. The reviewer also built nothing. The Lean
build in §2 is the build session's.

## 4. Duplicate search (to 2026-10-09)

openai/math: issues are disabled, and the only PR is #1 (manuscript update). google-deepmind/formal-conjectures
#6896 (issue), #6897 and #6898 (PRs) only link OpenAI's binary/real theorem. Community follow-up
repositories (about 30 cloned): no text applying the 179 argument to other alphabets. Web/arXiv: nothing
found that states Q, T3 or T1 as consequences of the 179 argument.

## 5. GAPs

- The difference-set ↔ perfect-array equivalence is not formalized; it is standard.
- Theorem F has no numerical cross-check beyond Lean. Its smallest non-vacuous case (pq = 65) is too
  large for exhaustive search.
- The computations are unit tests at toy sizes. The infinite families rest on the Lean proofs.
