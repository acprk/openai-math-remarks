# R7: OpenAI's height-two Hahn–Wilson counterexample works at every prime p ≥ 293

**One sentence.** OpenAI's paper *Counterexamples to the Hahn–Wilson conjecture at height two* states its theorem for "every sufficiently large prime". The only explicit sufficient condition it gives, (op:parameters), first holds at p = 17,073,437. If we choose the central window of the digit construction as large as the proof allows (R₀ = (M−3)/2 instead of ⌊M/8⌋) and restate the one displayed bound that depended on the old choice, the same argument works for **every prime p ≥ 293**.

> **The argument is OpenAI's; we only re-select parameters. No inequality had to be replaced by exact enumeration here: every change is a symbolic restatement in terms of the new parameter. We list every changed line below.**

## Precise statement (Before / After)

| | Before (OpenAI, as written) | After (this remark) |
|---|---|---|
| Theorem 1.1 range | "every sufficiently large prime p" | every prime p ≥ 293 |
| Parameter condition | (op:parameters), `operator.tex` L264–268: M odd, M ≥ 81, p > 16M², 2⌊M/8⌋ > p^{1/3} | M odd, M ≥ 11, R₀ = (M−3)/2, and (P1)–(P4) below |
| Smallest prime the condition admits | 17,073,437 (with M = 1033, R₀ = 129) | 293 (with M = 11, R₀ = 4, B₀ = 8, β = log₂₉₃ 8 ≈ 0.366) |
| Central window | R₀ = ⌊M/8⌋, B₀ = 2R₀ | R₀ = (M−3)/2, B₀ = M−3, I = {2,…,M−2} |

The new conditions, with q = p+1, are:

* **(P1)** (p+1)/(2M) − p(M+1)/(p−1) > 1. This is inequality (op:correction-bound) with the paper's own bound κ ≤ p/(p−1). It is the binding condition.
* **(P2)** (M−3)³ > p, i.e. β = log_p B₀ > 1/3 (`operator.tex` L503).
* **(P3)** p > 4M. This separates the first band from the second and gives b ≥ X/(2M) on the first band.
* **(P4)** ((M−3)p/M + 3)(p+1) ≤ p². This keeps same-cluster differences outside the stage and gives |lower portion| ≤ X/p.

Choice of M for a given p: take the least odd M ≥ 11 with (M−3)³ > p. That is M = 11 for 293 ≤ p < 512, M = 13 for 512 ≤ p < 1000, and so on. For M = 11, (P1) holds exactly when p ≥ 293. For M ≥ 13, minimality of M gives p ≥ (M−5)³ ≥ 2M(M+3), and this implies (P1). (P3) and (P4) are then immediate. The script checks every prime 293 ≤ p ≤ 2·10⁶ explicitly.

## Attribution

The construction, the theorem and its proof are OpenAI's work: family 319, preprint `preprints/Counterexamples-to-the-Hahn-Wilson-conjecture-at-height-two-September-26-2026/` in https://github.com/openai/math (local snapshot at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`; all line numbers below refer to `build/sections/*.tex` and `build/main.tex` at that commit). This remark only re-selects the parameters R₀, B₀, M of Proposition `op:operator` and checks that every step of the paper still goes through.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Every changed line

| File : lines | Original | Replacement | Why it still holds |
|---|---|---|---|
| `operator.tex` L264–268 (op:parameters) | M ≥ 81, p > 16M², 2⌊M/8⌋ > p^{1/3} | M odd, M ≥ 11, (P1)–(P4) | (P1)–(P4) are exactly the inequalities used below |
| `operator.tex` L269–273 | "hold at every sufficiently large prime … M ~ p^{2/5}" | "hold for every prime p ≥ 293 with M the least odd integer ≥ 11 such that (M−3)³ > p" | see "Choice of M" above |
| `operator.tex` L324–329 (op:correction-bound) | cites (op:parameters) | cites (P1) | the displayed inequality is (P1) verbatim |
| `operator.tex` L500 | R₀ = ⌊M/8⌋ | R₀ = (M−3)/2 | — |
| `operator.tex` L502–503 | B₀ = 2R₀, β = log_p B₀ > 1/3 | B₀ = M−3, β > 1/3 by (P2) | (P2) |
| `operator.tex` L507 | "I meets neither endpoint 0, M−1" | unchanged | I = {2,…,M−2} |
| `operator.tex` L592 | (p/4+3)X/p² | (2R₀p/M+3)X/p² = ((M−3)p/M+3)X/p² | the digit diameter of I∪(I−1) is ≤ 2R₀p/M+1; the old p/4 was this bound with 2R₀ ≤ M/4 |
| `operator.tex` L597–603 | "The bounds on p separate …" | cite (P3), (P4) | (P4) puts same-cluster differences below X/q and keeps the deviation ≤ 2/p; (P3) separates bands and gives b ≥ X/(2M) |
| `operator.tex` L694 | "inside I since M ≥ 81" | "inside I since R₀ ≥ 4" | the reached clusters are within three places of s₀ (in fact {s₀, s₀+1}, because p > 3M) |
| `consequences.tex` L14 | "a prime satisfying the operator bounds in (op:parameters)" | same, with the new (op:parameters) | — |
| `consequences.tex` L172–173, L187; `main.tex` L57–58, L115, L129 | "sufficiently large prime(s)" | "every prime p ≥ 293" | the only large-prime input is (op:parameters) |

No other line depends on the size of p. Outside `operator.tex` the paper requires only p > 5 (`geometry.tex` L14, `completion.tex` L16), p > 8 / D > 12 (`comparison.tex` L19, L25, L41; `transfer.tex` L21), D > 2 or D > 3 (`gluing.tex` L157, `completion.tex` L136, L171, L285), and q prime to p. All of these hold at p = 293.

## Reproduce

From `scripts/` (Python 3.10, only `numpy` beyond the standard library; no SageMath needed):

```
/usr/bin/python3 main/hw319_pmin.py              # ~10 s   -> expected_output/hw319_pmin.txt   (293 literal; 199 sharpened; original 17073437)
/usr/bin/python3 main/hw319_detail.py            # ~2 s    -> expected_output/hw319_detail.txt
/usr/bin/python3 -u independent/hw319_enumerate.py  # ~3-6 min; label-level exact enumeration, coverage to 2e6, controls -> expected_output/hw319_enumerate.txt (exit 0)
/usr/bin/python3 -u independent/indep_scan.py       # ~7 min, separate re-derivation with exact centres -> expected_output/indep_scan.txt
```

## What is NOT claimed

* We did not re-verify the rest of OpenAI's paper. The remark is only as correct as the original proof. We checked that no step outside `operator.tex` uses the size of p, but we did not check the correctness of those steps.
* The prime ranges of the external theorems the paper cites (Pstrągowski's algebraicity, p > 8 as quoted; Lee–Pstrągowski 2026, Theorems 1.3/5.4; Lee 2026, Corollary 1.9) were taken from the paper as quoted. We did not check them against the sources. The paper states them for p > 8 or "odd primes".
* We make no claim for p < 293. With R₀ ≥ 4, which the literal "three places" sentence at L693–694 needs, the mechanism fails for p < 293 at (P1). Two stronger variants appear in the scripts but are **not** part of this remark: (a) a sharper reading of the cut argument at L689–694 (R₀ ≥ 2 suffices) gives 199; (b) restating Lemma `op:coefficient-control` with the exact first-band lower end in place of X/(2M) gives 167/127/109. Both need extra rewriting of OpenAI's text that we have not had adversarially reviewed.
* Novelty: the earlier search session searched GitHub and arXiv (up to 2026-10-09) and found no explicit prime bound for this construction. To our knowledge, no one has stated one.
