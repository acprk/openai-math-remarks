# R12. Pixton-completeness counterexample: genus 10^60 → 2.2633·10^14

**What this remark shows.** OpenAI's family-053 preprint disproves Pixton's original completeness conjecture, in Chow and in rational cohomology. It does so with an explicit tautological class Y on M̄_{g,n}, at g = 10^60 and n = 3·C(g,3). In the preprint the genus enters **only** through the finite "colour budget" of Section 8 (`08-support.tex:326-510`). We transcribe each inequality of that budget exactly as an integer condition and optimize the free parameters (K, m₀). The same argument then closes at

- **g = 226,327,553,857,834 ≈ 2.2633·10^14** (K = 568, m₀ = 49,723,541,408), keeping the preprint's reservation of 4K+20 colour blocks, and
- g = 93,824,084,405,552 ≈ 9.38·10^13 (K = 567), if one reserves only the 3K+5 blocks that the preprint itself shows are needed (`08-support.tex:448-458`).

Both values are the exact least g admitted by the transcribed conditions. The number of markings drops from n ≈ 5.0·10^179 to n ≈ 5.8·10^42.

**The argument is OpenAI's; we only re-check the budget.** Everything below assumes that the preprint is correct (*premise: the original manuscript is correct*). We did not re-prove any lemma. We only redo the arithmetic that fixes g.

## Statement

| | g | K | m₀ | reserved blocks | n = 3·C(g,3) |
|---|---|---|---|---|---|
| **Before** (preprint, `08-support.tex:10, 331`) | 10^60 | 10^5 | 2000·g^{2/3} = 2·10^43 | 4K+20 | ≈ 5.0·10^179 |
| **After** (this remark) | **226,327,553,857,834** | 568 | 49,723,541,408 | 4K+20 = 2292 | ≈ 5.80·10^42 |
| After, reserving 3K+5 | 93,824,084,405,552 | 567 | 27,644,793,412 | 3K+5 = 1706 | ≈ 4.13·10^41 |

The genus shrinks by a factor of ≈ 4.4·10^45. The preprint says explicitly that its parameters "leave ample room … they are not proposed as minimal parameters" (`01-introduction.tex:138-140`). This remark measures how much room there is.

**Where the room was.** Two constraints compete as g decreases.
(iii) The large-vertex test needs d_e/2 ≥ loss + t_b with t_b = 25(b−1)² − 1 ≈ 25g². This forces m₀ ≈ 2·(300g²)^{1/3} ≈ 13.4·g^{2/3}; the preprint takes 2000·g^{2/3}.
(i) The (4K+20) reserved blocks must leave |S₀| > 280g/(K−4) colours for the large-deficit contradiction (L501-502). The preprint instead uses the cruder |S₀| > 0.9g together with 0.45(K−4) > 140.
Combining the two gives g^{1/3} ≳ 2·300^{1/3}(4K+20)/(1 − 280/(K−4)). This is minimized near K ≈ 570, where it gives g ≈ 2.26·10^14. Here |S₀| ≈ 0.5g, and C_K + 1 < 11g is far from binding. Because t_b ≈ 25g² cannot be reduced inside this framework, the framework cannot go much below ~10^13–10^14.

## Attribution

The class Y, the detector, the surface bound, the support induction and the colour-budget **inequalities themselves are OpenAI's work**. They come from family 053 of [github.com/openai/math](https://github.com/openai/math):

- [`preprints/A-high-arity-counterexample-to-Pixton-completeness-in-Chow-September-24-2026`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/A-high-arity-counterexample-to-Pixton-completeness-in-Chow-September-24-2026) (`build/sections/08-support.tex`, `06-coefficient.tex` Lemma `lem:large-vertex-test`)

The path is pinned at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The directory has not changed since the initial commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.

This remark adds the integer transcription of the budget, the optimized parameters, and the check that g enters nowhere else.

## What changes in the text

The only edits are numerical, all in `08-support.tex` (plus the displayed g in `01-introduction.tex:79`, `02-strata.tex:196`, `08-support.tex:412, 561`, `main.tex:65-66`):

| Location | Preprint | Replace by | Condition (exact integer form) |
|---|---|---|---|
| L10 | g = 10^60 | g = 226,327,553,857,834 | — |
| L331 | K = 10^5, m₀ = 2000g^{2/3} | K = 568, m₀ = 49,723,541,408 | — |
| L337-340 | (4K+20)m₀/g < 0.1 | \|S₀\| = g − 2292·m₀ = 112,361,196,950,698 | used only in L501-502 |
| L402-413 | 3gM_K + g(10K+2) + 2KC_K < 10^{17}g < 0.01D₀ | together with L425-435: loss + t_b ≤ d_e/2, e = 2⌊m₀/4⌋ | P5 |
| L421-422 | 3g−2 < (b−1)² − 1 | (g−K)² − 1 > 3g − 2 | P4 |
| L455 | 4(g−K−2) > 3g−3 | same | P3 |
| L470 | C_K + 1 < 11g | same | P6 |
| L501-502 | (0.45(K−4) − 140)g² > 0 | \|S₀\|(g(K−4)+2) > 280g² | P7 |
| L370 | b > g−K > g/2 | g > 2K | P2 |
| Lemma `lem:large-vertex-test` | b ≥ 5, e ∈ [m₀/2−2, m₀/2] | b > g−K ≥ 5; e = 2⌊m₀/4⌋ | — |

At the new point every slack is positive (`expected_output/verify_points.txt`).

## Why g enters nowhere else

We read Sections 1, 2, 7 and 8 in full and grepped every section for `sufficiently large`, `10^{60}`, `g\ge`, `g>`, `D_*`, `3g-3`. Outside Section 8, the hypotheses on g are only g ≥ 2 (`03-detector.tex:224, 498, 686`; `07-compactification.tex:132`), the stability and degree bound 2D_* < 3g−3+n (`02-strata.tex:233-236`, true for every g ≥ 2), and n = 3·C(g,3) by definition. The one "sufficiently large" (`06-coefficient.tex:474`) concerns multiples of a line bundle, not g. The genus-b hypotheses of the surface bound (b ≥ 5, dim S < (b−1)² − 1, `05-quot.tex`, `06-coefficient.tex:588-597`) are exactly P4 and P2 above.

## Reproduce

```bash
bash scripts/run_all.sh        # < 1 s: checks the three points (paper, 4K+20, 3K+5) inequality by inequality
bash scripts/run_all.sh full   # + exhaustive least-g search (≈ 6 min) + the primary session's script (≈ 10-20 min)
```

Python 3 stdlib only, single process. Outputs are diffed against `scripts/expected_output/`, and `scripts/SHA256SUMS` lists every file.

The note is `note.tex` / `note.pdf` (`latexmk -pdf note.tex`).

## Scope and limitations (what is NOT claimed)

- **No new mathematics.** All inequalities and all lemmas are OpenAI's. We change only the numbers that close them.
- **Premise: the original manuscript is correct.** We did not re-verify the detector, the Quot-scheme surface bound, the controlled lifts, or the support-advancement lemma.
- **Minimal only within the transcription.** "Least g" means least under the preprint's own Section 8 bounds (3g ordinary advancements, 11g cost, 140g², t_b = 25(b−1)²−1, …). Tightening those bounds, for example the 3g factor or the 99g² + 30g² allowance, is not attempted here, and neither is a different argument.
- **The 3K+5 variant** uses the count that the preprint itself proves (L448-458). The preprint reserves 4K+20 only as slack, so we headline the 4K+20 value.
- **Still not computable.** g ≈ 2·10^14 and n ≈ 6·10^42 remain far beyond any explicit computation in the strata algebra.
- **Novelty limited to our search.** To our knowledge (search of GitHub and arXiv up to 2026-10-09), no smaller genus for this counterexample has been published. `piiiico/openai-math-counterexamples` lists family 053 as "not-finite" and does not re-check the budget.

## Authorship

[AUTHOR NAME].

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
