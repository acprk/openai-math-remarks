# Verification record: R12 (colour budget of OpenAI family 053, genus 10^60 → 2.2633·10^14)

Date: 2026-10-09. Host: Linux x86_64, Python 3 (stdlib), single process. No Lean was run.

Pins: openai/math `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The family-053 directory has not changed since `adc7f124`, and the local checkout is clean.

## Who did what

| Role | Session | What |
|---|---|---|
| Primary computation | R7-skill rescan s0 (`reviews/R7-skill/rescan/s0/pixton_budget.py`) | Transcribed the Section 8 budget. Binary search over g for K ∈ {290, 300, …} with local K-refinement. Result: 2.263·10^14 (4K+20) and 9.38·10^13 (3K+5). |
| Main agent | — | Re-ran the integer test. |
| Adversarial review + independent re-implementation | this session (separate AI session) | Line-by-line reading of `08-support.tex` L326-510, Lemma `lem:large-vertex-test` / Lemma `lem:exterior-ideal` / Theorem `thm:surface` hypotheses, and Sections 1, 2, 7. New code `independent/indep_budget_053.py`, written from the tex. It uses an exhaustive K range and exact least-g certification. |

## (a) Adversarial reading

| # | Question | Finding |
|---|---|---|
| A1 | Does g enter outside Section 8? | Grep of all sections for `sufficiently large`, `10^{60}`, `g\ge`, `g>`, `D_*`, `3g-3`, followed by reading the hits. Outside §8 we find only g ≥ 2 (`03-detector.tex:224, 498, 686`; `07-compactification.tex:132`), deg Y = 2D_* < 3g−3+n (`02-strata.tex:235`, true for all g ≥ 2), and the definitions n = 3C(g,3) (`01-introduction.tex:79`, `02-strata.tex:196`). The single "sufficiently large" (`06-coefficient.tex:474`) concerns multiples of a line bundle K. The detector (§3) is valid for any g ≥ 3 (one needs triples of colours). **No hidden "g large".** |
| A2 | Is each inequality transcribed faithfully? | P2: uniqueness of the large vertex needs b > g/2 and is implied by b > g−K with g ≥ 2K (L370). P3: L455. P4: L421-422 with b ≥ g−K+1. P5: eq. small-loss (L402) plus eq. large-test-threshold (L433). The exact requirement is D₀ − loss ≥ D₀ − d_e/2 + t_b, i.e. loss + t_b ≤ d_e/2. The loss uses 3g (not 3g−2) advancements times M_K, exactly as written. t_b is bounded by 25(g−1)² − 1 since b ≤ g. P6: L470. P7: L498-502, exact difference \|S₀\|/2·(g(l−4)+2) − 140g², monotone in l, at l = K. P8: L448-458. e = 2⌊m₀/4⌋ satisfies m₀/2 − 2 ≤ e ≤ m₀/2 and is even, as Lemma `lem:large-vertex-test` requires. |
| A3 | Do the preprint's auxiliary numerics (m₀ ≥ 128, d_e > D₀/10, \|S₀\| > 0.9g, C_K < 1.334·10^15, "These are integers") impose anything extra? | No. They are used only to *derive* P5–P7 at the preprint's parameters (L336-340, L408-413, L436-440). We check P5–P7 directly. The new K and m₀ are integers. |
| A4 | Lemma `lem:large-vertex-test` hypotheses | b ≥ 5 (from b > g−K and P2/P3); dim S < (b−1)² − 1 (P4); t_b ≤ d_e/2 (implied by P5); e within [m/2−2, m/2] and even (by construction). Lemma `lem:exterior-ideal` needs 0 ≤ t ≤ r with 2r = d_e, which again is t_b ≤ d_e/2. |
| A5 | 3K+5 vs 4K+20 | The preprint proves that at most 3K+5 blocks are used (L448-458) and reserves 4K+20 as slack. The 3K+5 variant is therefore justified by the preprint's own count, but we headline the 4K+20 figure. |
| A6 | Do costs hold at small K? | The preprint states M_K = C_K + 1 (L389). This holds for K ≥ 4. P6 keeps "ordinary cost ≤ 11g" valid. The noise-test cost 10g+1 < 11g and h_v + 2 < 11g hold trivially. |

## (b) Independent computations

| Check | Command | Result |
|---|---|---|
| Positive control (paper point) | `independent/verify_points.py` | g = 10^60, K = 10^5, m₀ = 2·10^43, 4K+20 blocks: all of P1–P8 hold. The paper's sub-claims (C_K < 1.334·10^15, C(m₀,3) > 10^9 g², ratio 8.0004·10^{-12}) are reproduced. |
| New points, condition by condition | same | 4K+20: g = 226,327,553,857,834, K = 568, m₀ = 49,723,541,408, \|S₀\| = 0.496g. Every slack is positive: P5 slack 1.0·10^20, P7 slack 3.5·10^16. g−1 is infeasible with the optimal m₀. 3K+5: g = 93,824,084,405,552, K = 567; P8 is tight by design (blocks = 3K+5). |
| **Exhaustive least-g search** | `independent/indep_budget_053.py` (≈ 5.5 min) | For every K from 285 up to the point where P6 alone exceeds the optimum, the exact least g is found by integer bisection on j (m₀ = 4j) for the single sign change of G(j) − g_min(j), with a local certificate over 200 j below. Optimum: **226,327,553,857,834 at K = 568** (4K+20) and **93,824,084,405,552 at K = 567** (3K+5). |
| Cross-check against the primary script | `primary/pixton_budget.py` | 3K+5: identical value. 4K+20: the primary reports 226,327,553,876,041, which is 18,207 larger. Its plain bisection in g assumes monotonicity, but feasibility is not monotone in g because m₀ jumps by 4. Both points are feasible. The independent value is the least one. Relative difference 8·10^{-11}. |

The independent code shares no code with `primary/pixton_budget.py`. It shares only the reading of the tex.

## (c) Novelty / dedup (up to 2026-10-09)

- GitHub code search: "Pixton-completeness", "Pixton genus 10^60", "m_0=2000g". Hits are indexes and catalogues only (blisspixel/numinous, satyajitghana/ai, anghel4d/broadside-observer, jcraig949jfi/Prometheus, AmEskandari/ai4math). `piiiico/openai-math-counterexamples` marks 053 "not-finite" with no budget check. openai/math has issues disabled, and no related PR exists.
- arXiv/web: Canning–Larson–Schmitt (compact type, 2026) and the Gorenstein results concern small genus. We found no improvement of the 053 genus. A blog summary (kingy.ai) only lists the family.

## Verdict

**CONFIRMED** (premise: the preprint's proof is correct). The genus enters only through the §8 budget. Under the exact transcription the argument closes at g = 226,327,553,857,834 with the preprint's 4K+20 reservation, and at 93,824,084,405,552 with 3K+5.

## GAPs (not covered)

- No lemma of the preprint was re-proved. The detector, the surface bound (Quot schemes), the controlled lifts and support advancement are assumed.
- "Least g" is relative to the preprint's own bounds. The exhaustive search relies on the single-crossing property of G(j) − g_min(j), which is argued from growth rates (j^{3/2} vs. linear) and checked locally. It is not formally proved.
- `primary/pixton_budget.py` was re-run in this session; its output is in `expected_output/primary_pixton_budget.txt`.
- No Lean formalization.
