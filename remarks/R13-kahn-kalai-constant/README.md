# R13: family 176, the constant in the second Kahn–Kalai bound

**In short.** OpenAI proves p_c(n,H) ≤ min{1, 2048e^50 · p_E(n,H)(1+log₂h)} and p_c ≤ 6144e^50 · p_E log₂n. The same proof, with its numerical choices re-tuned and three elementary estimates tightened, gives **p_c ≤ min{1, 10⁵ · p_E(1+log₂h)}** and **p_c ≤ 2·10⁵ · p_E log₂n**. Both constants improve by about 10²⁰ (1.06·10²⁰ and 1.59·10²⁰).

## Attribution
The theorem, construction and proof are OpenAI's work: family 176, *The second Kahn–Kalai conjecture* (September 24, 2026), [openai/math `preprints/The-second-Kahn-Kalai-conjecture-September-24-2026`](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-second-Kahn-Kalai-conjecture-September-24-2026), commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. This remark only does the bookkeeping. It treats the manuscript's numbers as parameters, sharpens three elementary estimates, and certifies the resulting constants by computer.

## Before / After
| Site (file, lines) | Manuscript | This remark |
|---|---|---|
| margin a/ρ (`resampling.tex` L63) | e^{-50} | e^{-L}, L = 123/20 |
| denominator cut Z < e^{-αm} (L68) | α = 10 | α = 7/5 (exact sums kept, no rounding to e^{-9m}, e^{11m}) |
| pruning ε = e^{-βm} (`trees.tex` L98) | β = 1 | β = 2/5 |
| potential sum (L157–168) | Σ_{h<j} m_h ≤ m_j/15, sum < 1/4 | Σ_{h<j} m_h ≤ Σ_i ⌊m_j/6^i⌋; B ≤ 0.39591 |
| capacity growth / hierarchy (`trees.tex` L34, `graphs.tex` L67–72) | 16 / 17 | g = 6 / 7 |
| spread multiplier σ/p_E (`graphs.tex` L9, L64) | 128 | λ = 7⁷/6⁶ ≈ 17.651 (2·C(7r,r) < λ^r, Robbins) |
| spread inflation (`trees.tex` L244–252) | product over all m ≥ 1, < 4 | product over m = 2^v (Lemma 4): 7.207 at β = 0.4 |
| attempts (`trees.tex` L257–277) | 4s via Markov, probability 2/3 | exact min N with P(Bin(N, 0.604) ≥ s) ≥ 1/2; sup N(s)/s = 5/3 |
| Theorem 1.1 (`introduction.tex` L37, L41) | 2048e^50 ≈ 1.06·10²⁵; 6144e^50 ≈ 3.19·10²⁵ | **10⁵** (certified ≤ 99 379); **2·10⁵** (certified ≤ 197 985) |

**Why 1/2 instead of 2/3 is legitimate.** p_c is defined as inf{p : P(N(G(n,p),H) > 0) ≥ 1/2} (`introduction.tex` L17–19). The Lean comparator uses the same 1/2 (`criticalThreshold` in `lean/ComparatorChallenges/SecondKahnKalai.lean`). The manuscript's tree theorem proves 2/3, which is more than the main theorem needs. If one keeps 2/3 anyway, a rigorous Chernoff/Hoeffding tail gives C₁ ≤ 119 255 and C₂ ≤ 212 769.

### Where the gain comes from (C₁)
| factor | manuscript | here | ratio |
|---|---|---|---|
| margin e^L | e^50 | e^6.15 | 1.10·10¹⁹ |
| σ/p_E | 128 | 17.651 | 7.25 |
| inflation K | 4 | 7.207 | 0.555 |
| attempts per success | 4 | 5/3 | 2.40 |
| **total** | 1.062·10²⁵ | 9.94·10⁴ | **1.07·10²⁰** |

Almost all of the gain comes from the margin e^{-50} in the local resampling lemma. The proof only needs B < 1, and L = 6.15 is enough for that.

## Proof sketch
See `note.tex` for the full proof. Lemmas 1–2 are the manuscript's Lemma 2.1 and Lemma 3.2 with the numbers replaced by (L, α, β, g). The success probability of one reduction is at least 1 − B, where B = sup_m φ₀(m) + Σ_{j≥2} sup_{m≥g^{j−1}} φ₀(m)e^{cS(m)}. Lemma 3 (spread inflation) uses ⌊y⌋ ≥ 2^v⌊y/2^v⌋, so the capacities met by one level form a dyadic chain. Lemma 4 (attempts) uses induction to show that the successes dominate Bin(N, p), and bounds s > 200 by the Kaas–Buhrman median theorem. The hierarchy lemma needs 2·C(h,r)λ^{−r} < 1 for r ≥ h/7. This follows from C(7r,r) ≤ √(7/(12πr)) · λ^r (Robbins).

## Reproduce
Requirements: Python 3 with mpmath (primary script), SageMath 9.5 (independent script). Each run is single-threaded.
```
cd scripts
/usr/bin/python3 primary/kk176_verify.py 6.15 1.4 0.4 6               # primary: C1<=1.0018e5, C2<=2.033e5   (~15 s)
/usr/bin/python3 primary/kk176_verify.py 50 10 1 16 2/3               # primary positive control: B=1.234e-4 (~10 s)
sage -python independent/kk176_indep.py 6.15 1.4 0.4 6                # independent: C1<=99 379, C2<=197 985 (~20 s)
sage -python independent/kk176_indep.py 6.15 1.4 0.4 6 2/3            # independent, success prob 2/3       (~2 min)
sage -python independent/kk176_indep.py 50 10 1 16                    # independent control at paper params  (~20 s)
```
Expected outputs are in `scripts/expected_output/`.

## What is NOT claimed
- This is not a new method. Every change is a parameter choice or an elementary sharpening (dyadic product, Robbins' bound, exact binomial count). Following the project's rule, nothing was formalised in Lean, and the Lean comparator statement still hard-codes 2048e^50 and 6144e^50.
- The parameters were found by a coarse grid search. They are not optimal, and we make no claim about the best constant this method can give.
- We make no statement about the order of magnitude of the true constant in the second Kahn–Kalai conjecture.

## Premise
This remark assumes the manuscript is correct, in particular the posterior-resampling identity (`resampling.tex` L38–44), the outside-label independence (`trees.tex` L124–136) and the graph hierarchy (`graphs.tex` L50–112). We re-read each of these line by line and they carry over with our parameters (see VERIFICATION.md), but we did not re-prove them independently.

## Prior work
To our knowledge (search of GitHub and arXiv up to 2026-10-09; see VERIFICATION.md §4), no explicit constant for this bound has been published other than OpenAI's. Tran (arXiv:2609.20546) and Dubroff–Kahn–Park (arXiv:2508.14269, 2510.02005) prove related bounds with unspecified absolute constants.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
