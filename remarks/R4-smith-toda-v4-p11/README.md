# R4: OpenAI's Smith–Toda complex V(4) exists at every prime p ≥ 11, given their p = 1009 proof

> ## ⚠ Precondition
> **This remark is conditional on the correctness of the original p = 1009 manuscript, which we did not re-verify.** We show only that, *if* OpenAI's proof at p = 1009 is correct, the same proof gives V(4) at every prime p ≥ 11. Its homotopy-theoretic steps have not been independently checked by us or, to our knowledge, by anyone else. Those steps are the adjusted-filtration injectivity in §3, the Frobenius quotient-group calculation in §4, the extension of powers in §5, and the pairing in the appendix.

**One sentence.** OpenAI's preprint *A finite Smith–Toda complex V(4) at the prime 1009* constructs V(4) at the single prime p = 1009, a value chosen "to make the finite inequalities generous". Every step of that proof is either valid for all primes p ≥ 19 as written (with 1009 replaced by p) or, at p = 11, 13 and 17, can be replaced by an exact enumeration of the relevant May E₁ pages. So, granting the p = 1009 manuscript, **V(4) exists at every prime p ≥ 11.**

> **The argument is OpenAI's; we only re-select parameters and replace 3 specific inequalities by exact enumeration. We list every changed line.**

## Precise statement (Before / After)

| | Before (OpenAI, as written) | After (this remark) |
|---|---|---|
| Theorem 1.1 | at p = 1009 there is a finite p-local X with BP_*X ≅ BP_*/(p, v₁, v₂, v₃, v₄) as comodules | the same, for **every prime p ≥ 11** |
| How the finite claims are proved | inequalities instantiated at p = 1009 | p ≥ 19: the same inequalities, now in symbolic form; p = 11, 13, 17: three of them fail and are replaced by exact enumeration of the full May E₁ page in every relevant bidegree |
| p = 7 | — | negative control: the enumeration finds a non-empty product-obstruction group at e = (0,2,0,2), k = 3. This agrees with Toda (V(3) has no suitable unital product for p < 11) and with Nave (V((p+1)/2) = V(4) does not exist at p = 7). |

Context: the classical constructions give V(1), V(2), V(3) for p > 2, 4, 6. Nave proved that V((p+1)/2) does not exist for p ≥ 7. In the literature we found (Ravenel's green book 5.6.13; Ausoni–Rognes, arXiv 2204.05890; Culver–Zhang, Remark 3.29), V(4) is not constructed at any prime. OpenAI's companion preprint *Finite Smith–Toda Complexes at Varying Primes* proves existence of V(n) at *some* prime depending on n, with no explicit bound. p = 5 and p = 3 are not addressed here.

## Attribution

The construction and proof are OpenAI's work: family 308, preprint `preprints/A-finite-Smith-Toda-complex-V4-at-the-prime-1009-September-23-2026/` in https://github.com/openai/math (local snapshot at commit `08fcd618213aa43b0fcaff4526a8a48d66d02e55`; line numbers refer to `build/sections/*.tex`). This remark only (i) replaces the numerical value 1009 by a general prime p and (ii) replaces three finite inequalities by machine enumeration at p = 11, 13, 17.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Every changed line

### (a) The three inequalities that fail at p = 11, 13 or 17, replaced by enumeration

| Line | Inequality as written | Holds from | Used to prove | Replacement certificate at p = 11, 13, 17 |
|---|---|---|---|---|
| `02-lower-complexes.tex` L109–113 | p² − (27/2)p − 9 > 0 (⇒ k ≤ 2) | p ≥ 17 | Lemma `low:algebraic-vanishing`: every lower obstruction bidegree has empty May E₁ page | check (A) of `may308.py`: every bidegree (qk+1, u_n+k), n ≤ 3, and (qk+1−ρ, U+k) for all profiles e ∈ {0,1,2}^{n+1} with ρ ≥ 2, up to the normalized cutoff, has **no** E₁ monomial. Cross-checked by DP count (`may308_dp.py`) and by z3 + DP (`independent/sweep.py`) |
| `02-lower-complexes.tex` L154–155 | [−14, 3] contains no nonzero multiple of p | p ≥ 17 | the same lemma (coverage relations) | the same check (A) |
| `03-may-and-first-obstruction.tex` L161–163 | 2p − 19 > p − 2 (single high block) | p ≥ 19 | exhaustion of Table `may:n5-table` at w = u_n+1, n = 4, 5 | check (C): the full list at (q+1, u₅+1) is exactly the five entries of Table `may:n5-table`, with exponents p−2, p−3, p−3, p−3, p−4, and d₁ is injective on the length-3 and length-5 spans. The k = 1 list at (q+1, u₄+1) is {acb B^{p−2}, acybz B^{p−3}}. Cross-checked by `independent/sweep.py` and `independent/d1check.py` |

Checks (B), (D), (E), (F) of `may308.py` re-verify, at each p, the other finite tables the proof relies on: the v₄-target list (`may:first-list`, `may:remaining-list`: non-empty exactly at k = 1, 2, p−1, p, with the paper's monomials); the Frobenius quotient table (`frob:quotient-groups`); the small-stem lemma (`top:small-stems`, d₁-cohomology 0 at (6, W+1)); and the page bound (`top:early-disappearance`).

### (b) Lines that instantiate p = 1009 numerically (replaced by the symbolic expression; each inequality checked for all primes p ≥ 11 by `ineq308.py`, apart from the three in (a))

| File | Lines | Content → replacement |
|---|---|---|
| `00-introduction.tex` | L44, L64–65 | "At p = 1009" → "At every prime p ≥ 11"; drop "the choice 1009 makes the inequalities generous" |
| `01-framework.tex` | L6, L10–11 | p = 1009, q = 2016 → p ≥ 11 prime, q = 2p−2; drop the primality check of 1009 |
| `01-framework.tex` | L263–265 | numeric u₄, 2u₄ < p⁴ → 2u₄ < p⁴ < u₅ (all p ≥ 3) |
| `02-lower-complexes.tex` | L40, L45, L55–67 (Table `low:cutoffs`), L83–86 | numeric cutoffs → formulas; max w < p³ (p ≥ 5); at p = 11, 13, 17 the enumeration in (a) covers every k up to the formula cutoff |
| `02-lower-complexes.tex` | L111, L142 | numeric differences → p² − (27/2)p − 9 (see (a)); p² − 4p − 8 > 0 (p ≥ 7) |
| `02-lower-complexes.tex` | L155 | "multiple of p = 1009" → see (a) |
| `03-may-and-first-obstruction.tex` | L34, L55, L67 | "(p = 1009)" → (q/(q−1))u₄ < 4p³ < p⁴ (p ≥ 5); p² − 9p − 4 > 0 and p² − 9p − 3 > 0 (p ≥ 11) |
| `03-may-and-first-obstruction.tex` | L153, L163 | "p − 6 = 1003" → p − 6 > 4 (p ≥ 11); "2p − 19 = 1999" → see (a) |
| `03-may-and-first-obstruction.tex` | L273 | "3 ≠ 0 in 𝔽₁₀₀₉" → 3 ≠ 0 in 𝔽_p (p ≥ 5) |
| `04-frobenius-elimination.tex` | L4, L337, L502 | "p = 1009" → p ≥ 11; 3p+2 < p(p−1) (p ≥ 5); 6 ≠ 0 in 𝔽_p (p ≥ 5) |
| `05-final-obstruction.tex` | L4, L141–149 | numeric W, j ≤ 1012, w < p³ < u₄ → W = 2p²+4p+1, j ≤ ⌊(W−5)/(q−1)⌋, W+1+j_max < p³ < u₄ (p ≥ 5) |
| `05-final-obstruction.tex` | L168–169, L188 | 3(p(p−1)−1) > W+1 (p ≥ 11); [−3, 8] contains no nonzero multiple of p (p ≥ 11) |
| `05-final-obstruction.tex` | L843–844 | numeric K, N, S → qp+1−S = 2p−6 > 0, N−K = p−3 > 0 (p ≥ 5) |

All 33 p-size inequalities that the manuscript writes down are listed with their line numbers in `scripts/main/ineq308.py`. The literal threshold, if every one is used as written, is 19. The adversarial review read all sections, including the appendix, and found no other place where the size of p is used.

## Reproduce

From `scripts/`:

```
/usr/bin/python3 main/may308.py 11 13 17       # all May-page claims; ~2 s each        -> expected_output/may308_p{11,13,17}.txt
/usr/bin/python3 main/may308.py 7              # negative control: FAILS ['A']          -> expected_output/may308_p7.txt
/usr/bin/python3 main/may308_dp.py 7 11 13 17  # knapsack-DP cross-check (numpy), ~20 s -> expected_output/may308_dp.txt
/usr/bin/python3 main/ineq308.py               # the 33 inequalities as predicates in p (sympy) -> expected_output/ineq308.txt
(cd main && /usr/bin/python3 spot1009.py)      # positive control at p = 1009            -> expected_output/spot1009.txt
for p in 19 23 29 31 37 41 43 47; do /usr/bin/python3 main/may308.py $p | tail -1; done   # -> expected_output/may308_p19to47.txt
# independent implementation (needs z3-solver and numpy; pip install -r independent/requirements.txt):
(cd independent && for p in 7 11 13 17; do python sweep.py $p; python d1check.py $p; done; python t1009.py)
```

## What is NOT claimed

* **We do not claim that OpenAI's p = 1009 proof is correct.** We did not re-verify the homotopy theory (§3 adjusted-filtration injectivity, §4 Frobenius left-hand side, §5 extension of powers, appendix pairing). The remark transfers whatever that proof establishes to all p ≥ 11.
* We make no claim at p = 3 or 5. At p = 7 the mechanism fails, consistent with Nave's nonexistence theorem.
* We do not claim a ring structure on V(4). (The original does not either.)
* "Verbatim" means: the same sequence of lemmas, with the numerical instantiation replaced by the symbolic one. At p = 11, 13, 17 three sub-arguments are replaced by finite machine checks. Those checks are exact enumerations, cross-validated by two independent implementations (direct recursion and DP; z3 and DP).
* Novelty: the earlier rescan session searched GitHub and arXiv (up to 2026-10-09) and found no V(4) construction at any explicit prime other than OpenAI's p = 1009. To our knowledge, no one has published one.
