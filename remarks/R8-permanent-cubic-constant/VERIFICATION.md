# R8 verification record (2026-10-09)

## 1. Primary computation (lemma A/B attack session; `scripts/primary/`)
| Item | Command | Result |
|---|---|---|
| BOX unit test | `box_check.sing` (5 configurations) | len Z ≥ (q+1)², containment and Ann checks all 1. Note: `EX=1` with k=7 errors because block sizes exceed k. That is a script bug and was not used. |
| Constants | `constants_box.py 1408 5000` | BOX 125 317e / 120 334e; stated/BOX 44.1× / 46.0× |
| Ceiling | `constants_corrected.py 1408 5000` | A_ceil 30 872e / 31 857e (179× / 174×) |
| Lemma B | `bipartite_order.sing` K=16 | min order 14, local length 105 |
| Family dims | `family_exact.py`, `lemmaB_violation.py` | k=22: 253 < 256; k=31: 473 < 484 |

## 2. Independent re-implementation (adversarial reviewer; `scripts/independent/`, no primary code reused)
| Item | Script | Result |
|---|---|---|
| BOX steps 1–3 (k=7,8,9; includes q > c: (8,3), (9,3)) over ℚ, F_1000003, F_{2^10} | `box_indep.sage` | 18/18: X ∈ Γ, u regular, Hilbert polynomial = len Z ≥ (q+1)², I+(u) ⊆ box, low-degree GB ⊆ Ann |
| Gradient formula | `grad_check.sage` | Cayley–Hamilton adjugate formula agrees with brute-force derivative (k=5,6, ℚ and GF(2^20)) |
| Counterexamples | `grad_check.sage` | ∇G_D = 0 exactly for k=9,q=3 (codim 61 < 64) and k=22,q=7 (253 < 256), over ℚ and GF(2^20). Controls: random X has 81/484 nonzero entries; Prop F violated (ranks (6,0)) has 96 nonzero entries |
| Γ equations are combinations of almost-principal (r−1)-minors | `minor_span.sage` | True for (k,q) = (5,1), (6,1), (7,2), (7,1), over ℚ and GF(2^8) |
| Lemma B counterexample | `lemmaB_check.sage` | adj·M = det·I checked; min order 14; m^14/m^15 stabilisation gives local length = 105 exactly over ℚ |
| Constants | `constants_indep.py` | m=1408: 125 316.9e; m=5000: 120 334.0e; paper-tuned 262 835 / 245 793 |
| Uniform constant | `uniform_indep.py` | max over 1408 ≤ m ≤ 200 000 is 125 422.2e (m=1419). Tail: g′ numerator has largest real root 64.3; g(2·10⁵) = 118 893e → 118 800e |
| Ceiling | `ceiling_indep.py` | asymptotic 33 065.4e (β = 0.3919, 167.2×); m=1408: 30 871.8e; m=5000: 31 857.5e |

## 3. Adversarial review of the proofs (by hand, same session as §2)
- BOX: checked line by line against `adjugate.tex`. L139–248 (generation) and L261–270 (gradient ⇔ [s^q t^q u^c] adj = 0) carry no hypothesis on q; L291 is the only use of q ≤ c, and BOX bypasses it. The Ann step is coefficient extraction, so it is valid in char 2. For degrees d ≥ k the containment holds trivially, because B_d = S_d for d ≥ 2q+1 and 2q+1 ≤ k. u is regular because a diagonal adjugate entry is nonzero at each point P_i (L81–87). **Verdict: CONFIRMED.**
- Chain constraints: r ≥ 2, 2 ≤ d = Δ < k², d−1 < Δ (`transfer.tex` L144–145), δ = Δ ≤ v, q ≤ ⌊k/2⌋ (`selector.tex` L56), m ≥ 11k−2q. All are rechecked in `constants_indep.py` and `uniform_indep.py`.
- Prop F: proved by generalized Laplace expansion along the row blocks (note §4). Spot checks are in §2.
- **Corrections to the earlier report (REPORT.md 10/09):**
  - "174×" is the per-m ceiling at m=5000. The ceiling that applies to a *uniform* constant is the asymptotic one, 33 065e (167×). The per-m ceiling increases towards it (179× at m=1408).
  - The ceiling does not depend on char 2. The families lie in Γ_D over every field, so it binds every linear-section proof, including one that bounds the complex codimension directly.
  - The uniform statement for m ≥ 1408 is 1/(125 500 e), i.e. 44×. The figure 46× is the value at m=5000 (46.5× asymptotically).
  - REPORT's "½−γ²" uses γ = r/k; this remark writes β.
  - Conjecture A\* was not written down in REPORT. The formulation in the note is ours.

## 4. Paper Prop. 7.3 vs Lean slice size; follow-up literature (added 2026-10-09)
Sources: openai/math local checkout at `08fcd618`. The files `lean/OAI/Algebra/Determinantal/*`, `lean/ComparatorChallenges/PermanentCubic.lean` and the family-108 preprint are unchanged since `adc7f124` (`git log adc7f124..HEAD` on these paths is empty).

| # | Question | Finding |
|---|---|---|
| L1 | What does the paper use? | Prop. 6.1 (`dimension.tex`, `prop:gradient-codimension`): dim V(∇G_D) ≤ k² + k − T, T = C(q+2,2). Prop. 7.3 (`transfer.tex`, `prop:smooth-projection`): q = ⌊(k−1)/3⌋, d = ⌊k²/100⌋, with r−1 ≥ k/4 and d−1 ≥ k²/200. Thm 1.1 (`consequences.tex` L8–38): k = ⌊m/11⌋ ≥ m/12, so bdc ≥ k³/(3200e) ≥ m³/(5 529 600 e), where 5 529 600 = 12³·3200. |
| L2 | What does Lean prove? | `Main.lean` L214–245 `critical_dimension_effective`: krullDim ≤ (k+1)² − ⌊C(q+2,2)/18⌋ for size k+1, assuming 3q ≤ k. The /18 comes from `ConductorDimension.lean` L430–444 `half_curve_parameter_bound` (18n + 5ℓ + p ≤ 9k(k+1)) combined with `normalization_length_le_half` (defect ≤ 5ℓ + p). `Main.lean` L259 `sectionSize k := C(selectorIndex k + 2, 2)/18`. `sectionSize_bounds`: D−1 ≥ k²/330. `degree_effective_bound`: r−1 ≥ 0.33k. `explicit_constant` L313: m ≤ 11.1k, and (r−1)(D−1) ≥ k³/1000, so 4·1000·(11.1)³ = 5 470 524 ≤ 5 529 600. `permanent_cubic_lower_bounds` states 5 529 600, which is the comparator statement. |
| L3 | Is there an inconsistency? | **Yes, confirmed (Guo's remark is correct).** The Lean slice D = ⌊C(q+2,2)/18⌋ ≈ k²/324 is not Prop. 7.3's d = ⌊k²/100⌋. The Lean dimension bound (codim ≈ T/18) does not support d = ⌊k²/100⌋, which needs codim ≥ Δ ≈ k²/18. The Lean theorem is still a correct proof of the *stated constant* by another route. `ConditionalBound.lean` contains the paper's parameters (`sectionDimension = k²/100`, `codimension = T − 2k`) only together with the conditional `selector_lower_of_critical_dimension`, which takes the dimension bound as a hypothesis. |
| L4 | Which one does our accounting use? | **The paper.** The baseline 5 529 600 is taken from `introduction.tex` L48 / `consequences.tex` L36. The row "manuscript's own lemma, parameters re-chosen" uses Δ = C(q+2,2) − k from Prop. 6.1, and the BOX rows use (q+1)² − k. The number 5 529 600 is identical in paper and Lean, so **the ratios 44.1× (uniform m ≥ 1408), 46.0× (m=5000) and 46.5× (asymptotic) are unchanged.** |
| L5 | Does it affect the claim? | It affects the scope, not the number. The new constants 125 500e / 118 800e assume the paper's Prop. 6.1 (conductor-length + Severi count, codim T − k). Lean does not formalize this at that strength. It was already covered by the "Premise" paragraph and is now stated explicitly in README and note. The pipeline from a given dimension bound to bdc is Lean-formalized (conditional theorem above). |
| L6 | What would Lean-only inputs give? | `scripts/independent/lean_dimension_model.py` (expected output `scripts/expected_output/lean_dimension_model.txt`). With Lean's certified D = ⌊C(q+2,2)/18⌋, 3q ≤ k−1, and parameters re-chosen: c⁻¹ = 4 205 362e at m=1408 (1.315×), 4 268 907e at m=5000, asymptotically 4 289 904e (1.29×). Lean's own fixed choice gives 5 091 307e at m=1408, which is consistent with its 5 529 600. Hypothetical box lemma inside Lean's T/18 route (not formalized): ≈ 2.14·10⁶ e (2.59×). The remaining factor of about 17× between this and 125 500e is the /18 loss in the formalized curve-parameter count. |

### Follow-up literature (dedup)
- Guo, arXiv:2609.29568, checked on arxiv.org on 2026-10-09: v1 2026-08-26, v2 2026-10-08. The abstract states L_div(per_n) ≥ (n²/86400) log₂(n/48) for n ≥ 1408 *assuming Prop. 7.3*, and (n²/119790) log₂(n/44) for n ≥ 1936 from the Lean declarations plus a written argument. Also liminf L_div/(n² log₂log₂ n) ≥ 1/12. Its target is circuits with division, not the cubic bdc constant.
- Xie, Zenodo 23108131, checked on zenodo.org on 2026-10-09: 2026-10-02, bdc(per_n) ≥ n²(log₂n − 5)/(64e) for n ≥ 64, i.e. Ω(n² log n). This is weaker than cubic in n and does not touch 5 529 600.
- Neither changes 5 529 600, and neither overlaps with the box lemma or the ceiling.

## Not covered (GAPs)
- No human review. No Lean formalization (by project policy, small lemma + parameter choice).
- The premise that the manuscript is correct is not re-verified here.
- The optimality of the block families (Conjecture A\*) is open. Families with several ranked blocks were not explored exhaustively.
- `minor_span.sage` covers only k ≤ 7. The general statement rests on the expansion argument.
- Dedup was limited to the two follow-ups above (abstracts read; full texts not read). The paper's Prop. 6.1 at strength T − k is not Lean-verified, and our constants inherit this premise.
