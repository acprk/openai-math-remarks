# R8 — family 108: the cubic constant for the border determinantal complexity of the permanent

**One line.** A short "box lemma" (len Z ≥ (q+1)², with no condition q ≤ r−1) and re-chosen parameters improve OpenAI's bound bdc(per_m) ≥ m³/(5 529 600 e) to **m³/(125 500 e) for all m ≥ 1408** (44×), and to m³/((118 800+o(1)) e) asymptotically (46.5×). Explicit block-triangular matrices show that the gradient codimension is not (2q+2)². The manuscript's pipeline cannot do better than about 1/(33 065 e), which is 167× the stated constant.

## Attribution
The original construction and proof are OpenAI's work: family 108, *A cubic lower bound for border determinantal complexity of the permanent* (September 24, 2026), [openai/math `preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026`](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026), commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. This remark does three things:
- proves one short lemma that replaces Prop. `adjugate-length`;
- re-chooses parameters in OpenAI's chain;
- gives explicit families inside the gradient locus, which bound what the method can achieve.

## Before / After
| | Before (manuscript) | After (this remark) |
|---|---|---|
| Length bound (`adjugate.tex` L250–309) | len Z ≥ C(q+2,2), needs q ≤ c = r−1 (L291) | len Z ≥ (q+1)², any 0 ≤ q ≤ (k−1)/2 |
| `dimension.tex` L24–25 | dim Γ_D ≤ k²+k−C(q+2,2) | dim Γ_D ≤ k²+k−(q+1)² |
| Parameters (`transfer.tex` L219–222, `consequences.tex` L9) | q=⌊(k−1)/3⌋, d=⌊k²/100⌋, k=⌊m/11⌋ | d = Δ = (q+1)²−k, q ≈ 11k/32, k maximal with 11k−2q ≤ m |
| Theorem (`introduction.tex` L48, `consequences.tex` L36) | c = 1/(5 529 600 e), m₀ = 1408 | c = 1/(125 500 e), m₀ = 1408; asymptotically 1/((118 800+o(1)) e) |

### Breakdown of the gain (values of c⁻¹/e)
| | uniform m ≥ 1408 | m = 5000 | m → ∞ |
|---|---|---|---|
| stated | 5 529 600 | 5 529 600 | 5 529 600 |
| manuscript's own lemma, parameters re-chosen | 267 955 | 245 793 | 238 328 (= 8·31³) |
| box lemma + parameters | 125 422 | 120 334 | 118 800 (= 4·29 700) |
| **parameter choice only** | **20.6×** | 22.5× | **23.2×** |
| **box lemma (on top)** | **2.14×** | 2.04× | **2.006×** |
| total | 44.1× | 46.0× | 46.5× |

Most of the loss comes from crude estimates in the manuscript: d ≈ k²/100 when Δ ≈ k²/18 is available, k ≥ m/12, and r−1 ≥ k/4. Fixing these accounts for about 23×. Of the box lemma's factor 2.006, the change C(q+2,2) → (q+1)² gives 2.000, and dropping q ≤ c gives 1.003.

## Related work
Two follow-ups to family 108 bound other quantities and leave the cubic constant 5 529 600 unchanged:
- **Jizhou Guo**, *An n² log log n lower bound for permanent circuits with valid division*, [arXiv:2609.29568](https://arxiv.org/abs/2609.29568) (v1 2026-08-26, v2 2026-10-08). v2 moves the smooth slice of family 108 to circuits with division: assuming the manuscript's Prop. 7.3 it gets L_div(per_n) ≥ (n²/86400)·log₂(n/48) for n ≥ 1408; from the Lean declarations alone, plus a written argument, it gets (n²/119790)·log₂(n/44) for n ≥ 1936. Guo also notes that the Lean slice size D = ⌊C(q+2,2)/18⌋ differs from Prop. 7.3. We discuss this point below.
- **Ying Xie**, *Smooth slices and superquadratic border determinantal complexity of the permanent*, [Zenodo 23108131](https://zenodo.org/records/23108131) (2026-10-02): bdc(per_n) = Ω(n² log n), explicitly ≥ n²(log₂n − 5)/(64e) for n ≥ 64.

Neither work changes the cubic constant. This remark improves that constant only.

## Which source the baseline comes from: paper vs. Lean
- The baseline **5 529 600** is the constant of the paper's Theorem 1.1. It is derived in `consequences.tex` L8–38 from Prop. 7.3 (`transfer.tex`: d = ⌊k²/100⌋, r−1 ≥ k/4, d−1 ≥ k²/200, then k ≥ m/12, so 12³·3200). The Lean theorem `PermanentBorder.permanent_cubic_lower_bounds` states the same number.
- The Lean proof reaches that number by a **different and weaker** chain (`lean/OAI/Algebra/Determinantal/Main.lean` L214–388, unchanged since `adc7f124`). It proves the dimension bound only in the form dim ≤ k² − ⌊C(q+2,2)/18⌋, where the paper's Prop. 6.1 has k² + k − C(q+2,2). The factor 18 comes from `ConductorDimension.half_curve_parameter_bound`, 18n + 5ℓ + p ≤ 9k(k+1). Lean then takes the slice D = ⌊C(q+2,2)/18⌋ ≈ k²/324 and compensates with tighter estimates: r−1 ≥ 0.33k, D−1 ≥ k²/330, m ≤ 11.1k. Since 4·1000·11.1³ ≈ 5 470 524 ≤ 5 529 600, it lands on the same constant. Guo's observation is correct: the Lean slice is not the slice of Prop. 7.3, and Prop. 7.3's d = ⌊k²/100⌋ is not covered by the Lean dimension bound.
- **Our accounting follows the paper.** Every row of the table above, and our new constants, use the paper's Prop. 6.1 codimension (C(q+2,2) − k, or (q+1)² − k with the box lemma). The baseline number is the same in both sources, so **the reported ratios do not change**: 44.1× uniform, 46.5× asymptotic. What does change is the scope of the claim. The new constants depend on the paper's (unformalized at that strength) Prop. 6.1/conductor-length/Severi count. This was already part of our Premise, and it is **not** Lean-certified. The selector + smooth-section + smooth-border pipeline is formalized in Lean conditionally (`ConditionalBound.lean`, `selector_lower_of_critical_dimension`).
- For comparison, `scripts/independent/lean_dimension_model.py` computes what the Lean-certified dimension input would allow. With re-chosen parameters it gives c⁻¹ ≈ 4.21·10⁶ e at m = 1408 (1.31×) and 4.29·10⁶ e asymptotically (1.29×). Inserting the box lemma into Lean's T/18 route, which is hypothetical and not formalized, gives ≈ 2.14·10⁶ e (2.59×). The remaining gap to our 125 500 e (≈ 17×) is essentially the factor-18 loss in the formalized curve-parameter count.

## Proof of the box lemma
Let Λ = [s^q t^q u^c] and let I_Z be the saturated ideal of Z.
1. By Lemma `adjugate-generation` (`adjugate.tex` L139–248), every degree-(k−1) form in I_Z is a constant combination of adjugate entries, so Λ kills it. Hence for d ≤ k−1, (I_Z)_d ⊆ Ann Λ = (s^{q+1}, t^{q+1}, u^{c+1}). This is coefficient extraction, so it holds in every characteristic.
2. For d ≤ 2q ≤ k−1, step 1 gives (I_Z+(u))_d ⊆ B_d, where B = (s^{q+1}, t^{q+1}, u). For d ≥ 2q+1, B_d = S_d. So I_Z+(u) ⊆ B in every degree.
3. Z misses the line u = 0, because a diagonal adjugate entry is nonzero at each P_i (`adjugate.tex` L81–87). So u is a nonzerodivisor modulo I_Z, and len Z = dim S/(I_Z+u) ≥ dim S/B = (q+1)².

Everything after this step in the manuscript (the conductor bound, the Severi bound, the fibre bound and the transfer) does not depend on q.

## The gradient codimension is not (2q+2)²; limit of the method
- **Prop. F.** Each ∂G_D/∂X_ij is a combination of almost-principal (r−1)-minors. If X is coordinate block upper triangular with m′ blocks and diagonal-block ranks ρ_i satisfying Σρ_i + m′ ≤ r−1, then all these minors vanish, so X ∈ Γ_D. This holds for every field, characteristic and θ.
- **Counterexamples.** For k=9, q=3, the matrices [[0,Y],[0,0]] form a family of codimension 61 < 64 = (2q+2)². For k=22, q=7, blocks (16,6) with ranks (5,0) give codimension 253 < 256. Both were checked by exact computation over ℚ and GF(2^20). Asymptotically the families have codimension (½−β²)k², where β = r/k, compared with (2q+2)² ≈ (1−β)²k². The difference is 2(β−½)² k².
- **"Lemma B" is false.** For k=16 and X = [[0,Y],[0,0]] with 8×8 blocks, len Z ≥ 105 for every Y (the exact local length is 105). The family has dimension 64, but Lemma B would allow at most 62.
- **Ceiling.** Suppose a proof applies Theorem `smooth-border` to a linear section of G_D and uses the selector with m ≥ 11k−2q. A smooth section needs d ≤ codim Γ_D. So any uniform constant from this pipeline has c⁻¹ ≥ (33 065+o(1)) e, attained at β = (√601−1)/60 ≈ 0.392. That is at most **167×** better than the stated constant; the box lemma is 3.59× from this limit. Per-m ceilings are 30 872e at m=1408 (179×) and 31 858e at m=5000 (174×). These per-m values do not limit a uniform constant.
- **Open conjecture A\*.** For β ∈ (¼,½), codim Γ_D ≥ (½−β²−o(1))k², i.e. the block families are asymptotically extremal. A proof would give about 33 065e.

## Reproduce
Requirements: Python 3, SageMath 9.5 (with its bundled sympy), Singular 4.2.1. Each command is single-threaded.
```
cd scripts/independent
python3 constants_indep.py 1408 5000 20000     # per-m constants, paper-tuned vs BOX        (~3 s)
sage -python uniform_indep.py 1408 200000      # uniform constant + exact tail monotonicity (~2 min)
python3 ceiling_indep.py 1408 5000             # method ceiling                              (~3 min)
python3 lean_dimension_model.py                 # constants from the Lean-certified dimension bound (~2 s)
sage box_indep.sage both                       # BOX steps 1-3, QQ / F_1000003 / F_2^10       (~10 s)
sage grad_check.sage both                      # counterexamples k=9,q=3 and k=22,q=7         (~5 min)
sage lemmaB_check.sage 16                      # Lemma B counterexample, local length 105     (~15 s)
sage minor_span.sage                           # partials lie in span of almost-principal minors, k<=7 (~30 min)
cd ../primary
python3 constants_box.py 1408 5000; python3 constants_corrected.py 1408 5000
Singular -q -c 'int K=16; int CH=32003;' bipartite_order.sing
Singular -q -c 'int K=8; int Q=2; int CH=32003; int EX=4;' box_check.sing
```
Expected outputs are in `scripts/expected_output/`.

## What is NOT claimed
- The box lemma is a short observation in OpenAI's framework. We do not claim it as a new technique, and we did not formalize it in Lean.
- The ceiling applies only under the stated premises: linear sections of this G_D, Theorem `smooth-border`, and the 11k−2q selector. It does not cover initial forms at nonzero points, other coefficients, other selectors, or other lower-bound theorems.
- The true codimension of Γ_D is not determined. Conjecture A\* is open, and we have not checked whether families with several ranked blocks do better.
- The uniform constant for 1408 ≤ m ≤ 200 000 comes from an exhaustive machine search over a local window of (k,q). Any admissible choice gives a valid bound, so the search can only be pessimistic.

## Premise
This remark assumes the manuscript is correct, in particular Lemma `adjugate-generation`, Prop. `conductor-length`, Lemma `severi-dimension`, Prop. `matrix-fiber`, the transfer and smooth-slice propositions, Theorem `smooth-border`, and Prop. `selector`. These are the paper's statements. The Lean development proves the resulting dimension bound only in the weaker form k² − ⌊C(q+2,2)/18⌋ (see above), so our constants are not Lean-certified.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review. To our knowledge (dedup on GitHub/arXiv/Zenodo up to 2026-10-09; see Related work), the cubic constant has not been improved elsewhere. Guo (arXiv:2609.29568) and Xie (Zenodo 23108131) bound different quantities.
