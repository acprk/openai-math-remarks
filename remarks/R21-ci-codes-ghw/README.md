# R21 — Generalized Hamming weights of complete-intersection codes: planar and split cases

> **Premise and confidence (read first).**
> - This is a **direct extension of a known argument**. For r = 1 the Artinian-reduction / Cayley–Bacharach
>   argument of Eisenbud–Green–Harris (1996) is already known to give the result (Cooper–Seceleanu–Tohăneanu–Vaz Pinto–Villarreal 2020, Prop. 6.8).
>   We run the same argument with "≥ r" in place of "≥ 1", then check that the Cartesian code attains the resulting bound.
>   Every ingredient is in the literature. The new content is the assembly, which is short. We do not regard it as a new method.
> - The two unconditional cases rest on **published** EGH cases: two variables (Richert 2004; Caviglia–Maclagan 2008) and
>   the split case (Abedelfatah, PAMS 2015, Cor. 4.4). They do **not** use OpenAI's family-200 manuscripts.
>   Those manuscripts claim EGH only in characteristic 0, which does not apply to codes over F_q.
> - Under this project's gates this is a **parallel generalization** (r = 1 → all r with the same argument): it is recorded as a short remark only; no paper, no Lean.
> - Numerics (consistency checks, not part of the proof): primary 1050 + independent 1424 + adversarial 38 189 (d, r) pairs, 0 violations. The adversarial review (separate AI session) found no mathematical error: PASS-with-fixes, all fixes applied (see VERIFICATION.md).

**One-line result.** Camps Moreno–Salizzoni–San-José (arXiv:2608.19978, Aug 2026) conjecture that projective
Cartesian codes have the smallest generalized Hamming weights among evaluation codes on reduced complete
intersections of the same type (their Conjecture 3.1). This conjecture holds for **every r and every d, over every F_q**,
in two cases:
- every reduced complete intersection in **P²** (G1);
- every complete intersection whose generators f₃,…,fₙ (of the n−2 largest degrees) are **products of linear forms** (G3).

More generally, EGH for the type of X implies the conjecture for all r.

## Statement (Before / After)

Notation:
- X ⊂ Pⁿ(F_q) is a reduced complete intersection of F_q-points, with I_X = (f₁,…,fₙ), deg fᵢ = aᵢ, and 2 ≤ a₁ ≤ … ≤ aₙ.
- s = Σ(aᵢ − 1), and C_d(X) is the evaluation code of degree-d forms.
- μ_r(a, m) = min{ |𝒟| : 𝒟 a down-set of the box ∏[0, aᵢ−1] with #{β ∈ 𝒟 : |β| > m} ≥ r }.

| | Before | After (this remark) |
|---|---|---|
| Conjecture 3.1 of CMSS: d_r(C_d(X)) ≥ d_r(C_d(X_Car)) | r = 1: equivalent to CB12 (CSTPV Prop. 6.8); known in P² (Gold–Little–Schenck + CMSS), in P³ (CB12 in P³: Geramita–Kreuzer 2013 Cor. 4.4, Caviglia–De Stefani 2021 Cor. A), and for d < a₁ (CMSS Thm 2.3). All r: d = 1 (CMSS Thm 2.5; equigenerated earlier in CSTPV Cor. 6.6); d = s−1 (d_r = r+1, González-Sarabia–Camps–Sarmiento–Villarreal Prop. 9); d ≥ s (trivial); implicitly when in(I_X) = (x₁^{a₁},…,xₙ^{aₙ}) (footprint bound of GSMBVV + Beelen–Datta). | **Theorem.** If EGH holds for X, then d_r(C_d(X)) ≥ μ_r(a, s−1−d) = d_r(C_d(X_Car)) for all d ≥ 1 and all r. **Corollary.** Unconditional for (G1) every reduced CI in P² and (G3) n ≥ 3 with f₃,…,fₙ split, over every F_q (q ≥ aₙ so that X_Car exists; d_r ≥ μ_r holds for all q). |
| Role of EGH | EGH ⇒ CB12 (r = 1) | EGH ⇒ the full weight hierarchy bound (all r) |

## Proof sketch (note §2)

1. **Weight formula** (GSMBVV Thm 4.5): d_r(C_d(X)) = min{ |Γ| : Γ ⊆ X, h_X(d) − h_{X∖Γ}(d) ≥ r }.
2. **Cayley–Bacharach** (Davis–Geramita–Orecchia; EGH96 CB7): h_X(d) − h_{X∖Γ}(d) = |Γ| − h_Γ(s−1−d).
3. **Artinian reduction.** Pass to the algebraic closure and pick a linear form ℓ that vanishes at no point of X.
   - Then S/(I_Γ, ℓ) has Hilbert function Δh_Γ, and its ideal contains the regular sequence f̄₁,…,f̄ₙ.
   - EGH gives an ideal containing the pure powers yᵢ^{aᵢ} with the same Hilbert function.
   - The standard monomials of its initial ideal form a down-set 𝒟 of the box with |𝒟| = |Γ| and |Γ| − h_Γ(m) = #{β ∈ 𝒟 : |β| > m}.
   - Hence |Γ| ≥ μ_r(a, s−1−d).
4. **Cartesian code attains μ_r.** For γ in the down-set 𝒟, the products G_γ = ∏ᵢ ∏_{j>γᵢ}(xᵢ − α_{i,j}x₀) (padded with a power of x₀) are triangular on the grid.
   They have degree ≤ d exactly when |γ| > s−1−d.
   This recovers the Beelen–Datta weights of affine Cartesian codes.
5. **Unconditional cases.**
   - (G1): EGH in two variables holds over every field (Richert 2004; Caviglia–Maclagan 2008, Thm 2).
   - (G3): splitting into linear factors survives reduction mod ℓ, so Abedelfatah Cor. 4.4 (over an arbitrary field) applies.

## Reproduce

```bash
cd scripts && ./run_all.sh           # ~4-8 min, single thread; prints OK for 4 outputs
SKIP_SAGE=1 ./run_all.sh             # without Sage (independent checks only)
```
Requirements: Python 3.10 + numpy; SageMath 9.5 for the primary scripts.

| script | role |
|---|---|
| `primary/ghw_check.sage`, `primary/ghw_check2.sage` | **Primary computation** (Sage). Brute-force GHW over supports for random CIs of types (2,2), (2,3), (3,3), (2,4), (3,4) in P² and (2,2,2) in P³. 1050 (d, r) tests, 0 violations; Cartesian GHW = μ_r. |
| `independent/ghw_ci_indep.py` | **Independent re-implementation** (pure Python, no shared code). GHW from the column matroid, and μ_r by filtering all subsets of the box. CIs are certified by a Hilbert-function test plus an exact rational-point count. New types: (2,5), (3,5), (4,4) in P² and the split family G3 in P³: (2,2,2), (2,2,3), (2,3,3). **1424 (d, r) tests, 0 violations, 117 strict.** It includes a teeth control: the mis-shifted bound μ_r(a, s−d) is violated 988 times. |
| `independent/cb7_and_reduction_check.py` | Direct check of the two proof ingredients on random subsets: the CB7 identity (3280 checks) and the down-set profile of Δh_Γ (720 checks). 0 failures. |

## Attribution

- **This remark.** Theorem G is our assembly of known results; see the Premise above.
- **The question.** Posed by E. Camps Moreno, F. Salizzoni and R. San-José, arXiv:2608.19978, Conjecture 3.1.
- **The ingredients:**
  - Eisenbud–Green–Harris (CB7, CB12, and the EGH conjecture);
  - Davis–Geramita–Orecchia;
  - González-Sarabia–Martínez-Bernal–Villarreal–Vivares (generalized minimum distance functions);
  - Cooper–Seceleanu–Tohăneanu–Vaz Pinto–Villarreal (Prop. 6.8, the r = 1 equivalence);
  - Beelen–Datta (GHW of affine Cartesian codes);
  - Richert, Caviglia–Maclagan, Abedelfatah (EGH cases).
- **Related work we acknowledge.**
  - S. M. Cooper, "Subsets of complete intersections and the EGH conjecture" (2012), on Hilbert functions of subsets of CIs.
  - Geramita–Kreuzer, J. Algebra 391 (2013), on uniformity (r = 1) of complete intersections.
  - Tohăneanu's book *Commutative Algebra Methods for Coding Theory* (2024), §6.2.3.
- **Relation to OpenAI.** The topic comes from OpenAI family 200:
  - *Commuting Division-Coefficient Forms and the Artinian Eisenbud–Green–Harris Conjecture* and *The Artinian Lex-Plus-Powers Betti Theorem*.
  - Location: [github.com/openai/math](https://github.com/openai/math), `preprints/Commuting-Division-Coefficient-Forms-and-the-Artinian-Eisenbud-Green-Harris-Conjecture-September-23-2026`, commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.
  - Family 200 claims EGH in characteristic 0. Our theorem shows that a **characteristic-free** EGH would settle CMSS Conjecture 3.1 completely.
  - If family 200 is correct, the same argument gives the analogous bound for the generalized minimum distance function of reduced CIs of K-points with char K = 0. We have not verified family 200, and the corollary above does not depend on it.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction;
independent adversarial reviews were run by separate AI sessions and are not human peer review.

## What is NOT claimed

- **Not Conjecture 3.1 in general.** We prove it only in cases (G1) and (G3), and conditionally on EGH for the type of X.
  - Other known EGH cases very likely transfer the same way: Caviglia–Maclagan's fast-growing degrees, Caviglia–De Stefani (arXiv:2007.15467) Thm A (aᵢ ≥ Σ_{j<i}(a_j−1) for i ≥ 3, e.g. all CIs of types (2,2,2), (2,2,3), (2,3,3), (2,3,4) in P³), and Cooper's three-variable cases.
  - We did not check their field hypotheses and do not claim them.
- **Not new for r = 1.** That case is CB12 territory (CSTPV Prop. 6.8; P² via Gold–Little–Schenck; further r = 1 cases in Geramita–Kreuzer 2013 and Caviglia–De Stefani 2021).
- **Not a new technique.** It is the EGH96 residual argument with r in place of 1.
- **Nothing from OpenAI's characteristic-0 EGH claim** is used or verified here.
- **Limited numerics.** Our scripts use prime fields 𝔽₅–𝔽₁₇ and small types only.
  - The adversarial review added exhaustive searches over 𝔽₂, 𝔽₃, 𝔽₄ and random instances over 𝔽₄, 𝔽₅, 𝔽₈, 𝔽₉ and larger fields (its code is not shipped; see VERIFICATION.md).
  - The computations are consistency checks, not part of the proof.
- **Novelty.** To our knowledge (search of GitHub and arXiv up to 2026-10-09), no source states the r ≥ 2 bound via EGH, nor the unconditional P² and split cases for r ≥ 2.
  - arXiv: API abstract searches for generalized Hamming weight / generalized minimum distance + "complete intersection", Cayley–Bacharach + Hamming, EGH + code, and author Tohăneanu.
  - Full texts read: 2608.19978 (still v1), 1812.06529, 1707.03285, 1601.07604, 1611.06405, 1712.00387, 2006.14717, 2106.14759, 2103.14106, 1212.2653, math/0611488, math/0311129.
  - Citations followed on Semantic Scholar: of 2608.19978 and of Geramita–Kreuzer 2013.
  - Books: Google Books search inside Tohăneanu (2024).
  - Not accessed: the full texts of Geramita–Kreuzer 2013 and Cooper 2012 (publisher pages blocked). Our reading of them comes from their abstracts and from citing papers.
