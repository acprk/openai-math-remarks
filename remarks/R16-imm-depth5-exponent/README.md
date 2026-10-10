# R16 — Homogeneous depth-five IMM: exponent √n/400 → (1/(4√2) − ε)√n

**One-line result.** With the same method as OpenAI's *Homogeneous depth-five lower bounds for iterated
matrix multiplication*, every syntactically homogeneous ΣΠΣΠΣ circuit computing IMM_{n,n} over a field of
characteristic zero has size at least **n^{(1/(4√2) − ε)√n}** for n ≥ n₀(ε). The paper states
n^{√n/400}. The exponent constant improves by 100/√2 ≈ **70.7×**, and the thresholds are explicit.

> **The argument is OpenAI's; this remark re-tunes its parameters.** There is no new idea. The only step not
> already in the paper is expanding all high-degree factors of an upper product at once (a standard
> multilinear expansion). The IMM-side proof (Sections 4–5) is used unchanged.

## Statement (Before / After)

| | Before (OpenAI) | After (this remark) |
|---|---|---|
| Theorem 1.1 / Corollary 6.1 | S ≥ n^{√n/400}, n ≥ n₀ (n₀ not specified) | S ≥ n^{(1/(4√2)−ε)√n}, n ≥ n₀(ε) |
| Leading constant | 1/400 = 0.0025 | 1/(4√2) = 0.176777 (×70.7) |
| Explicit thresholds (natural log) | — | c = 1/400 for log n ≥ 17.6; c = 1/8 for log n ≥ 57.941; c = 1/6 for log n ≥ 296.735 |
| Upper bound (unchanged) | n^{√n+4} (Prop. 6.2) | n^{√n+4} |

Further thresholds: c = 3/320 for log n ≥ 18.24; c = 3/20 for log n ≥ 112.04 (`scripts/expected_output/explicit_bound.txt`).
At log n = 296.7 the bound is 0.166665 < 1/6, so the threshold for 1/6 is 296.735, not 296.7.

## Where the factor 70.7 comes from

| Step | Change | Factor |
|---|---|---|
| 1/400 → 3/320 | slack in the paper's own chain (its inequalities already give any c < 3/320) | ×3.75 |
| 3/320 → 1/8 | re-tuning: q ≤ n^{−2/5} → q ≈ n^{−1/2}; λ ≥ 3/8 → λ ≈ 1/2; expansion threshold and the scale of s | ×13.3 |
| 1/8 → 1/(4√2) | expand all |X| high-degree factors simultaneously, at cost S^{|X|} (standard multilinear expansion; the only new step) | ×√2 |

With one expansion and s = σ√n the bound is min(σ/4, 1/(16σ)), maximal (1/8) at σ = 1/2. Simultaneous
expansion gives min(σ/4, 1/(8σ)), maximal at σ = 1/√2.

## Attribution

The construction, the theorem and its proof are OpenAI's work: family 135 of
[github.com/openai/math](https://github.com/openai/math), commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/Homogeneous-depth-five-lower-bounds-for-iterated-matrix-multiplication-September-25-2026/`.
All line numbers below refer to `build/sections/*.tex` at that commit.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent checks were run by separate AI sessions and are not human peer review.

## What changes

| Id | File : lines | Original | Replacement |
|---|---|---|---|
| M1 | `02-rank.tex` L57–58 | s = largest integer ≤ √n with the parity of n | s = largest integer ≤ √(n/2) with the parity of n. Lemma 2.1 holds as written (√(n/2) − 2 ≥ √n/2 for n ≥ 95). |
| M2 | `03-circuits.tex` L123–205 (Lemma 3.3) | t₀ = n/(4s); H_e ≤ q^{d_e}(1+4q^{1−2d_e}); even corrections ≤ C₁√n | t₁ = n/(2s), x_e = se/(2n) < 1/4 for e < t₁. Even e: H_e ≤ q^x(1 + (q^{1−2x}+q)/(1−q)); odd e: H_e ≤ q^x(q^{1/2−2x}+q^{1/2})/(1−q); H₁ ≤ 2q^λ. Total corrections C_e = √n(1+o(1)), C_o = 4s(√q + q/(1−q)) = O(n^{1/4}), both explicit. |
| M3 | `03-circuits.tex` L250–299 (proof of Prop. 3.1) | either all factors low, or expand **one** factor of degree ≥ t₀; R(g) ≤ 2S²D e^{C√n} n^{−√n/100} | expand **every** factor occurrence of degree ≥ t₁ (set X): R(G) ≤ 2D S^{|X|} exp(−L_q s/2 − E_X c′ + C_e + C_o), with L_q = log(1/q), E_X = Σ_{j∈X} e_j, c′ = L_q(n−2s)/(2n) − log 2 |
| M4 | `06-completion.tex` L9–29 (proof of Thm 1.1) | D e^{−C₁√n} ≤ R(f) ≤ 2S²D e^{C₂√n} n^{−√n/100} ⇒ S ≥ n^{√n/400} | IMM side R(f) ≥ D e^{−C_I} with C_I = 2√n + (5/2 + 1/(2√2)) log n + 1.292 + o(1) (the paper's Sections 4–5 chain with all terms explicit); then (1+|X|) log S ≥ A + |X| t₁c′, so log S ≥ min(A, t₁c′), A = L_q s/2 − C_I − C_e − C_o − log 2. Both terms are (1/(4√2))√n log n (1 − o(1)). |

Consequential edits: the constant in Theorem 1.1 (`01-introduction.tex` L54–65) and Corollary 6.1
(`06-completion.tex` L35–42). The characteristic-zero transfer is unchanged.

## Limit of the framework

Within "Lemma 3.2 + arbitrary expansion, parameter family λ = 1/2 − σ/(2√n), q ≍ n^{−1/2}, the leading
constant cannot exceed 1/(4√2)." This is a first-order continuum argument, not a theorem: the best
provable constant in this family is τ*(σ) = min(σ/4, (√(σ²+4) − σ)/8), maximal at σ = 1/√2 with value
1/(4√2) (`scripts/main/framework_limit.py`; single- and two-degree adversaries in
`scripts/independent/v_framework.py`).

The real obstacle is on the IMM side, which forces q ≳ n^{−1/2}. With α = n^{−γ}, γ > 1/2, an exact
transfer-matrix evaluation of the word sum of Section 5 gives a loss of order n^γ/(√2 γ log n). This loss
eventually exceeds the circuit-side gain ((2γ−1)/(4√2))√n log n. At n = 10⁶ and γ = 0.75 it is already
larger: 0.096 vs 0.088, in units √n log n (`scripts/main/imm_wordsum.py`).

## Reproduce

```bash
cd scripts && ./run_all.sh     # python3, mpmath 1.3, numpy, scipy; single thread; about 4 minutes
```
`run_all.sh` prints `OK` for each of the eight scripts if its output matches `expected_output/`.

- `main/explicit_bound.py` — positive control (reproduces the paper's chain 1/20, 3/160 → 1/100 → 1/400
  and shows that the chain gives any c < 3/320). It then evaluates an explicit lower bound for
  log S/(√n log n), monotone in n, and bisects the thresholds.
- `main/hbound_check.py` — exact H_e against the M2 bounds with the real parameters at n = 10⁶, 4·10⁶, e^14.
- `main/framework_limit.py` — τ*(σ), exact finite-n framework value, constrained vs. unconstrained sum in
  Lemma 3.2, and the side conditions of Lemma 2.1 at n = e^14.
- `main/imm_wordsum.py` — exact word sum of Section 5 for α = n^{−γ}.
- `independent/` — separate re-implementation from the TeX source, with the exact integer parameters of (eq:parameters) (`02-rank.tex` L59–72):
  `v_params_H.py` (M2 inequalities and inequality (*) of M3 over many adversarial degree lists),
  `v_imm.py` (exact word sum vs. the paper's bound; explicit C_I), `v_threshold.py` (thresholds from
  exact parameters, also with s replaced by s − 2), `v_framework.py` (framework limit).

## Scope and limitations (what is NOT claimed)

- **Re-tuning, no new idea.** Every step is the paper's argument with different parameters, plus a
  standard simultaneous multilinear expansion.
- **IMM side unchanged.** Sections 4–5 are used as written; this remark relies on their correctness.
  C_I is the paper's own chain made explicit. It is loose (the exact word sum is ≈ 0.11√n at n = 6·10⁴),
  and tightening it would lower the finite thresholds but not the leading constant. This is not used.
- **Framework limit is not a theorem.** It is a first-order continuum computation within the stated
  parameter family.
- **No Lean.** The official Lean formalisation (`lean/OAI/Algebra/DepthFive/`) covers the original 1/400
  statement. Nothing here is formalised.
- **No peer review.** The checks were done with AI assistance and are not human peer review.
