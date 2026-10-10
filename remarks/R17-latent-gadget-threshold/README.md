# R17 — Unique Games latent gadget: rank threshold ≈ p^-16 → 1.5·p^-3, and a 1/(4p) lower bound

## Premise and confidence (read first)

- **Before → After.** Lemma 3.1 of *The Unique Games Theorem* needs a rank threshold r_*(p). The preprint only says it exists. Made explicit, its proof certifies **r_* ≈ p^-16** (≈ 2^450 at p = 10^-8), or **p^-23.1** (≈ 2^648) with its effective recipe. The same construction satisfies the lemma with **r_* ≤ ⌈1.5·p^-3⌉ + 1** (≈ 2^80.3 at p = 10^-8). We also prove **r_* ≥ 1/(4p)** for every gadget satisfying the lemma (≈ 2^24.6).
- **This is bookkeeping, not a new method.** Two steps change, and both re-analyse OpenAI's own construction:
  1. the final enlargement step, which squares the threshold, is dropped. The community presentation lgboim/unique-games-short-proof already did this for ℓ divisible by 3. We add a padding step for the other ℓ.
  2. the harmonic potential with Markov's inequality is replaced by the supermartingale e^{-3θj}/(r_j+1).

  The lower bound is a few lines of Fourier analysis. An expert could write all of it in an afternoon.
- **Confidence.** The proofs are short and were checked step by step against the preprint's source by an adversarial AI session. Verdict: confirmed with corrections, all applied; see `VERIFICATION.md`. We take the preprint's lemmas as given.

## Statement (Before / After)

| | Before (certified by the preprint's own proof) | After |
|---|---|---|
| r_*(p), upper | 3d+2 with d ≳ 3r_0², r_0 ≈ e^{8−γ}p^-8 ⇒ ~p^-16 (main text); ~p^-16/ln2 = p^-23.1 (effective recipe at the end of the proof of Lemma 3.1) | **⌈1.5 p^-3⌉ + 1** for every large ℓ; ⌈1.5 p^-3⌉ − 1 for ℓ ∈ 3ℤ |
| r_*(p), lower | none stated | **≥ 1/(4p)** (any gadget satisfying Lemma 3.1) |
| p = 10^-8 | 2^449.8 (main text), 2^647.8 (effective recipe) | 2^80.3; lower bound 2^24.6 |

The exponent of 1/p goes from 16 (or 23.1) to 3. The unknown truth for this gadget is between p^-1 and p^-3.

Intermediate point (for attribution): the paper's harmonic analysis without the enlargement step gives r_* = r_0 ≈ e^{8−γ}p^-8 (2^223.3 at p = 10^-8). This is the threshold of lgboim's presentation, which handles ℓ ∈ 3ℤ only.

## Attribution

The gadget, its construction and its analysis are OpenAI's work: family 102 of
[github.com/openai/math](https://github.com/openai/math), commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/The-Unique-Games-Theorem-September-23-2026/` (Section 3 of the PDF, source `build/sections/02-latent.tex`).
This remark only re-packages the last step, re-chooses the parameters of the rank-survival argument,
and adds an elementary lower bound.

Adjacent work, acknowledged: [lgboim/unique-games-short-proof](https://github.com/lgboim/unique-games-short-proof)
(2026-10-09), `papers/ugc_short_proof.tex`, lines 185–200 and 449–485. It takes K = B_blk with
ℓ = 3d, sets r_* = r_0, and starts the descent from an arbitrary S. It keeps the harmonic analysis
(r_0 ≈ e^{8−γ}p^-8). It does not contain the supermartingale, the p^-3 bound, the padding to
every ℓ, or a lower bound.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are
not human peer review.

## What changes (line numbers of `build/sections/02-latent.tex`)

| Where | Paper | Here | Effect on the exponent |
|---|---|---|---|
| l. 639–654, 674–745 (base detection, enlargement) | detects only families restricting onto all of B^*, then enlarges B to K via all injections, r_* = b+2 = 3d+2 ≈ 9r_0² | the descent works from any r_0-dimensional domain S (a lift may have any domain, l. 477–481); K = B for ℓ = 3d; pad with F_2^e, e ∈ {1,2}, otherwise; r_* = r_0 + 2 | ÷2 (16 → 8) |
| l. 534–612 (Prop. 3.9) | harmonic potential, E[drop] ≤ 4θ per level, Markov at ½ ⇒ sθ ≤ H_{r_0}/8 ⇒ r_0 ≈ e^{8−γ}p^-8 | supermartingale e^{-3θj}/(r̃_j+1) with the per-level loss rate 3rθ; losses at step 0 or at nongeneric steps are counted as failures (≤ 1/60 in total) | 8 → 3 |
| l. 537, 657–660 (s, r_0) | s = ⌊H_{r_0}/(8θ)⌋, exp(1 − H_{r_0}/8) ≤ p | s = ⌈ln(1/p)/θ⌉, r_0 = ⌈1.5p^-3⌉ − 1, ε_0 = 1/(120(ln(1/p)+1)2^{r_0}), θ ≤ 0.01, 2^{r_0}θ ≤ 1/120 | — |

Where the loss was: the **enlargement step squares the threshold** (×2 in the exponent). The
**Markov/harmonic step costs a factor 8/3** in the exponent of r_0. The factor 8 is 4 (drop rate,
3θ plus a θ allowance for nongeneric spaces) times 2 (Markov at ½).

Limit of this route (method only, not a lower bound on r_*): if the only input is the per-level
loss estimate 3rθ, the detection probability is about 1 − e^{−r_0 p³/2}. So r_0 must be at least
(2 ln(8/7) + o(1))·p^-3 ≈ 0.267·p^-3, and the constant 1.5 is within a factor 5.6 of this.

## Proof sketch

Full proofs are in `note.pdf`.
1. **Lower bound.** Let F_σ = (−1)^{σ∘C}. Equivariance forces the Fourier support of F_σ onto lifts of σ. The Fourier weights ν_σ are probability laws, and E_{σ≠0}(1 − s_σ) ≤ p' with s_σ = Σ_h ν_σ(h)μ̂(h). Sample σ_1..σ_r and h_i ~ ν_{σ_i}. By the union bound, E[detect] ≤ r·p'/2. With probability → 1 the σ_i are independent, and then detection is ≥ 1/8. Hence r_* ≥ 1/(4p).
2. **Upper bound.** Run the paper's descent from an arbitrary r_0-dimensional S. The bad events are a loss at step 0 (≤ (2^{r_0}−1)θ) or a loss at a nongeneric step (≤ sθ2^{r_0}ε_0); together they have probability ≤ 1/60. Before the first bad event, W_j = e^{-3θj}/(r_j+1) is a supermartingale. So Pr[extinct] ≤ e^{3θs}/(r_0+1) + 1/60 ≤ (2/3)e^{0.03} + 1/60 < 0.704 < 3/4. The leaf detects with probability ≥ ½, so the total is ≥ 1/8. The padding V ⊕ F_2^e handles ℓ ≢ 0 (mod 3).

## Reproduce

```bash
cd scripts && ./run_all.sh     # python3 + numpy + mpmath; single thread; ~5-6 min
```
It prints `OK` for both outputs.
- `primary/gadget_bounds.py` is the internal review's original computation, on the earlier O(p^-6) version of the claim. It checks the r = 1 identity and the union-bound inequality of the lower bound on random equivariant maps. It runs exact death-chain DPs against the supermartingale bound and against the paper's harmonic/Markov bound, and computes the paper's thresholds (main text and effective recipe). Its "new" column (2^166.7 at p = 10^-8) is the superseded O(p^-6) analysis. The final threshold is computed in the independent script. We changed one line: the printed lower bound is now log2(1/(4p)), matching the proposition (it was log2(3/(16p))).
- `independent/indep_gadget.py` is written separately:
  - **Part A:** exact rational checks of the lower-bound chain on 12 random maps.
  - **Part B:** exact rational rank chains, testing the supermartingale property and the extinction bound.
  - **Part C:** thresholds A/B/C/D/L in the table above.
  - **Part D:** worst-case chains at the remark's parameters, and the route limit.

## Scope and limitations (what is NOT claimed)

- The order of r_*(p) is open between p^-1 and p^-3.
- The only consequence for the Unique Games reduction is the size of its alphabet parameters. Later sections use η = 2^{-r_*}/16, so the alphabet stays astronomically large: about 2^{2^{80}} in the alphabet exponent at p = 10^-8, against 2^{2^{450}}. Nothing in the reduction's statement changes.
- "Before" means the threshold certified by the preprint's proof through the union bound of its Lemma 3.6. The preprint's effective algorithm chooses d by exact enumeration and might find a smaller d. The comparison applies the same union bound on both sides.
- We rely on the correctness of the preprint's Lemmas 3.2–3.8 and Proposition 3.9 (the conditional rank-loss estimates and fresh orientation). No Lean; the upstream Lean formalisation (`lean/docs/102.md`) covers the final reduction, not this constant.
- **Novelty.** To our knowledge (search of GitHub and arXiv up to 2026-10-09: GitHub repository searches for "unique games", "latent alphabet gadget" and "unique games theorem openai"; web search for explicit bounds on the gadget threshold), no explicit r_* has been published. The only related item is lgboim's presentation (p^-8 for ℓ ∈ 3ℤ), acknowledged above.
