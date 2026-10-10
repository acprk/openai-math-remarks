# R5 — Thorp shuffle: 1600d → 24d

*Using only the main paper (no companion results) the same bookkeeping gives 45d (route O, c ≈ 44.3006).*

**One-line result.** With the same method as OpenAI's *Optimal-order mixing of the Thorp shuffle*, the
Thorp shuffle on n = 2^d cards satisfies ‖q_d^{*⌊cd⌋} − U‖_TV → 0 for every
c > 2/log₂(8/7) + 2/log₂(10/9) = 23.5394…, so **24d** shuffles suffice for large d. The paper states
1600d. Using only the main paper's own operator-norm transfer, the same bookkeeping gives every
c > 44.3006… (**45d**).

> History: the selection table originally listed the 45d (paper-internal) version, and this directory was
> then called `R5-thorp-45d`. Combining the stage-1 improvements with the two-factor Hilbert–Schmidt
> transfer from the **same family's companion preprint** gives 24d; route H has since passed an
> independent adversarial review (see `VERIFICATION.md`). Both statements are proved in `note.pdf`.

## Statement (Before / After)

| | Before | After |
|---|---|---|
| Main paper, Thm 1.1 | ‖q^{*1600d} − U‖_TV → 0 | any c > 23.5394 (route H), or any c > 44.3006 (route O, main paper only) |
| Best explicit constant elsewhere in family 238 | 512d (*Random coordinate frames…*, Cor. 11.4, label `shear:cor:512`) | 24d (×21.75) |
| Lower bound (unchanged) | 2d − O(1) (support count) | 2d − O(1) |

General form (note, Theorem 1.1). Write L(b) = log₂(2b/(2b−1)), κ(b) = 2/L(b) and e(b) = 1 − 2/b.
For integers b_i ≥ 3, any c > Σκ(b_i) works if
- **(H)** Σe(b_i) > 1 (one factor is bounded with the companion's Hilbert–Schmidt bound), or
- **(O)** Σe(b_i) > 2 (all factors use the main paper's operator-norm bound).

The optima over integer block numbers are (4,5) → 23.5394 for (H) and (3,4,5,5) → 44.3006 for (O).
Both are certified with interval arithmetic.

## What is ours and what is not

The family's 512d (*Random coordinate frames…*, Cor. 11.4, `shear:cor:512`: 2 × 256d, eight equal blocks)
**already** combines stage 1 with the companion's two-factor Hilbert–Schmidt transfer. So the ×4 step
of the ledger below (eight factors → two) is the family's, not ours. What this remark adds is:
1. three improvements of stage 1: concentration loss θ → 0, the KL chain rule, and the local look-back;
2. different block numbers in the two factors, {4,5}.

Against the family's best explicit constant this is 512d → 23.54d, a factor **21.75**.
Against the main paper's 1600d it is ×67.97, but ×4 of that is the family's own two-factor transfer.

## Attribution

The original construction and proof are OpenAI's work: family 238 of
[github.com/openai/math](https://github.com/openai/math), commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, in particular
`preprints/Optimal-order-mixing-of-the-Thorp-shuffle-September-26-2026/` (main paper) and
`preprints/From-partial-permutation-information-to-Fourier-bounds-September-26-2026/` (companion; its
Theorem 1.1, estimate (1.2), is the Hilbert–Schmidt transfer used in route H, stated there for any
number of blocks of any sizes and every u > 0). This remark only re-does the bookkeeping of constants
inside these arguments. It takes the correctness of both preprints as given.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are
not human peer review.

## What changes (main paper `build/main.tex` line numbers)

| Where | Paper | Here | Gain |
|---|---|---|---|
| l. 325–339, eq. (4.2)–(4.3) | half-density ≥ 1/4 of m (β = 1/16, θ = 1/2) | any θ → 0: β → m/n (Hoeffding without replacement) | ×1.97 |
| l. 358–401, Lemma 4.1 / Prop. 4.2 | sequential TV + Cauchy–Schwarz: needs n³E_T → 0 | KL chain rule + χ² + Pinsker: needs n²E_T → 0 | ×1.46 |
| l. 343–356, eq. (4.4) | E_t ≤ γ^{t−2d} (offset 2d) | initial term in the energy identity + look-back to time 0 on subcubes: offset O(log d) | ×1.09 |
| l. 419–454, Lemma 5.1 | inverse-degree exponent 1 | any exponent σ > 0 (same three-case proof, or Liebeck–Shalev Thm 2.6) | (route O) |
| l. 456–539, Sec. 6 | 8 equal blocks | b blocks of sizes ⌊n/b⌋ or ⌈n/b⌉, b chosen per factor | ×1.82 (H) |
| l. 570–581, Sec. 7 | ‖B‖_HS ≤ √D‖B‖_op for all 8 factors | companion's HS bound for one factor: 2 factors suffice (H) | ×4 (H) — **already in the family's 512d, not ours** |

Full ledgers (both routes): `scripts/expected_output/constants_certificate.txt` and note Table 1.
The largest single losses in the original are (i) using eight factors where the family's own
companion shows two suffice (×4; already removed in the family's 512d), (ii) the slack choice T = 200d versus the 67.5d the paper's own
inequalities need (×2.96), and (iii) the concentration constant (×1.97).

## Proof sketch

1. **Stage 1** (note §2). Energy identity with the initial term for t < d (Lemma 2.1). Contraction
   b_s ≤ (1−β)E_s + βη_s for **all** s, with look-back to time max(0, s−d+1) on subcubes
   (Lemma 2.2). Closed-form solution of the recurrence (Lemma 2.3):
   E_t ≤ (2/(1−ρ))γ'^{t−s₀} + γ'^{t−d}/(m(1−ρ)) + η*, where s₀ = ⌈log₂(Kd)⌉. KL chain rule (Lemma 2.4).
   Together: a factor of ⌊c_b d⌋ shuffles, c_b > κ(b), makes average marginals of lists of length
   n − ⌊n/b⌋ o(1)-close to uniform.
2. **Stage 2** (note §3). Trimming as in the paper. Operator bound (paper Prop. 6.1 with D^{1/b}, C_σ),
   and the HS bound (companion Thm 1.1 (1.2); a short proof is included).
3. **Completion** (note §4). D‖∏ν̂_i‖²_HS ≤ D‖ν̂_1‖²_HS ∏_{i≥2}‖ν̂_i‖²_op ≤ A·D^{−ε}, and
   Σ_{ρ≠1,sgn} D^{−ε} → 0; the sign coefficient and the trimming are handled as in the paper.

**Frames.** The main paper returns to the physical frame only after blocks whose length is a multiple of d,
and ⌊c_i d⌋ is in general not a multiple of d. This is harmless. Since Y_T = R^T X_T, the physical law of a
block is the cyclic-frame law left-multiplied by the fixed position permutation R^T. Left multiplication by
a fixed permutation preserves every marginal TV distance to U_inj (it relabels positions bijectively) and
maps each coset gH_a to the coset R^T g H_a, so the coset caps are unchanged. Proof in note §4.

## Reproduce

```bash
cd scripts && PY=python3 ./run_all.sh     # python3 + mpmath>=1.3 + numpy; single thread; ~2-4 min
```
Prints `OK` for each of the five outputs:
- `constants_certificate.py` — **the only numerical input to the proof.** Exhaustive enumeration of
  block multisets (exact rationals) and interval evaluation (mpmath.iv) of all constants. It also
  certifies optimality, the σ-windows, the lower bounds of §"Scope" (integer points 20.7636 at b = 4;
  fractional relaxation infimum 20.738, b* ≈ 3.8647) and both ledgers.
- `energy_lemma_check.py` — illustration only. Iterates the extremal recurrence allowed by Lemmas 2.1–2.2,
  checks the closed form of Lemma 2.3 at every step, and shows the O(log d) offset (about +7 steps per
  doubling of d; without look-back the offset is about 1.1–2d).
- `small_d_exact.py` — unit test in exact rationals: t_mix(2) = 3, t_mix(3) = 5.
- `energy_mc.py 8 4 4000` — seeded Monte Carlo of the conditional flow: the identity of Lemma 2.1
  holds and the decay rate is ≈ the predicted root.
- `independent/indep_constants.py` — the adversarial reviewer's own re-derivation of the constants
  (separate code; does not import `constants_certificate.py`).

## Scope and limitations (what is NOT claimed)

- **Asymptotic only.** No explicit d₀: C_σ for the small σ needed (σ < 2/9 for (4,5)) is finite but not
  computed, and Σ D^{−ε} → 0 is slow for small ε.
- **Not a sharp constant.** The truth is between 2d − O(1) and 24d. With the present ingredients the constant is
  ≥ 20.73d (route H) and ≥ 41.47d (route O): every integer multiset costs more than 2κ(4) = 20.7636
  (resp. 41.527), and even the fractional relaxation over real b has infimum min κ/e = 20.738 (b* ≈ 3.8647),
  i.e. 41.476 for route O. These are strict lower bounds, not attained limits.
- **No Lean.** Nothing here is formalised. The official Lean covers the original 1600d statement and the
  companion's coset-to-Fourier estimate.
- **Premise.** We rely on the correctness of the main paper and of the companion's Theorem 1.1.
- **Novelty.** To our knowledge (search of GitHub and arXiv up to 2026-10-09: openai/math family-238
  preprints; GitHub repo and code search for "Thorp shuffle"; arXiv API query "Thorp shuffle"), no explicit
  constant below 512d has been published for this model. The only related community repository,
  yunjiangster/thorp-shuffle-mixing (2026-10-03), proves a 5.15·d^{5/3} bound, which is weaker in order. Most
  of the steps are standard tightenings; the remark is about the constant, not about a new method.

## Relation to earlier internal estimates

- **CONST-OPT (45d, adversarially reviewed).** This is route O. It kept the main paper's final step
  (‖B‖_HS ≤ √D‖B‖_op for every factor). That is why it needed four factors, and why its stated "framework
  bound" 44.30 / ≥ 41.47 applies only to that architecture. Its claim is reproduced and certified here
  (Corollary 1.2, route O). One fix: the earlier interval check certified (4,4,4,5) = 44.30299 but not
  the headline (3,4,5,5) = 44.30061. Both are now certified.
- **K3 (26.32d, proof not written out).** This estimate already used a two-factor HS transfer with equal
  blocks (b = 5, 2 × 13.16d). It also included non-asymptotic machinery: fixed-list low-irrep inputs and a
  modified, Plancherel-aggregated transfer whose proof was never written. Route H keeps only the
  asymptotic part, which needs nothing beyond the companion's published Theorem 1.1. Allowing
  *different* block numbers per factor ((4,5) instead of (5,5)) lowers 26.32 to 23.54. K3's
  finite-d numbers (e.g. TV ≤ 1/4 after 27–30d at d = 20–64) are **not** part of this remark.
