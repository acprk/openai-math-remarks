# R20 — Subset Sum: 2^{0.49n} → 2^{0.4889n} (dense instances: 2^{0.4741n})

## Premise and confidence (read first)

> **Everything here assumes that OpenAI's preprint *Subset Sum in Time O(2^{0.49n})* (family 138) is correct.**
> That preprint has no formal proof and, as far as we know, no published peer review. Nothing below re-proves it.
> If it is wrong, this remark says nothing.

- **Kind of result: parameter tuning plus one small variant.** The proof is the preprint's proof. We change parameters,
  one sampling rule in alias extraction (the "σ-variant"), and the filtering densities, and we redo the exponent
  bookkeeping. No new idea is involved, and by the project's own criteria this is **not** a new mathematical result.
- **Framework limit.** A heuristic calculation for the same framework (alias extraction + Fourier checksum) gives
  a limit of about **0.4798**. The preprint uses 49.5% of the headroom 1/2 − 0.4798; this remark uses 55%.
  Going much below 0.4889 would need a new ingredient (more selective filters, or a less wasteful representation),
  not more tuning.
- **What is certified.** The two filtering-moment inequalities, the only numerical input, are certified with
  interval arithmetic (mpmath) and again in ball arithmetic by an independent implementation (python-flint).
  All exponent inequalities are checked in interval / exact rational arithmetic.
- **Asymptotic only.** As in the preprint (ε = 10⁻⁸), the bound applies only for astronomically large n. No
  explicit N_c is given.

**One line.** Before → After: **O(2^{0.49n}) → O(2^{0.4889n})** in the preprint's word-RAM model, inputs of
polynomial bit length. The saving below meet-in-the-middle 2^{n/2} grows from 0.0100 to 0.0111 in the exponent;
the running-time ratio is 2^{0.0011n}. If in addition max{t, a_i} ≤ 2^{n/2}, the bound is **O(2^{0.4741n})**. With
the preprint's own densities the same dense-case argument gives about 2^{0.4748n}. The preprint does not state
a dense-case bound.

Appendix (small, elementary): **Proposition C.** If K ≤ G = (ℤ/p^e)^k has index M, then
|K ∩ G[p]| ≥ p^k · M^{−1/e}. Equality holds iff G/K ≅ (ℤ/p^e)^s. Hence any random homomorphic reduction of G to
a group of order ≤ M has collision probability ≳ M^{−1/e} on some fixed element of order p. So the preprint's
random-prime collision bound (Lemma 2.5, ≈ M^{−1}) has no homomorphic analogue over (ℤ/p^e)^k, e ≥ 2. This is
an obstruction to one proof step, not a lower bound for any problem, and we do not claim the fact is new.

## Statement (Before / After)

| | Before (preprint) | After (this remark) |
|---|---|---|
| Worst-case time, b ≤ n^c | O(2^{0.49n}) (all-history exponent 0.489995) | **O(2^{0.4889n})** (all-history exponent 0.4888802) |
| Dense instances, max{t, a_i} ≤ 2^{n/2} | not stated; its Thm 4.1 alone gives ≈ 2^{0.4748n} (numerical) | **O(2^{0.4741n})** (checksum run at τ = 0.47409, certified) |
| Alias extraction | shared block split m/2 : m/2, n^{10} trials, per-list caps | split σm, σ = 0.38; T = 2^{(1−h(σ))μn+O(log n)} trials; cumulative caps |
| Filtering densities | degree-6 r, s; one tilt κ = 1.0315 | degree-8 r, s; tilts κ_f = 1.034, κ_g = 1.0339 |
| (η, τ, λ₁, λ₂) | (0.51772, 0.48996, 8.5, 8.3) | (0.5176, 0.48888, 8.323, 8.132) |
| Filtering moments (continuous) | < 2^{0.5168} (our recheck: 2^{0.514974}) | < 2^{0.5170788} (needs < η − 48ε − slack) |

## Attribution

The algorithm, its analysis and every lemma used are OpenAI's work: family 138 of
[github.com/openai/math](https://github.com/openai/math), commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/`. This remark only re-tunes parameters, changes one
sampling rule and the filtering densities, and redoes the bookkeeping. It takes the preprint's correctness as
given.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction;
independent adversarial reviews were run by separate AI sessions and are not human peer review.

## What changes (preprint `build/sections/*.tex`, line numbers)

| Where | Preprint | Here |
|---|---|---|
| `setup.tex` l. 150–206, Prop. 2.4 | δ = 0.01004, cap 2^{⌈0.489985n⌉} | δ = 0.0112, cap 2^{⌈0.488825n⌉} (= 1/2 − δ + 1.5·10⁻⁵ + 10⁻⁵, the last term absorbs the +1/2 as in the preprint) |
| `filtering.tex` l. 24–29, 61–81 | κ = 1.0315; η = 0.51772, τ = 0.48996; degree-6 r, s | κ_f = 1.034, κ_g = 1.0339; η = 0.5176, τ₀ = 0.48888; degree-8 r, s (note §3) |
| `filtering.tex` l. 90–123, Prop. 3.1 / Lemma 3.2 | moments < 2^{0.5168} / 2^{0.5169}, weights 8.5, 8.3 | < 2^{0.5170788} / 2^{0.51718}, weights 8.323, 8.132 |
| `checksum.tex` l. 438–491, (4.22)–(4.25) | rounding 8.3(4a+3)ε < 30ε; error exponent −0.00082 | 8.132(4a+3) = 26.43 < 30; error exponent < −0.00042 |
| `aliases.tex` l. 28–30, 72–77 | η = .51772, m = 2⌊.0567n/2⌋, γ = .48925 | η = .5176, m = 50⌊.0877n/50⌋, γ = .3714 |
| `aliases.tex` l. 92–101, 120–143 | n^{10} trials, \|J∩M\| = m/2, masks k and m/2 − k | T = n^{10}⌈2^m / C(m, σm)⌉ trials, \|J∩M\| = σm, masks k and σm − k |
| `aliases.tex` l. 197–253, (5.8), (5.11) | lists 0.489909, pairs 0.489849 | lists 0.488777, pairs 0.488856 (**incl. guard +10⁻⁵**) |
| `aliases.tex` l. 260–273, (5.12) | each list and pair count capped at 2^{⌈τn⌉} | **one cumulative budget** 2^{⌈τ₀n⌉} over all trials |
| `aliases.tex` l. 327–352, (5.15) | 17 blocks of size m/2; Markov at 0.0106 | K = 15 blocks of size σm; Markov at θ′ = 0.0226 > θ = 0.022405 |
| `aliases.tex` l. 397, (5.16) | .4894 > γ | σ(1 − θ′) = 0.371412 > γ |
| `implementation.tex` l. 125–150, 274–290, (6.4) | polynomially many prime draws; all-history 0.489995 < 0.49 | 2^{ρn}·poly draws of P′ (total cost 2^{0.01997n}, coupling failure still o(1)); all-history 0.4888802 < 0.4889 |

**Repairs of the three holes found by the earlier adversarial review** (CONST-OPT round):
(A) with exponentially many trials, the caps must be cumulative (otherwise the cap is multiplied by T);
(B) the visited-pair exponent must carry the bit-guard term +10⁻⁵ (0.488846 → 0.488856);
(C) the checksum must be certified at the exponent at which it runs, τ₀ = 0.48888 (λ = 8.323 / 8.132; margin
5.2·10⁻⁴ in the continuous moment), not at the headline 0.4889. All three are in note §2–3 and in the scripts.

## Where the gain comes from

| Step (numerical, from the CONST-OPT notes; not certified) | Exponent |
|---|---|
| Preprint | 0.49 (0.489995 before slack) |
| (i) re-optimize λ, η, μ, δ with the preprint's densities | 0.48977 |
| (ii) σ-variant of alias extraction (overlap waste 0.1226 → ≈0.06 per unit μ) | 0.48911 |
| (iii) new degree-8 densities (checksum curve τ_cs(0.5176): 0.48955 → 0.48882) | 0.48881 |
| Certified operating point (this remark) | run at 0.48888, claim **0.4889** |

At the operating point the high-η regime is tight in three places at once: the preliminary search, the alias
lists and the alias pairs all meet the checksum curve τ_cs(η) near η ≈ 0.518. The preprint's λ-analysis was
already essentially tight (the exact Varadhan dual improves it by < 10⁻⁵). The heuristic framework limit
≈ 0.4798 is where an idealized checksum meets the alias-count lower bound 1 − η, near η ≈ 0.520. The gap to it
comes from filter selectivity (about 73% of the KL bound) and from representation waste.

## Reproduce

```bash
cd scripts && ./run_all.sh     # python3 + mpmath>=1.3 + numpy + python-flint; single thread; ~6 min
```
Prints `OK` for six outputs:
- `chain_check.py` — every exponent inequality of the modified algorithm (interval arithmetic), plus the
  preprint's chain for comparison.
- `hardcore_certificate.py main|dense|paper` — rigorous upper bounds for the two filtering moments (mpmath
  intervals, N = 128, analytic trapezoid bounds as in the preprint's Appendix A, κ split in two). `paper`
  rechecks the preprint's own Proposition 3.1 (2^{0.514974} < 2^{0.5168}).
- `theorem_c_exhaustive.py` — Proposition C over all subgroups of (ℤ/4)², (ℤ/4)³, (ℤ/8)², (ℤ/16)², (ℤ/9)²,
  (ℤ/27)², (ℤ/25)².
- `independent/indep_check.py` — independent re-implementation (python-flint balls, N = 192, different strips,
  scores from the definitions of D, A, E; exact rationals for the chain; subgroups as kernels of all
  endomorphisms).

## What is NOT claimed

- **Not a proof of the preprint.** We neither re-verify nor formalize the preprint's Sections 2–6. The premise is
  explicit.
- **Not a new method.** Parameter tuning plus a one-parameter variant. The exponent 0.4889 is not optimal even
  within the framework (independent re-optimization gives ≈ 0.48884), and the heuristic framework limit 0.4798
  is not a theorem.
- **No explicit constants or thresholds.** O(·) hides C_c and N_c, which are astronomically large.
- **Dense corollary** concerns max{t, a_i} ≤ 2^{n/2} only, and is a direct consequence of the preprint's
  Theorem 4.1. We state it because the preprint does not, not because it is deep.
- **Proposition C** is an elementary fact about finite abelian p-groups. It only shows that one step of the
  preprint (Lemma 2.5) has no homomorphic analogue over (ℤ/p^e)^k. It says nothing about the complexity of
  subset sum over groups.
- **No Lean.** Nothing here is formalized.

## Duplicate search (to 2026-10-09)

To our knowledge (search of GitHub and arXiv up to 2026-10-09), no exponent below 0.49 for worst-case Subset
Sum has been published. Searched: GitHub repository search for "subset sum" (sorted by update; nothing related
to the 0.49 algorithm), code search for "0.49n" / "2^{0.49n}" / "0.4889" (no hits); arXiv API, most recent 40
entries for "subset sum". Related and acknowledged: **L. Ye and B. Li, *Subset Sum via Partial Match*,
arXiv:2610.10446 (2026-10-07)**, which independently breaks 2^{n/2} with O(2^{0.499999n}) by a different method
(representation technique of Randolph–Węgrzycki + partial match). Also: Chen–Jin–Randolph–Servedio
(arXiv:2301.07134), 2^{n/2}/n^{0.5023}. Randolph–Węgrzycki (arXiv:2511.10823, 2604.04656) treat Subset Balancing.
Proposition C is folklore-level group theory; we did not search for it.

## Files

- `note.tex`, `note.pdf` — 4-page note: statements, the σ-variant with proofs, new densities, Proposition C,
  ledger.
- `scripts/` — see Reproduce; `params.json` holds every constant as an exact decimal.
- `VERIFICATION.md` — who checked what, the adversarial review, and the remaining gaps.
