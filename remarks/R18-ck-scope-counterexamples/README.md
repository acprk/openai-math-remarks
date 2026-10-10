# R18 — Courtade–Kumar / Hellinger (family 119): scope of the sharp bounds beyond the uniform binary cube

## Premise and confidence (read first)

- **This is not an error in OpenAI's work.** Both preprints of family 119 state their theorems only for the **uniform binary cube** {−1,1}^n with independent bit-flip noise. Every statement below is about extensions that 119 does **not** claim. This is an applicability boundary only.
- **Before → After.** Before: it is natural to expect that the refined Boolean bound (Cor. 1.2), the sharp soft-channel contraction (Thm 1.1), and the Hellinger bound "for every output bias" extend to other finite product spaces. After: **they do not** under the standard resampling noise (Bonami–Beckner operator; q-ary symmetric channel parametrised by its correlation). One coordinate of a uniform q-ary space (q ≥ 3) or of a p-biased bit (p ≠ ½) is already a counterexample, for all small noise levels. Example: q = 4, |A| = 1, λ = 0.1 gives Hellinger loss 0.005464 > 0.005013, and I(f;Y) = 0.004817 > ψ(λ r_m) ≈ 0.004446.
- **Trivial.** This is a second-order Taylor expansion, essentially a local χ² computation, that experts would regard as folklore. We record it only because these extensions are natural and a reader might expect them to hold.
- **Depends on the noise model.** Under independent bit flips applied to p-biased bits, the adversarial reviewer found no violation by exhaustive check for n ≤ 3. The plain Courtade–Kumar curve I ≤ ψ(λ) is not contradicted by any example here.
- **Confidence: high.** The counterexamples are certified with interval arithmetic and the expansions with exact symbolic series. An adversarial AI session reproduced all numbers and corrected the framing to name the noise model; see `VERIFICATION.md`.

## Statement

Let (Ω, μ) be finite and n = 1, and let A ⊂ Ω with μ(A) = w ∉ {0, ½, 1}. Put f = 2·1_A − 1 and m = 2w − 1 ≠ 0. Let R_λ be the resampling channel: Y = X with probability λ, otherwise Y is a fresh sample of μ. For uniform μ on two points this is the BSC with ρ = λ. Then, as λ → 0,

| quantity | expansion | 119's bound (extended) | expansion |
|---|---|---|---|
| I(f(X);Y) | λ²/2 + O(λ³) | ψ(λ r_m), r_m = ψ⁻¹(H(m)) < 1 | λ² r_m²/2 + O(λ⁴) |
| Hellinger loss √(1−m²) − E√(1−(R_λ f)²) | λ²/(2√(1−m²)) + O(λ³) | 1 − √(1−λ²) | λ²/2 + O(λ⁴) |

So the extensions of Cor. 1.2, of Thm 1.1 (with u = f, since I(f) = H(m)), and of the Hellinger theorem all fail for small λ.

| setting | λ | I(f;Y) | ψ(λ r_m) ≤ | Hellinger loss | 1−√(1−λ²) |
|---|---|---|---|---|---|
| q = 4, \|A\| = 1 (w = ¼) | 1/10 | 0.0048170 | 0.0044458 | 0.0054645 | 0.0050126 |
| q = 3, \|A\| = 1 (w = ⅓) | 1/10 | 0.0049341 | 0.0048070 | 0.0052003 | 0.0050126 |
| p-biased bit, p = 3/10 | 1/5 | 0.0193294 | 0.0188546 | 0.0207666 | 0.0202041 |

Violation windows: for w = ¼, (A) fails for 0 < λ < 0.43648 and Hellinger fails for 0 < λ < 0.32273. For w = ⅓ the windows are 0.34953 and 0.33746.

On the uniform binary cube no such witness exists. There, a function of one coordinate has m ∈ {0, ±1}. More generally, a biased Boolean function there cannot put all its variance at level one. Control: f = AND on {−1,1}², which has the same law as q = 4, |A| = 1, satisfies both bounds under bit flips (ρ = 0.1: I = 0.003345 < 0.004446; Hellinger loss 0.003864 < 0.005013).

## Attribution

The theorems are OpenAI's work: family 119 of [github.com/openai/math](https://github.com/openai/math), commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`:
- `preprints/Sharp-binary-information-contraction-on-the-discrete-cube-September-24-2026/` (`build/main.tex`: Thm 1.1 at lines 72–84, Cor. 1.2 at lines 92–104; the uniform measure is fixed at line 53, and line 113 says the input distribution "remains uniform");
- `preprints/Hellinger-contraction-with-arbitrary-Boolean-output-bias-September-24-2026/` (`build/main.tex`: Thm 1.1 at lines 54–61, "with respect to the uniform probability measure").

This remark changes no line of either preprint.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Proof

Write δ = λ(f − m), so Eδ = 0 and Eδ² = λ²(1 − m²). For smooth Φ, Φ(m) − EΦ(m+δ) = −½Φ''(m)·Eδ² + O(λ³). Then H''(m) = −1/(1−m²) gives λ²/2, and (√(1−t²))'' = −(1−m²)^{−3/2} gives λ²/(2√(1−m²)). Also ψ(x) = x²/2 + O(x⁴), and r_m < 1 because H(m) < log 2. In general, a concave profile Φ contributes the coefficient ½(1−m²)|Φ''(m)|. For the entropy this equals ½ for every m, which is why the plain CK curve is matched to second order but not exceeded.

## Reproduce

```bash
cd scripts && ./run_all.sh     # python3 + mpmath + sympy (+ numpy 2.2.6 for byte-identical primary output); < 1 min
```
It prints `OK` for each of three outputs.
- `primary/qary_refined_boolean.py` is the internal review's scan of q = 3..9, k = 1..q−1, λ = 0.1..0.9, single-coordinate Boolean f. It lists the violations of the extended Cor. 1.2. The q = 4, k = 2 line with gap 2.5e-32 is round-off at the equality case m = 0, not a violation.
- `primary/hellinger_counter.py` is the internal review's Hellinger scan for w ∈ {0.25, 0.3, 0.4, 0.45}.
- `independent/indep_scope.py` is written separately. It contains:
  - interval-arithmetic certificates of the three witnesses (exact rational inputs);
  - exact sympy series;
  - the violation windows;
  - the AND control.

## Scope and limitations (what is NOT claimed)

- No error in 119 is claimed. Both preprints restrict to the uniform binary cube.
- The failure is under the resampling noise operator, which is the standard one for q-ary and p-biased spaces. Under independent bit flips on p-biased bits we found no failure (n ≤ 3, adversarial reviewer's exhaustive check), and nothing is claimed for that model.
- Nothing is claimed about the plain Courtade–Kumar inequality I ≤ ψ(λ) on general product spaces. Our examples satisfy it.
- No Lean.
- **Novelty.** To our knowledge (search of GitHub and arXiv up to 2026-10-09), no prior statement of these counterexamples exists. The adversarial reviewer checked Mahdavifar–Beirami arXiv:2609.26444 (full text), Ky–Tran 2609.24184, Chen et al. 2609.24931, Kindler–O'Donnell–Witmer 1506.03167, Javanmard–Woodruff 2601.09679, Yu 2410.10147 and Weinberger–Shayevitz 1807.11289. We checked the GitHub repositories dpwoodru/general-courtade-kumar-lean, hessammahdavifar/most-informative-bit and deepinfra/courtade-kumar; all treat uniform inputs. The deepinfra notes explicitly exclude "an extension to biased inputs/asymmetric noise". The nearest statement is Mahdavifar–Beirami's remark that their argument needs the uniform input law. Given how elementary the computation is, it may well be folklore.
