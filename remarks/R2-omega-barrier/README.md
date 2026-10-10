# R2 — The 9/4 barrier for OpenAI's ω ≤ 9/4 persists under an enlarged system of inequalities

*(The basic barrier — the least profile P\* — was observed first in eliottcassidy2000/math THM-4568; see Prior work.)*

**One-line summary.** OpenAI's 9/4 argument (openai/math family 107) uses six properties of a symmetrized profile `P(a,b)`. These properties are satisfied by an explicit profile

```
P*(a,b) = φ(m) · (M + (m−1)/2),   m = min(a,b),  M = max(a,b),
φ(m)    = ∏_{j=1}^{m−1} (1 + 1/(3j)),
```

and that profile has exponent exactly 4/3. So these properties alone can never prove ω < 9/4.

A larger set of valid character inequalities, found with an OpenAI Codex session, also cannot force ω < 9/4. P* is the pointwise least profile, which was **already observed publicly**, with a proof (see Prior work). What this remark adds is an independent re-derivation written as one self-contained proof (its only new ingredient is the elementary bound `φ(m)³ ≤ 2m−1`), a precise statement of the scope, and above all the enlarged system, which is joint work with Codex.

## Precise statement (before / after)

**Before (OpenAI, Lemma 5.1 + proof of Thm 1.1).** Let `P : ℕ²₊ → ℝ>0` satisfy the following. Positivity, (Pos), is the sixth property, at l. 427 and l. 611. Line numbers refer to `build/paper.tex`.

| axiom | statement | paper.tex lines |
|---|---|---|
| (S) | `P(a,b) = P(b,a)` | 427, 611 |
| (B) | `P(1,b) = b` | 427, 611 |
| (C) | `2P(a,b) ≥ P(a,b+1) + P(a,b−1)` for `a ≥ 1, b ≥ 2` | 442–447, 612 |
| (T) | `P(a, 3h+a−1) ≥ 3P(a,h)` for `a,h ≥ 1` | 526–530, 613 |
| (U_t) | `P(a,b) ≤ (a+b−1)^{1/t}`; only the diagonal is used | 428, 657 |

Then `P(a,a) ≥ a^{4/3}`, hence `t ≤ 3/4` for every character. The detecting-character lemma then gives `ω ≤ 9/4`.

The axiom list is complete. Lemma 5.1 (l. 610–650) has only these hypotheses, and the final step (l. 652–660) adds only (U_t). The OpenAI Lean formalization has the same structure. `structure ScalarProfile` (`Growth/Profile.lean` l. 23–33) has exactly these six fields. `meanExponent_le_three_quarters` (`Polynomial/Inequalities.lean` l. 97–114) reaches t ≤ 3/4 only through `toScalarProfile`.

**After (this remark, note.tex Theorem 1).**
1. P* satisfies (S), (B), (C), (T) and (U_{3/4}), with equality in (U) only at (1,1).
2. Every P satisfying (S), (B), (C), (T) satisfies P ≥ P* pointwise. So P* is the least element.
3. Hence, for t > 0, the system Ax(t) is satisfiable **if and only if t ≤ 3/4**.
4. Every step of Lemma 5.1, up to `H_a ≥ ∏(1+1/(3m))`, is an equality on P*.
5. Only the two final simplifications are strict. They lose the constant `3/(2Γ(4/3)) ≈ 1.68` but not the exponent.

**Enlarged barrier (note.tex Theorem 2, joint with Codex).** Consider the balanced assignment

```
λ*(C(a,b)^π) = P*(a,b)^{3/4},   λ*(B_L(m)) = m^{3/4},   λ*(⟨e,h,l⟩) = (ehl)^{3/4}.
```

It satisfies, **for all sizes**, every family marked P in note.tex Table 1. The families are:

- **Nonsymmetrized entropy forms.** The nonsymmetrized entropy forms of the paper's lemmas.
- **Chains.** Alternating k-sector chains for odd k.
- **Determinant filtrations (new).** Three new determinant-filtration degenerations:
  - (A) transfer along `a+b = const`;
  - (B) `C ⊗ B_X(m)`;
  - (C) bivariate shared-output.
- **Full direct sum (D).** The bivariate full-direct-sum degeneration.
- **Fourier filtering** with its cost kept.
- **Outer-product and dot-product bounds.**
- **Wrap-around** to `ℂ[ℤ_n]`.
- **Scalar peeling.**
- **Small-CW laser.**
- **Scalar targets** from any tensor product of generators.
- **Single matrix block.** Monomial single-matrix-block extractions.
- **`T_m` relations.** `T_m ⇝ ⊕ C`, and `T_m ⊗ C(n,n)` for m ≥ 2.
- **Zero/rest laser.** The zero/rest laser partition.

## What is NOT claimed

- λ* is **not** claimed to be a tensor character. If some character had `λ(T_m) = m^{9/4}` for some m ≥ 2, then ω ≥ 9/4, and together with OpenAI's result ω = 9/4. Nothing here proves or refutes that.
- This is **not** a barrier for polynomial-multiplication tensors in general. It concerns the stated axiom systems only.
- Several routes remain **open**:
  - multi-block laser extraction from powers of `C(n,n)` after non-monomial basis changes;
  - higher-power repartitionings beyond the listed families;
  - `T ⊗ C` relations outside the proved sufficient condition;
  - relations that genuinely use asymmetry.
- Some families are checked only on finite ranges:
  - Kronecker restriction: exact checks;
  - interval laser recipes for n ≤ 30: interval-certified.

## Attribution

- **OpenAI.** The original construction and proof of ω ≤ 9/4 are OpenAI's work: family 107, [*An Upper Bound of 9/4 for the Matrix Multiplication Exponent*](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Matrix-Multiplication-Nine-Fourths-October-2-2026/build/paper.tex), commit `adc7f1241b42`, unchanged at HEAD `08fcd618213a`. This remark only analyses which properties the argument uses and shows that those properties stop at 9/4.
- **Prior work: the least profile is not new.** It was recorded publicly on **2026-10-07**, with an exact closed form, a proof, an LP check and an audit, in [eliottcassidy2000/math, THM-4568](https://github.com/eliottcassidy2000/math/blob/5966360a9c69b6bb0f3da39dc219f4a6c5dc3244/01-canon/theorems/THM-4568-the-nine-fourths-growth-lemma-is-exactly-sharp.md) (commit `450388a22f72`). It was also noted numerically on **2026-10-08** in [GasStationManager/SafeVerifyAgent `audits/mm94-intake/READ-F.md` §2d](https://github.com/GasStationManager/SafeVerifyAgent/blob/af3b9e1e8853b8bfc4033f18edbafdd8c1477e0f/audits/mm94-intake/READ-F.md). We found these after doing our own derivation. Our Theorem 1 should be read as an independent re-derivation with a written proof.
- **Complementary work: all fields.** [selanavot/matrix-multiplication-all-fields](https://github.com/selanavot/matrix-multiplication-all-fields) (Sela Navot, with GPT-6 Astra / GPT-6.1 Sol) extends OpenAI's Lean proof of ω ≤ 9/4 from ℂ to **every field**. That work shows the 9/4 argument is robust across fields. This remark shows that 9/4 is the limit of the argument's profile axioms. We do not repeat any of the field-port work.
- **OpenAI Codex collaboration.** The enlarged barrier (note.tex §4; Table 1 rows A–D, F, scalar, 1-block, T→C, T⊗C, wrap, peel, CW, 0/rest, interval) was produced jointly with an OpenAI Codex session: model `gpt-6-astra`, Codex CLI v0.157.0, three rounds on 2026-10-09.
  - Claude wrote the briefs and supplied the P* witness, which was found earlier and independently in THM-4568. Claude also supplied the unequal-h chain and Kronecker rows; the equal-h chain with its tangent-line argument is already in THM-4568 as its "gadget barrier". Codex proposed and proved the other families listed above and wrote their verification programs.
  - Claude independently re-implemented the checks and re-ran Codex's programs.
  - The Codex reports, which contain the proofs of the extension families, are in [`proofs/`](proofs/). Codex's verification programs are copied verbatim, with a header, as `scripts/codex_round{1,2,3}_verify.py`.

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction. Independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Proof sketch

The full proof is in note.tex §2.

- **(C).** The row `P*(a,·)` is linear for b ≥ a. For b < a its increments are `φ(b)(a/(3b) + 2/3)`. These are non-increasing because the comparison reduces to `2(a−b−1) ≥ 0`. At the diagonal the increments match, since `φ(a−1)(3a−2)/(3(a−1)) = φ(a)`.
- **(T).** Equality holds for h ≥ a−1. For h < a, view both sides as functions of a′ ≥ h. Their difference L − R is 0 at a′ = h and has non-negative increments `φ(a′)(2+h/a′) − 3φ(h)`.
- **(U_{3/4}).** By induction `φ(m)³ ≤ 2m−1`, using `(2m+1)(3m)³ − (2m−1)(3m+1)³ = 9m²+7m+1`. Then `P* ≤ (M+(m−1)/2)(2m−1)^{1/3} ≤ (a+b−1)^{4/3}`. Wendel's inequality also works, but it is not needed.
- **Minimality.** The paper's own chain gives `P(a,a) ≥ (3a−1)φ(a)/2 = P*(a,a)`. Off the diagonal, (C) and (T) give `P(a,h+1) ≥ P(a,h)(2h+a+1)/(2h+a−1)`, which telescopes to `P ≥ P*`.

## Reproduce

All scripts use Python 3 (`/usr/bin/python3`) and one thread. The core scripts use only the standard library, and all their decisions are exact rationals.

```sh
cd scripts
python3 r2_profile_axioms.py          # ~30 s; Theorem 1 checks (a,b ≤ 200)
python3 r2_extension_families.py      # ~50 s; Table 1 inequality families (independent re-implementation)
python3 r2_tensor_constructions.py    # ~25 s; coefficient-level checks of the degenerations
OMP_NUM_THREADS=1 python3 r2_least_profile_lp.py   # numpy+scipy; FLOAT LP cross-check (not a proof)
python3 codex_round1_verify.py        # sympy; Codex's own round-1 program
python3 codex_round2_verify.py        # stdlib; Codex round 2
OPENBLAS_NUM_THREADS=1 python3 codex_round3_verify.py [--all-cuts]   # numpy+scipy+mpmath; Codex round 3
```

Each script has an expected output in `scripts/expected_output/`. Every exact script ends with `ALL PASS` or the equivalent PASS lines.

## Files

- `note.tex`, `note.pdf`: a 6-page note with the statement, the full proof, the enlarged-barrier table with the status of each row, and references.
- `scripts/`: the checking programs and `expected_output/`.
- `proofs/`: the Codex reports for rounds 1–3, containing the validity proofs of the extension families.
- `VERIFICATION.md`: who verified what, with commands, results and gaps.

## Premises

- The paper's Lemmas 4.1 and 4.2 and Corollary 3.2 are taken as correct. They are formalized in Lean in openai/math.
- The Codex families additionally use standard facts: the tight-support laser method (Le Gall 2014, Thm 4.1), the small-CW laser value, and the border rank of truncated multiplication. These are stated in the Codex reports.
- To our knowledge, based on a search of GitHub code and repositories, arXiv (including arXiv:2610.10494) and the web up to 2026-10-09, the enlarged-barrier families (A–D, peel, 1-block, T→C, 0/rest) are not recorded elsewhere.
