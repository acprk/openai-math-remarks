# R24 — Family 179: the norm-m form of the alternating-product argument (parallel generalization)

## Premise and confidence (read first)

**This is a parallel generalization, and we label it as such.** OpenAI's circulant Hadamard preprint
(family 179) proves that for x ∈ ℤ[i][C_{u²}] with xx* = u², the alternating product Δ(x) of
character values is a root of unity of odd order. It then runs a first-order 2-adic comparison. The
preprint itself remarks that the proposition holds for *every* x with xx* = u². Its proof never uses
that u² is a square. Everything below replaces u² by an arbitrary norm m (odd, or prime to 3 over
ℤ[ω]). Each alphabet then gets a one-line normalization, and the comparison is run once. The
project lead's own verdict is that this is a parallel generalization that an expert could write in an
afternoon from the original paper, and that no new mathematical idea is involved. We publish it only as
a formally checked record of what the method gives as it stands.

- **Correctness:** high. Lean 4 proofs with 0 `sorry`, using only the three standard axioms, built on the
  unmodified upstream formalization. Small cases were cross-checked by independent exhaustive search.
- **Novelty:** low. The method, and almost every step of every proof, is OpenAI's.
- Related work by the same authors is in preparation.

**One line.** Before: the 179 argument is stated and formalized for binary sequences, i.e. ℤ[i], norm u²
and local part C₄. After: the same argument, formalized in Lean for norm m, gives
(Q) no perfect quaternary sequence of length 2m (m > 1 odd),
(T3) no perfect three-phase sequence of length 3m (m > 1, 3 ∤ m),
(T1) no perfect binary array / Menon–Hadamard difference set on C₂ × C₂ × C_{u²} (u > 1 odd), and
(F) a rank-one identity for x ∈ ℤ[i][C_{pq}] with xx* = pq.

## Statement (Before / After)

| | Before (OpenAI, family 179) | After (this remark) |
|---|---|---|
| Alternating product (Prop. 3.4, `prop:alternating`) | x ∈ ℤ[i][C_{u²}], xx* = u²: Δ(x) is a root of unity of odd order | **Prop′:** x ∈ ℤ[i][C_N], N > 1 odd, xx* = N; **Prop″:** x ∈ ℤ[ω][C_N], 3 ∤ N, xx* = N: root of unity of order prime to 2 (resp. 3); Δ = 1 if N has at least two distinct prime factors |
| First-order comparison (§4) | three specific elements, norm u², λ = 1+i | **Master Lemma (sum form):** any finite family x_j = x₀ + t_j y_j of norm-N elements with 0 ≠ t_j ∈ λ𝒪 and Σ y_j ≡ c·J (mod λ), ℓ ∤ c, is impossible; over ℤ[i] (λ = 1+i) and ℤ[ω] (λ = 1−ω) |
| Applications | no circulant Hadamard matrix of order > 4 | **Q**, **T3**, **T1**, **F** (statements in `note.pdf`, Theorem 3) |

Lean statements (from `scripts/lean/Checks.lean`):

```lean
GeneralNorm.QuaternaryQ.not_exists_perfect_quaternary (m : ℕ) [NeZero m] : 1 < m → Odd m →
  ¬∃ a, (∀ j : ZMod (2*m), a j = 1 ∨ a j = -1 ∨ a j = Complex.I ∨ a j = -Complex.I) ∧
        ∀ τ : ZMod (2*m), τ ≠ 0 → ∑ j, a j * conj (a (j + τ)) = 0
EisensteinNorm.TernaryT3.not_exists_perfect_ternary (m : ℕ) [NeZero m] : 1 < m → ¬ 3 ∣ m →
  ¬∃ a, (∀ j : ZMod (3*m), a j ^ 3 = 1) ∧ ∀ τ ≠ 0, ∑ j, a j * conj (a (j + τ)) = 0
GeneralNorm.MenonT1.not_exists_perfect_binary_array (N : ℕ) [NeZero N] : 1 < N → Odd N →
  ¬∃ f, (∀ g : ZMod 2 × ZMod 2 × ZMod N, f g = 1 ∨ f g = -1) ∧ ∀ τ ≠ 0, ∑ g, f g * f (g + τ) = 0
RankOne.theoremF_gaussian (p q prime, odd, p ≠ q) (x : ℤ[i][C_{pq}]) (hx : x * x^* = pq) (g) :
  ε(x) * x(g) = (pushforward p x)(g mod p) * (pushforward q x)(g mod q)
```

(T1 for non-square N is trivial, because the sum of f squares to 4N. The content is the case N = u².)

## Attribution

The original construction and proof are OpenAI's work: family 179 of
[github.com/openai/math](https://github.com/openai/math), commit
[`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb),
preprint `preprints/The-circulant-Hadamard-conjecture-September-23-2026/` and Lean
`lean/OAI/LinearAlgebra/CirculantHadamard/`. The remark that Prop. 3.4 holds for every x with xx* = u²
is in the preprint itself (`build/sections/introduction.tex` l. 191, `build/sections/contradiction.tex`
l. 3). This remark only (i) restates the proposition and the comparison for norm m and for ℤ[ω], and
(ii) supplies the normalizations x = (H₀ ± H₁)/(1+i) (quaternary), h(r)/(1−ω) (three-phase), and
h(±1, ±1)/2 (Klein four-group).

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are
not human peer review.

## What changes relative to the original proof

| Step | 179 preprint | Here |
|---|---|---|
| Norm | u² (after reduction to order 4u²) | any N (odd for ℤ[i]; 3 ∤ N for ℤ[ω]); the proof of Prop. 3.4 is unchanged |
| Local part / normalization | C₄, binary | C₂ with ÷(1+i) (Q); C₃ with ÷(1−ω) (T3); C₂×C₂ with ÷2 (T1) |
| Comparison | two differences | one (Q), two (T3), three (T1); abstracted as the sum-form Master Lemma |
| Coefficient ring | ℤ[i] | ℤ[i] and ℤ[ω] (Lean: `Eisenstein/`, a copy of the ℤ[i] development with λ = 1−ω, ℓ = 3) |
| Theorem F | — | Δ = 1 for two primes + relabelling by units + Fourier inversion on C_p × C_q |

## Reproduce

```bash
cd scripts && ./run_all.sh          # cc + python3/numpy; single thread; about 1 minute
```
This prints `OK perfect_search` and `OK klein_mitm`. These checks are independent of the Lean
development and of any earlier artifact (see `VERIFICATION.md`). To also re-check the Lean part, copy
`scripts/lean/OAI/` into the `lean/` directory of openai/math at commit `fd4aeeb2…` and run
`LEAN_PROJECT=/path/to/openai-math/lean ./run_all.sh`. This needs Lean v4.34.1 and the upstream Mathlib
cache. Expected: `9/9 printed theorems use only the three standard axioms`. See `scripts/lean/README.md`.

## Scope and limitations (what is NOT claimed)

- **No new method.** Prop′/Prop″ and the Master Lemma are restatements of OpenAI's Prop. 3.4 and §4
  argument for other norms and for ℤ[ω]. Q, T3, T1 and F follow from them with short normalizations.
- **No claim that Q or T3 newly excludes any particular length.** Many lengths in these families are
  already excluded by classical self-conjugacy arguments (e.g. Turyn's 2p^a for quaternary), by field
  descent, or by other work (e.g. Arasu–de Launey–Ma 2002; Ma–Ng 2009 for three-phase). This remark
  makes no claim about individual lengths. For T1, Turyn's exponent bound already excludes
  C₂×C₂×C_{u²} whenever some prime p ∣ u (p^a ∥ u) is self-conjugate modulo 2(u/p^a)². This covers
  every prime power u and every odd u < 39. The first u it misses are 39, 55, 117, 155, 165, 195, …
- **Theorem F** is proved only for two distinct odd primes and coefficients in ℤ[i]. It is a structural
  identity, not a nonexistence result. It is vacuous unless pq is a sum of two squares (ε(x) has norm pq),
  i.e. unless p, q ≡ 1 (mod 4).
- **Formalization scope.** The Lean statements are about ℂ-valued sequences (Q, T3) and ℤ-valued
  arrays (T1) with perfect periodic autocorrelation. The equivalence with Menon–Hadamard difference sets and with circulant Butson
  Hadamard matrices is standard and not formalized. Trust rests on the Lean kernel and on the
  statements printed in `Checks.lean`.
- **Computations are unit tests.** The exhaustive searches (lengths ≤ 15, groups of order ≤ 36) are
  sanity checks, not evidence for the infinite families.

## Duplicate search (to 2026-10-09)

We searched GitHub (openai/math has issues disabled; its only PR is #1, the 10-07 manuscript update), google-deepmind/formal-conjectures issue #6896 and PRs #6897/#6898, and community follow-up
repositories that track openai/math (about 30 cloned on 2026-10-09). We also ran web/arXiv searches
for perfect quaternary/three-phase sequences and Menon difference sets in C₂×C₂×C_N. To our knowledge
no prior or concurrent public text states Q, T3, T1 or F as consequences of the 179 argument. The
formal-conjectures items #6896–#6898 only link OpenAI's real (binary) circulant Hadamard theorem to an
existing statement and do not touch other alphabets. The recent BH(n,q) nonexistence criterion of
Czifra–Matolcsi–Szöllősi (J. Algebr. Comb. 2026, arXiv:2511.11230) uses positive-definite functions,
and we did not find an overlap with these families.

## Files

- `note.tex`, `note.pdf`: 3-page note (statements, why it is a parallel generalization, Lean, checks).
- `scripts/lean/`: the 17 Lean files (4085 lines), `Checks.lean`, `lean-toolchain`, `README.md`,
  `SHA256SUMS`, `expected_output/axioms.txt`.
- `scripts/independent/`: `perfect_search.c` (exhaustive search over μ_q-valued functions) and
  `klein_mitm.py` (meet-in-the-middle for C₂×C₂×P, |P| = 9).
- `scripts/expected_output/`, `scripts/run_all.sh`.
- `VERIFICATION.md`: who checked what.
