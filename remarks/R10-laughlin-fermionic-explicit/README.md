# R10. Fermionic Laughlin gap on the sphere: explicit for every N ≥ 2 (constant 0.05947 > 1/25)

**What this remark shows.** OpenAI's family-269 preprint proves the Fock-space inequality H_Q² ≥ γ H_Q for every γ < γ_* = 0.046167, but only for Q ≥ Q_γ with a non-explicit Q_γ. Its Laughlin gap H_{N,3(N−1)} ≥ (1/25)(I − P_L) therefore holds only for N ≥ N_0 with a non-explicit N_0. We keep the preprint's architecture. Its normal-ordered square, rotation averaging, three-body spectrum, tail lemma, Bessel step and kernel lemma are all reused. We replace its only non-effective step, the transfer Lemma 5.2, with a direct finite-flux positivity check, and we certify every flux rigorously. The result:

* **Theorem A.** For every integer Q ≥ 24, H_Q² ≥ 0.05947 · H_Q on the whole Fock space F_Q, for all particle numbers at once.
* **Theorem B.** For every N ≥ 2, H_{N,3(N−1)} ≥ 0.05947 · (I − P_{L,N}) on ∧^N U_{3(N−1)}.

## Statement (Before / After)

Notation follows the preprint. U_Q is spin Q/2. H_Q = Σ_{p=0}^{2Q−2} B_p†B_p, with coefficient one per unordered pair projector onto pair spin Q−1. This is the operator `LaughlinFock.hamiltonian`.

| | Fock inequality H_Q² ≥ γ H_Q | Laughlin gap at Q = 3(N−1) |
|---|---|---|
| **Before** (preprint Thm 1.1, Cor 1.2; Comparator `thm_fock`) | every γ < γ_* = 0.046167, for Q ≥ Q_γ. Q_γ is existential. The preprint says (01-introduction.tex l.75-77) that it "does not give an explicit flux threshold" | (1/25)(I − P_L) for N ≥ N_0, with N_0 existential |
| **After** (this remark) | γ = **0.05947** for **every integer Q ≥ 24**. Hence `thm_fock` holds with **Q_γ = 24 for every γ < γ_***, and in fact for every γ ≤ 0.05947 | **0.05947 (I − P_L) for every N ≥ 2**, i.e. N_0 = 2 |

Certified constants by range:

| flux range | Γ data file | certified γ_Q ≥ (worst Q) |
|---|---|---|
| 24 ≤ Q ≤ 48 (every integer) | `gamma_24_36_48.json` | 0.136071 (Q = 29) |
| 49 ≤ Q ≤ 72 | `gamma_48_60_72.json` | 0.166963 (Q = 66) |
| 73 ≤ Q ≤ 299 | `gamma_72_100_150_300.json` | 0.079313 (Q = 209) |
| Q ≥ 300 (168 boxes in x = 1/Q ∈ [0, 1/300]) | `gamma_300_1000_100000000.json` | 0.059470 (box x ∈ [0, 1/50400]) |
| N = 2..8 (Q = 3..21), Laughlin sector only | rigorous ED | g_N ≥ 0.9990, 0.7847, 0.6796, 0.6083, 0.5847, 0.5477, 0.5305 |

Where the constant goes at Q → ∞, compared with the preprint:

| | two-body price η | three-body tail δ_∞ | four-body allowance | γ |
|---|---|---|---|---|
| preprint | 0.935274 | 61/4096 = 0.014893 | 1222ε = 0.003666 | 0.046167 |
| here | 0.924439 | 0.014893 | Σ n_D ρ_D = 0.001159 (+ 3-body slack 4·10⁻⁵) | 0.059470 |

The limiting gain (×1.29) comes from re-solving the preprint's SDP with a larger Γ rank, and from paying the exact four-body deficit instead of a uniform ε. It is a re-tuning. We do not count it as a new idea. The substantive change is that the statement now holds at every finite flux, with N_0 = 2.

## Attribution

The model, the Fock-space comparison method, the normal-ordered square (Lemma 2.1), the spin formulas (Lemmas 3.1–3.3), the comparison operators (eq. (F)), the Bessel step (Prop. 5.3), the tail lemma (Lemma 6.1) and the kernel lemma (Lemma 6.2) are **OpenAI's work**. They come from family 269 of [github.com/openai/math](https://github.com/openai/math), commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a` (unchanged at `08fcd618213aa43b0fcaff4526a8a48d66d02e55`):

- [`preprints/A-Fock-space-inequality-and-the-Laughlin-spectral-gap-September-24-2026`](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-Fock-space-inequality-and-the-Laughlin-spectral-gap-September-24-2026) (`build/sections/0[1-8]-*.tex`)
- `lean/ComparatorChallenges/LaughlinFock.lean` (`thm_fock`, `sorry` in the challenge file) and `lean/docs/269.md`

This remark adds four things:

1. A finite-flux version of the comparison (note, Prop. 2). With an arbitrary rational PSD Γ = ⊕_t F_tF_tᵀ on the preprint's support, H_Q² − γH_Q ⪰ 𝒦_Q ⪰ 0 holds with γ = 1 − η − δ_Q − c₃(Q) − Σ_D n_D ρ_D. Every term is exact at the given Q.
2. **Lemma R:** rank W_D = hw(D) − zm(D) for every Q ≥ 24. This follows from the preprint's Lemma 6.2. It reduces the four-body condition to m × m matrices with m ≤ 9.
3. Four rational Γ's, ball-arithmetic certificates at every integer 24 ≤ Q ≤ 299, and order-6 Taylor-model certificates on 168 boxes covering x = 1/Q ∈ [0, 1/300].
4. Rigorous exact diagonalisation for N ≤ 8.

## What changed in the preprint's proof (tex line numbers)

| preprint | role | replaced by |
|---|---|---|
| 04-comparison.tex l.190-198 ((3certificate), (A3bound), "for all sufficiently large Q") | three-body limit comparison | exact finite-Q slacks s_z = q_z(Q) − (d_a/d_{3,z}) e_z^Q. A negative slack is paid as \|s_z\|(z+1)(d_{3,z}/d_a)H, by the argument of Lemma 6.1 (06-gap.tex l.250-256), which is valid at every Q |
| 04-comparison.tex l.262-278 ((4certificate) with ε = 3·10⁻⁶) and 05-transfer.tex l.46-166 (Lemma 5.2, ρ_Q → 0, non-explicit) | four-body limit certificate + transfer | positivity at the given Q of Ŵ_b(SᵀC_D^QS + ρ_D I)Ŵ_bᵀ on ∧⁴ (b-space), equivalently ρ_D G + G C_D^Q G ⪰ 0 (r-space, uses Lemma R). **No transfer and no limit.** |
| 05-transfer.tex l.173-232 (Prop. 5.3, 1222ε + o(1)) | four-body cost | Σ_D n_D ρ_D H, with the same Bessel step (l.199-229) |
| 06-gap.tex l.237-274 (δ_Q → 61/4096) | tail | exact finite sum (integer Q), or z ≤ 250 exact plus Σ_{z>250} 6z(z+1)(2/3)^z (boxes) |
| 06-gap.tex l.276-293, l.323-332 (proofs of Thm 1.1, Cor 1.2) | "for large Q" / N_0 | rigorous evaluation at every Q ≥ 24. N ≤ 8 by ED, which uses Lemma 6.2 (l.300-321) for every N ≥ 2 |
| `verification/rows.txt` (7 rows) | comparison data | 4 rational Γ files (`scripts/primary/data/`) |

Lemma 6.2 (06-gap.tex l.300-321) holds for every N ≥ 2. Its proof uses only Q = 3(N−1) ≥ 3, unique factorisation and degree counting, with no large-N assumption. The step Theorem A ⇒ Theorem B is: H² ≥ γH on ∧^N gives spec ⊂ {0} ∪ [γ, ∞), and the kernel is the Laughlin line.

## Relation to the bosonic V0 work of Shencongli

[Shencongli/bosonic-laughlin-spectral-gap](https://github.com/Shencongli/bosonic-laughlin-spectral-gap) (v0.1, 2026-10-07, AI-generated candidate proof curated by Xin Shen) adapts the same OpenAI comparison to the **bosonic V0** model. It proves H_Q² ≥ H_Q/3 for Q ≥ 3088, hence a gap for N ≥ 1545. Fluxes Q < 3088 and particle numbers N < 1545 are not covered. Its "stable rank" formula (paper.tex eq. (stable-rank), §"Transferring the finite certificate") is the bosonic analogue of our Lemma R: the same divisibility-of-the-kernel argument, applied to Δ_r² instead of ∏ b_ij³. Its interval certificates on 1/Q ∈ [0, 1/q_−] play the role of our Q ≥ 300 boxes.

Differences:

- **Model.** Here fermionic V1, ν = 1/3. There bosonic V0, ν = 1/2.
- **Coherent occupation bound.** Their n_β ≤ H + I is specific to bosons. The fermionic case does not need it, because ‖c_j‖ ≤ 1.
- **Mid-range flux.** Their threshold reflects a fixed target γ = 1/3 whose margin vanishes near Q = 3088 (A − 1/3 ≥ 4.2·10⁻⁷ there). They have the same "mid-range Q is not covered" situation. They leave it open, as they do N < 1545. Here the mid range 24 ≤ Q ≤ 299 is closed by re-solving the SDP on three flux windows and certifying each integer Q. The remaining N ≤ 8 are closed by rigorous ED.

We did not use their code or data.

**Gate assessment (project rule: stripping test).** Once the preprint's lemmas and the bosonic precedent are removed, the new mathematics here is Lemma R and the finite-flux bookkeeping of Prop. 2. Both are short consequences of the preprint's own lemmas. The rest is computation. This is an *effectivisation* of OpenAI's argument, which is the "explicit parameters" kind of remark in this repository. It is not a new method, and no Lean formalisation is attempted.

## Reproduce

```bash
cd scripts
bash run_all.sh               # verify mode (seconds): SHA-256 of code/data (SHA256SUMS) and shipped outputs (SHA256SUMS.outputs), re-aggregate, diff
bash run_all.sh full          # ~65 min wall, 4 single-threaded processes: recompute the primary certificate, then verify
bash run_all.sh independent   # ~45 min wall, <=4 processes: adversarial re-derivations (scripts/independent/)
```

Requirements: Python 3.12, python-flint 0.9.0 (Arb ball arithmetic), numpy 2.5, scipy 1.18, mpmath 1.4 (`PY=/path/to/python`). Expected outputs are in `scripts/expected_output/`. The rational Γ files are listed in `scripts/SHA256SUMS`. The float SDP solutions they were rounded from are `scripts/primary/data/sdp/*.npz` (regenerate with `rationalize.py`; untrusted input).

## Scope and limitations (what is NOT claimed)

- **Computer-assisted, not formalised.** There is no Lean proof. `thm_fock` remains `sorry` in OpenAI's challenge file. We only state that our Theorem A implies its statement with Q_γ = 24. We did not check the normalisation of `LaughlinGap.lean`'s `energy` against H_{N,Q}, so we make no claim about that Comparator file.
- **Residual trust (GAP list, see VERIFICATION.md):**
  - G1: Lemma R is a written argument. It was checked by exact rank computation at Q = 24, 300, 10³, 10⁴, 10⁸, and it is bypassed at integer Q ≤ 299 by the rank-free b-space check.
  - G2: transcription of the preprint's coupling formulas into code. Checked by an independent Clebsch–Gordan re-implementation and by direct sector tests of H² − γH − 𝒦_Q ⪰ 0.
  - G3: ED relies on Higham's Cholesky error bound. Re-checked by Lanczos with residuals and by LDLᵀ inertia.
  - G4: float conversions of ball endpoints, covered by a 10⁻¹⁴ margin.
  - G5: no Fock-form statement for Q < 24.
  - G7: the Γ's are not re-optimised for the box range.
- **Not optimal.** The true gaps are much larger. Sector tests give γ' ≈ 0.21 at Q = 24–29 for N ≤ 5, and the ED gaps are ≥ 0.53 for N ≤ 8. 0.05947 is what this certificate yields. The asymptotic constant is limited by η and by the tail 61/4096.
- **Fermionic V1 on the sphere only.** The planar statement of the preprint (Section 7) is not touched.
- **Premise.** We rely on the preprint's Lemmas 2.1, 3.1–3.3, 6.1 (bound part), 6.2 and the Schur/Bessel steps. Lemmas 2.1, 3.1, 3.2 and the coefficient-one normalisation were re-verified numerically from Clebsch–Gordan coefficients at Q = 24, 27. Lemma 6.2 was re-read line by line, and its kernel dimension was checked by ED for N ≤ 8.
- **Novelty limited to our search.** To our knowledge (search of GitHub and the web up to 2026-10-09), no explicit-threshold or all-N statement for the fermionic V1 spherical Laughlin gap has been published. The closest is the bosonic V0 work above.

## Authorship

[AUTHOR NAME].

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
