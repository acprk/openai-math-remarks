# Verification record: R10 (explicit fermionic Laughlin gap, OpenAI family 269)

Date: 2026-10-09. Host: Linux x86_64. Python 3.12 with python-flint 0.9.0 (Arb), numpy 2.5.3, scipy 1.18.1, mpmath 1.4.1. At most 4 single-threaded processes. No Lean was run.

Pins: openai/math `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. The preprint directory, `LaughlinFock.lean` and `docs/269.md` are unchanged at HEAD `fd4aeeb2…`, and the local checkout has no modifications.

## Who did what

| Role | Session | What |
|---|---|---|
| Primary computation | R7-skill 269 session (`reviews/R7-skill/269/cert/`, here `scripts/primary/`) | SDP for Γ, rationalisation, arb certificates for Q = 24..299 (r-space + rank-free b-space), 168 Taylor-model boxes for Q ≥ 300, rigorous ED for N ≤ 8, Jordan–Wigner sector cross-check (22 sectors). |
| Adversarial review + independent re-implementation | this session (separate AI session) | Instructed to refute Theorems A and B. Rebuilt the operators from Clebsch–Gordan coefficients with new code (`scripts/independent/`). The new code shares no code with `primary/` except: check1e and check3 read the certified Γ, γ_Q and ρ_D, which are the objects under test, and check3 calls the primary pipeline in order to compare against it. Results below. |

## Items checked

### 1. Definitions (GAP G2)

`independent/indep_fock.py` builds e_i as the Condon–Shortley basis and the pair-spin-(Q−1) states from the Racah formula (exact rational squares). It fixes v_p = σ|Q−1, Q−1−p⟩ with v_0 = −e_0∧e_1 (03-spin.tex l.10). Jordan–Wigner uses B(e_x∧e_y) = c_y c_x (02-operators l.211; same as `LaughlinFock.pairAnnihilator`).

| check (`check1_definitions.py`, Q = 24 and 27) | preprint | result |
|---|---|---|
| CG pair vectors vs eq. (pair) for all p, i<j | 03-spin l.12-16 | max diff 2.2e-16 |
| Σ_{i<j} P_ij^{(1)} on ∧³ ⊂ U^{⊗3} (first quantisation) vs Σ_p B_p†B_p | 01-intro l.28-35; 02-operators l.221-228 | max diff 4.4e-16 (coefficient one confirmed) |
| Normal-ordered square, N = 3 (matrix) and N = 4 (5 random vectors) | Lemma 2.1, 02-operators l.241-246 | 8.9e-16; relative 3.6e-16 |
| q_z(Q) = spectrum of W_3†W_3 − I on every spin block | Lemma 3.2, 03-spin l.48-53 | 4.4e-16 |

**Full operator inequality on sectors** (`check1e_sector.py`). This uses the certified rational Γ and γ_Q. K_Q = d_a·Haar average of K, computed by Schur's lemma on highest-weight spaces (04-comparison l.111-120). We use our own descendant frames, re-projected onto spin J at each level with the Casimir. A first version without re-projection showed spurious negatives of −2e-6 at N = 5: round-off was amplified by J_−^μ on higher-spin contamination. With re-projection, λ_min(H² − γ_Q H − K_Q) over all spins J:

| (Q, N) | γ_Q (file) | min λ | best γ' in sector |
|---|---|---|---|
| (24, 3/4/5) | 0.211779 (24_36_48) | −2.8e-16 / −6.8e-16 / −1.3e-15 | 0.2136 |
| (27, 3/4) | 0.145717 | −5.9e-16 / −1.2e-15 | 0.2136 |
| (29, 3/4) (worst Q of range) | 0.136071 | −3.2e-17 / −6.0e-16 | 0.2136 |
| (49, 3), (66, 3) | 0.175908, 0.166963 (48_60_72) | −1.6e-16, −2.5e-16 | 0.178 |
| (73, 3), (100, 3) | 0.085458, 0.089643 (72_100_150_300) | −4.9e-16, −5.4e-16 | 0.096 |
| negative control: v_p → (−1)^p v_p at (24, 3) | | **−778** (VIOLATION, as intended) | |

The test is therefore sensitive to the phase and normalisation conventions, and the certificate's conventions pass. Not covered: Γ_300 (Q ≥ 300) is too large for sector tests. It is tested by item 3 instead.

Correspondence between certificate and preprint formulas:

| certificate step | preprint |
|---|---|
| H² = H + L3(W3(W3†W3 − I)W3†) + L4(W4W4†) | Lemma 2.1, 02-operators l.241-246 |
| A_2 = ηH | 04-comparison l.170-174 |
| A_3 Schur form | (A3), l.177-180 |
| A_4 Schur form | (A4), l.228-240 |
| Schur averaging rule | l.111-120 |
| negative 3-body slack and tail | Lemma 6.1, 06-gap l.250-256 |
| four-body cost n_D ρ_D | Prop. 5.3 Bessel step, 05-transfer l.199-229 |
| compression equivalence | Lemma 5.1, l.24-43 |

We re-derived that the condition actually certified, Ŵ_b(SᵀCS + ρI)Ŵ_bᵀ ⪰ 0 on the weight space, is the one whose cost is n_D ρ_D H. We also showed that the r-space form ρG + GCG ⪰ 0 with G = SΠSᵀ is equivalent to it (minimum-norm argument). The `ρ G_J` term in `cert_core.py` is therefore correct, and `(G²)_J` is not needed.

### 2. Lemma R (GAP G1)

`check2_rank.py` computes rank W_D = rank U − rank(J_+U) directly. Here U holds the vectors v_p∧e_j∧e_k at index sum D+1. This does not use the kernel argument.

- (i) Exact rational arithmetic in the basis f_i = √C(Q,i) e_i, where U is an integer matrix and J_+ f_i = (Q−i+1) f_{i−1}. Run at Q = 24, 300, 1000, 10⁴, 10⁸.
- (ii) float64 SVD in the orthonormal basis at Q = 300, 1000, 10⁴. Null singular values are ≤ 6.6e-15 and the next singular value is ≥ 2.0.

All 23 values of D agree with hw(D) − zm(D) = 0,0,0,0,1,0,1,1,2,1,3,2,4,3,5,4,6,5,7,6,8,7,9. Mismatches: 0. We also read the written proof (the kernel of H_4 via Lemma 6.2's divisibility, and ran Ŵ = (ker H_4)^⊥) and found it correct. Only the upper bound is used in boxes: a smaller rank would make G_J singular and the check would fail.

### 3. Q ≥ 300 boxes

`check3_boxes.py` (own cover check; boxes #0, #1, #60, #120, #167):

- **Cover.** 168 boxes form an exact contiguous rational cover of [0, 1/300]. There are no FAIL records. Box #0 = [0, 1/50400] covers every Q ≥ 50400, and x = 0 is only a closed endpoint of an analytic function. No limit is taken.
- **Taylor coefficients.** For the tail, the three-body costs c_z (z = 2, 5, 7, 9) and entries of G_J and P for D = 13, 19, 21, 23, we re-derived the coefficients by degree-12 interpolation at 13 Chebyshev nodes of [x0 − 1.2h, x0 + 1.2h], using 256-bit point values and a 90-digit solve. They agree with the arb-series coefficients c_0..c_5 to a relative 2.6e-48 (scaled by h^i). In every case the order-6 enclosure |F_6| at the ball base contains both the interpolated c_6 and the point-wise c_6 at all nodes.
- **Point re-certification.** At the box endpoints and centres, and at integer fluxes in each box, we re-certified with the box's own ρ_D: r-space evaluate() and the rank-free b-space form. Integer fluxes tested: 50400, 50401, 10⁶, 10⁹, 10¹⁵, 25200, 25201, 827, 828, 840, 417, 418, 420, 300, 301. All pass, and γ(point) ≥ γ_lo(box) everywhere. For example box #0 gives γ = 0.059470 at Q = 10⁹ and 10¹⁵, against γ_lo = 0.059470.
- **Code reading.** The remainder logic (Lagrange form, affine-in-t convexity at t = ±h), the tail bound (ratio test for 6z(z+1)(2/3)^z, valid for Q ≥ 300 > 250) and the K-reduction (C rows vanish for r ≥ 17; w_{D,1} = 0 for even D) are correct.

### 4. Fock inequality ⇒ Laughlin gap

H² ≥ γH on ∧^N implies spec ⊂ {0} ∪ [γ, ∞) by the spectral theorem. Lemma 6.2 (06-gap.tex l.300-321) holds for **every N ≥ 2**. It needs Q = 3(N−1) ≥ 3, the factorisation of the pair-spin Q−r summands by b_ij^r, Gauss's lemma, and the degree count. There is no large-N or large-Q assumption. The only "for N ≥ N_0" is in the proof of Cor. 1.2 (l.323-332), where it comes from the threshold of Thm 1.1. Kernel dimension 1 was also confirmed numerically for N ≤ 8 (item 5).

### 5. ED for N ≤ 8 (GAP G3)

- **Primary.** `ed_cert.py` certifies X = H − gI + ψψᵀ ≻ 0 by float Cholesky with Higham Thm 10.3 (‖ΔA‖₂ ≤ γ_{n+1} tr(R̂ᵀR̂), ×10 safety) plus the entrywise rounding of fl(X). The bound is 3.3e-7 against τ = 2.7e-4 at N = 8. The Courant–Fischer step and the use of the L_z = 0 sector are valid: QN/2 ∈ ℤ, so every multiplet meets L_z = 0.
- **Independent.** `check5_ed.py` builds H from CG pair vectors with its own Slater rules. It runs (i) Lanczos (eigsh, k = 4) with residual bounds and (ii) LDLᵀ (Bunch–Kaufman) inertia of H − σI at σ = g_N, ±1e-6.
  - N = 8 (dim 8512): λ = 0, 0.5315767947, 0.5423831886, 0.5538054369, residuals ≤ 1.8e-14. #{λ < 0.5305} = 1, #{λ < 1e-6} = 1, #{λ < −1e-6} = 0.
  - N = 2..7 are identical in form, all OK, and the gaps agree with `ed_cert` to 10 digits.

### 6. Dedup / stripping test

- Shencongli/bosonic-laughlin-spectral-gap (clone of `fda5da1`, paper.tex l.20-60, 680-760, 860-905, 1334-1380, 1548-1600) treats bosonic V0 only: Q ≥ 3088, N ≥ 1545.
- It contains the bosonic analogue of Lemma R (eq. (stable-rank)) and interval certificates on 1/Q ∈ [0, 1/q_−]. Its boson-specific ingredient, the coherent occupation bound n_β ≤ H + I, is not needed for fermions, because ‖c_j‖ ≤ 1.
- It shares the mid-range problem: Q < 3088 and N < 1545 are left open. Its margin at the fixed target 1/3 vanishes at Q ≈ 3088.
- Web search (2026-10-09) found no fermionic explicit or all-N result.
- **Stripping test:** our new mathematics is Lemma R plus the finite-Q bookkeeping, both short consequences of the preprint. Everything else is certified computation. Classified as an effectivisation (remark-level). No Lean, per project rules.

## Residual GAPs

| ID | Description | Status |
|---|---|---|
| G1 | Lemma R | Written proof checked. Exact rank at 5 fluxes. Bypassed at integer Q ≤ 299 by the b-space check. Q ≥ 300 relies on the written proof for fluxes not sampled. |
| G2 | Transcription | Independent CG rebuild plus sector tests (13 sectors, 3 of the 4 Γ files). The coupling coefficients E_D^Q, S_r^Q are not separately re-implemented. They are tested only through the sector inequality. |
| G3 | Higham bound | Applies to LAPACK's blocked Cholesky. Independently corroborated by inertia and Lanczos. |
| G4 | Ball-endpoint → float conversions | Covered by a 1e-14 margin on γ. |
| G5 | Q < 24 | No Fock statement for Q < 24. Not needed for Theorem B. |
| G6 | Lean | Not formalised. |
| G7 | Γ's not re-optimised for boxes | Constant not optimal. |

## Verdict

**CONFIRMED.** We found no counterexample, definitional mismatch or certificate flaw. One false alarm, a numerical instability in our own first sector test, was diagnosed and fixed. Theorem A (Q ≥ 24, γ = 0.05947) and Theorem B (N ≥ 2, gap ≥ 0.05947) stand, subject to the residual GAPs above.
