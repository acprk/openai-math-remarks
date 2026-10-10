# R16 verification record

Claim: with s ≤ √(n/2) and simultaneous expansion of all factors of degree ≥ t₁ = n/(2s) (M1–M4 in the
README), the proof of OpenAI's Theorem 1.1 (family 135) gives S ≥ n^{(1/(4√2)−ε)√n}, with explicit
thresholds c = 1/8 for log n ≥ 57.941 and c = 1/6 for log n ≥ 296.735.

All runs: 2026-10-09, CPython 3.10.12, mpmath 1.3.0, NumPy 2.2.6, SciPy 1.15.3, one thread per script.
`scripts/run_all.sh` reproduces all eight expected outputs byte-for-byte (`OK` ×8).

## Layer 1: main computation (`scripts/main/`)

| Script | Check | Result |
|---|---|---|
| `explicit_bound.py` (2 s) | Positive control: the paper's chain | low exponent 1/20, high 3/160, Prop. 3.1 uses 1/100, Theorem 1/200 → 1/400; the chain as written gives any c < 3/320 |
| | Continuum optima | single expansion: σ = 1/2 → 1/8; simultaneous expansion: σ = 1/√2 → 1/(4√2) = 0.1767766953 |
| | Explicit lower bound B(n) = min(B₁, B₂) for log S/(√n log n). Every term is an explicit majorant/minorant, worst case over s ∈ (√(n/2)−2, √(n/2)] | monotone on the grid log n ∈ [14, 20009], step 5 |
| | Thresholds (bisection) | 1/400: 17.59; 3/320: 18.24; 1/8: 57.94; 3/20: 112.04; 1/6: 296.74; 0.17: 442.69; limit 0.1767767 |
| | Decomposition of K/√n at log n = 50 | C_I 2.000, C_e 1.000, C_o 0.0000 (units √n) |
| `hbound_check.py` (36 s) | Exact H_e = Σ q^{\|i−λe\|} against the M2 bounds for all e < t₁, real parameters | n = 10⁶, 4·10⁶, 1202604 (≈ e^14): max H_e/bound = 1.000000, H₁ ≤ 2q^λ, all OK |
| `framework_limit.py` (101 s) | (a) τ*(σ) on a grid of 300001 points | argmax σ = 0.70711, max 0.176776 |
| | (b) exact finite-n framework value (exact H_e, q = n^{−1/2}) | σ = 1/√2: 0.0640 / 0.0878 / 0.1038 / 0.1146 at n = 10⁴…10⁷ |
| | (c) constrained sum Σ_{Σi_j=k} vs. Π H_e (all degree 2) | difference/√n = 0.082 / 0.068 / 0.064 at n = 10⁴, 10⁵, 10⁶ (cost e^{O(√n)}) |
| | (d) Lemma 2.1 side conditions at n = 1202604 | 2n ≤ min(a,b)/2, a ≥ k, q = 9.16·10⁻⁴ ≤ 0.1, t₁ = 776.9 > 2, (1−1/(a+v))^{−n} ≤ 2: all true |
| `imm_wordsum.py` (55 s) | Exact word sum of Section 5 with α = n^{−γ} | at n = 10⁶: log Z/(√n log n) = 0.0057 (γ = 0.5), 0.0567 (0.7), 0.0960 (0.75 > gain 0.0884) |

## Layer 2: independent re-implementation (`scripts/independent/`)

Written from the TeX source in a separate session; does not import `main/`. Uses the exact integer
parameters of (eq:parameters) (k, m, v, u, a, b as integers; α, β, q in 80-digit arithmetic).

| Script | Check | Result |
|---|---|---|
| `v_params_H.py` (29 s) | M2 nearest-integer claims and H_e bounds for every e < t₁; inequality (*) of M3, −Σ_{j∉X} log H_{e_j} − E_X log H₁ ≥ L_q s/2 + E_X c′ − C_e − C_o, for 320+ adversarial degree lists per n (uniform degrees around t₁, 2t₁, 4t₁, random mixtures); per-S value ≥ min(A, B) | n = 10⁴, 100001, 10⁶, 2000001, 10⁷: all OK; min slack in (*) = 33.6 / 57.6 / 94.4 / 114.2 / 164.0; C_e/√n = 1.078 → 1.004 (exact all-degree-2 loss 1.055 → 1.003) |
| `v_imm.py` (1 s) | Exact word sum (transfer matrix) vs. the paper's bound, with the modified s | n = 2000 … 60001: exact 8.7 … 26.8 vs. bound 87.3 … 396.1 (exact/√n = 0.195 … 0.109); schedule ℓ_i ∈ {1,2} |
| | Explicit C_I − 2√n − (5/2 + 1/(2√2)) log n | 1.294 (log n = 20), 1.292 for log n = 40 … 296.7 |
| `v_threshold.py` (9 s) | Thresholds from exact parameters, exact C_I(n), also with s replaced by s − 2 | c(n) ≥ 0.124999 at log n = 57.94, 0.166665 at 296.7, 0.176177 at 5000; monotone on [20, 1520] step 0.5 (0 violations); thresholds 1/8: **57.941**, 1/6: **296.735** (1/400 and 3/320 already hold at the lower search limit log n = 20) |
| `v_framework.py` (15 s) | Framework limit, single-degree adversary and two-degree mixtures | max_σ = 0.176750 at σ = 0.7070 (grid); σ = 1/√2: 0.176777; two-degree mixtures at σ = 1/√2: 0.176777 |

The two layers agree on all thresholds (57.94 vs 57.941; 296.74 vs 296.735) and on the limit 1/(4√2).

## Layer 3: adversarial review

A separate AI session tried to refute the claim against the TeX source and the scripts. Verdict:
confirmed. This is not human peer review.

## What is numeric and what is exact

- The asymptotic statement (any c < 1/(4√2)) is analytic and uses no computation.
- The explicit thresholds rest on explicit majorants evaluated in 60-digit (main) and 80-digit
  (independent) arithmetic, with monotonicity checked on a grid. They are not interval-certified.
- The framework limit and the IMM-side obstacle (γ > 1/2) are first-order / numerical evidence, not theorems.

## Open gaps

- **GAP-1:** thresholds are high-precision floating point with a grid monotonicity check, not interval
  arithmetic.
- **GAP-2:** no Lean. The upstream Lean (`lean/OAI/Algebra/DepthFive/`) covers the original 1/400 statement.
- **GAP-3:** relies on the correctness of the IMM-side Sections 4–5 as written.
- **GAP-4:** the adversarial review checked modifications M1–M4 against the upstream source; the written proofs in `note.tex` were then composed from those checked inequalities and the scripts, and have not had a separate line-by-line review.
