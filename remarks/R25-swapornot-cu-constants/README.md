# R25 — Family 238: exact values of the reciprocal-degree constants C_u (C₁ = 473/144), and ε_n = Θ(n^{−e*}), e* = 3.03456…

## Premise and confidence (read first)

- **What kind of result this is: bookkeeping of constants. It is not a new theorem and not a new method.**
  Family 238 (Thorp shuffle) uses the constant C_u = sup_m Σ_{λ⊢m} f_λ^{−u} and proves only that it is
  finite. This remark computes the value of C_u for five values of u. It also points out that one O-bound in
  the main paper is loose. The arguments are routine: enumeration, a peeling recursion and an analytic tail.
- **Origin, and why it is trivial.** We computed these constants for an internal attempt at
  full-domain CCA security of swap-or-not (SN) in O(log N + λ) fixed rounds, via the family's transfer theorem.
  An internal audit found that this statement is **trivial**. The asymptotic statement follows in a few lines from the
  family's own corollary ("7n/8-label error o(1) ⇒ μ∗μ → U") combined with Hoang–Morris–Rogaway's Theorem 3
  and submultiplicativity of total variation. The explicit constants improve on that trivial route only by
  a factor of 3–5. **This remark therefore makes no claim about swap-or-not**, round counts, or the HMR
  conjecture. Only the constants themselves are kept.
- **Confidence.** Two separate routes certify every value: different code, different peeling inequality,
  different tail lemma, different arithmetic in the exact region. The overall architecture (exact region,
  recursion, tail) is the same elementary one, and route B's recursive region uses IEEE doubles with an
  explicit margin rather than formal intervals (see `VERIFICATION.md` §4). A separate adversarial review also checked them (see
  `VERIFICATION.md`). Both routes are computer-assisted proofs and are not formalised.

**One line.** Before → After:
- C_u: "C_u < ∞" (no value anywhere in family 238) → C₁ = 473/144 (attained at m = 6), plus four more exact values in the table below.
- Main paper constant A = √(32C₁): no value → √946/3 ≈ 10.2524.
- ε_n: "O(n^{−5/2})" → Θ(n^{−e*}), e* = 99·log₂(32/31) − 3/2 = 3.03456…

## Statement (Before / After)

| | Before (family 238) | After (this remark) |
|---|---|---|
| C₁ = sup_m Σ 1/f_λ | finite (main paper, Lemma `lem:reciprocal-degrees`) | **473/144 = 3.284722…**, attained at m = 6 |
| C_{1/2} | finite (companion, Thm `l:degree-all-powers`) | **5.9061465497…** (m = 8) |
| C_{3/10} | finite | **13.0732169915…** (m = 14) |
| C_{1/5} | finite | **34.5870288005…** (m = 22) |
| C_{3/20} | finite | **89.3331732895…** (m = 30) |
| C_{1/10} | finite | ≥ 562.0737 (Z_{1/10}(50)); upper bound **not claimed** |
| A = √(32C₁) in Prop. `prop:orbit-lift`, eq. (h-12) | symbolic | √946/3 = 10.2524… |
| ε_n, eq. (h-8) | O(n^{−5/2}) | ½n^{−e*} ≤ ε_n ≤ (√2/2)n^{−e*} for d ≥ 15, e* ∈ (3.03456, 3.03457) |

Each value of C_u is the exact Z_u(m*), enclosed in an interval of width < 10⁻¹⁰. The proof that no other m
does better is computer-assisted. Growth (a one-line remark in the note): f_λ ≤ √(m!) gives
log C_u ≥ (π²/6 + o(1))/(u log(1/u)). We claim no matching upper bound.

Relation to R5 (`../R5-thorp-24d/`, a separate entry that is not modified here): R5 needs C_σ for some
σ < 2/9 and lists it as "finite but not computed". C_{1/5} = 34.587… is such a value. This does **not**
give R5 an explicit d₀, because R5's other non-explicit step (Σ D^{−ε} → 0) remains.

## Attribution

The Thorp-shuffle results and the transfer theorem in which C_u appears are OpenAI's work: family 238 of
[github.com/openai/math](https://github.com/openai/math), commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.
The relevant files are
`preprints/Optimal-order-mixing-of-the-Thorp-shuffle-September-26-2026/` (main paper) and
`preprints/From-partial-permutation-information-to-Fourier-bounds-September-26-2026/` (companion).
The finiteness of C_u is classical (Liebeck–Shalev 2004; Gamburd 2006). This remark only evaluates the
constant and sharpens one stated exponent. It does not depend on the correctness of the preprints: the
statements about C_u are self-contained, and the ε_n statement is about the displayed formula.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction;
independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Where in the sources

| File, line | Text | Here |
|---|---|---|
| companion `sections/introduction.tex` l. 14; `sections/degrees.tex` l. 165 | C_u := sup Σ D_λ^{−u} < ∞ | values for five u |
| companion `sections/core.tex` l. 41, 133 | ‖f̂(ρ)‖²_HS ≤ bBC_u D^{−1+(u+2)/b} | C_u may be substituted |
| main `build/main.tex` l. 420–424 | C₁ := sup_s Σ (f^λ)^{−1} < ∞ | C₁ = 473/144 |
| main `build/main.tex` l. 495 | A = √(32C₁) | A = √946/3 |
| main `build/main.tex` l. 381–389, 401 | ε_n = … = O(n^{−5/2}), using 198 log₂(32/31) > 8 | exponent 99 log₂(32/31) − 3/2 = 3.03457 |

The ε_n sharpening changes nothing downstream: the paper uses only ε_n → 0. It is listed as a looseness, not
an error. (E1, the documentation-errata entry, deferred it to "a separate constants remark on family 238";
this is that remark.)

## Proof sketch

**Route A (primary).**
1. Exact enumeration for m ≤ 55, using the hook formula and MPFI intervals.
2. Peeling lemma (ballot words): f_λ ≥ f_μ·f_{(L,min(k,L))} for λ = (L, μ).
3. A 3-index recursion (size, first row, length) up to M_E, then a 1-D recursion up to M_1.
4. An analytic tail with a Catalan geometric series.

**Route B (independent).**
1. Exact enumeration for m ≤ 50. Here f is computed by the branching rule, and f^{−p/q} is enclosed by exact
   integer q-th roots.
2. Exact first-row hook identity f_λ = C(m,k) f_μ / Π with Π = ∏_{j≤μ₁}(1 + μ'_j/(L−j+1)). It is bounded by
   Chebyshev's sum inequality and AM–GM: Π ≤ (1 + k(H_L − H_{L−b})/b²)^b ≤ min{(m/L)^L, e^{k/(L−min(k,L)+1)}}.
3. A recursion like Route A's, using this bound instead of ballot words.
4. A different tail lemma (Lemma T: Robbins' binomial bound, three ranges), evaluated in `mpmath.iv`.

Both routes give the same enclosure of Z_u(m*). In both, every bound for m ≠ m* is strictly smaller than
that enclosure. Full statements and proofs are in `note.pdf` §2–4.

## Reproduce

```bash
cd scripts && ./run_all.sh        # ~8-11 min; SageMath (route A only) + python3 with numpy, mpmath
```
The script prints `OK` for each of six outputs:
- `exact_tables.py 55`: enumeration for route A. It uses NPROC=3 processes.
- `check_lemma.py`: unit test of the ballot lemma on all 4507 partitions with m ≤ 22. It finds 0 violations.
- `cu_bound.py u M0 ME M1`: **route A certificate** for the five values of u (run under `sage -python`).
- `eps_n_exponent.py`: interval certificate for ε_n.
- `independent/indep_Cu.py`: **route B certificate**. It shares no code with route A. Its output ends with
  `ALL VERIFIED`.
- `independent/indep_eps.py`: exponent bounds for ε_n using integer comparisons only.

## Scope and limitations (what is NOT claimed)

- **No swap-or-not / cryptographic claim.** That statement and its round counts (8132 rounds at N = 2³², λ = 128, and so
  on) are not part of this remark. The asymptotic statement is a few-line consequence of known results; see
  "Premise".
- **No new technique.** Every step is a standard peeling estimate. The content is the numbers and their
  certification.
- **C_{1/10}.** Route A gives an upper bound C_{1/10} ≤ 1267.27. Route B did not reproduce it with its
  parameters, so it is **not claimed**. Only C_{1/10} ≥ 562.0737 is confirmed twice.
- **Growth in u.** Only the lower bound log C_u ≳ π²/(6u log(1/u)) is stated. The internal notes also
  sketch an upper bound O(1/(u log 1/u)); it is not written out and is not claimed.
- **Not formalised.** The upstream Lean (`lean/docs/238.md`) defines C_u through `sSup` and proves only
  finiteness. Nothing here is in Lean.
- **ε_n.** This is a sharper reading of a displayed formula. It does not improve any theorem of the paper.
- **Small-d Thorp mixing times** (t_mix = 3, 5 for d = 2, 3) are in R5 (`small_d_exact.py`) and are not
  repeated here.
- **Novelty.** We searched GitHub code and repositories (queries "473/144", "Witten zeta symmetric group",
  C_u with Specht) and the web/arXiv (Witten zeta function of S_n; explicit sup_n ζ_{S_n}(s)) up to
  2026-10-09. We found only asymptotic and effective-tail results with unspecified constants: Liebeck–Shalev
  2004, Gamburd 2006, Magee–Naud–Puder arXiv:2003.10911 §3.3, and arXiv:2411.04347. We found no table of
  sup_m Z_u(m). To our knowledge the values are not recorded elsewhere, but they are easy to compute for small
  m, and the content of the remark is the certified supremum rather than the numbers.
