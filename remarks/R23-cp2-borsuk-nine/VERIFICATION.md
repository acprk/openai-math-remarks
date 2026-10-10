# Verification — R23 (family 156, CP² projector set, b = 9)

All runs used a single thread (OMP/OpenBLAS/MKL = 1). The headline claim is: for the nine points in
`scripts/P9.txt` (SHA-256 `07c6a2f1…dcab`) and T = 21/40, every nonzero x ∈ ℂ³ satisfies
max_i |⟨x,p_i⟩|² > T|p_i|²|x|². The same is certified at T = 263/500.

## Layer 1 — primary computation (`scripts/cert_interval.py`)

- Three charts x_j = 1 with free coordinates in [−1,1]⁴. These cover CP² after dividing by the coordinate of largest
  modulus.
- Interval arithmetic in numpy binary64, with `nextafter` outward rounding after every operation. Inputs are exact hex
  floats, and T|p_i|² is computed exactly in rationals and then enclosed.
- 16-way branch and bound.
- Results:

  | T | accepted leaves | boxes examined | min half-width |
  |---|---|---|---|
  | 21/40 | 1 897 683 | 2 024 192 | 2⁻¹² |
  | 263/500 | 2 282 103 | 2 434 240 | — |

  Both runs print `CERTIFIED`.
- Negative controls:
  - T = 53/100 fails with an explicit witness x = (1, 0.67578125+0.30859375i, −0.26171875+0.08984375i), whose
    normalized max overlap is 0.529505.
  - Deleting p₉ fails at the first level. This is expected: by the lower bound, no eight caps with T > 1/2 can cover.
- This replaces the earlier float-only certificate (internal audit, 2026-10-07). That certificate used a Lipschitz bound in a
  trigonometric parametrization and no directed rounding.

## Layer 2 — independent re-implementation (`scripts/independent/cert_exact.py`)

- Written separately; it shares only the input file. It uses no floating point: exact integers after converting the hex
  inputs exactly.
- Bound: an exact second-order Taylor expansion of the quadratic H_i (gradient term plus |Q| bound). Splitting is by bisection
  along the widest coordinate, depth-first.
- Results:

  | T | leaves | nodes | max depth |
  |---|---|---|---|
  | 21/40 | 16 876 | 33 749 | 36 |
  | 263/500 | 19 189 | 38 375 | 55 |

  Both print `CERTIFIED`.

## Layer 3 — adversarial review (separate AI session, instructed to refute)

**Verdict: UPHELD-WITH-FIXES. No Critical findings. The mathematics and both certificates are upheld.**
What the reviewer did:
- Re-ran `run_all.sh` on a copy (numpy 2.2.6); all lines were OK.
- Wrote a **third certificate** with exact `Fraction`s. It differs from both existing ones: it uses the opposite inner-product
  convention Σ conj(x_k)p_k, complex-rational evaluation, and a different lower bound H(c) − Σ|∂_kH|r_k − T|p|²Σr_k²
  (valid because Hessian/2 = td·Gram − tn·|p|²I). It also sums leaf volumes, which tile each chart box exactly (total 48 = 3·2⁴).
  Results: T = 21/40 CERTIFIED (16 879 leaves); T = 263/500 CERTIFIED (19 173 leaves); T = 0.53 no certificate (timeout,
  no false certification).
- Ran a global minimization with SLSQP in epigraph form from 4 000 random starts. It found min_x max_i |⟨x̂,p̂_i⟩|² =
  0.5260548764555, with four caps active, which matches our 0.526055. No counterexample was found.
- Audited chart coverage, the Re/Im formulas, the conjugation convention (irrelevant), outward rounding (including
  `sq()`, subnormals and the T|p|² enclosure), the integer scaling and bisection bookkeeping in `cert_exact.py`, the cap
  lemma, and every step of the complex lower-bound transfer (U(3)-invariant measure, ∫uu* = I/3, continuity of Q^{1/2} and
  A_±, disjoint supports on orthogonal spectral ranges, odd map S⁸ → S^{m−1}). No errors were found.

Findings and how they were handled:
1. **Major (framing).** The earlier wording "closes the route of complexifying the 156 witness" overclaimed. CP² is the
   complex analogue of RP² (k = 3), where the preprint does not exclude equality either. The analogue of the k = 4 witness is
   CP³ ⊂ ℝ¹⁵. **Fixed:** the README and note now claim only that no complex-projective projector set gives a
   counterexample below dimension 9, and they state the analogy caveat explicitly.
2. **Minor.** The note said the minimum half-width was 2⁻⁹; the actual value is 2⁻¹². **Fixed.**
3. **Minor (triviality).** The lower bound is a parallel transfer, and the overall result is elementary. The upper bound is not
   immediate from a classical configuration: the Hesse SIC gives exactly 1/2. **Fixed:** both points are stated in the README
   and note.
4. **Minor.** Add a data-mutation control. **Fixed:** `control_drop9`.

## Reproduction on this machine (2026-10-09)

`scripts/run_all.sh`: all 8 lines print `OK`, in about 68 s wall time on one core.

## GAPs (not covered)

- The primary certificate assumes IEEE-754 round-to-nearest binary64. The exact-integer certificate (and the reviewer's
  `Fraction` certificate) remove this assumption.
- The lower bound is checked by hand (by us and the reviewer), not formally. There is no Lean (per project policy for
  elementary/negative results).
- The literature search was by web search only (GitHub, arXiv, Kalai's survey and blog, up to 2026-10-09). No statement of
  b(X_CP²) was found. This is not exhaustive.
