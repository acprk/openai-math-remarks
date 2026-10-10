# VERIFICATION — R20 (Subset Sum exponent 0.49 → 0.4889; dense 0.4741; Proposition C)

All runs were on 2026-10-09 with CPython 3.10.12, mpmath 1.3.0, NumPy 2.2.6 and python-flint 0.9.0, one
thread each. `scripts/run_all.sh` reproduced all six expected outputs byte-for-byte (`OK` ×6, including
`independent/indep_check`).

**Premise for everything below: OpenAI's preprint (family 138) is correct.** None of the checks re-verify it.

## 1. Who checked what

| Item | Primary derivation | Independent re-implementation | Adversarial review | Status |
|---|---|---|---|---|
| σ-variant of alias extraction (coverage, sizes, cumulative caps) | internal session (2026-10-09); proofs written out in note §2 by this build | exponent inequalities recomputed in exact rationals + arb (`independent/indep_check.py`, Part 2) | round 1 (constant-optimisation reviewer): accepted, three holes A/B/C; round 2 (this build): no fault in the σ-variant | holds (premise) |
| Repairs A (cumulative caps), B (guard in pair exponent), C (checksum certified at τ₀ = 0.48888) | this build (note §2–3, `chain_check.py`, `hardcore_certificate.py main`) | Part 1–2 of the independent script | round 2 checked all three | closed |
| New densities, κ_f ≠ κ_g, Theorem 4.1 unchanged | constant-optimisation (differential evolution) | — | round 1: identities split correctly; round 2: re-derived u, v and checked that Thm 4.1 uses κ only through normalisation and relative entropies | holds |
| Filtering moments (the only numerical input) | `hardcore_certificate.py` (mpmath intervals, N = 128, analytic trapezoid + strip bounds) | `indep_check.py` Part 1 (python-flint balls, N = 192, strips 0.4 / 0.34, scores from the definitions of D, A, E, its own strip bounds) | round 2: own non-interval computation at N = 512, 1024; strip algebra re-derived | certified twice, same 12 digits |
| Dense corollary 0.4741 | this build | Part 1 (dense row) | round 2: **must-fix**, the checksum ran at τ = 0.4741, so its time 2^{(0.4741+20ε)n} exceeded the claim; repaired: run at τ = 0.47409, λ = (8.447, 8.236), recertified | closed |
| Proposition C | earlier internal session (short proof + exhaustive check on 3 groups) | `theorem_c_exhaustive.py` (closure BFS) and `indep_check.py` Part 3 (kernels of all endomorphisms) on 7 groups, same subgroup counts | round 2: proof correct; third enumeration method on (ℤ/4)², (ℤ/8)², (ℤ/49)², (ℤ/8)³ (802 subgroups), (ℤ/4)⁴ (1983 subgroups); 0 violations | holds |

## 2. Adversarial review, round 2 (separate AI session, 2026-10-09): **CONFIRMED WITH FIXES**

The reviewer was told to refute the note. It could not refute Theorem 1.1 (under the premise) or Proposition C.
Every number in note.tex matched the script outputs. Its own computation gave the same moments:
main 0.5170787970216 / 0.5169603832448, dense (old τ) 0.4935138 / 0.4996160, paper 0.5149740 / 0.5148515.
It also reproduced τ_cs(0.5176) = 0.488829 (new densities, with slack) vs 0.489561 (the preprint's).
Findings and what was done:

1. **must-fix — dense corollary ran at τ = 0.4741 exactly.** Fixed: τ = 0.47409, λ = (8.447, 8.236);
   log₂I₁ < 0.493598, log₂I₂ < 0.499710, MSE exponent < −0.000199, λ₂(4a+3) = 26.61 < 30; time exponent
   0.4740902 < 0.4741 (`hardcore_certificate.py dense` now prints this check).
2. **minor — the preliminary cap dropped the +1/2 rounding slack the preprint keeps.** Fixed: the cap is
   2^{⌈0.488825n⌉}.
3. **minor — the number of prime draws for P′ is now 2^{ρn}·poly, not polynomial.** Fixed: the note states the
   total cost 2^{(ρ+γμ/2)n} = 2^{0.01997n} and the coupling failure 2^{ρn}e^{−c₀n²} = o(1); `chain_check.py`
   lists the component.
4. **cosmetic — the generation bound mixed in poly(b).** Fixed: exponent ρ + μ + 1/4 + 10⁻⁵ = 0.34139.
5. **cosmetic — "same bounds to 12 digits" says little about independence.** Reworded: the agreement is
   expected; the independence lies in the code, the library and the error analysis. The reviewer could not
   open the Ye–Li abstract. This build fetched it through the arXiv API ("O(2^{0.499999n}) … representation
   technique framework of Randolph and Węgrzycki … partial match"), so the citation is accurate.

Round 1 (constant-optimisation reviewer, earlier on 2026-10-09) had accepted O(2^{0.4889n}) and reported holes A/B/C. It
also proved the quadrature gap with a rigorous banded bound (the basis of `hardcore_certificate.py`), and its
independent re-optimization found ≈ 0.48884 for the whole chain.

## 3. Commands and key outputs

```
python3 scripts/chain_check.py
  preliminary cap 0.488825 ; lists 0.4887769817 ; pairs (incl. guard) 0.4888560287 ; generation 0.3413897147
  P' sampling 0.01996560465 ; sparse 0.516455762 < eta ; checksum 0.4888802 ; tau_tot = 0.4888802 ; gap 1.98e-5
  coverage: K sigma mu = 0.49989 < 1/2 ; theta = 0.02240492908 < theta' = 0.0226 ; gamma < 0.371412
python3 scripts/hardcore_certificate.py main    -> log2 I1 < 0.517078797022 , log2 I2 < 0.516960383245 ; MSE < -0.000420723 ; CERTIFIED
python3 scripts/hardcore_certificate.py dense   -> log2 I1 < 0.493598280852 , log2 I2 < 0.499710378479 ; MSE < -0.00019914 ; CERTIFIED
python3 scripts/hardcore_certificate.py paper   -> log2 I1 < 0.514974048351 , log2 I2 < 0.514851451273 (< 0.5168, preprint's Prop. 3.1 holds)
python3 -u scripts/theorem_c_exhaustive.py      -> 15 / 129 / 37 / 83 / 23 / 76 / 45 subgroups, 0 violations, equality <=> free quotient
python3 scripts/independent/indep_check.py      -> same bounds (arb, N = 192), Part 2 all ok, Part 3 same subgroup counts, PASSES
```

## 4. What is exact and what is numerical

- Exact / rigorous: the two moment inequalities (interval and ball arithmetic with analytic quadrature
  bounds), every exponent inequality (interval arithmetic or rationals), and Proposition C (proof plus
  exhaustive checks).
- Numerical only (stated as such): the ledger of note §5 ((i) 0.48977, (ii) 0.48911, (iii) 0.48881), the
  τ_cs values, the dense value ≈ 0.4748 with the preprint's densities, and the heuristic framework limit ≈ 0.4798.

## 5. Remaining gaps

- **GAP-1 (premise).** The preprint (family 138) has no formal proof, and we know of no peer review. Everything
  here depends on it.
- **GAP-2.** No Lean. The note's proofs are prose modifications of the preprint's proofs.
- **GAP-3.** The framework limit 0.4798 is heuristic (idealized checksum vs. alias-count lower bound); it is
  not a theorem and is not used anywhere.
- **GAP-4.** Asymptotic only; with ε = 10⁻⁸ the bounds take effect only for astronomically large n (already so
  in the preprint).

## 6. Duplicate search (2026-10-09)

See README "Duplicate search". GitHub repository and code search found nothing on improving 0.49. The arXiv API
(40 most recent "subset sum" entries) found Ye–Li arXiv:2610.10446 (O(2^{0.499999n}), a different method), which
is cited.
