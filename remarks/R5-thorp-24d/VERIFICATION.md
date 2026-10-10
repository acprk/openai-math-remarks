# VERIFICATION — R5 (Thorp shuffle constant)

All runs were on 2026-10-09 with CPython 3.10.12, mpmath 1.3.0 and NumPy 2.2.6, each pinned to one
core. `scripts/run_all.sh` reproduced all five expected outputs byte-for-byte (`OK` ×5, including
`independent/indep_constants`).

## 1. Who checked what

| Item | Primary derivation | Independent check | Adversarial review | Status |
|---|---|---|---|---|
| Route O ingredients M1–M7 (θ→0, KL chain, local look-back, any σ, b blocks, mixed factors, cap) | internal session | separate session re-derived each step; reran `const_chain.py`; recomputed 2/L_b independently; checked search completeness | same separate session (did not refute; corrected the σ-window of (3,4,5,5) to σ < 2/59) | holds; reproduced here |
| Lemma 2.3 (closed form of the recurrence) | this build session (new, tighter than the internal session's: η* in place of η/(1−ρ)) | worst-case recurrence iterated at 40-digit precision; bound held at every step (`energy_lemma_check.py`; max ratio 2.6e-3) | none yet | proof written in note §2 |
| Route H (two factors, Σe(b_i) > 1) | this build session | HS bound re-derived by hand (coset Plancherel, note Prop. 3.2) and matched to companion Thm 1.1 (1.2), which allows any b and any u > 0 | **independent adversarial review (separate AI session), verdict CONFIRMED**: checked line by line against companion Thm 1.1 and upstream Lean `ComplementaryCosetsFiniteStatement`; separate interval certification (`scripts/independent/`) | holds (GAP-1 closed) |
| Physical vs cyclic frame (block lengths ⌊c_i d⌋ not multiples of d) | clarified after the adversarial review (raised by the reviewer) | physical law = R^T · (cyclic-frame law); fixed left translation preserves marginal TV and coset caps (note §4, "Frames") | — | proof written in note §4 |
| Constants 23.5394 / 44.3006, optimality, σ-windows, fractional lower bounds, ledgers | this build session | `constants_certificate.py`: exact rationals plus mpmath.iv outward rounding | adversarial reviewer: `independent/indep_constants.py` (separate code), same optima, runners-up, σ-windows and ledger | certified twice |
| Small-d unit test t_mix(2)=3, t_mix(3)=5 | earlier internal session (floating point) | `small_d_exact.py`, re-implemented in exact `Fraction` arithmetic | — | agree |
| Energy identity (Lemma 2.1) in simulation | constant-optimisation `energy_mc.py` | rerun d=8, b=4, 4000 runs: RHS = E_t within MC error | — | consistent |

## 2. Commands and key outputs

```
python3 scripts/constants_certificate.py
  route H optimum (4, 5):        23.53941309729 <= c <= 23.53941309729 ; runner-up (3,3,4) >= 25.5889
  route O optimum (3, 4, 5, 5):  44.300608089059 <= c <= 44.300608089059 ; runner-up (4,4,4,5) >= 44.30299
  sigma windows: (4,5) any sigma < 2/9 ; (3,4,5,5) any sigma < 2/59 ; (4,4,4,5) any sigma < 2/19
  fractional relaxation: inf_(real b>=3) kappa(b)/(1-2/b) in [20.738135, 20.738179] (b* ~ 3.8647)
  integer points: min = kappa(4)/(1/2) = 20.763572  => route H > 20.7636 d, route O > 41.5271 d
  (strict lower bounds for the present ingredients; the earlier figure 20.7326 was a valid but
   looser lower bound from coarser cells, not the infimum and not a limit)
  1600/c_H = 67.97 ; 512/c_H = 21.75 ; 1600/c_O = 36.12
python3 scripts/independent/indep_constants.py   (adversarial reviewer's code; expected output in
                                                  scripts/independent/expected_output/)
  route H optimum (4, 5): [23.53941309729, 23.53941309729] ; runner-up (3, 3, 4) >= 25.58892221
  route O optimum (3, 4, 5, 5): [44.3006080890591, 44.3006080890591] ; runner-up (4,4,4,5) >= 44.30298538
  sigma windows 2/9 and 2/59 ; real minimiser b* = 3.8646511, min kappa/e = 20.73817916 ;
  integer min kappa/e = 20.76357 ; route H ledger 1600 -> 539.98 -> 134.99 -> 68.44 -> 46.96 -> 42.96 -> 23.54
python3 scripts/energy_lemma_check.py   (b=4,5; theta=0.1; d=40,80,160)
  excess over 2d/L with look-back: 48 -> 55 -> 62 (b=4), 57 -> 65 -> 73 (b=5): O(log d)
  without look-back: 1.76d -> 1.45d -> 1.27d ; closed form never violated
python3 scripts/small_d_exact.py
  d=2: TV(t=3) = 1/6 <= 1/4, TV(t=2) = 1/3 ; d=3: TV(t=5) = 1171/5760, TV(t=4) = 169/315
```

## 3. What is numeric and what is exact

- The proof of Theorem 1.1 is analytic and uses no computation.
- Corollary 1.2 (the specific constants and their optimality) is a finite computation. It is fully
  certified: rationals for admissibility, interval arithmetic for logarithms, and a proved elementary
  lower bound κ(b) ≥ 2(2b−1)ln2 for completeness.
- An earlier internal script (`const_chain.py`) used double precision for the search and certified only
  (4,4,4,5) with intervals. That gap is closed here.
- `energy_lemma_check.py` (40-digit floating point) and `energy_mc.py` (Monte Carlo) are illustrations.
  Nothing depends on them.

## 4. Duplicate search (2026-10-09)

- **openai/math (HEAD fd4aeeb):** every family-238 preprint was searched for explicit constants. Found:
  1600d (main paper), 2048d (*Conditional information…*), 32800d and 512d (*Random coordinate frames…*,
  Cor. 11.4), and 16,040,400d and larger (*Conditional permutations…*). The remaining siblings state only
  "an absolute number of sweeps".
- **GitHub:** a repo search for "thorp", code searches for "Thorp shuffle", "Thorp 1600d" and
  "Thorp mixing constant", and yunjiangster/thorp-shuffle-mixing (last push 2026-10-03; 5.1537·d^{5/3},
  weaker in order).
- **arXiv:** an API query for "Thorp shuffle", sorted by date. The newest hit is 0912.2759 (Morris,
  O(log³ n)). No 2026-10 entries were found.

## 5. Independent adversarial review of route H (separate AI session, 2026-10-09)

A separate AI session (not human peer review) tried to refute route H. Verdict: **CONFIRMED**.
- It checked note Prop. 3.2(H) and §4 line by line against the companion's Theorem 1.1, estimate (1.2)
  (disjoint nonempty blocks of arbitrary sizes, every u > 0), and against the upstream Lean statement
  `ComplementaryCosetsFiniteStatement`, which formalises the companion's finite coset-to-Fourier estimate.
- It certified the constants independently with its own interval-arithmetic code
  (`scripts/independent/indep_constants.py`, mpmath.iv, exact `Fraction` admissibility, completeness via
  κ(b) ≥ 2(2b−1)ln2): route H optimum (4,5) = 23.53941309729, runner-up ≥ 25.5889; route O optimum
  (3,4,5,5) = 44.3006080890591. The script and its expected output are included and are run by
  `run_all.sh`.
- Clarifications added to the note after the review: (a) the physical/cyclic frame for block lengths that are
  not multiples of d (harmless, proof in §4); (b) the fractional relaxation infimum is 20.7382, not the
  20.7326 previously printed (that number is only a looser lower bound; "limit" wording removed);
  (c) attribution: the ×4 two-factor step is already in the family's 512d (Random-coordinate-frames
  Cor. 11.4, `shear:cor:512`), so the remark's own gain over the family's best is 512/23.54 = 21.75.

## 6. Open gaps (GAP)

- **GAP-1 (closed 2026-10-09):** route H (24d) has now had an independent adversarial review (§5).
  Route O (45d) had one before.
- **GAP-2:** no explicit d₀.
- **GAP-3:** no Lean.
- **GAP-4:** the result relies on the correctness of the main paper (formalised upstream for 1600d) and
  of companion Thm 1.1. Upstream Lean covers the companion's finite estimate, but not its two-factor
  corollary.
