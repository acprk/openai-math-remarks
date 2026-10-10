# VERIFICATION — R21 (complete-intersection codes: GHW ≥ Cartesian, cases G1/G3)

All runs on 2026-10-09: CPython 3.10.12, numpy 2.2.6, SageMath 9.5, one core each.
`scripts/run_all.sh` reproduced all 4 expected outputs (`OK` ×4, about 3 min wall time).

## 1. Who checked what

| Item | Primary (review session, 10-07) | Independent re-implementation (this build) | Adversarial review (separate session) | Status |
|---|---|---|---|---|
| Proof of Theorem G (d_r formula, CB7 shift, Artinian reduction, EGH ⇒ down-set, final inequality) | REPORT §I4 sketch | written out in full in note §2 | checked line by line; edge cases Γ = X, r = dim C_d, d ≥ s | correct |
| CB7 normalisation h_X(t) − h_{X∖Γ}(t) = \|Γ\| − h_Γ(s−1−t) | fixed numerically (first script used s−d, corrected) | `cb7_and_reduction_check.py`: 3280 identity checks, 0 failures | matched against CSTPV Prop 6.8 proof; 5090 checks incl. F₄, F₈, F₉, 0 failures | holds |
| Down-set profile of Δh_Γ (what EGH delivers) | — | 720 random subsets, every profile realised by a box down-set | — | holds on tests |
| Cartesian code attains μ_r (Prop. 3) | Sage brute force | matroid GHW, all tested types incl. (4,4), (3,5), (2,3,3) | proof checked (deg G_γ, triangularity); 3482 Cartesian pairs, all equal | proved + checked |
| d_r(C_d(X)) ≥ μ_r on non-Cartesian CIs | `primary/ghw_check*.sage`: 1050 (d,r), 0 violations | `ghw_ci_indep.py`: 1424 (d,r), 0 violations, 117 strict; teeth control (shift s−d) violated 988 times | own Sage code: 3011 CIs, 38 189 (d,r), 0 violations, 135 strict; exhaustive over F₂, F₃, F₄ for small types; reducible curves, unions of lines; non-prime fields F₄, F₈, F₉ | no counterexample |
| Hypotheses of cited EGH cases | — | Abedelfatah Cor 4.4 text read (arbitrary field, a₁ ≤ … ≤ aₙ, f₃..fₙ split) | Caviglia–Maclagan (any field via flat extension; n = 2 automatic) and Richert for n = 2; Abedelfatah: arbitrary field, induction for n ≥ 4 only needs split-case EGH in n−1 variables | holds |
| Literature / novelty | — | wide dedup session (see §3) | own search | no prior statement for r ≥ 2 found |

The independent program shares no code with the primary Sage scripts. GHW come from the column matroid: d_r = N − max{|S| : rank G_S ≤ k − r}, with all 2^N subsets explored by DFS with an incremental echelon basis. μ_r is obtained by filtering all 2^{|box|} subsets. CIs are certified by a Hilbert-function test plus an exact count of distinct rational points; there are no Gröbner bases.

## 2. Adversarial review (separate AI session): PASS-with-fixes

The reviewer found no mathematical error and no counterexample. Its issues and our responses:

1. **Major, attribution.** r = 1 is already known in P³: CB12 in P³ (GK13 Cor 4.4; Caviglia–De Stefani 2021 Cor A) combined with CSTPV Prop 6.8. **Fixed**: note §1 and the README "Before" column.
2. **Major/Minor, attribution.** Some all-r special cases were already known: d = s−1 (arXiv:1611.06405 Prop 9), d ≥ s (trivial), and d = 1 in the equigenerated case (CSTPV Cor 6.6). **Fixed**: all listed.
3. **Minor.** Caviglia–De Stefani arXiv:2007.15467 Thm A (aᵢ ≥ Σ_{j<i}(a_j−1) for i ≥ 3) would give more unconditional types in P³. **Mentioned** in note Remark (iv) and in README "What is NOT claimed", but **not claimed**: per the task scope, only G1 and G3 are claimed, and we did not audit that paper's field hypotheses.
4. **Minor, rigour.** Abedelfatah's Thm 4.2 assumes full EGH in n−1 variables. For n ≥ 4 the induction behind Cor 4.4 only needs the split case. **Fixed**: added a sentence to the proof of the Corollary, and stated that the split generators are those of the largest degrees.
5. **Minor, reproducibility.** expected_output was empty at review time. **Fixed** (4 files). The field range in the computations paragraph was also inaccurate: the primary run includes F₅, and its (2,2,2) instances are not necessarily split. **Fixed.**
6. **Minor.** "for every q" should be "for every q ≥ aₙ" so that X_Car exists. **Fixed**; the note now also says that d_r ≥ μ_r holds for all q.
7. **Minor.** Caviglia–De Stefani 2021 was described too narrowly. **Fixed**: the note now says Cor A covers P³ and Thm B is a partial result in Pⁿ.
8. **Gate, not a math error.** Under the project rules this is a parallel generalization: the same EGH96 argument with r in place of 1. **Accepted.** It is recorded as a short remark only, with no paper and no Lean.

The reviewer's own scripts were kept outside this repository and are not shipped. Its numerical results are quoted in the table above as a third, independent check.

## 3. Duplicate search (2026-10-09)

A separate dedup session did the main search. arXiv sources:
- API abstract searches: GHW / generalized minimum distance + "complete intersection", Cayley–Bacharach + Hamming, EGH + code, and author Tohăneanu, over dates up to 2026-10-10.
- Full texts read: 2608.19978 (only v1 exists), 1812.06529, 1707.03285, 1601.07604, 1611.06405, 1712.00387, 2006.14717, 2106.14759, 2103.14106, 1212.2653, math/0611488, math/0311129, 2308.13728, 1905.12136.

Other sources:
- Semantic Scholar: citations of 2608.19978 (one, unrelated) and of GK13.
- Google Books: search inside Tohăneanu (2024). §6.2.3 contains only r = 1 statements.
- The Kreuzer publication page.

Conclusion: no source states the r ≥ 2 bound via EGH, or the unconditional P² / split cases for r ≥ 2.

## 4. Open gaps

- **GAP-1.** The full texts of Geramita–Kreuzer 2013, Cooper 2012 and EGH96 were not accessible: publisher pages were blocked, and the stored EGH96 PDF is an HTML error page. Our use of EGH96 CB7 is cross-checked against the CSTPV Prop 6.8 proof and numerically. Our reading of GK13 and Cooper comes from their abstracts and from citing papers.
- **GAP-2.** Richert 2004, the n = 2 EGH reference, was not read directly. The n = 2 case is also covered by Caviglia–Maclagan Thm 2 (read; any field via flat extension).
- **GAP-3.** The independent script certifies split CIs in P³ through plane restrictions. Instances of type (2,2,2,2) in P⁴ (G3 with n = 4) were not tested by any program.
- **GAP-4.** The numerics are consistency checks on small types and small fields only. The proof does not depend on them.
