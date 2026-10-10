# VERIFICATION — R2 (the profile axioms of ω ≤ 9/4 stop at 9/4)

Date: 2026-10-09.

Machine: Linux x86-64. All runs used `/usr/bin/python3` (Python 3, sympy 1.14, scipy 1.15.3, mpmath 1.3.0) with one thread.

No Lean build was run for this remark.

## 0. Falsification step: is the axiom list complete?

The barrier would be void if the paper used some property of P beyond (Pos), (S), (B), (C), (T), (U). We checked line by line in `openai/math` `preprints/Matrix-Multiplication-Nine-Fourths-October-2-2026/build/paper.tex`. The file is identical at `adc7f1241b42` and at HEAD `08fcd618213a` (`git log adc7f12..HEAD -- <dir>` is empty).

| use of P in the paper | lines | property |
|---|---|---|
| definition of the profile | 420–422 | — |
| basic facts (symmetry, P(1,b)=b, rank bound) | 426–429 | (S), (B), (U) |
| Lemma 4.1 (concavity) | 442–447 | (C) |
| Lemma 4.2 (tripling) | 526–530 | (T) |
| Lemma 5.1 hypotheses | 610–614 | exactly (Pos), (S), (B), (C), (T) |
| Lemma 5.1 proof | 617–650 | uses only those hypotheses |
| final step | 652–660 | adds only the diagonal case of (U) |
| detecting characters | 198–202, 662–670 | used only after t ≤ 3/4; constrains nothing about P |

Lean cross-check, under `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/`:
- `structure ScalarProfile t` (`Growth/Profile.lean` l. 23–33) has exactly the fields `positive`, `symmetric`, `boundary`, `concave`, `tripling` and `rank_bound`.
- `meanExponent_le_three_quarters` (`Polynomial/Inequalities.lean` l. 110–114) proves t ≤ 3/4 only via `toScalarProfile` (l. 97–107) and `ScalarProfile.exponent_le_three_quarters` (`Profile.lean` l. 108–116).

**Verdict:** the list is complete. No property of P was missed, so the barrier stands for the paper's axiom system.

**Tightness note (exact checks).** (T) holds with equality exactly when h ≥ a−1, and (C) holds with equality across b = a−1, a, a+1. Every step of Lemma 5.1 up to `H_a ≥ ∏(1+1/(3m))` is an equality on P*.

## 1. Main computation (Claude, this session): independent re-implementation

None of these scripts shares code with `barrier_check.py`, `barrier_extra.py` or the Codex programs.

| script | what it checks | size | result |
|---|---|---|---|
| `r2_profile_axioms.py` | (Pos), (S), (U) as `P*³ ≤ (a+b−1)⁴` (equality only at (1,1)), (T) with its exact equality set; tightness of every step of Lemma 5.1; `m ≤ φ³ ≤ 2m−1` and both induction identities; the row recursion of the minimality proof; the (T) increment for h < a; Wendel's bound with an interval-certified `Γ(4/3)^{-3} ≤ 1.40435051` | a,b ≤ 200 (b < 800 for (C)) | all checks PASS (exact rationals) |
| `r2_extension_families.py` | Table 1 inequality families (A), (B), (C), (D), chains, Kronecker, (F), (O), (L), the sandwich, wrap, peel, CW, T⊗C, the zero/rest head and tail, and the (A2) slack closed form | about 0.9 M comparisons | ALL PASS (exact; 3/4-powers bracketed by integer 4th roots) |
| `r2_tensor_constructions.py` | coefficient- and support-level degenerations: OpenAI tripling (a,h ≤ 12), alternating k-sector chains (k = 3,5,7), peeling, truncated → CW, the zero/rest support, the D-adapted bivariate products (7600 basis-pair products; covers Lemma 4.1 and families B, C, D), and the two-stage weights for D | small sizes | ALL PASS |
| mutation test (ad hoc) | deliberately corrupted variants of the above | — | 7600/7600, 1784/7600 and 132/144 failures detected, so the checks are not vacuous |
| `r2_least_profile_lp.py` | FLOAT LP: min P(n,n) under (S), (B), (C), (T) on a (4n−1)² grid | n ≤ 12 | agrees with P*(n,n) to 1.4e-14 (numerical only) |

Commands and expected outputs are in `scripts/` and `scripts/expected_output/`.

## 2. Codex programs (OpenAI Codex, gpt-6-astra), re-run by Claude

| program | result on re-run |
|---|---|
| `codex_round1_verify.py` | 239,009 exact comparisons plus 780 apolar and 1,296 bivariate coefficient checks: PASS |
| `codex_round2_verify.py` | (C) 700,416, D-certificate 65,536 and the other families: PASS |
| `codex_round3_verify.py` (default) | table for n ≤ 30 reproduced to 4 decimals; zero/rest head and tail PASS; truncated recipes q = 2..29 below q |
| `codex_round3_verify.py --all-cuts` | 8,712 root cut triples for n ≤ 12; all upper endpoints below `F(n,n)^{3/4}` |

## 3. Status of each "balanced-assignment" inequality class

**P** = proved for all sizes. **E** = exact finite check only. **IC** = interval-certified finite recipe search. **N** = numerical only.

| class | status | proof location | finite evidence |
|---|---|---|---|
| paper (C), (T), any q per permutation | P | note Thm 1; balanced-lift reduction (Codex R1 §7) | ours |
| alternating odd-k chains | P | tangent line `φ(a)(b+(a−1)/2) ≥ P*` (this note; also the auditor in EC THM-4568) | E (k = 3, 5); tensor (k ≤ 7) |
| (A) transfer `C⊗B_Z(m)` | P | Codex R1 §2 (Jensen; we re-checked the slack closed form exactly) | E; tensor: Codex only (780 coefficients) |
| (B) `C⊗B_X(m)` | P | Codex R1 §3 (Jensen) | E; tensor: ours + Codex |
| (C) bivariate shared output | P | Codex R2 §2 | E; tensor: ours + Codex |
| (D) bivariate full direct sum | P | Codex R2 §3 (rational certificate) | E (true 3/4 power); tensor: ours |
| (F) Fourier with cost | P | monotonicity | E |
| (O), (L) outer and dot bounds | P | `φ(a) ≤ 2a/3` | E |
| wrap `C(n,n) ≥ ℂ[ℤ_n]` | P | `φ³ ≥ n` | E |
| peel `C(n,n) ⇝ 1 ⊕ C(n−1,n−1)` | P | Codex R3 (5)–(6) | E; tensor: ours |
| small-CW laser | P (assuming the standard small-CW laser value) | Codex R2 §4 | E; support: ours |
| scalar targets from any product | P | flattening sandwich (Codex R2 (5), (12)) | E (sandwich) |
| one matrix block, monomial | P | abelian-group lemma (Codex R2 §5) | — |
| `T_m ⇝ ⊕ C` | P | Codex R2 §6.1 | constant check |
| `T_m ⊗ C(n,n)`, m ≥ 2 | P (sufficient condition only) | Codex R2 §6.2 | E |
| zero/rest laser, children ≤ border rank | P (computer-assisted head q ≤ 15) | Codex R3 §5 | E/IC; re-run |
| Kronecker restriction | E | — | a, b, a′, b′ ≤ 9 |
| consecutive-interval recipes | IC | — | n ≤ 30; all root cuts n ≤ 12 |
| LP minimality of P* | N (proved separately in Thm 1(ii)) | — | n ≤ 12 |

## 4. Prior-art search (done 2026-10-09)

We searched GitHub code (`gh search code` for "3h+a-1", "shifted tripling" and "symmetrized profile"), GitHub repositories, arXiv (2610.10494 and the October 2026 listings found by web search), the web (VibeMathed, press coverage) and `<local path>`.

- **Found prior work.** The least profile and the "9/4 is the ceiling of the argument" statement already appear in:
  - `eliottcassidy2000/math` THM-4568 (commit `450388a22f72`, 2026-10-07T16:03Z), with a proof, an LP check and an audit;
  - `GasStationManager/SafeVerifyAgent` `audits/mm94-intake/READ-F.md` §2d (2026-10-08), numerically.

  Copies of the eliottcassidy2000 files are kept for provenance in `0940project/reviews/A-fastarith/allfields/prior_art_eliottcassidy2000/`.
- **Their "gadget barrier" overlaps our chain row.** It covers shared-leg gadgets `P(a, mh+c) ≥ mP(a,h)` and the minimal shift `c ≥ (m−1)(a−1)/2`.
- **Not found elsewhere.** We found the enlarged families nowhere else: the A–D determinant/bivariate filtrations, peel, CW, 1-block, T→C, T⊗C and the zero/rest laser, together with the balanced lift.
- **Other sources.** `selanavot/matrix-multiplication-all-fields`, `jzuiddam/omega-nine-quarters-all-fields`, arXiv:2610.10494 (Uznański, rectangular) and VibeMathed contain no barrier statement.

## 5. Adversarial review

A separate AI session (a Claude subagent) reviewed note.tex and README.md with read-only access and the instruction to refute. It re-checked the axiom list, every step of Theorem 1 (exact Fractions, a,h ≤ 40, b < 120) and three rows of Table 1 (chain, wrap, D).

It found **nothing fatal**. Its should-fix items were all applied:
- Theorem 1(iii) now says "for t > 0".
- Table 1 now separates "λ* satisfies" from "the relation is valid", and adds the odd-k chain weighting.
- The attribution of P* and of the equal-h chain is credited to THM-4568.
- The README's claims of contribution were toned down.

Its cosmetic items (counting six properties, defining r in row D, the a′ typo, the existence-lemma step) were fixed as well.

Not human peer review.

The Codex rounds also acted as adversarial checks on Claude's original BARRIER.md. They corrected an overstatement: BARRIER.md said "only one route remains" and "λ* is not a true character". The present text says instead that several routes remain open and that λ* is not *claimed* to be a character.

## 6. Gaps (not covered)

- No Lean formalization of Theorem 1. Instantiating `ScalarProfile (3/4)` with P* would be straightforward. It was not run, per the build rules and the gate in project CLAUDE.md.
- Tensor-level check of family (A) (apolar `B_Z(m)` filtration): Codex's 780-coefficient check only, for a, b, m ≤ 4.
- Tensor checks at small sizes support the general proofs but do not replace them.
- The small-CW laser value, tight-support laser theorem (Le Gall 2014) and border rank of truncated multiplication are used as standard facts.
- Open routes (multi-block laser after basis change, higher-power repartitionings, T⊗C outside the sufficient condition, asymmetric relations) are not addressed.
