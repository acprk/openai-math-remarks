# R7 verification record

Claim: with R₀ = (M−3)/2 and conditions (P1)–(P4) (see README), the proof of Theorem 1.1 in OpenAI's height-two Hahn–Wilson preprint goes through verbatim for every prime p ≥ 293. The original explicit condition (op:parameters) first holds at p = 17,073,437.

All runs: Python 3.10.12 (`/usr/bin/python3`), numpy 2.2.6, single thread, on 2026-10-09. Every script uses exact `int`/`Fraction` arithmetic.

## Layer 1: main computation (rescan session, re-run by the coordinating agent)

`scripts/main/hw319_pmin.py` encodes the conditions of `operator.tex` with line references: C1 (L322–326), CB (L503), CS (L590–597, with the L592 bound rewritten), CE (L507), first-band bounds (L583–603), J (L692–694), and the base stage (L546–562). It searches primes 11 ≤ p < 10⁵ and all odd M, R₀.

```
literal R0>=4 p_min = 293 witnesses [(11, 4)]
sharpened R0>=2 p_min = 199 witnesses [(9, 3)]
paper p_min = 17073437  count in window 5960
paper p_min over all M: 17073437
```

`scripts/main/hw319_detail.py` gives the per-M table: the C1 threshold is 127/199/293/397/521 for M = 7/9/11/13/15. At (p, M, R₀) = (293, 11, 4), the C1 margin is 259/803 and B₀³ − p = 219.

## Layer 2: independent re-implementations

1. `scripts/independent/indep_scan.py` was written in a separate session and does not import the `hw319_*` scripts. It uses exact centres d_i = ⌊(2i+1)p/(2M)⌋ in place of the paper's p/M ± 1 approximations, and band bounds computed from the actual centre differences. Output (`expected_output/indep_scan.txt`):
   `A paper-C1, literal J (R0>=4) -> (293, 11, 4)`, `B ... R0>=2 -> (199, 9, 3)`. It also lists the failing condition for every prime p < 293.
2. `scripts/independent/hw319_enumerate.py` was written in this adversarial session. It does not re-encode the inequalities. It builds the actual surviving lower-label sets of the induction table (op:induction) at p = 293, M = 11, R₀ = 4 for stages r = 2, 3, 4, 5 (X up to 2.16·10¹²; 11, 88, 704 and 5632 labels per side). It then computes the exact depth of **every** column/row pair and checks each claim the proof makes about them (claim → L-number):
   * same-cluster pairs never have an allowed depth (L597–598);
   * every k₁ = 1 pair has an allowed depth ≥ X/(2M) (hypothesis of Lemma op:coefficient-control);
   * max band-1 depth < min depth of all allowed k₁ ≥ 2 pairs (L601–603);
   * the band-1 depths lie in [X(1/M−2/p), X(1/M+2/p)], as printed in (op:first-band);
   * no cyclic wrap for column clusters in I (L613–615);
   * the two sides of each cluster pair have equal size B₀^{r−2} (L628–629);
   * the clusters reached by crossing columns lie in I and within three places of s₀ (L689–694);
   * d(I) contains a digit < δ and d(I−1) a digit > δ (L705–708);
   * the initial depth d = ⌊p/M⌋ (L547–551), C1, and CB.

   Observed: band-1 depth ratios b/X ∈ [0.086867, 0.094651] (stable in r); the first k₁ = 2 depth is at 0.179 X; the reached clusters are {5, 6} = {s₀, s₀+1}.
   The script also does the following:
   * It verifies (P1)–(P4) and the limiting cut condition for the chosen M at **every prime 293 ≤ p ≤ 2·10⁶** (0 uncovered; M ranges over 11..129). It checks the analytic tail inequality (M−5)³ ≥ 2M(M+3) for 13 ≤ M < 10⁵.
   * It re-runs the stage enumeration at (509, 11), (521, 13), (1009, 15) and (10007, 25).
   * It recomputes both p_min values (293 and 17,073,437) by its own search.
   * It runs four negative controls, all of which fail as expected: p = 283 fails C1; R₀ = (M−1)/2 makes I meet the endpoint M−1; p = 37 fails the band separation; (M, R₀) = (9, 3) fails CB and the literal three-places condition.
   * Exit code 0, `VERDICT: ALL CHECKS PASS`.

## Layer 3: adversarial review (this session, 2026-10-09). Verdict: **CONFIRMED**

Reading done:
* Read in full: `main.tex`, `operator.tex` (721 lines), `gluing.tex` (357 lines), and `obstruction.tex` L1–30, L290–345, L495–1199. That covers the criterion's statement, the coefficient extension, Lemma obs:principal, the growth lemma obs:growth, Lemma obs:right-annihilator, Proposition obs:dual-ghost and Theorem obs:nonmembership. Also read: `consequences.tex` L1–40, L160–200 and `transfer.tex` L15–40.
* A separate read-only sub-session read `geometry.tex`, `completion.tex`, `comparison.tex`, `transfer.tex` and `consequences.tex` line by line, and listed every size condition on p. It found: p > 5 (geometry L14, completion L16); p > 3 (geometry L629); p > 4 (geometry L702); D > 3 (completion L136, L171); D > 2 (completion L286, gluing L157); p > 8 / D > 12 (comparison L19, L25, L41, L345, L686; transfer L21); D > 6 (comparison L708); and q prime to p. The only references to the operator parameters are `consequences.tex` L14, L172–173 and L187. The remaining lines of `obstruction.tex` (L30–290 and L345–495: the abstract finite-tower criterion and the filtered microlocalization) were checked by grep for p, M, R₀, β and large-prime wording; no hit was found. They were not read line by line.

Every use of p, M, R₀, B₀, β in `operator.tex`, and how it fares after the change:

| Line(s) | Use | After the change |
|---|---|---|
| L85–96, L182–195 | displacement bound, κ ≤ p/(p−1), N_{p^s} ≥ p^{3(s−1)} | no size condition |
| L264–273 | (op:parameters) | replaced by (P1)–(P4) |
| L305–307 | h_g^{pX} ≡ 1 | general |
| L311 | hypothesis b ≥ X/(2M) | holds on band 1 by (P3); checked by enumeration |
| L322–329 | (op:correction-bound) | exactly (P1) |
| L346 | M < p | holds |
| L369, L383–388 | initialization, N₂ = p² | general |
| L469–493 | pivot rule, ≤ M protected digits, hyperplane count ≤ p · p^{−2} | general |
| L497–509 | centres, s₀, R₀, I, B₀, β, "I meets neither endpoint" | R₀ = (M−3)/2 gives I = {2,…,M−2}; β > 1/3 by (P2) |
| L547–564 | initial depth d = ⌊p/M⌋: qd < p², d ≥ p/(2M), no borrow, no second meeting depth | general for M < p; checked |
| L584–603 | first band, the (p/4+3) bound, separation | bound becomes (2R₀p/M+3); (P3), (P4); checked by enumeration r ≤ 5 |
| L605–637 | full/nonfull matching, no wrap for i ∈ I, equal side sizes | min I = 2; checked |
| L680–698 | finite crossings, "within three places … since M ≥ 81" | reached clusters {s₀, s₀+1} when p > 3M; R₀ ≥ 4 |
| L700–720 | unmatched counts B₀^{r−1}, needs digits < δ in d(I) and > δ in d(I−1) | R₀ ≥ 2 |
| obstruction L327, L1082, L1123 | β > 1/3 vs the O(N^{1/3}) of obs:growth | (P2) |

Things checked and found not to be issues:
* The bound (p/4+3)X/p² at L592 is tied to R₀ = ⌊M/8⌋. With the larger R₀ it must be rewritten; the sound replacement is (2R₀p/M+3)X/p², because the digit diameter of I∪(I−1) is ≤ 2R₀p/M+1 and the lower positions contribute < X/p². This rewrite is the only non-mechanical change, and the enumeration confirms it.
* The phrase "cyclic positive gaps p/M ± 1" (L506) is loose (the gaps are ⌊p/M⌋ or ⌈p/M⌉). The enumeration uses exact centres, so the claim is not needed in that form.
* Lemma obs:growth has constants depending on p, but no size hypothesis.
* Lemma op:coefficient-control at the X = p stage (initial depth): b = d = 26 ≥ 293/22 ≈ 13.3.

## GAPs (not covered)

1. **The correctness of OpenAI's proof itself was not re-verified.** That includes the operation calculus, the gluing, the finite-tower criterion, the transfer and the comparison. The remark is conditional on it in the same way the original theorem is.
2. External prime ranges are taken from the paper as quoted. Pstrągowski's algebraicity is quoted as p > n²+n+2, i.e. p > 8. Lee–Pstrągowski 2026, Theorems 1.3/5.4, and Lee 2026, Corollary 1.9 are quoted for odd primes. We did not check these against the sources. They sit on the main proof chain (comparison, exact fp-type, uniqueness of the BP⟨2⟩ form).
3. The label enumeration covers stages r ≤ 5 at p = 293 and r ≤ 3 (or r = 2) at the sample primes. Stages r > 5 rest on the analytic inequalities (P3)–(P4), which are proved for every r.
4. Parts of `obstruction.tex` (L30–290, L345–495) were checked by grep only.
5. The sharper variants (199, and 109–167) are not claimed. They change more of OpenAI's text and have not been adversarially reviewed.
