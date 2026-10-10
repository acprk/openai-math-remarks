# VERIFICATION — R19 (Jacobsthal's function in F_q[x], first two values)

All runs: 2026-10-09, CPython 3.10.12, OR-tools 9.15, python-sat, SageMath (system install), single core.
`scripts/run_all.sh` reproduces all three outputs byte-for-byte (`OK` ×3).

## 1. Who checked what

| Item | Primary | Independent re-implementation | Adversarial review | Status |
|---|---|---|---|---|
| Covering model (Lemma 1.1 of the note) | internal review (D-numtheory §4, "Lemma 1") | proof written out in the note (CRT for low- and high-degree primes together) | **correct**. The reviewer also checked: f can be taken of large degree; the coset V_m itself contains 1 and is never covered; monotonicity of h_q in m | holds |
| K_q(1) = q, K_q(2) = 2q − 1 (all prime powers q) | internal review proof; CP-SAT for q ∈ {2,3,5,7} | branch-and-bound with its own GF(q), q ∈ {2,3,4,5,7}; Sage factorisation of the construction for q ∈ {2,3,4,5,7,8,9} | **correct**. The reviewer recomputed with its own Sage code: explicit polynomials and irreducibles, its own DFS **and** a GLPK MILP (a third solver), including q = 8 (15) and q = 9 (17) | holds |
| Equality characterisation (q ≥ 3) | internal review proof | Sage: all residue vectors for q ∈ {3,4,5}; optimum 2q−1 attained exactly by the q² affine vectors | **true**. The reviewer enumerated q = 3, 5, 7: each has q² optimal configurations, all concurrent. Wording fixed: the q − 1 extra primes each divide a distinct uncovered element, and M has no other factors | holds |
| Table K_2(3..5) = 4, 6, 9; K_3(3) = 8 | CP-SAT (`primary/run_table.py`) | branch-and-bound | reviewer's DFS + GLPK agree | holds |

## 2. Commands and key outputs

```
python3 primary/run_table.py              (CP-SAT, single worker)
  q=2 m=1..5: 2 3 4 6 9 ; q=3 m=1..3: 3 5 8 ; q=5: 5 9 ; q=7: 7 13   (all OPTIMAL)
python3 independent/indep_Kqm.py           (branch-and-bound, own GF(q))
  same values, plus q=4: 4 7 ; every m<=2 line ends "[claim ...: OK]"
sage independent/construction_check.sage
  q=2,3,4,5,7,8,9: omega(M) = 2q-1, all q^2 elements of f+V_2 non-coprime: True
  q=3,4,5: residue vectors attaining 2q-1 with all q lines: q^2, of which non-affine: 0
```

## 3. Adversarial review (separate AI session; not human peer review)

Verdict: **CONFIRMED-WITH-CORRECTIONS**. The corrections concern wording only.
- (low; applied) The equality statement now says explicitly that each of the q − 1 extra primes divides a different uncovered element and that M has no other prime factors.
- (low; applied) The note now says that the incidence count uses only intersections with the first line, so it is independent of the order of the lines.
- (gate; acknowledged) The reviewer flags that the result is a one-paragraph incidence count, i.e. a free function-field analogue of trivial small cases. It fails the project's "expert afternoon" gate and does not justify Lean. The README states this plainly. This repository collects such remarks, and no Lean was written.
- Duplicate search by the reviewer:
  - jdb19937/jacobsthal (full TeX) is integer-only;
  - Kalmynin–Konyagin 2302.00459 is about integer polynomial values;
  - Kumchev–McNew–Park 2310.02495 concerns k-free short intervals;
  - covering systems of F_q[x] (Azlin 2011; Li–Wang–Wang–Yi 2024; 2402.03810; 2408.10460) are a different problem;
  - Jamison and Brouwer–Schrijver are close in spirit (now cited);
  - no prior statement of K_q(2) = 2q − 1 was found.

## 4. Duplicate search (2026-10-09)

- **GitHub:**
  - repository searches: "jacobsthal", "jacobsthal finite field", "jacobsthal polynomial", "jacobsthal Fq";
  - jdb19937/jacobsthal cloned (HEAD a91ffe3, 2026-10-08) and its TeX grepped with `command grep -P` for F_q / function field / polynomial ring: no hits.
- **arXiv/web:** searches for Jacobsthal function with F_q[x], F_q[t], polynomial ring, function field, short interval and covering. The arXiv API itself was unreachable from this machine, so the arXiv check went through web search and the reviewer's checks.

## 5. Open gaps

- **GAP-1:** no formula for m ≥ 3, and no asymptotics. The values 2, 3, 4, 6, 9 for q = 2 are unexplained.
- **GAP-2:** no Lean (deliberately, as the result is elementary).
