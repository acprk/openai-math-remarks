# R15 verification record

All computations were run on the shared host with ≤ 3 threads under `nice`; no Lean/lake. Raw outputs are in `expected_output/`.

## 1. Main computation (this session, code written from scratch for R15)

| Item | Tool / script | Result |
|---|---|---|
| Prop. 1 (square of a split character) vs direct decomposition of χ₊² in GAP's A_n tables, **every** self-conjugate λ ⊢ n, 5 ≤ n ≤ 20 | `an_square_check.g` (GAP generic tables; A_n constituents identified by restriction, not by GAP's labels) | `TOTAL formula mismatches: 0` (`an_square_check.out`) |
| (b) S_n-universal but A_n-non-universal, ε = +1 | same run | exactly (3,1,1) ⊢ 5 and (4,4,4,3) ⊢ 15 among n ≤ 20; for (4,4,4,3) the only missing A_15 constituent is Res χ^(14,1) (g = 1, a = −1) |
| (c) g(ρ_m,ρ_m,λ) mod p = 2³¹−1 for all λ ⊢ N_m, exact a(λ) = χ^λ(ĥ) | `kron_stair.cpp` (MN on β-sets, head/tail split, Horner) | m = 3,4,5,7,8,9,11,12: `eps=1`, 0 collisions g ≡ \|a\| (a<0), 0 cases g ≡ 0 (a=0), g(ρ,ρ,ρ) ≢ 1 ⇒ χ₊ universal. m = 6, 10: ε = −1 (collisions are exactly the trivial and sign characters). Wall time m = 12: 1434 s, 11 GB |
| internal checks | built into `kron_stair` | Σ g·dim λ ≡ dim(ρ)² and Σ g·a ≡ 1 (mod p) "OK" for every m ≤ 12; g(ρ,ρ,(N)) = g(ρ,ρ,(1^N)) = 1 |
| threshold independence | T = 5 vs 7 (m = 10), T = 7 vs 9 (m = 11) | bit-identical vectors (FNV hashes equal; `threshold_consistency.txt`) |
| exact integers, m ≤ 9 | 4 primes + `crt_exact.py` | exact identity Σ g·dim = dim(ρ)² over ℤ: True for m = 3,4,5,7,8,9; g(ρ,ρ,ρ) ≡ 1 mod 4 (as Prop. 1(3) predicts); minimal margins g − \|a\| = 2, 2, 14, 34, 6, 62 |
| full-vector cross-check | `gap_stair_g.g` (GAP S_n tables) vs kron_stair | identical for m = 3,4,5,6 (`gap_vs_kron_stair.txt`) |
| earlier independent value | earlier separate hook computation | g(ρ₈,ρ₈,(8,1^28)) = 84465 reproduced |
| 205 lower bound w_m | `w205.sage` (Sage 9.5 lrcalc, own Kostka via Pieri) | w ≤ g for every λ, m ≤ 8; two-sided test max(w(λ),w(λ′)) ≥ \|a\|+1 passes for all a<0, m ≤ 8; agrees **exactly** with the earlier separate-code tables (not shipped) (0 differences out of 17977 at m = 8) |
| hook obstruction | `w205.sage hooks 16` | w_m((m,1^{N−m})) = 1 for all 4 ≤ m ≤ 16, conjugate hook 6 … 816504 (`w205_hooks_m1-16.txt`); proof for all m ≥ 5 in note §4 |

## 2. Independent adversarial review (separate AI session, read-only, own code in /tmp)

Verdict: no mathematical or certificate error. Details:
* Prop. 1 re-derived by hand, including (3) for ε = ±1. An independent Sage check of the universality claims for n ≤ 20 agrees.
* (4,4,4,3): independent Sage `kronecker_product` confirms S_15-universality (g_λ = 9, g(λ,λ,(14,1)) = 1). The n ≤ 20 uniqueness of {(3,1,1),(4,4,4,3)} was reproduced.
* `kron_stair.cpp` audit: MN sign, rank/unrank, next_part, z_μ, head/tail uniqueness, Horner, COLLIDE logic and int64 safety are all correct. Re-runs of m = 3..9, 11 reproduce the hashes.
* Independent methods agree with `kron_stair`:
  * Sage `kronecker_product` for all 3718 λ at m = 7;
  * the hook generating function with a separate bitmask MN code, all hooks for m = 7..12 (incl. m = 12);
  * a per-λ odd-class character sum: 25 random λ with a < 0 at m = 9, and 28 λ with a < 0 at m = 11 (the smallest g with a < 0 at m = 11 is 99, at (64,2)).
* §4: w_m((m,1^{N−m})) = 1 is confirmed, and the two-sided hook test passes for all m ≤ 16.
* Required fixes, **all applied**:
  1. cite Bessenrodt–Bowman arXiv:2202.03066 (Thm 6.1, Cor 6.2) for the m ≡ 2 (mod 4) fact and the (n−1,1) count, and stop presenting it as new;
  2. say "the *proof* of Prop. 3.2" and state that w bounds the multiplicity in W_m;
  3. soften "cannot be carried to all m" to "we see no way…", since the obstructions are heuristic.
* Novelty: arXiv:2609.17319 (Sept 2026) still calls the A_n case open. No A_n result at n = 28..78 and no (4,4,4,3)-type example was found. Lee arXiv:2512.15035 was withdrawn.

## 3. Duplicate search (this session, 2026-10-09)

* arXiv/web: HSTZ, Saxl + alternating, A_n tensor square, 2609/2610 listings. Hits: Bessenrodt 2018 (GAP to n ≤ 25), HSTZ 2013 (O'Brien n ≤ 17), Bessenrodt–Bowman 2022/23 (S²/A² refinements, realness; no A_n universality), arXiv:2609.17319 (A_n still open), 2601.01379 (mentions only).
* GitHub: local clones of the openai/math commentary aggregators in /tmp/dedup (eliottcassidy2000, vibemathed, x_satyajitghana_ai, supsup, mahrud, piiiico, Rogerhu12, aster2024, zjukop3, SafeVerifyAgent, blisspixel, jkbennitt, eumemic, mxym, …) grepped for Heide/HSTZ/Saxl/alternating-square. The only hits were bibliography entries and copies of 205 itself. `gh search` repos/code: keithadler/saxl-lean (S_n only, no A_n), d3v1nmb/TensorProductAnalysis (S_n staircase k = 8 data only).
* Conclusion: to our knowledge, HSTZ for A_28, A_36, A_45, A_66, A_78 and the failure of S_n → A_n transfer for (4,4,4,3) are not in the literature or the community repositories.

## 4. Not covered (gaps)

* n = 66, 78 rest on one modular run each (p = 2³¹−1) plus mod-p checksums, not exact integers; no second prime at m = 12.
* m ≥ 13 not attempted (m = 13: p(91) ≈ 7·10⁷ partitions, beyond the memory budget of this implementation).
* The conditional 205 route was independently re-run only for m ≤ 8. The earlier inductive extension to m = 12 (not shipped) was **not** re-run and is not needed for the theorem.
* GAP's A_n tables are trusted for n ≤ 20; Sage/GAP generic character tables are trusted in the cross-checks.
