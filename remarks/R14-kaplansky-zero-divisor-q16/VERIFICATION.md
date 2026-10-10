# R14 verification record (2026-10-09)

## 1. Primary computation (topt-lean session; `scripts/primary/`)
| Item | Command | Result |
|---|---|---|
| q = 16 certificate with the paper's f | `kz196_exact.py` | max (Mf)/f = 7138933/7179627 = 0.994332 < 1, Σf = 288 |
| q = 32 | same | 398926565/405561387 = 0.98364 |
| q = 8 Collatz–Wielandt | same | μ = 1.0117488 > 1 |
| q = 128 control, ρ(M) for q = 2…32 | `kz196_spectral.py` | both paper rationals reproduced; ρ = 1.6238, 1.1825, 1.0117, 0.9750, 0.9784 (independent of the pairing over 40 random pairings) |

Note: the primary builds the weights from the closed-form table and does not use the plane geometry. The `lines` argument of `build_W` is unused.

## 2. Independent re-implementation (adversarial reviewer; `scripts/independent/kz196_indep.py`, `kz196_delta.py`)
This was written from `types.tex` and shares no code with the primary script. The reviewer had read the primary script beforehand.
| Item | Result |
|---|---|
| PG(2,q) from Sage GF(q), q = 2, 4, 8, 16, 32 | q+1 points per line; every two lines meet in exactly 1 point |
| Vertex types at n = 2v² | integral counts; letter frequency exactly p on both sides; odd \|S_x∩S_y\| for **all** type pairs when q ≤ 16 (200 000 random pairs when q = 32) |
| Weights from the geometry (B-side frequencies) | only the values 1/(q+1), p, 1 occur; A-side ≤ B-side checked for every turn |
| Equitable 3-class partition | row sums into each class are constant and equal the quotient Q, checked on the full matrix for every q tested |
| q = 16 certificate | 7138933/7179627 (agrees with primary) |
| q = 8 Collatz–Wielandt | 70905723398671/70082342865504 = 1.011749 (agrees) |
| ρ(M), exact (cubic) vs. numpy | q = 2: 1.623834, 4: 1.182472, 8: 1.011749, 16: 0.975008 (both agree) |
| q = 128 positive control | the class-row formulas equal 500731261911/504183783481 and 820350863/1363395845 exactly |
| δ for q = 16 | 288·λ^{h−1} ≤ e^{−h/300} for h ≥ 2412, so δ = 1/600 stays valid |
| Σ_W P(W)² via the quotient | q = 8: 85, 244, 9.0·10⁶ at h = 10, 100, 1000 (grows); q = 16: 220, 22.5, 2.9·10⁻⁹ (decays) |
| ρ for q = 10, 12, 14 (formula only) | 0.9912, 0.9816, 0.9770 |

## 3. Adversarial review (by hand)
**(a) Does q enter only through parity, the squared-word radius, degree ≥ 3, and r/k > 2? CONFIRMED, with one addition.**
We read every section and `grep`-checked every occurrence of q, 128, 129, 132, 16513, 16516, 8258, 64.5, 62.5, d_* and |T| across the five sections. The q-dependent sites are listed in the note (§1). The addition is that the existence-only constants (c in (2.5), and therefore D₀, d, D, ε) depend on q. These need no check, but we verified that c = 1/100 still satisfies (2.5). Several numerical illustrations specific to q = 128 also change: the counts at n = 2v², 16528 < e^{10}, 149/150, h ≥ 3002, and the range [1/129, 1]. None of them is load-bearing. `topology.tex` does not involve q, and `factors.tex` uses q only through the parity condition (L139).
Specific checks:
- The A-side conditional frequencies are ≤ w (extra–extra: 1/2 ≤ 1).
- The E-row case (t extra) is covered by the manuscript's "t̄ ordinary" case.
- The diameter lemma's (k/n)^{r−2k} needs (q+1)/2 > 2.
- The short closure needs 3·2^{k−1} non-backtracking walks, i.e. degree ≥ 3 after deleting 2 edges.

**(b) Independent certificates: CONFIRMED** (see §2).

**(c) Significance and limits.** We confirm the README statement. The generating set shrinks from 8 258 to 138 and the degrees drop to 17–20. The group is still non-explicit, because n → ∞ is unspecified and there is one relator family per component. Minimality holds within this type design and among plane orders known to exist. We added a nuance not present in SELECTION: ρ < 1 already at q = 10, 12, 14, so what blocks q < 16 is the existence of the plane (order 10 does not exist; 12 and 14 are open), not the spectral condition. The statement "smallest power of two" is correct as written.

**Overall R14 verdict: CONFIRMED.**

## 4. Duplicate search (2026-10-09)
- GitHub: we searched community repos about openai/math (35 repos updated 10/06–10/09) and checked tree listings of aster2024/openai-math-followups, Rogerhu12/openai-math-extensions, piiiico/openai-math-counterexamples, mahrud/openai-math-comments, fortnow/AI-Papers and Kohze/…review. None mentions 196, Kaplansky or the zero-divisor construction. `gh search repos "kaplansky zero divisor"` found nothing relevant.
- arXiv/web: the latest related work is Mian–Siddique 2609.22380 (Gardam's Ã₂ lattice, non-unique-product, no zero divisors), Garg–Mineyev 2501.07646 and Shin 2026. None of them modifies OpenAI's construction. Press summaries of the release give no parameters.
- Conclusion: we found no prior record of the q = 16 reduction.

## GAPs
- No human review, and no Lean (the formalisation is fixed at q = 128).
- We did not re-verify the manuscript's correctness for q = 128; this is the premise.
- Alternative type designs for q = 8 were not explored.
