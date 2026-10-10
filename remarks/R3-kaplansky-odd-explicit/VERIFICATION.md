# R3 verification record

Claim: OpenAI's odd-characteristic Kaplansky counterexample (family 197) goes through with (N, k, s) = (86, 41, 13), m = C(86,41), and the explicit prime p = 6054252965430758779982953 (the least prime > m with p ≡ 1 mod 4). The required changes are exactly the lines listed in the README.

Environment: Python 3.10.12 (`/usr/bin/python3`, sympy 1.14.0) and PARI/GP 2.13.3. All arithmetic is exact (`int`, `Fraction`, GP integers). Runs on 2026-10-09.

## Layer 1: main computation (primary session; re-run by the coordinating agent and again here)

* `scripts/main/tk_exact.py` computes m, t, t/m = 0.0102593 < 1/96, and p (25 digits). It checks the fit identities F(a) ∈ {0, 1} for a = 0..41 and the maximum degree Δ.
* `scripts/main/tk_opt.py` returns the first admissible (N, k, s) for each constant: (86,42,14) for 1/100, (86,41,13) for 1/96, (74,35,11) for 1/48. The original (1200,600,200) has m ≈ 10^359.6.
* `scripts/main/alt_keep_1over100.py` covers the (86,42,14) variant: p = m + 313 = 6486699605818670121410413.

## Layer 2: independent re-implementation (separate session, `scripts/independent/`)

* `check197.py` computes the multilinear coefficients c_r of Q(Σ X) through Stirling numbers, in exact rationals. It checks F(x) = Σ_r C(x,r) c_r against ∏(x−k)/(41−k) for every x = 0..41 (diagonal = 1, non-edges = 0). It also checks that the denominators lie in [1, 28], that p > m, p ≡ 1 mod 4, and that no c_r has a denominator divisible by p. As an end-to-end test of the multilinearisation, it expands Q symbolically for the toy instance (10,4,2) and gets 0 errors over 210 vertices. Output: `expected_output/check197.txt`.
* `minN.py` searches exhaustively over N ≤ 120, all k, s with 3k−3(s−1) > N. The least m with t/m < 1/96 is (86,41,13); with 1/48 it is (74,35,11). The README gives the argument that N > 120 cannot do better.
* `pcheck.gp` runs APRCL `isprime(p,2) = 1`. It factors p−1 = 2³·3·19·60101·220909311629180017, proves each factor prime with APRCL, and finds Pocklington witnesses [[2,5],[3,2],[19,2],[60101,2],[220909311629180017,2]]. It confirms that p is the least prime > m with p ≡ 1 mod 4, and that p ∤ (m!)²+1 (so the definition of p must change). It also checks the side conditions h/p³ > 1/48, h/M < 1/4, (p³−1)/M < 1/4 and p > 600, and gives (t/m)/(1/96) = 0.98489. The original script called `primecertify`, which GP 2.13 lacks, so the certificate is now built explicitly.

## Layer 3: adversarial review (independent session, 2026-10-09). Verdict: **CONFIRMED**

* All seven sections (`01`–`06`, `04a`) were read in full; `build/` has no appendix. With the new parameters (86, 41, 13) and p = m + 193, every step holds. The review checked for hidden conditions (vertex transitivity, regularity, maximum degree, Johnson-scheme eigenvalues, p ∤ binomials) and found that none is used anywhere.
* The review found that the first "five places" list was incomplete. It added: `01` L16–17 (the definition of p), `01` L66 and `03` L4 (1/100 → 1/96), `02` L16–20 (−1 is a square via p ≡ 1 mod 4), and `02` L79–87 (the estimate replaced by an exact value). None of these breaks the proof. The README table is the completed list.
* d/w > 1/96 depends only on p ≥ 5 (`03` L207–209). The margin is (t/m)/(1/96) = 0.98489.
* This session also re-checked by grep that `m!`, `1200`, `600`, `400`, `200` and `1/100` occur only on the listed lines.

## GAPs

1. The correctness of OpenAI's proof itself (the gluing, the sector construction and the compression) was not re-verified. The remark is conditional on it.
2. The Lean formalization (`lean/ComparatorChallenges/OddKaplansky.lean`, `lean/OAI/Algebra/OddKaplansky/*`) hard-codes the original source parameters. It was not redone. Redoing it needs the new parameters and a Lean-checkable certificate for p.
3. The (86,42,14) variant and the (74,35,11) direction were not adversarially reviewed. They are not claimed.
