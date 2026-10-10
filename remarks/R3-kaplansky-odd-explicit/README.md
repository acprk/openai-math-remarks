# R3: an explicit 25-digit prime for OpenAI's odd-characteristic Kaplansky counterexample

**One sentence.** OpenAI's counterexample to Kaplansky's direct-finiteness conjecture in odd characteristic is stated over a field of order p⁴, where p is the least prime factor of (m!)² + 1 with m = C(1200,600) ≈ 4·10³⁵⁹. That prime is not known and is out of reach of any computation. The same proof goes through with the smaller threshold-Kneser graph m = C(86,41) ≈ 6.05·10²⁴ and the **explicit, certified prime p = 6054252965430758779982953**.

> **The argument is OpenAI's; we only re-select parameters and replace 1 specific inequality (the estimate t/m ≤ 401(5/7)¹⁰⁰ < 1/100) by an exact computation. We list every changed line.**

## Precise statement (Before / After)

| | Before (OpenAI, as written) | After (this remark) |
|---|---|---|
| Graph Γ | 600-subsets of [1200], adjacent iff \|I∩J\| < 200 | 41-subsets of [86], adjacent iff \|I∩J\| < 13 |
| m = #vertices | C(1200,600) ≈ 3.97·10³⁵⁹ | C(86,41) = 6054252965430758779982760 ≈ 6.05·10²⁴ |
| Fit family 𝒯 | \|T\| ≤ 400, t = Σ_{j≤400} C(1200,j) | \|T\| ≤ 28, t = Σ_{j≤28} C(86,j) = 62112284749374537743384 |
| Rank comparison | t/m < 1/100 < 1/96 < d/w | t/m = 0.0102593… < 1/96 < d/w (margin (t/m)/(1/96) = 0.98489) |
| Prime p | least prime factor of (m!)²+1 (not computable in practice) | least prime > m with p ≡ 1 (mod 4): **6054252965430758779982953** = m + 193 |
| Field K | order p⁴ | order p⁴, now with an explicit p |

Theorem (after): there are an explicit field K of order p⁴ with p = 6054252965430758779982953, a finitely generated group G containing torsion, and finite sums a, b ∈ K[G] with ab = 1 and ba ≠ 1. The proof is OpenAI's, with the parameter changes below. As in the original, this does not make the group small or the elements practically computable: the local group H still has order (p⁴−1)^{p³}·p³.

Within the threshold-Kneser mechanism, (86,41,13) is the smallest choice. `minN.py` searches every (N, k, s) with N ≤ 120 such that 3k − 3(s−1) > N and t/m < 1/96, and finds no smaller m. For N > 120 no instance can be smaller: triangle-freeness gives k − s ≥ ⌊N/3⌋, so t ≥ C(N, ⌊N/3⌋), and t/m < 1/96 forces m > 96·C(N, ⌊N/3⌋) ≥ 96·C(121, 40) > 10³⁴.

## Attribution

The construction and proof are OpenAI's work: family 197, preprint `preprints/A-Counterexample-to-Kaplanskys-Direct-Finiteness-Conjecture-in-Odd-Characteristic-September-26-2026/` in https://github.com/openai/math (local snapshot at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`; line numbers refer to `build/sections/*.tex` there). The graph and fitting polynomial are the threshold-Kneser specialization of Haviv and Golovnev–Haviv that the paper itself cites. This remark only changes the parameters (N, k, s) and the definition of p.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Every changed line

| File : lines | Original | Replacement |
|---|---|---|
| `01-introduction.tex` L16–17 | m = C(1200,600), p = min{h > 1 : h \| (m!)²+1} | m = C(86,41), p = least prime > m with p ≡ 1 (mod 4) = 6054252965430758779982953 |
| `01-introduction.tex` L66 | t/m < 1/100 | t/m < 1/96 |
| `02-finite-data.tex` L8–9 (statement of Lemma op:field) | unchanged | unchanged (p > m odd prime, −1 a square) |
| `02-finite-data.tex` L16–20 (its proof) | least-divisor argument and (m!)² = −1 | p > m by definition; −1 is a square because p ≡ 1 (mod 4) |
| `02-finite-data.tex` L44–45 | 600-subsets of [1200]; \|I_i∩I_j\| < 200 | 41-subsets of [86]; \|I_i∩I_j\| < 13 |
| `02-finite-data.tex` L47–48 | \|T\| ≤ 400; t = Σ_{j=0}^{400} C(1200,j) | \|T\| ≤ 28; t = Σ_{j=0}^{28} C(86,j) |
| `02-finite-data.tex` L54 | t/m < 1/100 | t/m < 1/96 |
| `02-finite-data.tex` L57–58 | 3·600 − 3·199 = 1203 > 1200 | 3·41 − 3·12 = 87 > 86 |
| `02-finite-data.tex` L63–64, L67 | ∏_{k=200}^{599} (ΣX − k)/(600−k); denominators in {1,…,400} | ∏_{k=13}^{40} (ΣX − k)/(41−k); denominators in {1,…,28} |
| `02-finite-data.tex` L71, L73–74 | ∏_{k=200}^{599}; nonneighbour intersections in [200, 599] | ∏_{k=13}^{40}; nonneighbour intersections in [13, 40] |
| `02-finite-data.tex` L79–87 | t/m ≤ 401(5/7)¹⁰⁰ < 1/100 | **exact**: 96t = 5962779335939955623364864 < m = 6054252965430758779982760 |
| `03-modular.tex` L4 | t/m < 1/100 | t/m < 1/96 |
| `05-compression.tex` L93 | t/m < 1/100 < 1/96 < d/w | t/m < 1/96 < d/w |

Lines that use p or m but need **no** change: `03-modular.tex` L58–61 (h/p³ > 1/48 for p ≥ 3); L139–146 (uses p > 600, still true); L207–209 (d/w > h/(2p³) > 1/96, needs only p ≥ 5); `04-gluing.tex` L391 ("residues are distinct since p > m", still true by definition).

The Lean formalization in OpenAI's repository (`lean/ComparatorChallenges/OddKaplansky.lean` and `lean/OAI/Algebra/OddKaplansky/*`) hard-codes the original parameters. We did **not** redo it.

### Variant that keeps the paper's constant 1/100
(N, k, s) = (86, 42, 14) gives m = C(86,42) = 6486699605818670121410100 and t/m = 0.009575 < 1/100. The least prime p > m with p ≡ 1 (mod 4) is 6486699605818670121410413 (= m + 313; APRCL-prime). With this choice, the lines `01` L66, `02` L54, `03` L4 and `05` L93 stay unchanged. This variant was computed but not separately adversarially reviewed. We record it only as an option.

## Reproduce

From `scripts/` (Python 3.10 with sympy; PARI/GP 2.13):

```
/usr/bin/python3 main/tk_exact.py            # m, t, t/m, p, fit identities         -> expected_output/tk_exact.txt
/usr/bin/python3 main/tk_opt.py              # first admissible (N,k,s) per constant -> expected_output/tk_opt.txt
/usr/bin/python3 main/alt_keep_1over100.py   # the (86,42,14) variant               -> expected_output/alt_keep_1over100.txt
/usr/bin/python3 independent/check197.py     # exact Stirling multilinearisation + explicit (10,4,2) expansion -> expected_output/check197.txt
/usr/bin/python3 independent/minN.py         # exhaustive minimality, ~75 s         -> expected_output/minN.txt
gp -q independent/pcheck.gp < /dev/null      # APRCL + Pocklington certificate, side conditions -> expected_output/pcheck_gp.txt
```

## What is NOT claimed

* We did not re-verify the correctness of OpenAI's proof. The remark inherits it. The adversarial review checked that no step uses the size of m or p, or the special form of p, beyond the lines listed.
* "Explicit" refers to the prime and the field. The group G and the elements a, b remain astronomically large, exactly as in the original.
* The Lean certificate was not redone. It would need the new parameters and a Pratt or Pocklington certificate for the 25-digit prime.
* A smaller instance, (N, k, s) = (74, 35, 11) with p = 1569396701548619508553 (22 digits), would follow from replacing the crude bound \|𝓡\| ≥ \|A₀\|/(2p³) (`03-modular.tex` L137–154) by the true proportion ≈ 1/(48p). That needs a new estimate and is **not** claimed here.
* Characteristic two is not addressed (that is a separate OpenAI preprint).
* Novelty: the earlier search session searched GitHub and arXiv (up to 2026-10-09) and found no explicit prime for this construction. The only related community item found, formal-conjectures #6919, concerns Gottschalk's conjecture. To our knowledge, no one has given one.
