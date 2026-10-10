# R14: family 196, the Kaplansky zero-divisor construction works with q = 16

**In short.** OpenAI builds its torsion-free group with zero divisors from PG(2,128): 16 513 ordinary letters, an alphabet of 16 516 letters, and a rose with 8 258 edges. **The same construction and proof work with q = 16**: 273 points, 276 letters, and a rose with **138** edges. For q = 8 the key lemma is false (spectral radius 1.0117 > 1), so 16 is the smallest power of two for which this construction works.

## Attribution
The construction, theorem and proof are OpenAI's work: family 196, *A Torsion-Free Group Algebra with Zero Divisors* (September 23, 2026), [openai/math `preprints/A-Torsion-Free-Group-Algebra-with-Zero-Divisors-September-23-2026`](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-Torsion-Free-Group-Algebra-with-Zero-Divisors-September-23-2026), commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. This remark audits every place where q enters, gives exact certificates for q = 16 (works) and q = 8 (fails), and lists the values that change.

## Before / After
| Item (`types.tex` line) | Manuscript (q = 128) | This remark (q = 16) |
|---|---|---|
| v, \|T\|, rose edges (L12–18) | 16 513, 16 516, 8 258 | 273, 276, **138** |
| degrees, d_* (L45–46) | 129…132, 132 | 17…20, 20 |
| counts per line class at n = 2v² (L26) | A: 129 per double-extra type / 32 639 empty; B: 258 / 32 768 | A: 17 / 495; B: 34 / 512 |
| Lemma 2.1 row bounds (L83, L89) | 500731261911/504183783481, 820350863/1363395845 < 149/150 | **7138933/7179627 = 0.99433**, 77141/124215 |
| Σf, explicit δ (L93–100) | 16 528; δ = 1/600 for h ≥ 3002 | 288; δ = 1/600 for h ≥ 2412 |
| exact ρ(M) | 0.99280 | 0.97501 |
| expansion exponent r/k − 2 (L167) | 62.5 | 6.5 |
| closure degree (q+1)−2 (`planar.tex` L63) | 127 | 15 |
| c = 1/100 in (2.5) (L111) | valid | valid (0.091 < 1, 0.060 < 1) |

**q = 8 fails:** ρ(M) = 1.011749. The Collatz–Wielandt certificate is min (Mg)/g = 70905723398671/70082342865504 > 1, so Σ_W P(W)² does not decay. q = 2 and q = 4 also fail (ρ = 1.62 and 1.18).

## Why only these places matter
After a line-by-line read of all five sections, q enters only through:
1. **parity**: the line intersections (1 or q+1) and the extra parts (0 or 2) give odd |S_x∩S_y| if and only if q is even;
2. **squared words**: ρ(w²) < 1 (Lemma 2.1). The entropy section uses only δ > 0, the fixed |T| and the least weight (a₀);
3. **degree ≥ 3** in the short-closure lemma, and r/k > 2 in the expansion lemma;
4. **existence-only constants** (c, D₀, d, D, ε).

Letters fall into three classes: O₀ (t and t̄ ordinary), O₁ (t̄ extra) and E (t extra). This partition is equitable for M, so ρ(M) is the Perron root of an explicit 3×3 matrix. The manuscript's weight f is constant on the classes, so its two row bounds are exact row values.

## Reproduce
Requirements: SageMath 9.5 (Python 3, numpy). Each run is single-threaded.
```
cd scripts
sage -python independent/kz196_indep.py     # PG(2,q) from GF(q), empirical weights, certificates, q=128 control (~30 s)
sage -python independent/kz196_delta.py     # explicit delta for q=16; growth/decay of sum P(W)^2           (~3 s)
cd primary
sage -python kz196_exact.py                 # primary certificates (q=16, 32 and Collatz-Wielandt q=8)       (~10 s)
sage -python kz196_spectral.py              # primary: q=128 control rationals, rho(M) for q=2..32             (~1 min)
```
Expected outputs are in `scripts/expected_output/`.

## What is NOT claimed
- **Not an explicit group.** G still comes from random matchings with n → ∞ unspecified, with one relator family per graph component. A smaller alphabet does not make the counterexample explicit or computer-checkable, and it does not bound the number of relators.
- **Minimality only within this type design.** For q = 8 we show that *this* choice of types and extra letters fails. Other designs (more extras, different extra distributions) were not explored. We call q = 16 minimal among powers of two, and among projective-plane orders known to exist. The formula gives ρ < 1 for q = 10, 12 and 14, but no plane of order 10 exists, and planes of order 12 and 14 are open problems.
- No new method. The Lean formalisation hard-codes q = 128 (at least 7 files under `lean/OAI/Algebra/GroupRing/`, e.g. `TypedDiameter.lean`, `WordEstimate.lean`) and was not changed. Following the project's rule, no Lean was written.

## Premise
This remark assumes the manuscript is correct for q = 128, in particular the bounded-pattern estimate, the planar extraction and the cone-picture topology. Our check covers only the q-dependence of these arguments.

## Prior work
To our knowledge (search of GitHub and arXiv up to 2026-10-09; see VERIFICATION.md §4), this reduction of the alphabet has not been recorded elsewhere.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
