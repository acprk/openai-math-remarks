# R11. Zariski multiplicity counterexample (multiplicities 2 and 3): ambient dimension ≈ 4.6·10^143 → 1,802,976

**What this remark shows.** OpenAI's family-059 preprint builds two reduced, isolated, weighted-homogeneous hypersurface germs in ℂ^N. They are ambiently homeomorphic but have multiplicities 2 and 3. N is fixed by an integer "product relation" (Lemma `seed:product`) through N = 4s₀ + 16. The relation shipped with the preprint has 140-digit coefficients, which gives N ≈ 4.58·10^143. We feed the **same construction** a much shorter relation, found by lattice reduction. Its coefficients satisfy |e_b| ≤ 29 and Σe_b = 206. The result is **N = 1,802,976** (= 8 · 225,372). The preprint's own verifier accepts this relation.

**The argument is OpenAI's; we only find a shorter relation.** Everything below assumes that the preprint is correct (*premise: the original manuscript is correct*). We did not re-check its lattice-theoretic or topological parts.

## Statement

Let d = 3^11, ℬ = {b : 0 < b < d/2, 3 ∤ b}, A_b = (d−b)/b, R_b = (2d+b)(2d−b)/((3d−b)(d+b)). Lemma `seed:product` (seeds.tex:7-20) asks for integers e_b with E = Σe_b > 0 and ∏A_b^{e_b} = ∏R_b^{e_b} = 1.

| | relation | support | max \|e_b\| | E = Σe_b | Σ\|e_b\| | s₀ | N = 4s₀+16 |
|---|---|---|---|---|---|---|---|
| **Before** (preprint's supplied `explicit-relation.json`) | adjugate-type kernel vector | 849 | 140 digits | ≈ 3.10·10^139 | ≈ 1.89·10^141 | ≈ 1.14·10^143 | **≈ 4.576·10^143** |
| **After** (this remark, `scripts/data/relation_d3pow11_N1802976.json`) | short vector from weighted LLL | 1160 | 29 | 206 | 7414 | 450,740 | **1,802,976** |

The ratio is ≈ 2.5·10^137. Hessian ranks become 16E = 3296 for f₁ and 0 for f₂, so the coranks are N−3296 and N (Cor. `topo:arnold-corank`). The preprint never prints a number for N. Theorem 1.1 only asserts that a finite N divisible by 8 exists, and introduction.tex:54-56 says the construction "is not intended to minimize it". The "Before" value is what the preprint's own verifier reports for the supplied relation. The rank-certificate/adjugate relation in the proof (seeds.tex:58-70) is never computed explicitly. Its coefficients are of comparable size.

## Attribution

The construction and every proof step are **OpenAI's work**. They come from family 059 of [github.com/openai/math](https://github.com/openai/math):

- [`preprints/Ambiently-homeomorphic-isolated-hypersurfaces-of-multiplicities-two-and-three-September-24-2026`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/Ambiently-homeomorphic-isolated-hypersurfaces-of-multiplicities-two-and-three-September-24-2026) (`build/sections/{seeds,local,global,conclusion,certificate}.tex`, `verification/`)

The paths are pinned at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The directory has not changed since the initial commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. `scripts/primary/openai_explicit_relation.py` and `scripts/data/openai_explicit-relation.json` are verbatim copies (Apache-2.0, © OpenAI). The only change is a two-line provenance header added to the .py file.

This remark adds one thing: a new integer vector (e_b) satisfying Lemma `seed:product`, together with the search that found it.

## Why the substitution is legitimate (what the proof uses from the relation)

We read the whole preprint. The relation enters only through Lemma `seed:product`, which is an existence statement. We list every property of (e_b) that is used afterwards:

| Where | Property used | Holds for the new relation? |
|---|---|---|
| seeds.tex:131-146 (weight table) | e_b ∈ ℤ, b ∈ ℬ (so the weights b/d, b/3d, (d∓b)/3d lie in (0,½) with the stated reduced denominators) | yes (checked) |
| seeds.tex:156-166 | E > 0, so 4E > 0 balancing squares lie on side 1 and give ord₀p₁ = 2. All signed counts 2e_b are even, so s₀ is even | yes, E = 206 |
| seeds.tex:197-212 (Case 1) | ∏A_b^{2e_b} = 1 and ∏(R_b/A_b)^{2e_b} = 1 | yes (both products are exactly 1) |
| seeds.tex:214-248 (Case 2) | the triple-angle identity, row by row. It holds for every integer e_b | independent of e |
| seeds.tex:250-258, 269-273 | two common cubes; monodromy order divides D = 3^12 | independent of e (same d) |
| local.tex, global.tex | only Cor. `seed:prepared` (equal real Seifert forms, equal integral monodromy modules, odd order dividing 3^12·5·7, s even) | inherited |
| conclusion.tex:19-24, 284-287, 391-394 | N = 4s₀+16 ≡ 0 (mod 8); Hessian rank 16E with E > 0 | yes |

The adjugate prescription (seeds.tex:58-76) only proves existence and makes the construction finite. Nothing downstream uses the specific shape of that vector. The preprint itself already treats a *different* witness (its 849-row explicit relation) as an acceptable alternative (certificate.tex:11-20).

## How the relation was found

We keep the columns b ∈ ℬ for which all six of d−b, b, 2d±b, 3d−b, d+b are 6000-smooth (the preprint uses 4450). We prune singly-supported columns iteratively, as the preprint does, which leaves 1,405 columns and 1,150 nonzero valuation rows. We take the ℚ-kernel (dimension 255) and saturate it to the integer kernel lattice. On that lattice we run LLL with column weights c_b = 4 + 2Σ_{c∈{3b,b,d−b,d+b}}(ℓ(c)−1), where ℓ is the chain length, and then search short combinations. The exact identity

  N = 4·Σ_b |e_b|·c_b + 8E + 24

makes the weighted ℓ¹-norm the right target. LLL is used here only as a heuristic. Nothing about it is used in the proof.

## Reproduce

```bash
bash scripts/run_all.sh     # Python 3 stdlib only, single process, ≈ 15 s
```

The script checks `SHA256SUMS`. It then runs the preprint's own verifier on the new relation and, as a positive control, on the original relation. Next it runs our independent checker (`independent/indep_check_059.py`, which shares no code with the verifier) with a sampled mod-2 spectral test. Finally it runs two negative controls (perturbed relations must be rejected). Every output is diffed against `scripts/expected_output/`.

Optional search (SageMath 9.5, ≈ 55 min, 1 core): `cd scripts/primary && python3 prune.py 11 6000 && sage ktime.sage 11 6000 && sage fromK.sage 11 6000`. This produces `relK_11_6000.json`, identical in content to `data/relation_d3pow11_N1802976.json` (sha256 `e3b747b1…`); see `kt_11_6000.log`.

The note is `note.tex` / `note.pdf` (`latexmk -pdf note.tex`).

## Scope and limitations (what is NOT claimed)

- **No new mathematics in the topology or lattice theory.** The whole proof is OpenAI's. We supply a different witness for one finite lemma.
- **Premise: the original manuscript is correct.** The checks here cover only the finite arithmetic (the relation, the chains, the variable counts, and a sampled check of the mod-2 spectral equality). The preprint's classification inputs (Steenbrink/Némethi, Orlik–Randell/Hertling–Mase, strong approximation, Levine) and its local–global argument were read, not re-proved.
- **Not minimal.** 1,802,976 is an upper bound obtained with this construction. We did not optimize over d, the smoothness bound, or the weight table, and the weight table does not cancel common chain tails across rows. Partial runs found no kernel vector with E ≠ 0 for d ≤ 3^8, nor for d = 3^9 with bound 10^5. Larger cases were not completed, so we make **no** lower-bound claim.
- **Not a different counterexample.** The polynomials f₁, f₂ have the same shape as in the preprint (Thom–Sebastiani sums of chains, squares, cubes, x⁵+y⁷, six cubes, x¹¹+y¹³). Only the number of copies changes. They remain far from the preprint's separate ℂ⁴ multiplicity-4/5 example.
- **Novelty limited to our search.** To our knowledge (search of GitHub and arXiv up to 2026-10-09), no smaller explicit dimension for this multiplicity-2/3 pair has been published. The only third-party computation we found (`piiiico/openai-math-counterexamples`, row 059) re-verifies the original 140-digit relation.

## Authorship

[AUTHOR NAME].

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
