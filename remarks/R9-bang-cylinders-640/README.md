# R9. Bang's half-area cylinder question: a counterexample with 640 cylinders instead of 16,000,000

**What this remark shows.** OpenAI's family-100 preprint covers a regular tetrahedron by finitely many cylinders whose total perpendicular-base area is below half the minimum projection area. This answers Bang's half-area question negatively. The preprint's explicit parameters need 16,000,000 cylinders. We keep its construction and its parameter-free sector-selection lemma. We drop the global radial estimate, give each angular sector its own cutoff, and certify coverage with 102,400 exact rational linear programs. The result is a cover by **640 cylinders** with Σ|B_i| / A_min(K) ≤ 0.49997733980102355 < 1/2.

## Statement

K = {(x, y, √2·t) : 0 ≤ t ≤ 1, |x| ≤ 1−t, |y| ≤ t} is a regular tetrahedron of edge 2 with A_min(K) = √2. A cylinder is B + ℝu, where u is a unit vector and B ⊂ u^⊥. Bang's question, as recorded by Bezdek–Litvak (J. Geom. Anal. 19 (2009), Introduction) and Bezdek (arXiv:0903.4637, Problem 3.1), asks whether every finite cylinder cover satisfies Σ|B_i| ≥ ½·A_min(K).

| | tilt ε | sectors per family n | cylinders | Σ\|B_i\| / A_min(K) | how coverage is proved |
|---|---|---|---|---|---|
| **Before** (preprint Thm 1.1 / Cor 1.2) | 1/2000 | ⌈2/ε²⌉ = 8,000,000 | 16,000,000 | ½ − (13/6000)ε² + O(ε⁴) ≈ ½ − 5.4·10⁻¹⁰ | uniform margin η = 1/1000; one inequality 2η − (2+2η)ε − ε²/4 > 0 |
| **After** (this remark) | 2/5 | 320 | **640** | ∈ [0.49997733980102354, 0.49997733980102355] | per-sector cutoffs T_j; Z_jk < 0 checked exactly for all 320² sector pairs |

The cylinder count drops by a factor of 25,000 against the preprint's displayed parameters, and by a factor of 1,562.5 against the best count obtainable inside the preprint's own certificate (10⁶, see Related work). The directionwise ratio Σ|B_i| / |π_{u_i^⊥}K| is ≈ 0.4511, so the 640-cylinder cover also violates the directionwise "1-codimensional cylinder covering conjecture" (Bezdek–Khan).

(An internal planning table quoted the old count as 16,000,002. The preprint's own formula 2⌈2/ε²⌉ at ε = 1/2000 gives exactly 16,000,000.)

## Attribution

The tetrahedron, the angular construction (axes v_j = (1, εα_j, √2εβ_j), w_j = (−εα_j, 1, √2εβ_j), intercept triangles), the shared-boundary identity, the first-crossing selection lemma, and the area formula are all **OpenAI's work**. They come from family 100 of [github.com/openai/math](https://github.com/openai/math):

- [`preprints/Finite-angular-cylinder-covers-below-the-half-area-bound-September-27-2026`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/Finite-angular-cylinder-covers-below-the-half-area-bound-September-27-2026) (`build/sections/{introduction,construction,coverage,area}.tex`)
- `lean/ComparatorChallenges/TriangularCovering.lean`, `CylinderCovering.lean`, `lean/docs/100.md`

The paths are pinned at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. These files have not changed since the initial commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.

This remark adds three things: (i) the per-sector cutoffs T_j and the parameters ε = 2/5, n = 320; (ii) the observation that the selection lemma holds for every ε, so that coverage reduces to finitely many exact LPs; (iii) the exact certificate.

## Related work and acknowledgement

An independent parameter-window optimisation within the original certificate gives 10^6 cylinders (supsup/new_math_review, 2026-10-08). Our improvement to 640 replaces the coverage argument itself (per-sector truncations T_j and exact pairwise LPs over 320² sector pairs), which allows ε = 2/5 instead of ε ≤ 1/2000.

Details, which we gratefully acknowledge: [supsup/new_math_review](https://github.com/supsup/new_math_review), commit [`cd6549b2`](https://github.com/supsup/new_math_review/commit/cd6549b24054e6869f480f598f314f263b27e6cf) (2026-10-08T19:11Z, "Add family 100 (cylinder covering of the tetrahedron) card by Lattice"); `showcase.md`, card "Covering a tetrahedron with a million thin cylinders" (≈ L662 in the current version), and `possible_further_findings.md`, item 1. They keep the preprint's two corollaries, the coverage condition η ≥ (2τ+τ²/4)/(2(1−τ)) and the area condition η < 1/480 − τ², and take τ = 1/500, η = 1/498 with 2/n = τ², i.e. 2n = 1,000,000 cylinders (16× fewer than 16,000,000), with exact margins: coverage budget 1751/249000000 > 0 and area saving ≥ 2959/5187500000000 ≈ 5.7·10⁻¹⁰ of A_min. They also note that the window of that certificate first opens near τ ≈ 1/482 (929,752 cylinders), and describe their result as "a parameter choice in the paper's own construction, not a new cover design". Their limit is the preprint's first-order coverage inequality; ours removes that inequality.

## Why it works (what changed)

- `coverage.tex:13-45` (first crossing) and `construction.tex:66-84` (shared boundaries) use only the exact identity α_j − qβ_j = (1−q²)/4 at q = q_j, q_{j+1}, together with φ(±1) = 0. **No smallness of ε is used.** So every point of K lies in some wedge W_j of family 1 and some wedge V_k of family 2.
- The preprint then proves that the two radial cutoffs cannot both fail (`coverage.tex:47-124`). That proof uses the rough-angle bounds |y − q/2|, |x − p/2| ≤ ε and |x² − p²/4| ≤ (3/2)ε, which are first order in ε. Against a gain of only 2η, these produce the condition 2η − (2+2η)ε − ε²/4 > 0, i.e. ε ≤ 1/2000. **This is where the constant is lost.**
- We replace `coverage.tex:47-124` with: for each pair (j,k), Z_jk = max{min(t0 − T_j, s0 − T_k) : p ∈ W_j ∩ V_k}. This is an LP in (x, y, t, z) with 12 constraints. If every Z_jk < 0, the cover is complete (note, Thm 3). The condition is sufficient, not necessary.
- The needed margin varies across sectors by two orders of magnitude, from E_j ≈ 2.4·10⁻⁵ near q = 0 to 2.7·10⁻³ at |q| ≈ 1. With one uniform margin max E_j, the area is 0.50056 > ½. Per-sector cutoffs are therefore what makes ε = 2/5 work.

## Reproduce

Requirements: SageMath with `sage -python` (numpy, scipy, PPL backend). Uses at most 3 processes.

```bash
bash scripts/run_all.sh        # ~12 min wall: params -> 102,400 exact LPs (3 procs) -> summary, area, A_min, 1.1M-point geometric test
```

Each step is `diff`ed against `scripts/expected_output/`. `params_T.txt` (exact cutoffs, md5 `0f4291bd4ee5cd4686db5d688b8d0f1e`) is regenerated from `scripts/primary/persector_0.4_320.npy` and checked. The per-pair LP output `pairs_all.txt` (26 MB) is not shipped. It is regenerated by the script; md5 `64d83b2e619c2840a6a53bfd59df86f1` with SageMath 9.5/PPL on our host. The maximizer points may differ across PPL versions, but the optimal values Z_jk cannot. `scripts/MD5SUMS` lists every shipped file.

Optional: the original certificate via a different method (exact dual vectors) is `cd scripts/primary && sage -python cyl_cert.py 2/5 320 persector_0.4_320.npy` (~19 min; compare with `cert_04_320_ps.log`). Margin generation is `sage -python cyl_persector.py 0.4 320 3 0.002` (~50 min, floating point; its output only defines the cover and is not trusted).

The note is `note.tex` / `note.pdf` (3 pp., `latexmk -pdf note.tex`).

## Scope and limitations (what is NOT claimed)

- **No new construction or idea.** Everything geometric is OpenAI's. The contribution is a parameter choice (per-sector cutoffs, ε = 2/5, n = 320) and a computer-verified certificate.
- **Computer proof, not formalized.** Coverage rests on exact rational LP (PPL) and, independently, on exact dual vectors (`primary/cyl_cert.py`). There is no Lean proof. The OpenAI Lean statement `TriangularCovering.main` is unaffected and was not touched.
- **Not minimal.** We do not claim that 640 is optimal, even within this family. With n = 128 the heuristic margins gave area 0.50045 > ½ (not certified either way). An earlier extrapolation of "≈ 600 minimum" is not claimed.
- **The margin is thin.** The worst exact Z_jk is −2.4·10⁻⁸. Most of this comes from the +10⁻⁷ rounding pad on the margins. That is fine for an exact certificate, but the cover is tight.
- **Regular tetrahedron only.** As in the preprint, the A_min statement concerns the regular tetrahedron. The affine extension (preprint Cor. 1.2) covers the directionwise ratio only.
- **Novelty limited to our search.** To our knowledge (search of GitHub and arXiv up to 2026-10-09), no explicit cover with fewer than 10⁶ cylinders below the half-area bound has been published. The 10⁶ cover of supsup/new_math_review (2026-10-08, see Related work) predates this remark and stays within the preprint's coverage certificate; another third-party treatment of family 100 that we found re-checks the preprint by sampling.
- **Premise.** The two preprint lemmas used (membership and first crossing) were re-proved in the note. Nothing else from the preprint is relied on except A_min(K) = √2, which is proved there and checked numerically here.

## Authorship

[AUTHOR NAME].

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
