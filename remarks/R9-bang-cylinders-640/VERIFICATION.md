# Verification record: R9 (640-cylinder cover below Bang's half-area bound, OpenAI family 100)

Date: 2026-10-09. Host: Linux x86_64, SageMath 9.5 (`sage -python`, PPL backend), at most 3 processes. No Lean was run.

Pins: openai/math `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The family-100 preprint directory and `TriangularCovering.lean` were last changed in `adc7f124`. The local checkout had no modifications to them (`git status` clean).

## Who did what

| Role | Session | What |
|---|---|---|
| Primary computation | R7-skill rescan s1 (`reviews/R7-skill/rescan/s1/`) | Chose per-sector margins (`cyl_persector.py`, floating point). Wrote the certificate `cyl_cert.py`: exact rational dual vectors, one per (j,k), then the area via rational bounds r_j ≥ 1/√(·). Log `cert_04_320_ps.log`: max U_jk = −1.5·10⁻⁷, 0 uncertified, S/A_min ≤ 0.4999773398. |
| Adversarial review + independent reimplementation | this session (separate AI session) | Re-derived the geometry from the tex. Wrote new code (`scripts/independent/`) that reuses only the margin data file. Its exact primal LP, interval area, direct geometric test and A_min check are described below. |

## (a) Logic of the reduction (adversarial reading)

| # | Question | Finding |
|---|---|---|
| A1 | Does the first-crossing lemma (`coverage.tex:13-45`) hold for every ε? | **Yes.** It uses only F_0 = −t, F_n = t (from φ(±1) = 0) and the exact identity α_j − q_iβ_j = (1−q_i²)/4 (i = j, j+1; `construction.tex:80-84`). This identity makes "F_j ≤ y ≤ F_{j+1}" literally equivalent to q_j t0 ≤ y0 ≤ q_{j+1} t0. Monotonicity of F_i is not used, as the text says at line 34. The second family is the same with G_i. t0 ≥ 0 follows from Δ > 0. Re-proved in the note (Lemma 2). |
| A2 | Is "Z_jk < 0 for all (j,k)" a valid substitute for the preprint's coverage condition? | **Sufficient, and that is all that is needed.** Given p ∈ K, A1 gives (j,k) with p ∈ W_j ∩ V_k. Then Z_jk < 0 forces t0 < T_j or s0 < T_k, and membership (Lemma 1: unique intercept, since v_j has x-component 1 and w_k has y-component 1) gives coverage. The condition is **not necessary**: in a control run with the margins removed, 222 pairs had Z > 0, yet all 222 maximizer points were covered by some other cylinder. So the LP is conservative, and its failure would not disprove coverage. |
| A3 | Are there hidden "ε small" uses in the parts we still rely on? | We rely on: the shared-boundary identity (exact algebra), first crossing (A1), the intercept/membership formula (exact for any ε), the area formula |B_j| = |P_j|/|v_j| (the normal e1 of the intercept plane has e1·v_j = 1 for any ε), and A_min(K) = √2 (independent of ε). Everything the preprint does with ε small is in `coverage.tex:47-124`: the bound t0 + s0 ≤ 1 + ε/2, rough angles, linearization of x², the Δ/(4ε) term, and the requirement ε ≤ η/2. The LP certificate **replaces** all of these, so none is used. Cutoffs T_j ∈ [0.5000042, 0.5204]: these are positive, so the triangles are nondegenerate. Intercept triangles may stick out of K, which is harmless. q ranges over [−1,1], and the partition is the preprint's. |
| A4 | Same cutoff for both families? | Yes. The construction is symmetric under x↔y, t↔1−t. `cyl_cert.py` uses the same T_j for the failure rows of both families and for the area (factor 2), and the independent code does the same. |
| A5 | Is the primary's encoding of the LP correct? | Checked row by row against the geometry: the K rows, the wedge rows in F/G form, and the failure rows z·ε² ≤ t0 − ½ − ε²(d_j+E_j), likewise for s0. Correct. Its "infeasible" branch can never certify (it would need y ≥ 0 with Σy = −1), but it is never reached because W_j ∩ V_k ≠ ∅ for all pairs. That is harmless, and the log reports 0 uncertified. |

## (b) Independent computations

| Check | Command | Result |
|---|---|---|
| Exact cutoffs | `sage -python make_params.py ../primary/persector_0.4_320.npy` | 320 rationals, same rounding rule as `cyl_cert.py`; md5 `0f4291bd…`. E_j ∈ [2.39·10⁻⁵, 2.67·10⁻³]. |
| **Exact primal LP, all pairs** | `pair_lp_exact.py` (PPL, exact ℚ optimum; LP written directly in intercept form t0, y0, s0, x0, not the F/G form) | 102,400/102,400 pairs; **0 with Z_jk ≥ 0**; exact max Z = −29290260249784175320141/1220427626492079325376000000000 ≈ **−2.39999977·10⁻⁸** at (78,84), (84,241), (235,78), (241,235). Divided by ε² this is −1.4999999·10⁻⁷, which equals the primary's max U_jk (the dual bound is tight). ≈ 6 min with 3 processes. |
| **Rigorous area** | `area_interval.py` (256-bit RealIntervalField) | S/A_min ∈ [0.499977339801023544…, 0.499977339801023544…], upper < ½; ½ − S ≥ 2.2660·10⁻⁵. Formula sanity: ε = 0, T = ½ gives exactly ½ (Bang's cover). Agrees with the primary's 0.4999773398. |
| **Direct geometric test** | `geom_sample.py params_T.txt 20261009 pairs_all.txt` | 1,102,405 points with exact rational coordinates: 600k uniform in K, 250k in \|t−½\| ≤ 1/20, 100k on the 4 faces plus the t = 0, 1 edges, 50k on the 4 slanted edges and within 10⁻³ of the vertices, the 4 vertices and the centre, and all 102,400 exact LP maximizers ("hardest" points). For each point, membership in one of the 640 cylinders is verified exactly (floats are used only to order candidates). **Uncovered: 0.** Wall time 5 min. |
| Negative control (detector works) | `R9_SCALE=0.04 … geom_sample.py params_control_half.txt` (same axes, all T_j = ½) | 959 of 40,005 points uncovered. The test is not vacuous. |
| Negative control (LP works) | `pair_lp_exact.py params_control_nomargin.txt` for rows 0, 78, 160 | 222 pairs with Z > 0, max 4.07·10⁻⁴. The LP detects removed margins. See A2 for why these points are nevertheless covered. |
| A_min and directionwise ratio | `amin_check.py` | Grid of 2·10⁵ directions plus Nelder–Mead: min shadow 1.414213562373 = √2. Shadows in the 640 axis directions lie in [1.5479, 1.6033]. Directionwise relative cost ≈ 0.451085. |
| Uniform-margin comparison | `area_controls.py` | Same axes with uniform margin max E_j: area 0.5005642 > ½. Without margins: 0.49972, but coverage is not certified. |

The independent code shares **no code** with `primary/`. The only shared input is the float margin file `persector_0.4_320.npy`, which defines the cover and is not trusted. Shared assumptions: the reading of `construction.tex:111-123` (cylinder definitions) and the first-crossing lemma, which was re-proved.

## (c) Problem statement

| # | Question | Finding |
|---|---|---|
| C1 | What exactly did Bang ask? | Bezdek–Litvak, arXiv:0903.4842v1, Introduction (text extracted and read): "Bang raises also the important related question whether the sum of the base areas of finitely many cylinders covering a 3-dimensional convex body is at least half of the minimum area of a 2-dimensional projection of the body … sharp due to a covering of a regular tetrahedron by two cylinders." Their §3 definition: a cylinder is C = ℓ + B, ℓ a line through 0, B measurable in ℓ^⊥, base area vol₂(B). |
| C2 | Do the remark's bases match "base area"? | Yes. Each cylinder is written as π_{v^⊥}P_j + ℝv with base π_{v^⊥}P_j ⊂ v^⊥ (perpendicular), and |B_j| = |P_j|·(e1·v̂_j) = |P_j|/|v_j|. These are the same objects as in the preprint's `construction.tex:124-128` and the Lean `TriCylinder` (unit axis, base ⊥ axis). |
| C3 | Normalization | A_min(K) = min over planes of the orthogonal-projection area = √2 (preprint Lemma; numerically confirmed). The needed inequality is Σ|B_i| < ½·A_min(K), i.e. S/A_min < ½, which is exactly what is certified. The width/plank version (Bang's theorem proper) is a different statement and is not involved. |
| C4 | Old count | The preprint's explicit parameters are ε ≤ 1/2000 (Thm 1.1, Cor 1.2), m = 2⌈2/ε²⌉ = 16,000,000. Other constructions in the preprint (endpoint sectors, cubic buffer with τ ≤ 1/4000, appendix table) are stated "for sufficiently small" parameters, or have larger counts. The "16,000,002" in SELECTION.md has no source in the preprint. |

## (d) Prior third-party improvement (dedup, added 2026-10-09)

| # | Question | Finding |
|---|---|---|
| D1 | Source and date | [supsup/new_math_review](https://github.com/supsup/new_math_review), commit `cd6549b24054e6869f480f598f314f263b27e6cf`, author/committer date 2026-10-08T19:11:59Z, message "Add family 100 (cylinder covering of the tetrahedron) card by Lattice" (verified via `gh api repos/supsup/new_math_review/commits/cd6549b2`; a shallow clone at `97c8b75` (2026-10-09) contains the same card at `showcase.md` L662–686 and `possible_further_findings.md` item 1). |
| D2 | What they claim | Within the preprint's two corollaries — coverage η ≥ (2τ+τ²/4)/(2(1−τ)) and area η < 1/480 − τ² — the choice τ = 1/500, η = 1/498, 2/n = τ² gives 2n = 1,000,000 cylinders (16× fewer than 16,000,000). Exact margins stated: coverage budget 1751/249000000, area saving ≥ 2959/5187500000000. Window edge: 1/N tilts open at N = 483 (933,156 cylinders); rational tilts just below 1/482.1 give 929,752. Self-described as "a parameter choice in the paper's own construction, not a new cover design". |
| D3 | Recheck | Exact arithmetic here: 2η − (2+2η)τ − τ²/4 = 1751/249000000 at (τ,η) = (1/500,1/498), and 1/480 − τ² − η = 2959/41500000; times the preprint's area factor 2τ² this is 2959/5187500000000. Both agree with their numbers. |
| D4 | Relation to R9 | Their result stays inside the preprint's coverage certificate, whose first-order inequality forces τ ≲ 1/482 (> 9·10⁵ cylinders). R9 replaces that certificate (per-sector T_j + exact LPs over 320² pairs), so ε = 2/5 is admissible. R9 does not use their result. Improvement over the best in-certificate count: 10⁶/640 = 1,562.5×. README, note §1 and the note's bibliography now credit them. The earlier README sentence "no explicit cover with fewer cylinders ... has been published" was corrected to "fewer than 10⁶". |

## Verdict

**CONFIRMED.** The 640-cylinder cover (ε = 2/5, n = 320, cutoffs `params_T.txt`) covers K, and its total perpendicular-base area satisfies S/A_min ≤ 0.49997733980102355 < 1/2. This was established by two exact computations by different methods (dual vectors in the primary; exact primal LP optimum here) that agree to all printed digits, plus an exact 1.1M-point direct membership test and rigorous interval area.

## GAPs (not covered)

- No Lean formalization of the 640-cylinder certificate.
- Trust in PPL/Sage exact arithmetic and in Python `fractions`. Two independent exact pipelines mitigate this.
- `cyl_cert.py` was not re-run in this session (19 min). Its log is shipped and its numbers coincide with the new exact LP.
- The n = 128 comparison (0.50045) and the "≈600 lower extrapolation" come from the primary session and were not re-verified.
- The novelty search was not redone in this session (the primary's search: GitHub and arXiv, up to 2026-10-09). A later dedup pass found supsup/new_math_review (§(d)); its 10⁶-cylinder cover predates ours, is weaker (10⁶ vs 640 cylinders), and is now credited.
