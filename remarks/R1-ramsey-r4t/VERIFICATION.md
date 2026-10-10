# Verification record: R1 (r(4,t) from OpenAI family 170)

Date: 2026-10-09. Host: Linux x86_64. Toolchain: `leanprover/lean4:v4.34.1`. Mathlib: `d13f23b723b8a846827a245b89c10fc7d3f11612` (pinned by `lake-manifest.json`). All Lean runs used `flock` and `taskset -c 0-7`.

Pins:

- openai/math: `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.
- `lean/OAI/Combinatorics/SharpRamsey` git tree: `1a2568639e7e4fae939625dbef5de2aa60d6ce2c`.
  - The local checkout's tree hash is identical to upstream `origin/main` and to the GitHub API value at that commit.
  - The directory has not changed since the initial commit `adc7f124`.
- No local modifications in `SharpRamsey`, `ComparatorChallenges`, `lakefile`, `lean-toolchain`, `lake-manifest.json`, or the Mathlib checkout (`git diff`/`git status` were empty).

## Falsification checks (done before writing)

| # | Question | Method | Result |
|---|---|---|---|
| F1 | Is `ramsey` the standard r(s,t)? | Read `Model.lean`. `RamseyProperty s t N := ∀ G : SimpleGraph (Fin N), (∃ A, G.IsNClique s A) ∨ (∃ B, Gᶜ.IsNClique t B)` and `ramsey s t := sInf {N \| RamseyProperty s t N}`. | Yes. It is the least N forcing K_s or an independent t-set (a t-clique of the complement). There is no off-by-one. `ramsey (3+1) t` is r(4,t). A lower bound cannot hold vacuously: if the set were empty then sInf = 0, which contradicts the positive bound. In any case, `order_lt_ramsey` is proved via `ramsey_has_property`. |
| F2 | Do our definitions agree with the comparator's? | `diff` of `RamseyProperty`/`ramsey`/`MainBounds`/`MainLimit` between `ComparatorChallenges/SharpLogRamsey.lean` and `OAI/.../Model.lean`. | Textually identical. |
| F3 | Where do the restrictions s ≥ 6 and s = 5 come from? | Read `Main.lean`, `FromConstruction.lean`, and `PrimeConstruction.lean`. | `SharpLogRamsey.main (hs : 6 ≤ s)` calls `prime_construction (s-1) (_ : 5 ≤ s-1)`. That function is the one-line wrapper `actual_prime_construction (d-3)`. The case s = 5 is a separate paper with its own comparator (`RamseyFive.lean`). |
| F4 | Is d = 3 covered by the formal proof, or bypassed somehow? | Read `actual_prime_construction`, `Marking.eventually_graph_witness`, and `InitialExperiment.eventually_source_sizes`. Grep the dependency closure. | (1) The proofs are uniform in `i : ℕ`, with vector dimension `i+4`. The only arithmetic side condition is `2 ≤ i+3` in `eventually_source_sizes`. (2) `command grep -P` over all of `SharpRamsey/` found no `sorry`, `admit`, `axiom`, `native_decide`, `ofReduceBool`, `implemented_by`, `extern`, `unsafe`, `opaque`, `debug.*`/`skipKernelTC`, or `addDecl`. There is no `set_option` anywhere. (3) The local instances in `GraphWitness.lean` are `Finite`/`Fintype` instances built from `Finite.of_injective`/`Fintype.ofFinite`, which are harmless. (4) `SharpRamsey` imports nothing from `OAI` outside itself; everything else comes from Mathlib. |
| F5 | Do the lower-bound transfer and the upper bound need s ≥ 6? | Read the hypotheses. | `all_epsilon_lower` needs `1 ≤ d`. `Upper.ramsey_upper` needs `3 ≤ s`. |
| F6 | Does the manuscript give a mathematical reason to exclude d = 3? | Searched both manuscripts for the d ≥ 5 restriction, d = 3, s = 4, and r(4,·). Read §2's intro, the "Repeated projection" proposition, the high-rank lemma, and the final transfer. | No reason is given. The geometric core is done in projective dimensions 2 and 3, and dimensions j ≥ 4 are reduced to these. The r(5,t) paper notes that one counting bound "also covers dimension three". The r(5,t) paper also proves r(4,t) ≤ C t³/log² t as a step. This is a scan, not a full re-audit at d = 3 (see GAP-1). |
| F7 | Is the upper bound classical? | Checked the citation. | Ajtai, Komlós, Szemerédi, "A note on Ramsey numbers", J. Combin. Theory Ser. A 29(3) (1980) 354–360, doi:10.1016/0097-3165(80)90030-8. This gives r(k,n) ≤ c_k n^{k−1}/(log n)^{k−2}. The 1+o(1) constant is in Li, Rousseau, Zang, Graphs Combin. 17 (2001) 123–128. The 170 library proves it in Lean for s ≥ 3. |
| F8 | Has this already been done? | Web search for "r(4,t) lower bound t^3/(log t)^2 2026", Bradač 2605.28793, and "openai math r(4,t)". `gh search issues/prs/code` on openai/math for r(4,t), SharpRamsey, actual_prime_construction, and Ramsey. | No source states a log power below 4 for r(4,t). Mattheus–Verstraëte (Ann. Math. 199 (2024) 919–941) is still cited as state of the art in arXiv:2605.25843 (2026). Bradač's arXiv:2605.28793v3 gives log^{2s−4} = log⁴ at s = 4. openai/math has issues and discussions disabled. Its only PR (#1, 2026-10-08) is a bulk update, and the upstream `SharpLogRamsey.lean` and `docs/170.md` still say s ≥ 6. A code search finds `actual_prime_construction` only in openai/math itself. Search scope: web, arXiv, and GitHub up to 2026-10-09. |

## Mechanical checks

| Check | Command | Result |
|---|---|---|
| Build the library entry point | `lake build OAI.Combinatorics.SharpRamsey.Main` (in `math/lean`) | `Build completed successfully (9331 jobs)`. The SharpRamsey closure had previously been built from scratch: 9329 jobs in 1069 s on 8 cores. |
| Elaborate the remark | `lake env lean lean/R4Lower.lean` | Both theorems elaborate, about 37 s. |
| Axioms | `#print axioms` in the file | `r4_lower` and `main_four` each depend on `[propext, Classical.choice, Quot.sound]`, matching the comparator's `permitted_axioms`. |
| Fresh kernel replay | `LEAN_NUM_THREADS=1 lake env leanchecker --fresh OAI.Combinatorics.SharpRamsey.Main` | Exit 0, no errors; 1:24:22 wall clock (3232 s user CPU, single thread, max RSS 8.7 GB). This replays all of Mathlib plus the SharpRamsey closure in a fresh environment. |
| Script syntax | `bash -n scripts/check.sh` | ok |
| Note | `latexmk -lualatex note.tex` | 3 pages, no errors, no overfull boxes, no undefined references. |

An earlier session had independently checked a 5-line version of `r4_lower` (importing `PrimeConstruction` and `Lower` directly) with the same axiom output.

## GAPs (not covered)

- **GAP-1.** We did not re-audit the 4000-line manuscript at d = 3. The d = 3 claim rests on the Lean proof, which in turn relies on the Lean kernel, Mathlib at the pinned commit, and the definitions checked in F1/F2.
- **GAP-2.** The leanchecker replay covers the `SharpRamsey.Main` closure. Our own file was checked by elaboration plus `#print axioms`, not by a separate replay. The file adds two short theorems on top of that closure.
- **GAP-3.** No human expert has reviewed this. The adversarial checks above were done by AI sessions.
- **GAP-4.** The novelty search is limited to the sources listed in F8. A preprint from the last few days could have been missed.
