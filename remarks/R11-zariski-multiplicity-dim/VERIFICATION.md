# Verification record: R11 (shorter product relation for OpenAI family 059, N = 1,802,976)

Date: 2026-10-09. Host: Linux x86_64, Python 3 (stdlib), SageMath 9.5 (search only). At most 3 threads were used. No Lean was run.

Pins: openai/math `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The family-059 directory has not changed since `adc7f124`. In the local checkout the only extra file is the untracked `verification/support/explicit-relation-verification.json`, a report written by an earlier local run of the preprint's verifier. It is not a modification of shipped files.

## Who did what

| Role | Session | What |
|---|---|---|
| Primary computation | R7-skill rescan s0 (`reviews/R7-skill/rescan/s0/`) | Ran `prune.py` (6000-smooth columns, iterative pruning), `ktime.sage` (ℚ-kernel, saturation; 3250 s), and `fromK.sage` (weighted LLL, best E≠0 vector). Output `relK_11_6000.json`, sha256 `e3b747b1…`. The preprint's verifier accepts it with N = 1,802,976. |
| Main agent | — | Re-ran the preprint's verifier on the new relation: `verification_passed: true`. |
| Adversarial review + independent re-implementation | this session (separate AI session) | Full read of the preprint to list the properties of the relation that the proof uses (table (a)). New checker `independent/indep_check_059.py`, written from the tex without the verifier. Negative controls. Dedup search. |

## (a) Adversarial reading: what does the proof need from (e_b)?

Read in full: introduction, seeds, local, global, conclusion, certificate. Background was read by statement.

| # | Question | Finding |
|---|---|---|
| A1 | Role of E ≠ 0 | E > 0 places 4E squares on side 1 (seeds.tex:156). That gives ord₀p₁ = 2 and Hessian rank 16E > 0 (conclusion.tex:284-287, 391-394). If E < 0, the preprint flips the sign (seeds.tex:66). New relation: E = 206 > 0. |
| A2 | s₀ vs N | N = 4(s₀+2)+6+2 = 4s₀+16 (conclusion.tex:19-24). s₀ comes from the weight table, chains, tails, 4E squares and 2 cubes. Exact identity: N = 4Σ\|e_b\|c_b + 8E + 24. The paper's verifier, our checker and the search script all give 1,802,976. |
| A3 | Parity of chain lengths / variable counts | Signed counts are 2e_b, the tail copies have even multiplicity and 4E is even, so s₀ is even for **any** integer e_b (seeds.tex:156-161). The verifier picks 2 cubes, as the paper does. N ≡ 0 (mod 8): 1,802,976 = 8·225,372. |
| A4 | Monodromy order | Divides D = 3^12 because all chain weights have denominator dividing D = 3d. This depends only on d, which is unchanged. The added x⁵+y⁷ and the later 11/13 summands are independent of e. |
| A5 | Mod-2 spectra on both sides | Case 1 uses exactly the two product identities (seeds.tex:197-212). Case 2 (triple angle) holds row by row for any e_b. Case 3 uses only the common cubes. Verified both symbolically (products = 1) and by a sampled Fourier test on the full multisets. |
| A6 | Use of the specific adjugate structure? | **None.** The adjugate appears only in the existence proof and in the finiteness remark (seeds.tex:58-76). The preprint itself offers the 849-row explicit relation, which is not the adjugate vector, as an alternative witness (certificate.tex:11-20). local/global/conclusion reference only Cor. `seed:prepared` and E > 0 (`grep seed:` confirms this). |
| A7 | Requirements outside the verifier's scope | The weight-table range 0<b<d/2 and 3∤b (checked); reduced denominators d and 3d (checked by both codes); isolatedness and degree ≥ 3 of chains (Lemma `seed:chain`, valid for any start 0<c<D/2, re-checked for all 4,640 chains). No other numerical requirement on e_b exists. |
| A8 | Is "N" in the theorem a stated number that we change? | No. Thm 1.1 asserts existence of some N ≡ 0 mod 8, and the preprint says it does not minimize N (introduction.tex:54-56). The "before" value ≈ 4.576·10^143 is the verifier's output on the preprint's own supplied relation. |

## (b) Independent computations

| Check | Command | Result |
|---|---|---|
| Preprint's verifier, new relation | `python3 primary/openai_explicit_relation.py data/relation_d3pow11_N1802976.json` | `verification_passed: true`; 6960 factorizations; s₀ = 450,740; N = 1,802,976; Hessian ranks 3296 / 0 |
| Preprint's verifier, original relation (positive control) | same on `data/openai_explicit-relation.json` | passes; N = 457582925…034696 (144 digits, ≈ 4.576·10^143) |
| **Direct big-integer products** (no factoring) | `independent/indep_check_059.py … --spectral` | ∏A_b^{e_b}: numerator = denominator (118,503 bits); ∏R_b^{e_b}: equal (272,384 bits) |
| Sieve valuations (own code) | same | all prime valuations 0 for both products; on the original relation as well (positive control) |
| Seed construction rebuilt as explicit multisets | same | side counts 29,244 / 30,068 initial; difference 824 = 4E; s₀ = 450,740; N = 1,802,976; original relation reproduces the verifier's 144-digit N exactly |
| Sampled mod-2 spectral test | same | 398 characters k mod 4D (all k<60 plus 400 random; odd 3\|k skipped as trivial); max discrepancy 9.6·10⁻¹⁰ after exact angle reduction |
| Negative controls | `independent/negative_control.py` | +1 on one coefficient, or one sign flipped: both product tests fail, and the spectral discrepancy is 0.57 / 9.2 → rejected |
| Whole suite | `bash scripts/run_all.sh` | all 5 outputs match `expected_output/`; ≈ 15 s |

The independent checker shares no code with the preprint's verifier or with the search scripts. It shares only the reading of the weight table (seeds.tex:138-146) and of the chain algorithm (seeds.tex:94-101).

## (c) Novelty / dedup (up to 2026-10-09)

- GitHub code/repo search: "explicit-relation.json", "Ambiently-homeomorphic", "zariski multiplicity", "1802976". The only third-party computation found is `piiiico/openai-math-counterexamples` row 059, which re-verifies the original 849-row relation (sieve plus mod-p) and does not change N. Catalogue/atlas repos (blisspixel/numinous, satyajitghana/ai, anghel4d/broadside-observer, jcraig949jfi/Prometheus) only index the preprint. openai/math has issues disabled, and no PR mentions this family.
- arXiv/web search: the 2025 papers by Fernandes–Jelonek–Sampaio and others on metric and quasihomogeneous versions. No paper gives an explicit or smaller dimension for this pair.

## Verdict

**CONFIRMED** (premise: the preprint's proof is correct). The new relation satisfies Lemma `seed:product`. Every property of the relation that the preprint uses holds for it. The unchanged construction gives N = 1,802,976.

## GAPs (not covered)

- The topological and lattice-theoretic parts (Sections 4–6 of the preprint) and the cited classification theorems were read, not re-proved. The remark inherits them.
- The relation search (`ktime.sage`, ≈ 55 min) was not re-run in this session. Only its output is relied on, and that output is verified independently.
- No minimality or lower bound. The partial negative runs (d ≤ 3^8; 3^9 with B = 10^5) come from the primary session and are not certified.
- The spectral test is sampled (398 of 4D ≈ 2.1·10⁶ characters). The full equality follows symbolically from the verified product identities.
- No Lean formalization.
