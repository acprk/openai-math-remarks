# Remarks on proofs in OpenAI's mathematics repository

This repository collects short, independently checkable **remarks** on manuscripts in
[openai/math](https://github.com/openai/math). Each remark either

- **observes a consequence** that OpenAI's existing argument or Lean formalization already implies, but that is not stated in the repository;
- **re-selects parameters** inside OpenAI's own argument, and replaces a few specific inequalities with exact computations, so that a non-effective or astronomically large parameter becomes explicit and much smaller; or
- **certifies a limit** of the published method: no choice of parameters within the method's stated hypotheses can do better.

All constructions and main arguments discussed here are **OpenAI's work**. Each remark states exactly what it adds, lists every line of the original manuscript it changes, and states its premises. In particular, a remark that re-selects parameters is only as correct as the manuscript it builds on.

> Status: **local draft, not yet published.** Author attribution is a placeholder.

## Remarks

| ID | OpenAI family | Remark | Kind | Verification |
|---|---|---|---|---|
| R1 | 170 | [r(4,t) = t³/(log t)^{2+o(1)}](remarks/R1-ramsey-r4t/): OpenAI's existing Lean proof already covers s = 4 | consequence + Lean (`main_four`) | Lean kernel replay; adversarial review |
| R3 | 197 | [Odd-characteristic Kaplansky counterexample made explicit](remarks/R3-kaplansky-odd-explicit/): m ≈ 10^360 → 6·10^24, explicit 25-digit prime | explicit parameters | independent re-implementation; adversarial review |
| R4 | 308 | [Smith–Toda V(4) exists if and only if p ≥ 11](remarks/R4-smith-toda-v4-p11/) ("only if" unconditional; "if" conditional on the p = 1009 manuscript) | explicit parameters | two independent enumerations; adversarial review; actual Ext at p = 5, 7 |
| R5 | 238 | [Thorp shuffle: 1600d → 24d](remarks/R5-thorp-24d/) (45d from the main paper alone) | constant | interval certificates (twice); adversarial review |
| R6 | 161 | [A Sidorenko counterexample with 17 vertices and 30 edges](remarks/R6-sidorenko-17-30/) (conditional; premise audited line by line; original 35/66) | smaller counterexample | exact pinning enumeration; adversarial review |
| R7 | 319 | [Hahn–Wilson counterexamples for every prime p ≥ 293](remarks/R7-hahn-wilson-p293/) (was 17,073,437) | explicit parameters | stage enumeration at p = 293; adversarial review |
| R8 | 108 | [Permanent vs determinant cubic constant: 44× uniformly](remarks/R8-permanent-cubic-constant/), and a 167× ceiling for the method | constant + method limit | independent Sage checks; adversarial review |
| R9 | 100 | [Bang half-area cylinder counterexample: 16,000,000 → 640 cylinders](remarks/R9-bang-cylinders-640/) | explicit parameters | exact LP + 1.1M exact points + interval area; adversarial review |
| R10 | 269 | [Fermionic Laughlin gap on the sphere, explicit for every N ≥ 2](remarks/R10-laughlin-fermionic-explicit/): H_Q² ≥ 0.05947 H_Q for all Q ≥ 24 (was: γ < 0.046167 for Q ≥ Q_γ, non-explicit) | explicit parameters | Arb ball/Taylor-model certificates + rigorous ED; independent CG re-implementation; adversarial review |
| R11 | 059 | [Zariski multiplicity counterexample (2 vs 3) in ambient dimension 1,802,976](remarks/R11-zariski-multiplicity-dim/) (was ≈ 4.6·10^143) | explicit parameters | OpenAI's own validator + independent big-integer check; adversarial review |
| R12 | 053 | [Pixton completeness counterexample in genus ≈ 2.26·10^14](remarks/R12-pixton-genus/) (was 10^60) | explicit parameters | independent exact budget re-implementation; adversarial review |
| R13 | 176 | [Second Kahn–Kalai conjecture: constant 2048e^50 ≈ 10^25 → 10^5](remarks/R13-kahn-kalai-constant/) | constant | Arb ball arithmetic re-implementation; adversarial review |
| R14 | 196 | [Kaplansky zero-divisor counterexample over PG(2,16)](remarks/R14-kaplansky-zero-divisor-q16/) (was q = 128; 138 instead of 8,258 generators) | explicit parameters | independent PG(2,q) construction in Sage; adversarial review |


## How these remarks were produced and checked

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction. Each result was checked in three layers:

1. the primary computation;
2. an independent re-implementation, which did not reuse the primary code;
3. an adversarial review by a separate AI session, instructed to refute the claim.

These checks are **not human peer review**. Each remark's `VERIFICATION.md` records the exact scope and the remaining gaps.

## License

Apache License 2.0, the same license as openai/math. Third-party tools (for example plantri, PARI/GP, SageMath, Singular) are not redistributed; download instructions are given where they are needed.
