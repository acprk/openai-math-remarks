# R1. The existing OpenAI family-170 Lean proof already gives r(4,t) = t³/(log t)^{2+o(1)}

**What this remark shows.** OpenAI's family-170 manuscripts state the sharp logarithmic exponent of r(s,t) only for s ≥ 5. Their Lean library already proves the prime-field construction for every projective dimension d ≥ 3. Applying it at d = 3 gives the s = 4 case, r(4,t) = t³/(log t)^{2+o(1)}, through a ten-line Lean file whose axioms are the three standard ones.

## Statement

r(s,t) is the least N such that every simple graph on N vertices contains K_s or an independent set of size t.

**Before.** These are the stated results of family 170 and the literature:

| s | stated in openai/math 170 | best lower bound for r(4,t) we know of |
|---|---|---|
| 4 | not stated | c·t³/(log t)⁴ (Mattheus–Verstraëte, Ann. Math. 199 (2024) 919–941; also Bradač, arXiv:2605.28793, log^{2s−4}) |
| 5 | t⁴/(log t)^{3+o(1)} (companion paper; `ComparatorChallenges/RamseyFive.lean`) | |
| ≥ 6 | t^{s−1}/(log t)^{s−2+o(1)} (`ComparatorChallenges/SharpLogRamsey.lean`, `theorem main (s) (hs : 6 ≤ s)`) | |

**After.** There is a C > 0 such that for every ε > 0 and all sufficiently large t,

    t³/(log t)^{2+ε} ≤ r(4,t) ≤ C·t³/(log t)²,   so   (3 log t − log r(4,t)) / log log t → 2.

In Lean this is `MainBounds 4 ∧ MainLimit 4`, the comparator statement with s := 4. The upper bound is classical (Ajtai–Komlós–Szemerédi, J. Combin. Theory Ser. A 29 (1980) 354–360; leading constant 1+o(1) in Li–Rousseau–Zang, Graphs Combin. 17 (2001) 123–128). The 170 library also proves it for all s ≥ 3.

## Attribution

The construction, every proof, and the formalization are OpenAI's work: family 170 of [github.com/openai/math](https://github.com/openai/math), specifically

- [`preprints/Sharp-Logarithmic-Exponents-for-Fixed-Off-Diagonal-Ramsey-Numbers-September-24-2026`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/Sharp-Logarithmic-Exponents-for-Fixed-Off-Diagonal-Ramsey-Numbers-September-24-2026)
- [`preprints/The-Sharp-Logarithmic-Exponent-of-r-5-t-September-24-2026`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/The-Sharp-Logarithmic-Exponent-of-r-5-t-September-24-2026)
- [`lean/OAI/Combinatorics/SharpRamsey/`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/lean/OAI/Combinatorics/SharpRamsey) and `lean/docs/170.md`

These files are pinned at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (main, 2026-10-08). The `SharpRamsey` directory is unchanged since the initial commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`; its git tree is `1a2568639e7e4fae939625dbef5de2aa60d6ce2c`.

This remark does one thing. It points out that the existing formal proof covers d = 3, and it supplies the instantiation. **No new mathematics is claimed.**

## Why it works

- `Construction/PrimeConstruction.lean` has two relevant declarations:
  - `actual_prime_construction (i : ℕ) : PrimeConstruction (i+3)` holds for every i ≥ 0. Its proof uses vector dimension i+4, so it covers every projective dimension d ≥ 3.
  - `prime_construction (d) (hd : 5 ≤ d)` is only a wrapper that calls `actual_prime_construction (d-3)`. It is the only place where d ≥ 5 appears.
- `Bounds/Lower.lean` contains `all_epsilon_lower (d) (hd : 1 ≤ d)`, which transfers any `PrimeConstruction d` to t^d/(log t)^{d−1+ε} ≤ ramsey (d+1) t.
- `Bounds/SamplingParameter.lean` contains `Upper.ramsey_upper (hs : 3 ≤ s)`, which proves the AKS-type upper bound.
- `Bounds/FromConstruction.lean` contains `main_of_primeConstruction (s) (hs : 6 ≤ s)`. Its proof only uses 3 ≤ s and 1 ≤ s−1. `lean/R4Lower.lean` repeats that proof with s = 4.
- The manuscript is consistent with d = 3. Its sparse-pair description is proved directly in projective dimensions 2 and 3, and the "Repeated projection" proposition (§2) reduces higher dimensions to those two. Neither manuscript gives a reason for excluding s = 4.

## Reproduce

```bash
git clone https://github.com/openai/math && cd math && git checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
bash /path/to/R1-ramsey-r4t/scripts/check.sh "$PWD"
```

Expected output is in `scripts/expected_output/lean_check.txt`: both theorems depend only on `[propext, Classical.choice, Quot.sound]`.

If the SharpRamsey oleans are already built, `lake env lean R4Lower.lean` takes about 40 s. A cold build of the SharpRamsey closure took 1069 s on 8 cores.

The note is `note.tex` / `note.pdf` (3 pp., compile with `latexmk -lualatex note.tex`).

## Scope and limitations (what is NOT claimed)

- **No new mathematics.** The construction, its proofs, and their formalization are OpenAI's.
- **No re-audit of the manuscript at d = 3.** The evidence for d = 3 is the machine-checked Lean proof. This relies on trusting the Lean kernel, Mathlib `d13f23b`, and the definitions in `Model.lean`. Those definitions are textually identical to the comparator's and match the standard r(s,t).
- **No constant-factor result.** We do not claim r(4,t) = Θ(t³/(log t)²). The lower bound loses (log t)^ε.
- **The d = 4 analogue is not pursued.** It would recover r(5,t) from the general library.
- **Novelty is limited to our search.** To our knowledge (search of arXiv, the web, and the openai/math GitHub repository up to 2026-10-09; the repository has issues and discussions disabled), no earlier source states a lower bound for r(4,t) with a log power below 4, or points out that the 170 formalization covers s = 4.

## Authorship

[AUTHOR NAME].

Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.
