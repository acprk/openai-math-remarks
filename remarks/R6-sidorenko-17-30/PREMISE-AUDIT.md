# R6 — Line-by-line audit of premise (P) for K₇

Premise (P) says that Lemma 3.1 (rank-layer limit) and its proof in Sections 4–7 of the original paper
(`transverse.tex`, `orbits.tex`, `planted.tex`, `tail.tex`) hold for the 7-point complex K₇ once the
complex-dependent constants are re-instantiated. This file records an audit of that premise, made on 2026-10-09
by a separate AI session that had not written the construction.

**What this audit is, and is not.**
- It is a reading of every line of Sections 4–7 that uses information about the complex. Each such line was
  classified, and every complex-dependent number was recomputed by a script written from scratch:
  `scripts/premise_audit/rank_layer_audit.py`.
- It is **not** a re-derivation of the proof, **not** human review, and **not** a Lean formalisation. R6
  therefore stays **conditional**. The audit lowers the risk of premise (P); it does not remove the premise.
- Calibration: on the paper's 13-point complex the same script reproduces every number the paper states:
  degrees (4,6,4,6,5,4,5,5,5,6,5,6,5), Σ C(2deg,2) = 618, m² − 4m + 618 ≤ 1575, remainder 1641, budget 33M,
  threshold M* + 44C*, β_D^{-35}, and exposure orders for all 6 bit classes.

## Result

The audit found 21 places that use information about the complex. Each was put into one of three classes:
- **(a) constant only.** The argument is unchanged after substituting the K₇ value.
- **(b) uses only (1)–(3)** of Prop. `complex`.
- **(c) uses some other structural property** of the 13-point complex, for example minimum degree, absence of
  short cycles, a complex-specific case enumeration, condition (4), or a margin that needs V or F to be large.

**No place falls into class (c).** All margins come from "M large" and "D large", independently of V and F.
Smaller complexes give smaller constants.

| # | Location | Use in the original | K₇ value | Class |
|---|---|---|---|---|
| 1 | T L22 | "all 33 point pairs" | 15 | a |
| 2 | T L299–301, L583–588 | union bound over 𝒯, 2^{−\|E\|} expansion, "E fixed" | \|E\| = 15 | a |
| 3 | T L454–457 | point degree ≤ 6 from Table 1; degree = number of incident faces | degrees (6,4,5,3,4,5,3), max 6; degree = #faces holds | a + b(1) |
| 4 | T L461, L508–509 | D ≥ 12 from s ≤ 6, b(N−s) ≥ 15 > s | max degree 6 ⇒ still D ≥ 12 (recomputed) | a |
| 5 | T L611–622 | C_D^{\|J\|−1}, 8^{\|J\|} | \|J\| = 10 | a |
| 6 | T L593–630 | every activation set B ⊆ 𝒞 | 2³⁰ sets | a |
| 7 | O L193, L201, L221–222 | 35 Lagrangians, β_D^{-35} | β_D^{-17} | a |
| 8 | all of O and Pl | geometry of a single face / triple (orbit masses, planted gain, g₂(h)) | — | a (local) |
| 9 | Ta L39–41 | a_e ≥ 1/3, some a_e > 1/3 | 3a_e histogram {1:11, 2:2, 3:2} | b(3) |
| 10 | Ta L53–59 | six classes of eleven faces, thirteen points, ten faces each add one point | six classes of 5 faces, 7 points, 4 faces each add one point | b(2) |
| 11 | Ta L84–93 | P(A)⁶ ≤ ∏ P(A_t); each face in 3 classes; e internal in 3(1−a_e) classes | identical (3 bits) | b(2,3) |
| 12 | Ta L129–131 | "each pair occurs twice" ⇒ 𝒞 − Σg = ΣS_j + ℰ; the 1/6 baseline | identical | b(1,3) |
| 13 | Ta L140–345 | lower-active lemma, slack identity, coercivity (1/342) | — | a (single face) |
| 14 | Ta L363–373 | budget ≤ 33M = 22 · (3M/2) | ≤ 15M = 10 · (3M/2) | a |
| 15 | Ta L409–422 | face graph connected, propagation along shortest paths, diameter ℓ* | connected, ℓ* = 3 (original 5) | b(1) |
| 16 | Ta L424–433 | a pair e* with a_{e*} > 1/3 gives ℰ ≥ (a − 1/3)M²/16 | e.g. e* = {2,3}, a = 1 | b(3) |
| 17 | Ta L462–473 | D > M* + 44C*, "22 faces" | D > M* + 20C*, 10 faces | a |
| 18 | Ta L475–509 | integer equality cases c ≥ 2, c = 1, c = 0 with (p,t) ∈ {(0,0),(1,1)}; propagation; Σg = (2/3)Σb(h_e) | enumeration is inside one face; propagation uses (1), contradiction uses (3). Residual cases: all c_j = 0, h_e ∈ {0,2}, h_e = 0 on the 4 pairs with a_e > 1/3 | a + b(1,3) |
| 19 | Ta L566–616, L673–696 | non-direct case: Σ C(2deg,2) = 618; m ≤ 33; m² − 4m + 618 ≤ 1575; gain ≤ 66 ⇒ 1641 | **242**; m ≤ 15; ≤ **407**; gain ≤ 30 ⇒ **437**, so D > 437 suffices | a |
| 20 | Ta L619–655 | direct case: two distinct pairs of a simple graph share at most one endpoint | (I,E) is a simple graph | a |
| 21 | Ta L682–686 | 6𝒞 is an integer because 3a_e is an integer | identical (3 bits) | b(3) |

T = `transverse.tex`, O = `orbits.tex`, Pl = `planted.tex`, Ta = `tail.tex`, line numbers at openai/math commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.

**Substitutions for a written K₇ version.** 33 → 15, 35 → 17, 66 → 30, 22 → 10, "eleven / thirteen / ten" →
"five / seven / four", 618 → 242, 1575 → 407, 1641 → 437, 33M → 15M, 44C* → 20C*, β_D^{-35} → β_D^{-17}.
D ≥ 12 is unchanged, and ℓ* = 3. Table 1 and Table 2 are replaced by the bit table and exposure orders in
`README.md`.

## Structural remark

H₇ is a "twisted double" along the separating triangle {0,2,5}. The two halves each have 2 points and 5 faces,
and the automorphism φ = (1 4)(2 5)(3 6) swaps them. However, φ also swaps the points 2 and 5 of the separator,
so t(H₇) = E[g(x₀,x₂,x₅) g(x₀,x₅,x₂)] is not a square, and the reflection (Cauchy–Schwarz) argument misses by
exactly this step. The paper's 13-point complex has no non-trivial separating set of size 3; larger separators
were not checked.

## Limits

Lemma 3.1 is an asymptotic statement as q → ∞, with D > max(437, M* + 20C*). That is far beyond any size that
can be computed, so no numerical test of (P) was attempted. A toy-scale test would only be a unit test.
