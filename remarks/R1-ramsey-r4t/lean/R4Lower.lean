/-
  R1 remark: the s = 4 case of OpenAI family 170 (openai/math).

  Check inside the openai/math Lean project (directory `lean/`), pinned to
  openai/math commit fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb (main, 2026-10-08).
  lean/OAI/Combinatorics/SharpRamsey is unchanged since the initial commit
  adc7f1241b42e322a6451854ab7e4b4c146bf78a (git tree 1a2568639e7e4fae939625dbef5de2aa60d6ce2c).
  Toolchain leanprover/lean4:v4.34.1, Mathlib d13f23b723b8a846827a245b89c10fc7d3f11612.

      cd math/lean
      lake build OAI.Combinatorics.SharpRamsey.Main
      lake env lean /path/to/R4Lower.lean

  Expected output: both `#print axioms` lines list exactly
  [propext, Classical.choice, Quot.sound].

  Nothing below is new mathematics: every step is a library theorem of
  family 170, instantiated at d = 3 (s = 4).  `ramsey` is the library's
  `OAI.SharpLogRamsey.ramsey`, the least N such that every simple graph on
  Fin N has an s-clique or a t-clique in its complement (= r(s,t)).
-/
import OAI.Combinatorics.SharpRamsey.Main

open OAI.SharpLogRamsey

namespace OAIRemarks.R1

/-- Lower bound, s = 4 (d = 3): for every ε > 0, eventually
`t^3 / (log t)^(2+ε) ≤ r(4,t)`.  Uses the library construction
`actual_prime_construction 0 : PrimeConstruction (0+3)`. -/
theorem r4_lower (ε : ℝ) (hε : 0 < ε) :
    ∀ᶠ t : ℕ in Filter.atTop,
      (t:ℝ)^3/(Real.log t)^(((3-1:ℕ):ℝ)+ε) ≤ (ramsey (3+1) t : ℝ) :=
  LowerTransfer.all_epsilon_lower 3 (by norm_num)
    (LowerTransfer.actual_prime_construction 0) ε hε

/-- The comparator-form statement of family 170 at s = 4:
`MainBounds 4 ∧ MainLimit 4`, i.e. the exact statement of
`ComparatorChallenges/SharpLogRamsey.lean` with `s := 4`.
The proof is `main_of_primeConstruction` with the hypothesis `6 ≤ s`
replaced by what its proof actually uses (`3 ≤ s` for the upper bound,
`1 ≤ s-1` for the lower-bound transfer). -/
theorem main_four : MainBounds 4 ∧ MainLimit 4 := by
  have hb : MainBounds 4 := by
    obtain ⟨C, hC, t₀, ht₀⟩ := Upper.ramsey_upper (s := 4) (by norm_num)
    refine ⟨C, hC, fun ε hε => ?_⟩
    have hl := LowerTransfer.all_epsilon_lower 3 (by norm_num)
      (LowerTransfer.actual_prime_construction 0) ε hε
    obtain ⟨t₁, ht₁⟩ := Filter.eventually_atTop.mp hl
    exact ⟨max t₀ t₁, fun t ht => ⟨ht₁ t ((le_max_right _ _).trans ht),
      ht₀ t ((le_max_left _ _).trans ht)⟩⟩
  exact ⟨hb, mainLimit_of_mainBounds 4 hb⟩

end OAIRemarks.R1

#print axioms OAIRemarks.R1.r4_lower
#print axioms OAIRemarks.R1.main_four
