import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremQ
import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremT1
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.TheoremT3
import OAI.LinearAlgebra.CirculantHadamard.RankOne.TheoremF

/-! Statement and axiom check for this remark. Run from the root of the upstream Lean
project (`lean/` of openai/math) after copying `OAI/` from this directory into it:
`lake env lean Checks.lean`. -/

open OAI.CirculantHadamard

-- Prop' (alternating product, Z[i], any odd N, norm N) and Prop'' (Z[omega], 3 ∤ N)
#check @GeneralNorm.alternatingProduct_odd_torsion
#check @GeneralNorm.alternatingProduct_eq_one_of_two_primes
#check @EisensteinNorm.alternatingProduct_coprime_three_torsion
-- Master Lemma, sum form, Z[i] and Z[omega]
#check @GeneralNorm.master_lemma
#check @EisensteinNorm.master_lemma
-- Theorems Q, T1, T3, F
#check @GeneralNorm.QuaternaryQ.not_exists_perfect_quaternary
#check @GeneralNorm.MenonT1.not_exists_perfect_binary_array
#check @GeneralNorm.MenonT1.not_exists_menon_C2_C2_square
#check @EisensteinNorm.TernaryT3.not_exists_perfect_ternary
#check @RankOne.theoremF_gaussian

#print axioms GeneralNorm.alternatingProduct_odd_torsion
#print axioms EisensteinNorm.alternatingProduct_coprime_three_torsion
#print axioms GeneralNorm.master_lemma
#print axioms EisensteinNorm.master_lemma
#print axioms GeneralNorm.QuaternaryQ.not_exists_perfect_quaternary
#print axioms GeneralNorm.MenonT1.not_exists_perfect_binary_array
#print axioms GeneralNorm.MenonT1.not_exists_menon_C2_C2_square
#print axioms EisensteinNorm.TernaryT3.not_exists_perfect_ternary
#print axioms RankOne.theoremF_gaussian
