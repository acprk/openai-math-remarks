import OAI.LinearAlgebra.CirculantHadamard.AlternatingLocalization
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.CyclotomicAwayLocal
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.PrimeCharacterEvaluations

namespace OAI

/-!
# Local representatives of the alternating product over ℤ[ω], arbitrary norm N with 3 ∤ N

Generalization of the section `ActualAlternatingProduct` of
`OAI.LinearAlgebra.CirculantHadamard.AlternatingLocalization`, which treats only
`x : Elem EisensteinRing (u ^ 2)` with `x * x^* = u ^ 2`. Here the cyclic group is
`ZMod N` for an arbitrary odd `N` and the norm is `N` itself. All ring-generic
lemmas are imported unchanged from the original file.
-/

noncomputable section

namespace CirculantHadamard.AlternatingLocalizationW

open scoped BigOperators
open CirculantHadamard.AlternatingLocalization

section ActualAlternatingProduct

open CyclotomicRingsW PrimeComponentsN PrimeCharacterEvaluationsW CyclicRing
open CyclotomicRings (rho rho_isPrimitiveRoot rho_pow)

instance BW_local_algebra (N : ℕ) (m : MaximalSpectrum (B N)) :
    Algebra (B N) (AtMaximal (B N) m) :=
  OreLocalization.instAlgebra (R := B N) (R₀ := B N) (S := m.asIdeal.primeCompl)

instance BW_fraction_algebra (N : ℕ) : Algebra (B N) (FractionRing (B N)) :=
  OreLocalization.instAlgebra (R := B N) (R₀ := B N) (S := nonZeroDivisors (B N))

/-- The selected element of B satisfies the actual cyclotomic relation. -/
theorem rhoB_cyclotomic_root (N : ℕ) (p : PrimeIndex N) :
    (Polynomial.cyclotomic p.val (B N)).eval (rhoB N p) = 0 := by
  have hp := prime_of_mem p.property
  apply Subtype.val_injective
  change (B N).val ((Polynomial.cyclotomic p.val (B N)).eval (rhoB N p)) =
    (B N).val (0 : B N)
  rw [map_zero]
  change (B N).val.toRingHom
    ((Polynomial.cyclotomic p.val (B N)).eval (rhoB N p)) = 0
  rw [← Polynomial.cyclotomic.eval_apply]
  exact (rho_isPrimitiveRoot p.val hp.ne_zero).isRoot_cyclotomic hp.pos

theorem subset_evaluations_local_associated (N : ℕ) [NeZero N]
    (hu : 0 < N) (huodd : ¬ 3 ∣ N)
    (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (p : PrimeIndex N) (m : MaximalSpectrum (B N))
    (hpm : (p.val : B N) ∈ m.asIdeal)
    (S : Finset (PrimeIndex N)) (hpS : p ∉ S) :
    LocalComparison.ResidueOneAssociated (IsLocalRing.residue (AtMaximal (B N) m))
      (algebraMap (B N) (AtMaximal (B N) m) (subsetEvaluation N hu (insert p S) x))
      (algebraMap (B N) (AtMaximal (B N) m) (subsetEvaluation N hu S x)) := by
  let f := (awayInclusion N p.val).toRingHom
  let A := AtComap f m
  let : NeZero p.val := ⟨(prime_of_mem p.property).ne_zero⟩
  let : CharZero A := CyclotomicAwayLocalW.awayLocal_charZero N p.val m.asIdeal
  let : IsDiscreteValuationRing A :=
    CyclotomicAwayLocalW.awayLocal_isDiscreteValuationRing N p.val huodd p.property m.asIdeal hpm
  have hmax : IsLocalRing.maximalIdeal A = Ideal.span {(p.val : A)} :=
    CyclotomicAwayLocalW.awayLocal_maximalIdeal_eq_span N p.val huodd p.property m.asIdeal hpm
  let e := exponent N p.val - 1
  have he : exponent N p.val = e + 1 := by
    dsimp [e]
    have hpos := exponent_pos hu p.property
    omega
  have hn : order N p.val = p.val ^ (e + 1) := by rw [order, he]
  have hc : IsUnit (cofactor N p.val : A) :=
    coprime_natCast_isUnit (m.asIdeal.comap f) p.val (cofactor N p.val)
      (cofactor_coprime hu p.property)
      (CyclotomicAwayLocalW.prime_mem_awayPrime N p.val m.asIdeal hpm)
  let F := partialEvaluation N hu p S x
  let K := partialEvaluation N hu p S (ringStar x)
  have hnum : N = p.val ^ (e + 1) * cofactor N p.val := by
    rw [← hn]
    exact (order_mul_cofactor N p.val).symm
  have hFK : F * K = scalar (order N p.val)
      ((p.val : awayRing N p.val) ^ (e + 1) * (cofactor N p.val : awayRing N p.val)) := by
    dsimp [F, K]
    rw [partialEvaluation_product N hu p S x (ringStar x) hx]
    apply congrArg (scalar (order N p.val))
    exact_mod_cast hnum
  have h := actual_local_character_association f m p.val (prime_of_mem p.property)
    hpm hmax (rhoB N p) (rhoB_pow N p) (rhoB_cyclotomic_root N p)
    (order N p.val) e (cofactor N p.val) hn
    (root_order_pow N hu p (rhoB N p) (rhoB_pow N p)) hc F K hFK
  change LocalComparison.ResidueOneAssociated _
    (algebraMap (B N) (AtMaximal (B N) m) (partialAtRoot N hu p F))
    (algebraMap (B N) (AtMaximal (B N) m) (partialAtOne N p F)) at h
  simpa only [F, partialEvaluation_at_root_apply,
    partialEvaluation_at_one_apply N hu p S hpS x] using h

/-- The denominator is nonzero because it is an actual evaluated norm factor. -/
theorem exists_local_edge_unit (N : ℕ) [NeZero N]
    (hu : 0 < N) (huodd : ¬ 3 ∣ N)
    (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (p : PrimeIndex N) (m : MaximalSpectrum (B N))
    (hpm : (p.val : B N) ∈ m.asIdeal)
    (S : Finset (PrimeIndex N)) (hpS : p ∉ S) :
    ∃ v : (AtMaximal (B N) m)ˣ,
      toFraction (FractionRing (B N)) m (v : AtMaximal (B N) m) =
        algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu (insert p S) x) /
          algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu S x) ∧
      IsLocalRing.residue (AtMaximal (B N) m) (v : AtMaximal (B N) m) = 1 := by
  apply local_edge_of_associated m _ _
  · intro hz
    apply subsetEvaluation_ne_zero N hu S x hx
    apply IsFractionRing.injective (B N) (FractionRing (B N))
    simpa only [map_zero] using hz
  · exact subset_evaluations_local_associated N hu huodd x hx p m hpm S hpS

/-- At every maximal ideal over p, the one fixed global alternating product
is the image of a unit with residue one. No edge comparison is an input. -/
theorem exists_local_unit_alternatingProduct (N : ℕ) [NeZero N]
    (hu : 0 < N) (huodd : ¬ 3 ∣ N)
    (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (p : PrimeIndex N) (m : MaximalSpectrum (B N))
    (hpm : (p.val : B N) ∈ m.asIdeal) :
    ∃ v : (AtMaximal (B N) m)ˣ,
      toFraction (FractionRing (B N)) m (v : AtMaximal (B N) m) =
        alternatingProduct (Finset.univ : Finset (PrimeIndex N))
          (fun S => algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu S x)) ∧
      IsLocalRing.residue (AtMaximal (B N) m) (v : AtMaximal (B N) m) = 1 := by
  classical
  apply exists_unit_of_edges (toFraction (FractionRing (B N)) m).toRingHom
    (IsLocalRing.residue (AtMaximal (B N) m)) Finset.univ p (Finset.mem_univ p)
  intro S hS
  have hpS : p ∉ S := fun hp =>
    (Finset.mem_erase.mp ((Finset.mem_powerset.mp hS) hp)).1 rfl
  exact exists_local_edge_unit N hu huodd x hx p m hpm S hpS

/-- Integrality of the actual subset-character alternating product. The
maximal ideals over prime divisors use the proved local comparison; every
other maximal ideal uses the actual norm identity. -/
theorem exists_integral_alternatingProduct (N : ℕ) [NeZero N]
    (hu : 0 < N) (huodd : ¬ 3 ∣ N)
    (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing)) :
    ∃ z : B N, algebraMap (B N) (FractionRing (B N)) z =
      alternatingProduct (Finset.univ : Finset (PrimeIndex N))
        (fun S => algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu S x)) := by
  classical
  refine exists_base_of_edges_and_norm N hu.ne' Finset.univ
    (fun p : PrimeIndex N => p.val) ?_
    (fun S => subsetEvaluation N hu S x)
    (fun S => star (subsetEvaluation N hu S x) * subsetEvaluation N hu S x *
      star (subsetEvaluation N hu S x)) ?_ ?_
  · intro p hp
    exact ⟨⟨p, hp⟩, Finset.mem_univ _, rfl⟩
  · intro S _
    have h := subsetEvaluation_norm N hu S x hx
    calc subsetEvaluation N hu S x * (star (subsetEvaluation N hu S x) *
          subsetEvaluation N hu S x * star (subsetEvaluation N hu S x))
        = (subsetEvaluation N hu S x * star (subsetEvaluation N hu S x)) ^ 2 := by ring
      _ = (N : B N) ^ 2 := by rw [h]
  · intro m p _ hpm S hS
    have hpS : p ∉ S := fun hp =>
      (Finset.mem_erase.mp ((Finset.mem_powerset.mp hS) hp)).1 rfl
    exact exists_local_edge_unit N hu huodd x hx p m hpm S hpS

end ActualAlternatingProduct

end CirculantHadamard.AlternatingLocalizationW

end

end OAI
