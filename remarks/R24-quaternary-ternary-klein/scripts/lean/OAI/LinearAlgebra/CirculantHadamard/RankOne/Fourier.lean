import Mathlib.Analysis.Fourier.FiniteAbelian.PontryaginDuality

/-!
# Rank-one decomposition from the character identity (generic Fourier step)

Let `G` be a finite abelian group with two homomorphisms `πA : G →+ A`, `πB : G →+ B`
such that `g ↦ (πA g, πB g)` is a bijection `G ≃ A × B`. For `f : G → ℂ` write

  `F(χ, ψ) = ∑ g, f g * χ (πA g) * ψ (πB g)`

for `χ : AddChar A ℂ`, `ψ : AddChar B ℂ`. If for **all** pairs of characters

  `F(χ, ψ) * F(1, 1) = F(χ, 1) * F(1, ψ)`,

then `f` is "rank one" after weighting by its total mass:

  `(∑ h, f h) * f g = (∑ h with πA h = πA g, f h) * (∑ h with πB h = πB g, f h)`.

The proof is Fourier inversion on `A × B` via Mathlib's orthogonality relation
`AddChar.sum_apply_eq_ite`; only injectivity of `(πA, πB)` is used for uniqueness and
only surjectivity for the transform of the product of the two marginals.
-/

namespace OAI

namespace CirculantHadamard.RankOne

open Finset

variable {G A B : Type*} [AddCommGroup G] [Fintype G]
  [AddCommGroup A] [Fintype A] [DecidableEq A]
  [AddCommGroup B] [Fintype B] [DecidableEq B]

/-- The two-variable Fourier transform of `f : G → ℂ` along `(πA, πB)`. -/
def transform (πA : G →+ A) (πB : G →+ B) (f : G → ℂ)
    (χ : AddChar A ℂ) (ψ : AddChar B ℂ) : ℂ :=
  ∑ g, f g * χ (πA g) * ψ (πB g)

/-- The `A`-marginal of `f`, viewed again as a function on `G`. -/
def marginalA (πA : G →+ A) (f : G → ℂ) (g : G) : ℂ :=
  ∑ h with πA h = πA g, f h

/-- The `B`-marginal of `f`, viewed again as a function on `G`. -/
def marginalB (πB : G →+ B) (f : G → ℂ) (g : G) : ℂ :=
  ∑ h with πB h = πB g, f h

omit [Fintype A] [DecidableEq A] [Fintype B] [DecidableEq B] in
theorem transform_one_one (πA : G →+ A) (πB : G →+ B) (f : G → ℂ) :
    transform πA πB f 1 1 = ∑ g, f g := by
  simp [transform]

omit [Fintype A] [DecidableEq A] [Fintype B] [DecidableEq B] in
theorem transform_sub (πA : G →+ A) (πB : G →+ B) (f₁ f₂ : G → ℂ)
    (χ : AddChar A ℂ) (ψ : AddChar B ℂ) :
    transform πA πB (fun g => f₁ g - f₂ g) χ ψ =
      transform πA πB f₁ χ ψ - transform πA πB f₂ χ ψ := by
  simp only [transform, ← Finset.sum_sub_distrib, sub_mul]

omit [Fintype A] [DecidableEq A] [Fintype B] [DecidableEq B] in
theorem transform_const_mul (πA : G →+ A) (πB : G →+ B) (c : ℂ) (f : G → ℂ)
    (χ : AddChar A ℂ) (ψ : AddChar B ℂ) :
    transform πA πB (fun g => c * f g) χ ψ = c * transform πA πB f χ ψ := by
  simp only [transform, Finset.mul_sum, mul_assoc]

/-- Uniqueness of the Fourier transform along an injective pair of homomorphisms. -/
theorem eq_zero_of_transform_eq_zero (πA : G →+ A) (πB : G →+ B)
    (hinj : Function.Injective (fun g => (πA g, πB g))) (f : G → ℂ)
    (hf : ∀ χ ψ, transform πA πB f χ ψ = 0) : f = 0 := by
  classical
  funext g₀
  have key : ∑ χ : AddChar A ℂ, ∑ ψ : AddChar B ℂ,
      χ (-πA g₀) * ψ (-πB g₀) * transform πA πB f χ ψ =
        ((Fintype.card A : ℂ) * (Fintype.card B : ℂ)) * f g₀ := by
    have hterm : ∀ (χ : AddChar A ℂ) (ψ : AddChar B ℂ),
        χ (-πA g₀) * ψ (-πB g₀) * transform πA πB f χ ψ =
          ∑ g, f g * (χ (πA g - πA g₀) * ψ (πB g - πB g₀)) := by
      intro χ ψ
      rw [transform, Finset.mul_sum]
      refine Finset.sum_congr rfl fun g _ => ?_
      rw [sub_eq_add_neg, sub_eq_add_neg, AddChar.map_add_eq_mul, AddChar.map_add_eq_mul]
      ring
    simp_rw [hterm]
    have hswap : ∑ χ : AddChar A ℂ, ∑ ψ : AddChar B ℂ, ∑ g : G,
        f g * (χ (πA g - πA g₀) * ψ (πB g - πB g₀)) =
        ∑ g : G, ∑ χ : AddChar A ℂ, ∑ ψ : AddChar B ℂ,
          f g * (χ (πA g - πA g₀) * ψ (πB g - πB g₀)) := by
      calc
        _ = ∑ χ : AddChar A ℂ, ∑ g : G, ∑ ψ : AddChar B ℂ,
            f g * (χ (πA g - πA g₀) * ψ (πB g - πB g₀)) :=
          Finset.sum_congr rfl fun χ _ => Finset.sum_comm
        _ = _ := Finset.sum_comm
    rw [hswap]
    have hinner : ∀ g : G, ∑ χ : AddChar A ℂ, ∑ ψ : AddChar B ℂ,
        f g * (χ (πA g - πA g₀) * ψ (πB g - πB g₀)) =
          f g * ((∑ χ : AddChar A ℂ, χ (πA g - πA g₀)) *
            (∑ ψ : AddChar B ℂ, ψ (πB g - πB g₀))) := by
      intro g
      rw [Finset.sum_mul_sum, Finset.mul_sum]
      refine Finset.sum_congr rfl fun χ _ => ?_
      rw [Finset.mul_sum]
    simp_rw [hinner, AddChar.sum_apply_eq_ite, sub_eq_zero]
    rw [Finset.sum_eq_single g₀]
    · simp
      ring
    · intro g _ hg
      by_cases hA : πA g = πA g₀
      · by_cases hB : πB g = πB g₀
        · exact absurd (hinj (Prod.ext hA hB)) hg
        · simp [hB]
      · simp [hA]
    · intro h
      exact absurd (Finset.mem_univ g₀) h
  have hzero : ∑ χ : AddChar A ℂ, ∑ ψ : AddChar B ℂ,
      χ (-πA g₀) * ψ (-πB g₀) * transform πA πB f χ ψ = 0 := by
    simp [hf]
  rw [hzero] at key
  have hcard : ((Fintype.card A : ℂ) * (Fintype.card B : ℂ)) ≠ 0 :=
    mul_ne_zero (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)
      (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)
  have := (mul_eq_zero.mp key.symm).resolve_left hcard
  simpa using this

omit [Fintype A] [Fintype B] in
/-- The transform of the product of the two marginals factors. -/
theorem transform_marginal_mul (πA : G →+ A) (πB : G →+ B)
    (hbij : Function.Bijective (fun g => (πA g, πB g))) (f : G → ℂ)
    (χ : AddChar A ℂ) (ψ : AddChar B ℂ) :
    transform πA πB (fun g => marginalA πA f g * marginalB πB f g) χ ψ =
      transform πA πB f χ 1 * transform πA πB f 1 ψ := by
  classical
  -- Every pair of fibres meets in exactly one point.
  have hcount : ∀ h₁ h₂ : G, ∀ c : ℂ,
      (∑ g : G, if πA h₁ = πA g ∧ πB h₂ = πB g then c else 0) = c := by
    intro h₁ h₂ c
    obtain ⟨g₀, hg₀⟩ := hbij.2 (πA h₁, πB h₂)
    simp only [Prod.mk.injEq] at hg₀
    rw [Finset.sum_eq_single g₀]
    · simp [hg₀.1, hg₀.2]
    · intro g _ hg
      split_ifs with hc
      · exact absurd (hbij.1 (Prod.ext (hc.1.symm.trans hg₀.1.symm)
          (hc.2.symm.trans hg₀.2.symm))) hg
      · rfl
    · intro h
      exact absurd (Finset.mem_univ g₀) h
  have hlhs : transform πA πB (fun g => marginalA πA f g * marginalB πB f g) χ ψ =
      ∑ g : G, ∑ h₁ : G, ∑ h₂ : G,
        if πA h₁ = πA g ∧ πB h₂ = πB g then
          (f h₁ * χ (πA h₁)) * (f h₂ * ψ (πB h₂)) else 0 := by
    unfold transform marginalA marginalB
    refine Finset.sum_congr rfl fun g _ => ?_
    dsimp only
    rw [Finset.sum_filter, Finset.sum_filter, Finset.sum_mul_sum, Finset.sum_mul,
      Finset.sum_mul]
    refine Finset.sum_congr rfl fun h₁ _ => ?_
    rw [Finset.sum_mul, Finset.sum_mul]
    refine Finset.sum_congr rfl fun h₂ _ => ?_
    by_cases hA : πA h₁ = πA g
    · by_cases hB : πB h₂ = πB g
      · simp only [hA, hB, and_self, ite_true]
        ring
      · simp [hB]
    · simp [hA]
  rw [hlhs]
  have hswap : (∑ g : G, ∑ h₁ : G, ∑ h₂ : G,
        if πA h₁ = πA g ∧ πB h₂ = πB g then
          (f h₁ * χ (πA h₁)) * (f h₂ * ψ (πB h₂)) else 0) =
      ∑ h₁ : G, ∑ h₂ : G, ∑ g : G,
        if πA h₁ = πA g ∧ πB h₂ = πB g then
          (f h₁ * χ (πA h₁)) * (f h₂ * ψ (πB h₂)) else 0 := by
    calc
      _ = ∑ h₁ : G, ∑ g : G, ∑ h₂ : G,
          if πA h₁ = πA g ∧ πB h₂ = πB g then
            (f h₁ * χ (πA h₁)) * (f h₂ * ψ (πB h₂)) else 0 := Finset.sum_comm
      _ = _ := Finset.sum_congr rfl fun h₁ _ => Finset.sum_comm
  rw [hswap]
  simp_rw [hcount]
  rw [transform, transform, Finset.sum_mul_sum]
  simp

/-- **Generic rank-one theorem.** The full family of character identities
`F(χ,ψ) F(1,1) = F(χ,1) F(1,ψ)` forces the mass-weighted rank-one identity. -/
theorem rankOne_of_transform (πA : G →+ A) (πB : G →+ B)
    (hbij : Function.Bijective (fun g => (πA g, πB g))) (f : G → ℂ)
    (hF : ∀ (χ : AddChar A ℂ) (ψ : AddChar B ℂ),
      transform πA πB f χ ψ * transform πA πB f 1 1 =
        transform πA πB f χ 1 * transform πA πB f 1 ψ)
    (g : G) :
    (∑ h, f h) * f g = (∑ h with πA h = πA g, f h) * (∑ h with πB h = πB g, f h) := by
  have hzero := eq_zero_of_transform_eq_zero πA πB hbij.1
    (fun g => (∑ h, f h) * f g - marginalA πA f g * marginalB πB f g) (by
      intro χ ψ
      rw [transform_sub, transform_const_mul, transform_marginal_mul πA πB hbij,
        ← transform_one_one πA πB f, mul_comm, hF, sub_self])
  have := congrFun hzero g
  simpa [marginalA, marginalB, sub_eq_zero] using this

end CirculantHadamard.RankOne

end OAI
