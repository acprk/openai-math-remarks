import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.AlternatingTorsion
import OAI.LinearAlgebra.CirculantHadamard.BinaryLocalRing
import OAI.LinearAlgebra.CirculantHadamard.BinaryProductTransport
import OAI.LinearAlgebra.CirculantHadamard.BinaryCharacterSum
import OAI.LinearAlgebra.CirculantHadamard.BinaryResidue

namespace OAI

/-!
# The λ-adic Master Lemma over `ℤ[i]` (λ = 1 + i)

Let `N > 1` be odd and let `x₀, x_j (j ∈ s)` be elements of `ℤ[i][C_N]` with
`x * x^* = N`. Suppose `x_j - x₀ = t_j y_j` with `0 ≠ t_j ∈ (1+i)ℤ[i]`, and
`∑_j y_j = e J + (1+i) W`, where `J` is the all-ones element of the group ring and
`e ≡ 1 (mod 1+i)`. Then we have a contradiction.

Proof: in the localization `O` of `B = ℤ[i, ρ_p : p ∣ N]` at a maximal ideal over
`2`, every character value of every `x` is a unit (odd norm). By Prop′
(`GeneralNorm.alternatingProduct_odd_torsion`) the quotient of alternating products
`Δ(x_j)/Δ(x₀)` is a root of unity of odd order with residue `1`, hence equals `1`.
The first-order homomorphism along `t_j` then gives
`∑_S (-1)^{|S|} y_{j,S}/x_{0,S} ≡ 0`; summing over `j` gives
`∑_S (-1)^{|S|} J_S / x_{0,S} ≡ 0`, but only `S = ∅` survives (`J_∅ = N` odd).

The original circulant Hadamard proof (`BinaryLocalFamily`) is the special case
`x₀ = c`, `x₁ = d`, `x₂ = g` with square `N = u ^ 2`.
-/

noncomputable section

namespace CirculantHadamard.GeneralNorm

open scoped BigOperators
open CyclotomicRings PrimeComponentsN PrimeCharacterEvaluationsN CyclicRing
open BinaryLocalRing BinaryProductTransport

/-- The all-ones element `J = ∑_g g` of `ℤ[i][ZMod N]`. -/
abbrev Jelem (N : ℕ) [NeZero N] : Elem GaussianRing N :=
  AddMonoidAlgebra.ofCoeff (BinaryCoefficients.JOver GaussianRing)

section JSum

variable (N : ℕ) [NeZero N] (hN : 0 < N)

omit [NeZero N] in
theorem subsetCharacter_nontrivial (S : Finset (PrimeIndex N)) (hS : S.Nonempty) :
    ∃ g : Multiplicative (ZMod N), subsetCharacter N hN S g ≠ 1 := by
  obtain ⟨a, ha⟩ := subsetCharacter_nonempty N hN S hS
  exact ⟨Multiplicative.ofAdd a, ha⟩

/-- `∑_S (-1)^{|S|} ψ(χ_S(J)) d_S = N · d_∅`. -/
theorem subset_J_sum {T : Type*} [CommRing T] (ψ : B N →+* T)
    (d : Finset (PrimeIndex N) → T) :
    (∑ S : Finset (PrimeIndex N), (subsetSign S : T) *
      ψ (subsetEvaluation N hN S (Jelem N)) * d S) = ((N : ℕ) : T) * d ∅ := by
  exact BinaryCharacterSum.subset_map_evaluate_J_mul (gaussianToB N).toRingHom ψ
    (subsetCharacter N hN) (fun S => (subsetSign S : T)) d (by simp)
    (subsetCharacter_empty N hN) (subsetCharacter_nontrivial N hN)

theorem subset_J_residue_ne_zero
    {T k : Type*} [CommRing T] [CommRing k] [Nontrivial k] [CharP k 2]
    (ψ : B N →+* T) (residue : T →+* k)
    (d : Finset (PrimeIndex N) → T) (c₀ : Tˣ)
    (hd : d ∅ = ↑(c₀⁻¹)) (hodd : Odd N) :
    residue (∑ S : Finset (PrimeIndex N), (subsetSign S : T) *
      ψ (subsetEvaluation N hN S (Jelem N)) * d S) ≠ 0 := by
  rw [subset_J_sum N hN ψ d, hd]
  exact BinaryCharacterSum.residue_odd_mul_unit_inv_ne_zero residue hodd c₀

end JSum

section LocalUnits

/-- An element of odd norm `N` becomes a unit in the binary local ring. -/
theorem norm_isUnit_odd (N : ℕ) (hodd : Odd N) (x : B N)
    (hx : x * star x = (N : B N)) : IsUnit (toLocal N x) := by
  apply norm_isUnit N hodd x (star x * x * star x)
  calc x * (star x * x * star x) = (x * star x) ^ 2 := by ring
    _ = (N : B N) ^ 2 := by rw [hx]

/-- The unit attached to one character value of a norm-`N` element. -/
def evalUnit (N : ℕ) [NeZero N] (hN : 0 < N) (hodd : Odd N)
    (x : Elem GaussianRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : GaussianRing))
    (S : Finset (PrimeIndex N)) : (O N)ˣ :=
  (norm_isUnit_odd N hodd (subsetEvaluation N hN S x)
    (subsetEvaluation_norm N hN S x hx)).unit

@[simp] theorem evalUnit_val (N : ℕ) [NeZero N] (hN : 0 < N) (hodd : Odd N)
    (x : Elem GaussianRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : GaussianRing))
    (S : Finset (PrimeIndex N)) :
    (evalUnit N hN hodd x hx S : O N) = toLocal N (subsetEvaluation N hN S x) :=
  IsUnit.unit_spec _

end LocalUnits

theorem gaussianToLocal_injective (N : ℕ) : Function.Injective (gaussianToLocal N) := by
  intro a b h
  apply Subtype.ext
  have h' : gaussianToB N a = gaussianToB N b := toLocal_injective N h
  have h2 := congrArg (fun z : B N => (z : ℂ)) h'
  simp only [coe_gaussianToB] at h2
  exact h2

theorem residue_gaussian_one_add_i (N : ℕ) :
    residue N (gaussianToLocal N (1 + gaussianI)) = 0 := by
  rw [map_add, map_one, map_add, map_one]
  change 1 + residue N (imaginary N) = 0
  rw [residue_imaginary]
  exact CharTwo.add_self_eq_zero 1

/-- Every multiple of `1 + i` reduces to zero in the binary residue field. -/
theorem residue_gaussian_of_mem (N : ℕ) (t r : GaussianRing)
    (ht : t = (1 + gaussianI) * r) : residue N (gaussianToLocal N t) = 0 := by
  rw [ht, map_mul, map_mul, residue_gaussian_one_add_i, zero_mul]

/-- **Master Lemma** (ℤ[i] version, λ = 1 + i). -/
theorem master_lemma (N : ℕ) [NeZero N] (hN : 1 < N) (hodd : Odd N)
    {ι : Type*} (s : Finset ι)
    (x₀ : Elem GaussianRing N)
    (hx₀ : x₀ * ringStar x₀ = scalar N ((N : ℕ) : GaussianRing))
    (x : ι → Elem GaussianRing N)
    (hx : ∀ j ∈ s, x j * ringStar (x j) = scalar N ((N : ℕ) : GaussianRing))
    (t : ι → GaussianRing) (ht0 : ∀ j ∈ s, t j ≠ 0)
    (htl : ∀ j ∈ s, ∃ r : GaussianRing, t j = (1 + gaussianI) * r)
    (y : ι → Elem GaussianRing N)
    (hxy : ∀ j ∈ s, x j = x₀ + scalar N (t j) * y j)
    (e : GaussianRing) (he : ∃ r : GaussianRing, e = 1 + (1 + gaussianI) * r)
    (W : Elem GaussianRing N)
    (hsum : ∑ j ∈ s, y j = scalar N e * Jelem N + scalar N (1 + gaussianI) * W) :
    False := by
  classical
  have hN0 : 0 < N := lt_trans Nat.zero_lt_one hN
  let φ := residue N
  let ev : Finset (PrimeIndex N) → Elem GaussianRing N →+* B N :=
    fun S => subsetEvaluation N hN0 S
  let c₀ : Finset (PrimeIndex N) → (O N)ˣ := evalUnit N hN0 hodd x₀ hx₀
  let L : ι → Finset (PrimeIndex N) → O N := fun j S =>
    toLocal N (ev S (y j)) * (↑((c₀ S)⁻¹) : O N)
  -- evaluation of a scalar
  have ev_scalar (S) (r : GaussianRing) :
      ev S (scalar N r) = gaussianToB N r := by
    simp only [ev, subsetEvaluation, evaluate_scalar]
    rfl
  -- first-order identity for each j
  have hfirst : ∀ j ∈ s, ∑ S : Finset (PrimeIndex N),
      (subsetSign S : k N) * φ (L j S) = 0 := by
    intro j hj
    let cj : Finset (PrimeIndex N) → (O N)ˣ := evalUnit N hN0 hodd (x j) (hx j hj)
    let T : O N := gaussianToLocal N (t j)
    have hT0 : T ≠ 0 := by
      intro h
      apply ht0 j hj
      apply gaussianToLocal_injective N
      simpa only [map_zero] using h
    obtain ⟨r, hr⟩ := htl j hj
    have hφT : φ T = 0 := residue_gaussian_of_mem N (t j) r hr
    have hrep : ∀ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        ((cj S * (c₀ S)⁻¹ : (O N)ˣ) : O N) = 1 + T * L j S := by
      intro S _
      have hval : (cj S : O N) = (c₀ S : O N) + T * toLocal N (ev S (y j)) := by
        simp only [cj, c₀, evalUnit_val]
        change toLocal N (ev S (x j)) = toLocal N (ev S x₀) + T * toLocal N (ev S (y j))
        rw [hxy j hj, map_add, map_mul, ev_scalar, map_add, map_mul]
        rfl
      rw [Units.val_mul, hval, add_mul, Units.mul_inv]
      simp only [L]
      ring
    have hodd' := local_alternatingProduct_ratio_has_odd_power N Finset.univ
      (fun S => ev S (x j)) (fun S => ev S x₀) cj c₀
      (fun S _ => evalUnit_val N hN0 hodd (x j) (hx j hj) S)
      (fun S _ => evalUnit_val N hN0 hodd x₀ hx₀ S)
      (alternatingProduct_odd_torsion N hN hodd (x j) (hx j hj))
      (alternatingProduct_odd_torsion N hN hodd x₀ hx₀)
    have hmapOne : ∀ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        φ ((cj S * (c₀ S)⁻¹ : (O N)ˣ) : O N) = 1 := by
      intro S hS
      rw [hrep S hS, map_add, map_one, map_mul, hφT, zero_mul, add_zero]
    have hOne : (∏ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        (cj S * (c₀ S)⁻¹) ^ (subsetSign S)) = 1 :=
      units_eq_one_of_residue_one_of_odd φ _
        (binary_signed_product_residue_one φ Finset.univ
          (fun S => cj S * (c₀ S)⁻¹) subsetSign hmapOne)
        (by simpa only [alternatingProduct, Finset.powerset_univ] using hodd')
    exact first_order_sum_eq_zero_of_prod_eq_one φ T hT0 hφT Finset.univ
      (fun S => cj S * (c₀ S)⁻¹) (L j) subsetSign hrep hOne
  -- sum over j
  have htotal : ∑ S : Finset (PrimeIndex N),
      (subsetSign S : k N) * φ (∑ j ∈ s, L j S) = 0 := by
    simp only [map_sum, Finset.mul_sum]
    rw [Finset.sum_comm]
    exact Finset.sum_eq_zero hfirst
  -- rewrite ∑_j L_j S using the hypothesis on ∑ y_j
  obtain ⟨re, hre⟩ := he
  have hE : φ (toLocal N (gaussianToB N e)) = 1 := by
    change residue N (gaussianToLocal N e) = 1
    rw [hre, map_add, map_one, map_add, map_one, residue_gaussian_of_mem N _ re rfl, add_zero]
  have hI : φ (toLocal N (gaussianToB N (1 + gaussianI))) = 0 :=
    residue_gaussian_one_add_i N
  have hsumL (S : Finset (PrimeIndex N)) : φ (∑ j ∈ s, L j S) =
      φ (toLocal N (ev S (Jelem N)) * (↑((c₀ S)⁻¹) : O N)) := by
    have h1 : (∑ j ∈ s, L j S) = toLocal N (ev S (∑ j ∈ s, y j)) * (↑((c₀ S)⁻¹) : O N) := by
      simp only [L, map_sum, Finset.sum_mul]
    have hev : ev S (∑ j ∈ s, y j) = gaussianToB N e * ev S (Jelem N) +
        gaussianToB N (1 + gaussianI) * ev S W := by
      rw [hsum, map_add, map_mul, map_mul, ev_scalar, ev_scalar]
    rw [h1, hev]
    generalize gaussianToB N (1 + gaussianI) = gI at hI ⊢
    generalize gaussianToB N e = gE at hE ⊢
    simp only [map_add, map_mul, hE, hI, one_mul, zero_mul, add_zero]
  have hzero : φ (∑ S : Finset (PrimeIndex N), (subsetSign S : O N) *
      toLocal N (ev S (Jelem N)) * (↑((c₀ S)⁻¹) : O N)) = 0 := by
    rw [← htotal, map_sum]
    apply Finset.sum_congr rfl
    intro S _
    rw [hsumL, map_mul, map_mul, map_intCast, mul_assoc, ← map_mul]
  exact subset_J_residue_ne_zero N hN0 (toLocal N) φ
    (fun S => (↑((c₀ S)⁻¹) : O N)) (c₀ ∅) rfl hodd hzero

end CirculantHadamard.GeneralNorm

end

end OAI
