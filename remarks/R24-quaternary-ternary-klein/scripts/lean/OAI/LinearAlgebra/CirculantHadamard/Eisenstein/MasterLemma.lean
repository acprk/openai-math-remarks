import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.AlternatingTorsion
import OAI.LinearAlgebra.CirculantHadamard.BinaryCharacterSum
import OAI.LinearAlgebra.CirculantHadamard.FirstOrderResidue

namespace OAI

/-!
# The λ-adic Master Lemma over `ℤ[ω]` (λ = 1 - ω, ℓ = 3)

Let `N > 1` with `3 ∤ N` and let `x₀, x_j (j ∈ s)` be elements of `ℤ[ω][C_N]` with
`x * x^* = N`. Suppose `x_j - x₀ = t_j y_j` with `0 ≠ t_j ∈ (1-ω)ℤ[ω]`, and
`∑_j y_j = e J + (1-ω) W`, where `e ≡ c (mod 1-ω)` for an integer `c` with `3 ∤ c`.
Then we have a contradiction.

The proof is the ternary analogue of `GeneralNorm.master_lemma`: work in the
localization `O` of `B = ℤ[ω, ρ_p : p ∣ N]` at a maximal ideal over `3` (residue
characteristic `3`, `ω ↦ 1`); by Prop″ (`EisensteinNorm.alternatingProduct_coprime_three_torsion`)
each `Δ(x_j)/Δ(x₀)` has torsion order prime to `3` and residue `1`, hence equals `1`;
the first-order homomorphism along `t_j` and the character sum of `J` finish.
-/

noncomputable section

namespace CirculantHadamard.EisensteinNorm

open scoped BigOperators
open CyclotomicRingsW PrimeComponentsN PrimeCharacterEvaluationsW CyclicRing
open AlternatingLocalization AlternatingLocalizationW

/-- The all-ones element `J` of `ℤ[ω][ZMod N]`. -/
abbrev Jelem (N : ℕ) [NeZero N] : Elem EisensteinRing N :=
  AddMonoidAlgebra.ofCoeff (BinaryCoefficients.JOver EisensteinRing)

/-- `λ = 1 - ω`. -/
abbrev lam : EisensteinRing := 1 - omegaE

section Generic

universe uO uk uI

/-- A signed product of residue-one units has residue one. -/
theorem signed_product_residue_one {O : Type uO} {k : Type uk} {ι : Type uI}
    [CommRing O] [Field k] (φ : O →+* k) (s : Finset ι) (u : ι → Oˣ) (ε : ι → ℤ)
    (hu : ∀ j ∈ s, φ (u j : O) = 1) :
    φ ((∏ j ∈ s, (u j) ^ (ε j) : Oˣ) : O) = 1 := by
  let f : Oˣ →* kˣ := Units.map φ.toMonoidHom
  have hp : f (∏ j ∈ s, (u j) ^ (ε j)) = 1 := by
    rw [map_prod]
    apply Finset.prod_eq_one
    intro j hj
    rw [map_zpow]
    have : f (u j) = 1 := Units.ext (hu j hj)
    rw [this, one_zpow]
  exact congrArg (fun z : kˣ => (z : k)) hp

/-- A unit with residue one and torsion order prime to the residue
characteristic `3` is one. -/
theorem units_eq_one_of_residue_one_of_coprime_three {O : Type uO} {k : Type uk}
    [CommRing O] [IsDomain O] [Field k] [CharP k 3]
    (f : O →+* k) (u : Oˣ) (hmap : f (u : O) = 1)
    (hfinite : ∃ m : ℕ, 0 < m ∧ m.Coprime 3 ∧ u ^ m = 1) : u = 1 := by
  obtain ⟨m, _, hm, hu⟩ := hfinite
  have hpow : (u : O) ^ m = 1 := by
    simpa only [Units.val_pow_eq_pow_val, Units.val_one] using
      congrArg (fun v : Oˣ => (v : O)) hu
  apply Units.ext
  exact eq_one_of_residue_one_of_coprime f (u : O) Nat.prime_three hmap hpow hm

end Generic

variable (N : ℕ) [NeZero N]

/-- A maximal ideal of `B N` over `3`. -/
def threeIdeal : MaximalSpectrum (B N) :=
  ⟨Classical.choose (exists_maximal_ideal_over_nat_prime_of_charZero (B N) Nat.prime_three),
    (Classical.choose_spec
      (exists_maximal_ideal_over_nat_prime_of_charZero (B N) Nat.prime_three)).1⟩

omit [NeZero N] in
theorem three_mem : ((3 : ℕ) : B N) ∈ (threeIdeal N).asIdeal :=
  (Classical.choose_spec
    (exists_maximal_ideal_over_nat_prime_of_charZero (B N) Nat.prime_three)).2.2

/-- The ternary local ring. -/
abbrev O3 := AtMaximal (B N) (threeIdeal N)

abbrev k3 := IsLocalRing.ResidueField (O3 N)

instance k3_charThree : CharP (k3 N) 3 :=
  residue_charP (threeIdeal N) 3 Nat.prime_three (three_mem N)

def toLocal : B N →+* O3 N := algebraMap (B N) (O3 N)

def res : O3 N →+* k3 N := IsLocalRing.residue (O3 N)

omit [NeZero N] in
theorem toLocal_injective : Function.Injective (toLocal N) :=
  IsLocalization.injective (O3 N) (threeIdeal N).asIdeal.primeCompl_le_nonZeroDivisors

def eisToLocal : EisensteinRing →+* O3 N := (toLocal N).comp (eisensteinToB N).toRingHom

theorem eisToLocal_injective : Function.Injective (eisToLocal N) := by
  intro a b h
  apply Subtype.ext
  have h' : eisensteinToB N a = eisensteinToB N b := toLocal_injective N h
  have h2 := congrArg (fun z : B N => (z : ℂ)) h'
  simp only [coe_eisensteinToB] at h2
  exact h2

theorem res_omega : res N (eisToLocal N omegaE) = 1 := by
  have h := residue_root_eq_one (threeIdeal N) 3 Nat.prime_three (three_mem N)
    (eisensteinToB N omegaE)
    (by rw [← map_pow]; exact Subtype.ext (by simpa using omega_pow_three))
  exact h

theorem res_lam : res N (eisToLocal N lam) = 0 := by
  rw [map_sub, map_sub, map_one, map_one, res_omega, sub_self]

theorem res_of_mem (t r : EisensteinRing) (ht : t = lam * r) :
    res N (eisToLocal N t) = 0 := by
  rw [ht, map_mul, map_mul, res_lam, zero_mul]

theorem res_intCast_ne_zero (c : ℤ) (hc : ¬ (3 : ℤ) ∣ c) :
    res N (eisToLocal N (c : EisensteinRing)) ≠ 0 := by
  rw [map_intCast, map_intCast]
  intro h
  exact hc ((CharP.intCast_eq_zero_iff (k3 N) 3 c).mp h)

theorem res_natCast_ne_zero (n : ℕ) (hn : ¬ 3 ∣ n) : res N ((n : ℕ) : O3 N) ≠ 0 := by
  rw [map_natCast]
  intro h
  exact hn ((CharP.cast_eq_zero_iff (k3 N) 3 n).mp h)

/-- A character value of a norm-`N` element is a local unit when `3 ∤ N`. -/
theorem norm_isUnit (h3 : ¬ 3 ∣ N) (x : B N) (hx : x * star x = (N : B N)) :
    IsUnit (toLocal N x) := by
  apply (IsLocalRing.residue_ne_zero_iff_isUnit _).mp
  intro h0
  have h := congrArg (fun z => res N (toLocal N z)) hx
  simp only [map_mul, map_natCast] at h
  change IsLocalRing.residue (O3 N) (toLocal N x) * _ = _ at h
  rw [h0, zero_mul] at h
  exact res_natCast_ne_zero N N h3 (by rw [map_natCast]; exact h.symm)

def evalUnit (hN : 0 < N) (h3 : ¬ 3 ∣ N) (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (S : Finset (PrimeIndex N)) : (O3 N)ˣ :=
  (norm_isUnit N h3 (subsetEvaluation N hN S x) (subsetEvaluation_norm N hN S x hx)).unit

@[simp] theorem evalUnit_val (hN : 0 < N) (h3 : ¬ 3 ∣ N) (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (S : Finset (PrimeIndex N)) :
    (evalUnit N hN h3 x hx S : O3 N) = toLocal N (subsetEvaluation N hN S x) :=
  IsUnit.unit_spec _

/-- Transport of torsion from the fraction field to the local units. -/
theorem local_alternatingProduct_torsion (hN : 1 < N) (h3 : ¬ 3 ∣ N)
    (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing)) :
    ∃ m : ℕ, 0 < m ∧ m.Coprime 3 ∧
      (alternatingProduct (Finset.univ : Finset (PrimeIndex N))
        (evalUnit N (lt_trans Nat.zero_lt_one hN) h3 x hx)) ^ m = 1 := by
  obtain ⟨m, hm, hmc, hpow⟩ := alternatingProduct_coprime_three_torsion N hN h3 x hx
  refine ⟨m, hm, hmc, ?_⟩
  let F : (O3 N)ˣ →* FractionRing (B N) :=
    (toFraction (FractionRing (B N)) (threeIdeal N)).toRingHom.toMonoidHom.comp
      (Units.coeHom (O3 N))
  have hF : Function.Injective F := by
    intro X Y h
    apply Units.ext
    exact toFraction_injective (FractionRing (B N)) (threeIdeal N) h
  apply hF
  rw [map_pow, map_one, map_alternatingProduct]
  convert hpow using 2
  apply alternatingProduct_congr
  intro S _
  change toFraction (FractionRing (B N)) (threeIdeal N)
    ((evalUnit N _ h3 x hx S : O3 N)) = _
  rw [evalUnit_val]
  exact AlgHom.commutes _ _

omit [NeZero N] in
theorem subsetCharacter_nontrivial (hN : 0 < N) (S : Finset (PrimeIndex N)) (hS : S.Nonempty) :
    ∃ g : Multiplicative (ZMod N), subsetCharacter N hN S g ≠ 1 := by
  obtain ⟨a, ha⟩ := subsetCharacter_nonempty N hN S hS
  exact ⟨Multiplicative.ofAdd a, ha⟩

set_option maxHeartbeats 1600000 in
/-- **Master Lemma** (ℤ[ω] version, λ = 1 - ω). -/
theorem master_lemma (hN : 1 < N) (h3 : ¬ 3 ∣ N)
    {ι : Type*} (s : Finset ι)
    (x₀ : Elem EisensteinRing N)
    (hx₀ : x₀ * ringStar x₀ = scalar N ((N : ℕ) : EisensteinRing))
    (x : ι → Elem EisensteinRing N)
    (hx : ∀ j ∈ s, x j * ringStar (x j) = scalar N ((N : ℕ) : EisensteinRing))
    (t : ι → EisensteinRing) (ht0 : ∀ j ∈ s, t j ≠ 0)
    (htl : ∀ j ∈ s, ∃ r : EisensteinRing, t j = lam * r)
    (y : ι → Elem EisensteinRing N)
    (hxy : ∀ j ∈ s, x j = x₀ + scalar N (t j) * y j)
    (e : EisensteinRing)
    (he : ∃ (c : ℤ) (r : EisensteinRing), e = (c : EisensteinRing) + lam * r ∧ ¬ (3 : ℤ) ∣ c)
    (W : Elem EisensteinRing N)
    (hsum : ∑ j ∈ s, y j = scalar N e * Jelem N + scalar N lam * W) :
    False := by
  classical
  have hN0 : 0 < N := lt_trans Nat.zero_lt_one hN
  let φ := res N
  let ev : Finset (PrimeIndex N) → Elem EisensteinRing N →+* B N :=
    fun S => subsetEvaluation N hN0 S
  let c₀ : Finset (PrimeIndex N) → (O3 N)ˣ := evalUnit N hN0 h3 x₀ hx₀
  let L : ι → Finset (PrimeIndex N) → O3 N := fun j S =>
    toLocal N (ev S (y j)) * (↑((c₀ S)⁻¹) : O3 N)
  have ev_scalar (S) (r : EisensteinRing) : ev S (scalar N r) = eisensteinToB N r := by
    simp only [ev, subsetEvaluation, evaluate_scalar]
    rfl
  have hfirst : ∀ j ∈ s, ∑ S : Finset (PrimeIndex N),
      (subsetSign S : k3 N) * φ (L j S) = 0 := by
    intro j hj
    let cj : Finset (PrimeIndex N) → (O3 N)ˣ := evalUnit N hN0 h3 (x j) (hx j hj)
    let T : O3 N := eisToLocal N (t j)
    have hT0 : T ≠ 0 := by
      intro h
      apply ht0 j hj
      apply eisToLocal_injective N
      simpa only [map_zero] using h
    obtain ⟨r, hr⟩ := htl j hj
    have hφT : φ T = 0 := res_of_mem N (t j) r hr
    have hrep : ∀ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        ((cj S * (c₀ S)⁻¹ : (O3 N)ˣ) : O3 N) = 1 + T * L j S := by
      intro S _
      have hval : (cj S : O3 N) = (c₀ S : O3 N) + T * toLocal N (ev S (y j)) := by
        simp only [cj, c₀, evalUnit_val]
        change toLocal N (ev S (x j)) = toLocal N (ev S x₀) + T * toLocal N (ev S (y j))
        rw [hxy j hj, map_add, map_mul, ev_scalar, map_add, map_mul]
        rfl
      rw [Units.val_mul, hval, add_mul, Units.mul_inv]
      simp only [L]
      ring
    have htor : ∃ m : ℕ, 0 < m ∧ m.Coprime 3 ∧
        (alternatingProduct (Finset.univ : Finset (PrimeIndex N))
          (fun S => cj S * (c₀ S)⁻¹)) ^ m = 1 := by
      obtain ⟨a, ha, hac, hapow⟩ := local_alternatingProduct_torsion N hN h3 (x j) (hx j hj)
      obtain ⟨b, hb, hbc, hbpow⟩ := local_alternatingProduct_torsion N hN h3 x₀ hx₀
      refine ⟨a * b, Nat.mul_pos ha hb, Nat.coprime_mul_iff_left.mpr ⟨hac, hbc⟩, ?_⟩
      rw [alternatingProduct_mul, alternatingProduct_inv, mul_pow, inv_pow]
      change (alternatingProduct Finset.univ cj) ^ (a * b) *
        ((alternatingProduct Finset.univ c₀) ^ (a * b))⁻¹ = 1
      rw [pow_mul, hapow, one_pow, Nat.mul_comm, pow_mul, hbpow, one_pow, inv_one, mul_one]
    have hmapOne : ∀ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        φ ((cj S * (c₀ S)⁻¹ : (O3 N)ˣ) : O3 N) = 1 := by
      intro S hS
      rw [hrep S hS, map_add, map_one, map_mul, hφT, zero_mul, add_zero]
    have hOne : (∏ S ∈ (Finset.univ : Finset (Finset (PrimeIndex N))),
        (cj S * (c₀ S)⁻¹) ^ (subsetSign S)) = 1 :=
      units_eq_one_of_residue_one_of_coprime_three φ _
        (signed_product_residue_one φ Finset.univ
          (fun S => cj S * (c₀ S)⁻¹) subsetSign hmapOne)
        (by simpa only [alternatingProduct, Finset.powerset_univ] using htor)
    exact first_order_sum_eq_zero_of_prod_eq_one φ T hT0 hφT Finset.univ
      (fun S => cj S * (c₀ S)⁻¹) (L j) subsetSign hrep hOne
  have htotal : ∑ S : Finset (PrimeIndex N),
      (subsetSign S : k3 N) * φ (∑ j ∈ s, L j S) = 0 := by
    simp only [map_sum, Finset.mul_sum]
    rw [Finset.sum_comm]
    exact Finset.sum_eq_zero hfirst
  obtain ⟨ce, re, hre, hce⟩ := he
  have hE : φ (toLocal N (eisensteinToB N e)) = φ (eisToLocal N (ce : EisensteinRing)) := by
    change φ (eisToLocal N e) = _
    rw [hre, map_add, map_add, res_of_mem N (lam * re) re rfl, add_zero]
  have hI : φ (toLocal N (eisensteinToB N lam)) = 0 := res_lam N
  have hsumL (S : Finset (PrimeIndex N)) : φ (∑ j ∈ s, L j S) =
      φ (eisToLocal N (ce : EisensteinRing)) *
        φ (toLocal N (ev S (Jelem N)) * (↑((c₀ S)⁻¹) : O3 N)) := by
    have h1 : (∑ j ∈ s, L j S) =
        toLocal N (ev S (∑ j ∈ s, y j)) * (↑((c₀ S)⁻¹) : O3 N) := by
      simp only [L, map_sum, Finset.sum_mul]
    have hev : ev S (∑ j ∈ s, y j) = eisensteinToB N e * ev S (Jelem N) +
        eisensteinToB N lam * ev S W := by
      rw [hsum, map_add, map_mul, map_mul, ev_scalar, ev_scalar]
    rw [h1, hev]
    generalize eisensteinToB N lam = gI at hI ⊢
    generalize eisensteinToB N e = gE at hE ⊢
    rw [map_add (toLocal N), map_mul (toLocal N), map_mul (toLocal N), add_mul, map_add φ,
      map_mul φ, map_mul φ, map_mul φ, map_mul φ, hE, hI, zero_mul, zero_mul, add_zero,
      map_mul φ (toLocal N (ev S (Jelem N))), mul_assoc]
  have hzero : φ (eisToLocal N (ce : EisensteinRing)) *
      φ (∑ S : Finset (PrimeIndex N), (subsetSign S : O3 N) *
        toLocal N (ev S (Jelem N)) * (↑((c₀ S)⁻¹) : O3 N)) = 0 := by
    rw [← htotal, map_sum, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro S _
    rw [hsumL, map_mul φ, map_mul φ, map_intCast φ (subsetSign S),
      map_mul φ (toLocal N (ev S (Jelem N)))]
    ring
  have hJ : (∑ S : Finset (PrimeIndex N), (subsetSign S : O3 N) *
      toLocal N (ev S (Jelem N)) * (↑((c₀ S)⁻¹) : O3 N)) =
      ((N : ℕ) : O3 N) * (↑((c₀ ∅)⁻¹) : O3 N) :=
    BinaryCharacterSum.subset_map_evaluate_J_mul (eisensteinToB N).toRingHom (toLocal N)
      (subsetCharacter N hN0) (fun S => (subsetSign S : O3 N))
      (fun S => (↑((c₀ S)⁻¹) : O3 N)) (by simp)
      (subsetCharacter_empty N hN0) (subsetCharacter_nontrivial N hN0)
  rw [hJ, map_mul] at hzero
  rcases mul_eq_zero.mp hzero with h | h
  · exact res_intCast_ne_zero N ce hce h
  rcases mul_eq_zero.mp h with h | h
  · exact res_natCast_ne_zero N N h3 h
  · exact (IsUnit.map φ (c₀ ∅)⁻¹.isUnit).ne_zero h

end CirculantHadamard.EisensteinNorm

end

end OAI
