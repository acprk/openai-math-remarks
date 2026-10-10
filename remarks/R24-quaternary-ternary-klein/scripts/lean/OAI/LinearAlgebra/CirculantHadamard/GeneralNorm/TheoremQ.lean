import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.MasterLemma
import Mathlib.Data.ZMod.QuotientRing

namespace OAI

/-!
# Theorem Q: no perfect quaternary sequence of length `2m`, `m > 1` odd

A sequence `a : ZMod (2m) → {1, -1, i, -i} ⊂ ℂ` is *perfect* if all its
nontrivial periodic autocorrelations vanish:
`∑_j a_j * conj (a_{j+τ}) = 0` for every `τ ≠ 0`.

Proof. Split `ZMod (2m) ≅ ZMod 2 × ZMod m` (CRT) and write `h = H₀ + Z H₁`.
Every quaternary unit is `≡ 1 (mod 1+i)`, so `c = (H₀+H₁)/(1+i)` and
`d = (H₀-H₁)/(1+i)` lie in `ℤ[i][C_m]`, satisfy `c c^* = d d^* = m`, and
`d = c + (i-1) H₁` with `H₁ = J + (1+i) W₁`. This is the Master Lemma with one
index, `t = i - 1 = (1+i) i` and `e = 1`.
-/

noncomputable section

namespace CirculantHadamard.GeneralNorm.QuaternaryQ

open scoped BigOperators
open CyclotomicRings CyclicRing

/-- The quaternary alphabet `{1, -1, i, -i}` inside `ℂ`. -/
def IsQuaternary (z : ℂ) : Prop :=
  z = 1 ∨ z = -1 ∨ z = Complex.I ∨ z = -Complex.I

theorem quaternary_lift {z : ℂ} (hz : IsQuaternary z) :
    ∃ w : GaussianRing, z = 1 + (1 + Complex.I) * (w : ℂ) := by
  rcases hz with rfl | rfl | rfl | rfl
  · exact ⟨0, by simp⟩
  · refine ⟨-1 + gaussianI, ?_⟩
    simp only [Subalgebra.coe_add, Subalgebra.coe_neg, Subalgebra.coe_one, coe_gaussianI]
    linear_combination (-1 : ℂ) * Complex.I_sq
  · refine ⟨gaussianI, ?_⟩
    simp only [coe_gaussianI]
    linear_combination (-1 : ℂ) * Complex.I_sq
  · exact ⟨-1, by simp⟩

theorem quaternary_mul_conj {z : ℂ} (hz : IsQuaternary z) :
    z * (starRingEnd ℂ) z = 1 := by
  rcases hz with rfl | rfl | rfl | rfl <;> simp

/-- Coefficientwise description of a scalar multiple in the cyclic group ring. -/
theorem scalar_mul_coeff {R : Type*} [Semiring R] {n : ℕ} [NeZero n]
    (r : R) (f : Elem R n) (a : ZMod n) :
    (scalar n r * f).coeff a = r * f.coeff a := by
  classical
  rw [mul_apply]
  simp only [scalar_apply, ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ,
    ite_true, sub_zero]

/-- If `(1+i) X = P` coefficientwise and `P` has autocorrelation `2m δ`, then
`X X^* = m`. -/
theorem norm_of_scaled (m : ℕ) [NeZero m] (X : Elem GaussianRing m) (P : ZMod m → ℂ)
    (hP : ∀ i, (1 + Complex.I) * ((X.coeff i : GaussianRing) : ℂ) = P i)
    (hauto : ∀ k : ZMod m, ∑ i, P i * (starRingEnd ℂ) (P (i - k)) =
      if k = 0 then (2 * m : ℂ) else 0) :
    X * ringStar X = scalar m ((m : ℕ) : GaussianRing) := by
  ext k
  rw [mul_ringStar_apply, scalar_apply]
  have hsum : (2 : ℂ) * ∑ i, ((X.coeff i : GaussianRing) : ℂ) *
      (starRingEnd ℂ) ((X.coeff (i - k) : GaussianRing) : ℂ) =
      ∑ i, P i * (starRingEnd ℂ) (P (i - k)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    rw [← hP i, ← hP (i - k), map_mul]
    have h2 : (1 + Complex.I) * (starRingEnd ℂ) (1 + Complex.I) = 2 := by
      rw [map_add, map_one, Complex.conj_I]
      linear_combination (-1 : ℂ) * Complex.I_sq
    linear_combination (((X.coeff i : GaussianRing) : ℂ) *
      (starRingEnd ℂ) ((X.coeff (i - k) : GaussianRing) : ℂ)) * h2.symm
  rw [hauto k] at hsum
  have hcoe : ((∑ j, X.coeff j * star (X.coeff (j - k)) : GaussianRing) : ℂ) =
      ∑ i, ((X.coeff i : GaussianRing) : ℂ) *
        (starRingEnd ℂ) ((X.coeff (i - k) : GaussianRing) : ℂ) := by
    push_cast
    rfl
  rw [hcoe]
  by_cases hk : k = 0
  · simp only [hk, ite_true] at hsum ⊢
    simp only [SubringClass.coe_natCast]
    linear_combination hsum / 2
  · simp only [hk, ite_false] at hsum ⊢
    simp only [ZeroMemClass.coe_zero]
    linear_combination hsum / 2

section Sequence

variable (m : ℕ) (hodd : Odd m)

theorem coprime_two (hodd : Odd m) : Nat.Coprime 2 m :=
  (Nat.coprime_two_left).mpr hodd

/-- The CRT splitting `ZMod (2m) ≃+* ZMod 2 × ZMod m`. -/
def crt2 : ZMod (2 * m) ≃+* ZMod 2 × ZMod m :=
  ZMod.chineseRemainder (coprime_two m hodd)

/-- Autocorrelation of a split sequence. -/
theorem autocorr_split [NeZero m] (a : ZMod (2 * m) → ℂ) (s : ZMod 2) (k : ZMod m) :
    ∑ j, a j * (starRingEnd ℂ) (a (j + (crt2 m hodd).symm (s, -k))) =
      ∑ i, (a ((crt2 m hodd).symm (0, i)) *
          (starRingEnd ℂ) (a ((crt2 m hodd).symm (0 + s, i - k))) +
        a ((crt2 m hodd).symm (1, i)) *
          (starRingEnd ℂ) (a ((crt2 m hodd).symm (1 + s, i - k)))) := by
  let e := crt2 m hodd
  rw [← e.symm.toEquiv.sum_comp, Fintype.sum_prod_type]
  rw [Finset.sum_add_distrib]
  change (∑ f : Fin 2, _) = _
  rw [Fin.sum_univ_two]
  congr 1
  · apply Finset.sum_congr rfl
    intro i _
    congr 3
    change e.symm (0, i) + e.symm (s, -k) = e.symm (0 + s, i - k)
    rw [← map_add]
    congr 1
    ext <;> simp [sub_eq_add_neg]
  · apply Finset.sum_congr rfl
    intro i _
    congr 3
    change e.symm (1, i) + e.symm (s, -k) = e.symm (1 + s, i - k)
    rw [← map_add]
    congr 1
    ext <;> simp [sub_eq_add_neg]

end Sequence

/-- **Theorem Q.** For odd `m > 1` there is no perfect quaternary sequence of
length `2m`: no `a : ZMod (2m) → {1, -1, i, -i}` has all nontrivial periodic
autocorrelations equal to zero. -/
theorem not_exists_perfect_quaternary (m : ℕ) [NeZero m] (hm : 1 < m) (hodd : Odd m) :
    ¬ ∃ a : ZMod (2 * m) → ℂ,
      (∀ j, a j = 1 ∨ a j = -1 ∨ a j = Complex.I ∨ a j = -Complex.I) ∧
      ∀ τ : ZMod (2 * m), τ ≠ 0 →
        ∑ j, a j * (starRingEnd ℂ) (a (j + τ)) = 0 := by
  classical
  rintro ⟨a, hq, hauto⟩
  let e := crt2 m hodd
  -- lifts of the entries
  let w : ZMod (2 * m) → GaussianRing := fun j => Classical.choose (quaternary_lift (hq j))
  have hw (j) : a j = 1 + (1 + Complex.I) * (w j : ℂ) :=
    Classical.choose_spec (quaternary_lift (hq j))
  let W : ZMod 2 → ZMod m → GaussianRing := fun f i => w (e.symm (f, i))
  let A : ZMod 2 → ZMod m → ℂ := fun f i => a (e.symm (f, i))
  have hA (f i) : A f i = 1 + (1 + Complex.I) * (W f i : ℂ) := hw _
  -- the values of the split autocorrelations
  have hτ (s : ZMod 2) (k : ZMod m) (hsk : ¬ (s = 0 ∧ k = 0)) :
      e.symm (s, -k) ≠ 0 := by
    intro h
    apply hsk
    have h' := congrArg e h
    rw [RingEquiv.apply_symm_apply, map_zero] at h'
    have h1 := congrArg Prod.fst h'
    have h2 := congrArg Prod.snd h'
    simp only [Prod.fst_zero, Prod.snd_zero, neg_eq_zero] at h1 h2
    exact ⟨h1, h2⟩
  have hA0 (k : ZMod m) : ∑ i, (A 0 i * (starRingEnd ℂ) (A 0 (i - k)) +
      A 1 i * (starRingEnd ℂ) (A 1 (i - k))) = if k = 0 then (2 * m : ℂ) else 0 := by
    have h := autocorr_split m hodd a 0 k
    simp only [add_zero] at h
    change _ = ∑ i, (A 0 i * (starRingEnd ℂ) (A 0 (i - k)) +
      A 1 i * (starRingEnd ℂ) (A 1 (i - k))) at h
    rw [← h]
    by_cases hk : k = 0
    · subst hk
      simp only [neg_zero, ite_true]
      have he0 : e.symm (0, 0) = 0 := by
        rw [← Prod.zero_eq_mk, map_zero]
      rw [he0]
      simp only [add_zero, quaternary_mul_conj (hq _), Finset.sum_const, Finset.card_univ,
        ZMod.card, nsmul_eq_mul, mul_one]
      push_cast
      ring
    · simp only [hk, ↓reduceIte]
      exact hauto _ (hτ 0 k (fun h => hk h.2))
  have hA1 (k : ZMod m) : ∑ i, (A 0 i * (starRingEnd ℂ) (A 1 (i - k)) +
      A 1 i * (starRingEnd ℂ) (A 0 (i - k))) = 0 := by
    have h := autocorr_split m hodd a 1 k
    have h11 : (1 : ZMod 2) + 1 = 0 := by decide
    simp only [zero_add, h11] at h
    change _ = ∑ i, (A 0 i * (starRingEnd ℂ) (A 1 (i - k)) +
      A 1 i * (starRingEnd ℂ) (A 0 (i - k))) at h
    rw [← h]
    exact hauto _ (hτ 1 k (fun h => absurd h.1 (by decide)))
  -- the group-ring elements
  let H₁ : Elem GaussianRing m := ofCoeffs fun i => 1 + (1 + gaussianI) * W 1 i
  let W₁ : Elem GaussianRing m := ofCoeffs fun i => W 1 i
  let c : Elem GaussianRing m := ofCoeffs fun i => W 0 i + W 1 i + (1 - gaussianI)
  let d : Elem GaussianRing m := ofCoeffs fun i => W 0 i - W 1 i
  have hc : c * ringStar c = scalar m ((m : ℕ) : GaussianRing) := by
    apply norm_of_scaled m c (fun i => A 0 i + A 1 i)
    · intro i
      simp only [c, ofCoeffs_apply, hA, Subalgebra.coe_add, Subalgebra.coe_sub,
        Subalgebra.coe_one, coe_gaussianI]
      linear_combination (-1 : ℂ) * Complex.I_sq
    · intro k
      rw [← hA0 k, ← add_zero (∑ i, (A 0 i * (starRingEnd ℂ) (A 0 (i - k)) +
        A 1 i * (starRingEnd ℂ) (A 1 (i - k)))), ← hA1 k, ← Finset.sum_add_distrib]
      apply Finset.sum_congr rfl
      intro i _
      rw [map_add]
      ring
  have hd : d * ringStar d = scalar m ((m : ℕ) : GaussianRing) := by
    apply norm_of_scaled m d (fun i => A 0 i - A 1 i)
    · intro i
      simp only [d, ofCoeffs_apply, hA, Subalgebra.coe_sub]
      ring
    · intro k
      rw [← hA0 k, ← sub_zero (∑ i, (A 0 i * (starRingEnd ℂ) (A 0 (i - k)) +
        A 1 i * (starRingEnd ℂ) (A 1 (i - k)))), ← hA1 k, ← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro i _
      rw [map_sub]
      ring
  have hcd : d = c + scalar m (gaussianI - 1) * H₁ := by
    ext i
    rw [AddMonoidAlgebra.coeff_add, Finsupp.add_apply, scalar_mul_coeff]
    simp only [c, d, H₁, ofCoeffs_apply]
    simp only [Subalgebra.coe_add, Subalgebra.coe_sub, Subalgebra.coe_mul,
      Subalgebra.coe_one, coe_gaussianI]
    linear_combination (-((W 1 i : GaussianRing) : ℂ)) * Complex.I_sq
  have hsum : ∑ _j ∈ (Finset.univ : Finset Unit), H₁ =
      scalar m (1 : GaussianRing) * Jelem m + scalar m (1 + gaussianI) * W₁ := by
    ext i
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_unit, one_smul,
      AddMonoidAlgebra.coeff_add, Finsupp.add_apply, scalar_mul_coeff, scalar_mul_coeff]
    simp only [H₁, W₁, ofCoeffs_apply, BinaryCoefficients.JOver_apply, one_mul]
  exact master_lemma m hm hodd (Finset.univ : Finset Unit) c hc (fun _ => d)
    (fun _ _ => hd) (fun _ => gaussianI - 1)
    (fun _ _ h => by
      have h' := congrArg (fun z : GaussianRing => (z : ℂ)) h
      simp only [Subalgebra.coe_sub, coe_gaussianI, Subalgebra.coe_one,
        ZeroMemClass.coe_zero] at h'
      have : Complex.I = 1 := sub_eq_zero.mp h'
      simpa using congrArg Complex.re this)
    (fun _ _ => ⟨gaussianI, by
      apply Subtype.ext
      simp only [Subalgebra.coe_sub, Subalgebra.coe_mul, Subalgebra.coe_add,
        Subalgebra.coe_one, coe_gaussianI]
      linear_combination (-1 : ℂ) * Complex.I_sq⟩)
    (fun _ => H₁) (fun _ _ => hcd) 1 ⟨0, by simp⟩ W₁ hsum

/-- Sanity check: the hypothesis `1 < m` cannot be dropped. For `m = 1`
(`ZMod (2 * 1) = ZMod 2`) the sequence `(1, i)` of length `2` is perfect. -/
theorem exists_perfect_quaternary_length_two :
    ∃ a : ZMod 2 → ℂ,
      (∀ j, a j = 1 ∨ a j = -1 ∨ a j = Complex.I ∨ a j = -Complex.I) ∧
      ∀ τ : ZMod 2, τ ≠ 0 →
        ∑ j, a j * (starRingEnd ℂ) (a (j + τ)) = 0 := by
  refine ⟨![1, Complex.I], ?_, ?_⟩
  · intro j
    fin_cases j <;> simp
  · intro τ hτ
    fin_cases τ
    · exact absurd rfl hτ
    · change ∑ x : Fin 2, ![1, Complex.I] x * (starRingEnd ℂ) (![1, Complex.I] (x + 1)) = 0
      simp [Fin.sum_univ_two]

end CirculantHadamard.GeneralNorm.QuaternaryQ

end

end OAI
