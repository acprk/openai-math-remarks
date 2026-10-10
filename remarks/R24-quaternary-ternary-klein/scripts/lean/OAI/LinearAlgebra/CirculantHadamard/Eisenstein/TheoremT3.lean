import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.MasterLemma
import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremQ
import Mathlib.Data.ZMod.QuotientRing

namespace OAI

/-!
# Theorem T3: no perfect ternary sequence of length `3m`, `m > 1`, `3 ∤ m`

A sequence `a : ZMod (3m) → ℂ` with `a_j ^ 3 = 1` is *perfect* if
`∑_j a_j * conj (a_{j+τ}) = 0` for every `τ ≠ 0`.

Proof. Split `ZMod (3m) ≅ ZMod 3 × ZMod m` and write `h = H₀ + Z H₁ + Z² H₂`.
Every cube root of unity is `≡ 1 (mod λ)`, `λ = 1 - ω`, so `c_r = h(r)/λ`
(`r = 1, ω, ω²`) lie in `ℤ[ω][C_m]` and have norm `m` (`λ λ̄ = 3`). With `x₀ = c_1`,
`c_ω - c_1 = λ y₁`, `c_{ω²} - c_1 = λ y₂` and `y₁ + y₂ = 2ω² J + λ ω² (W₁ + W₂)`;
the ℤ[ω] Master Lemma applies with `e = 2ω² ≡ 2 (mod λ)`.
-/

noncomputable section

namespace CirculantHadamard.EisensteinNorm.TernaryT3

open scoped BigOperators
open CyclotomicRingsW CyclicRing

theorem omega_ne_one : omega ≠ 1 := omega_isPrimitiveRoot.ne_one (by norm_num)

theorem conj_omega : (starRingEnd ℂ) omega = -1 - omega := by
  change star omega = _
  rw [star_omega, omega_sq]

theorem cube_root_cases {z : ℂ} (hz : z ^ 3 = 1) : z = 1 ∨ z = omega ∨ z = omega ^ 2 := by
  obtain ⟨i, hi, rfl⟩ := omega_isPrimitiveRoot.eq_pow_of_pow_eq_one hz
  interval_cases i <;> simp

theorem cube_root_lift {z : ℂ} (hz : z ^ 3 = 1) :
    ∃ w : EisensteinRing, z = 1 + (1 - omega) * (w : ℂ) := by
  rcases cube_root_cases hz with rfl | rfl | rfl
  · exact ⟨0, by simp⟩
  · exact ⟨-1, by simp⟩
  · refine ⟨-1 - omegaE, ?_⟩
    simp only [Subalgebra.coe_sub, Subalgebra.coe_neg, Subalgebra.coe_one, coe_omegaE]
    ring

theorem cube_root_mul_conj {z : ℂ} (hz : z ^ 3 = 1) : z * (starRingEnd ℂ) z = 1 := by
  rcases cube_root_cases hz with rfl | rfl | rfl
  · simp
  · rw [conj_omega]; linear_combination (-1 : ℂ) * omega_sq
  · rw [map_pow, conj_omega]
    linear_combination (omega ^ 2 + omega - 1) * omega_sq

theorem lam_mul_conj : (1 - omega) * (starRingEnd ℂ) (1 - omega) = 3 := by
  rw [map_sub, map_one, conj_omega]
  linear_combination (-1 : ℂ) * omega_sq

theorem sum_zmod_three {M : Type*} [AddCommMonoid M] (F : ZMod 3 → M) :
    ∑ a, F a = F 0 + F 1 + F 2 :=
  Fin.sum_univ_three F

/-- If `λ X = P` coefficientwise and `P` has autocorrelation `3m δ`, then `X X^* = m`. -/
theorem norm_of_lam (m : ℕ) [NeZero m] (X : Elem EisensteinRing m) (P : ZMod m → ℂ)
    (hP : ∀ i, (1 - omega) * ((X.coeff i : EisensteinRing) : ℂ) = P i)
    (hauto : ∀ k : ZMod m, ∑ i, P i * (starRingEnd ℂ) (P (i - k)) =
      if k = 0 then (3 * m : ℂ) else 0) :
    X * ringStar X = scalar m ((m : ℕ) : EisensteinRing) := by
  ext k
  rw [mul_ringStar_apply, scalar_apply]
  have hsum : (3 : ℂ) * ∑ i, ((X.coeff i : EisensteinRing) : ℂ) *
      (starRingEnd ℂ) ((X.coeff (i - k) : EisensteinRing) : ℂ) =
      ∑ i, P i * (starRingEnd ℂ) (P (i - k)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    rw [← hP i, ← hP (i - k), map_mul]
    linear_combination (((X.coeff i : EisensteinRing) : ℂ) *
      (starRingEnd ℂ) ((X.coeff (i - k) : EisensteinRing) : ℂ)) * lam_mul_conj.symm
  rw [hauto k] at hsum
  have hcoe : ((∑ j, X.coeff j * star (X.coeff (j - k)) : EisensteinRing) : ℂ) =
      ∑ i, ((X.coeff i : EisensteinRing) : ℂ) *
        (starRingEnd ℂ) ((X.coeff (i - k) : EisensteinRing) : ℂ) := by
    push_cast
    rfl
  rw [hcoe]
  by_cases hk : k = 0
  · simp only [hk, ↓reduceIte] at hsum ⊢
    simp only [SubringClass.coe_natCast]
    linear_combination hsum / 3
  · simp only [hk, ↓reduceIte] at hsum ⊢
    simp only [ZeroMemClass.coe_zero]
    linear_combination hsum / 3

/-- Character twist of a block autocorrelation. -/
theorem char_autocorr {n : ℕ} [NeZero n] (χ : Multiplicative (ZMod 3) →* ℂ)
    (hχ : ∀ f, χ (Multiplicative.ofAdd f) * (starRingEnd ℂ) (χ (Multiplicative.ofAdd f)) = 1)
    (A : ZMod 3 → ZMod n → ℂ) (k : ZMod n) :
    ∑ i, (∑ f, χ (Multiplicative.ofAdd f) * A f i) *
        (starRingEnd ℂ) (∑ f, χ (Multiplicative.ofAdd f) * A f (i - k)) =
      ∑ s, (starRingEnd ℂ) (χ (Multiplicative.ofAdd s)) *
        ∑ f, ∑ i, A f i * (starRingEnd ℂ) (A (f + s) (i - k)) := by
  have hstep : ∀ i, (∑ f, χ (Multiplicative.ofAdd f) * A f i) *
      (starRingEnd ℂ) (∑ f, χ (Multiplicative.ofAdd f) * A f (i - k)) =
      ∑ f, ∑ s, (starRingEnd ℂ) (χ (Multiplicative.ofAdd s)) *
        (A f i * (starRingEnd ℂ) (A (f + s) (i - k))) := by
    intro i
    rw [map_sum, Finset.sum_mul_sum]
    apply Finset.sum_congr rfl
    intro f _
    rw [← Equiv.sum_comp (Equiv.addLeft f)]
    apply Finset.sum_congr rfl
    intro s _
    simp only [Equiv.coe_addLeft]
    have hmul : χ (Multiplicative.ofAdd (f + s)) =
        χ (Multiplicative.ofAdd f) * χ (Multiplicative.ofAdd s) := by
      rw [← map_mul]; rfl
    rw [map_mul, hmul, map_mul]
    linear_combination ((starRingEnd ℂ) (χ (Multiplicative.ofAdd s)) * A f i *
      (starRingEnd ℂ) (A (f + s) (i - k))) * hχ f
  rw [Finset.sum_congr rfl (fun i _ => hstep i)]
  simp only [Finset.mul_sum]
  conv_lhs => rw [Finset.sum_comm]
  conv_lhs => arg 2; ext f; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]

section Sequence

variable (m : ℕ)

theorem coprime_three (h3 : ¬ 3 ∣ m) : Nat.Coprime 3 m :=
  (Nat.Prime.coprime_iff_not_dvd Nat.prime_three).mpr h3

/-- The CRT splitting `ZMod (3m) ≃+* ZMod 3 × ZMod m`. -/
def crt3 (h3 : ¬ 3 ∣ m) : ZMod (3 * m) ≃+* ZMod 3 × ZMod m :=
  ZMod.chineseRemainder (coprime_three m h3)

theorem autocorr_split3 [NeZero m] (h3 : ¬ 3 ∣ m) (a : ZMod (3 * m) → ℂ)
    (s : ZMod 3) (k : ZMod m) :
    ∑ j, a j * (starRingEnd ℂ) (a (j + (crt3 m h3).symm (s, -k))) =
      ∑ f, ∑ i, a ((crt3 m h3).symm (f, i)) *
        (starRingEnd ℂ) (a ((crt3 m h3).symm (f + s, i - k))) := by
  let e := crt3 m h3
  rw [← e.symm.toEquiv.sum_comp, Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro f _
  apply Finset.sum_congr rfl
  intro i _
  congr 3
  change e.symm (f, i) + e.symm (s, -k) = e.symm (f + s, i - k)
  rw [← map_add]
  congr 1
  ext <;> simp [sub_eq_add_neg]

end Sequence

/-- The character `f ↦ r ^ f` of `ZMod 3` for a cube root of unity `r`. -/
def chi (r : ℂ) (hr : r ^ 3 = 1) : Multiplicative (ZMod 3) →* ℂ :=
  CyclicPolynomial.powerCharacter 3 r hr

theorem chi_mul_conj (r : ℂ) (hr : r ^ 3 = 1) (f : ZMod 3) :
    chi r hr (Multiplicative.ofAdd f) * (starRingEnd ℂ) (chi r hr (Multiplicative.ofAdd f)) = 1 := by
  simp only [chi, CyclicPolynomial.powerCharacter_apply, map_pow, ← mul_pow,
    cube_root_mul_conj hr, one_pow]

theorem chi_sum (r : ℂ) (hr : r ^ 3 = 1) (F : ZMod 3 → ℂ) :
    ∑ f, chi r hr (Multiplicative.ofAdd f) * F f = F 0 + r * F 1 + r ^ 2 * F 2 := by
  rw [sum_zmod_three]
  simp only [chi, CyclicPolynomial.powerCharacter_apply]
  have h0 : (0 : ZMod 3).val = 0 := rfl
  have h1 : (1 : ZMod 3).val = 1 := rfl
  have h2 : (2 : ZMod 3).val = 2 := rfl
  rw [h0, h1, h2, pow_zero, pow_one, one_mul]

/-- **Theorem T3.** For `m > 1` with `3 ∤ m` there is no perfect ternary sequence
of length `3m`: no `a : ZMod (3m) → ℂ` with `a_j ^ 3 = 1` has all nontrivial
periodic autocorrelations equal to zero. -/
theorem not_exists_perfect_ternary (m : ℕ) [NeZero m] (hm : 1 < m) (h3 : ¬ 3 ∣ m) :
    ¬ ∃ a : ZMod (3 * m) → ℂ,
      (∀ j, a j ^ 3 = 1) ∧
      ∀ τ : ZMod (3 * m), τ ≠ 0 →
        ∑ j, a j * (starRingEnd ℂ) (a (j + τ)) = 0 := by
  classical
  rintro ⟨a, hq, hauto⟩
  let e := crt3 m h3
  let w : ZMod (3 * m) → EisensteinRing := fun j => Classical.choose (cube_root_lift (hq j))
  have hw (j) : a j = 1 + (1 - omega) * (w j : ℂ) :=
    Classical.choose_spec (cube_root_lift (hq j))
  let W : ZMod 3 → ZMod m → EisensteinRing := fun f i => w (e.symm (f, i))
  let A : ZMod 3 → ZMod m → ℂ := fun f i => a (e.symm (f, i))
  have hA (f i) : A f i = 1 + (1 - omega) * (W f i : ℂ) := hw _
  let Aut : ZMod 3 → ZMod m → ℂ := fun s k =>
    ∑ f, ∑ i, A f i * (starRingEnd ℂ) (A (f + s) (i - k))
  have hAut (s k) : Aut s k = ∑ j, a j * (starRingEnd ℂ) (a (j + e.symm (s, -k))) :=
    (autocorr_split3 m h3 a s k).symm
  have hτ (s : ZMod 3) (k : ZMod m) (hsk : ¬ (s = 0 ∧ k = 0)) : e.symm (s, -k) ≠ 0 := by
    intro h
    apply hsk
    have h' := congrArg e h
    rw [RingEquiv.apply_symm_apply, map_zero] at h'
    have h1 := congrArg Prod.fst h'
    have h2 := congrArg Prod.snd h'
    simp only [Prod.fst_zero, Prod.snd_zero, neg_eq_zero] at h1 h2
    exact ⟨h1, h2⟩
  have hAut0 (k : ZMod m) : Aut 0 k = if k = 0 then (3 * m : ℂ) else 0 := by
    rw [hAut]
    by_cases hk : k = 0
    · subst hk
      simp only [neg_zero, ↓reduceIte]
      have he0 : e.symm (0, 0) = 0 := by rw [← Prod.zero_eq_mk, map_zero]
      rw [he0]
      simp only [add_zero, cube_root_mul_conj (hq _), Finset.sum_const, Finset.card_univ,
        ZMod.card, nsmul_eq_mul, mul_one]
      push_cast
      ring
    · simp only [hk, ↓reduceIte]
      exact hauto _ (hτ 0 k (fun h => hk h.2))
  have hAutne (s : ZMod 3) (hs : s ≠ 0) (k : ZMod m) : Aut s k = 0 := by
    rw [hAut]
    exact hauto _ (hτ s k (fun h => hs h.1))
  have hP (r : ℂ) (hr : r ^ 3 = 1) (k : ZMod m) :
      ∑ i, (A 0 i + r * A 1 i + r ^ 2 * A 2 i) *
        (starRingEnd ℂ) (A 0 (i - k) + r * A 1 (i - k) + r ^ 2 * A 2 (i - k)) =
        if k = 0 then (3 * m : ℂ) else 0 := by
    have h := char_autocorr (chi r hr) (chi_mul_conj r hr) A k
    simp only [chi_sum] at h
    rw [h, sum_zmod_three]
    change (starRingEnd ℂ) (chi r hr (Multiplicative.ofAdd 0)) * Aut 0 k +
      (starRingEnd ℂ) (chi r hr (Multiplicative.ofAdd 1)) * Aut 1 k +
      (starRingEnd ℂ) (chi r hr (Multiplicative.ofAdd 2)) * Aut 2 k = _
    rw [hAutne 1 (by decide) k, hAutne 2 (by decide) k, hAut0 k]
    simp [chi]
  -- the three group-ring elements
  let ω := omegaE
  let c1 : Elem EisensteinRing m := ofCoeffs fun i => (1 + 1 + ω) + W 0 i + W 1 i + W 2 i
  let cw : Elem EisensteinRing m := ofCoeffs fun i => W 0 i + ω * W 1 i + ω ^ 2 * W 2 i
  let cw2 : Elem EisensteinRing m := ofCoeffs fun i => W 0 i + ω ^ 2 * W 1 i + ω * W 2 i
  let y1 : Elem EisensteinRing m := ofCoeffs fun i => ω ^ 2 - W 1 i + ω ^ 2 * W 2 i
  let y2 : Elem EisensteinRing m := ofCoeffs fun i => ω ^ 2 + ω ^ 2 * W 1 i - W 2 i
  let WW : Elem EisensteinRing m := ofCoeffs fun i => ω ^ 2 * (W 1 i + W 2 i)
  have hω3 : omega ^ 3 = 1 := omega_pow_three
  have hω23 : (omega ^ 2) ^ 3 = 1 := by rw [← pow_mul, mul_comm, pow_mul, hω3, one_pow]
  have h1 : c1 * ringStar c1 = scalar m ((m : ℕ) : EisensteinRing) := by
    apply norm_of_lam m c1 (fun i => A 0 i + 1 * A 1 i + 1 ^ 2 * A 2 i) _ (hP 1 (by norm_num))
    intro i
    simp only [c1, ω, ofCoeffs_apply, Subalgebra.coe_add, Subalgebra.coe_one, coe_omegaE, hA]
    linear_combination (-1 : ℂ) * omega_sq
  have hw1 : cw * ringStar cw = scalar m ((m : ℕ) : EisensteinRing) := by
    apply norm_of_lam m cw _ _ (hP omega hω3)
    intro i
    simp only [cw, ω, ofCoeffs_apply, Subalgebra.coe_add, Subalgebra.coe_mul,
      Subalgebra.coe_pow, coe_omegaE, hA]
    linear_combination (-1 : ℂ) * omega_sq
  have hw2 : cw2 * ringStar cw2 = scalar m ((m : ℕ) : EisensteinRing) := by
    apply norm_of_lam m cw2 _ _ (hP (omega ^ 2) hω23)
    intro i
    simp only [cw2, ω, ofCoeffs_apply, Subalgebra.coe_add, Subalgebra.coe_mul,
      Subalgebra.coe_pow, coe_omegaE, hA]
    linear_combination (-(omega ^ 2 - omega + 1) + omega * (1 - omega) ^ 2 * (W 2 i : ℂ)) *
      omega_sq
  have hdiff : ∀ (z yy : Elem EisensteinRing m),
      (∀ i, ((z.coeff i : EisensteinRing) : ℂ) =
        ((c1.coeff i : EisensteinRing) : ℂ) + (1 - omega) * ((yy.coeff i : EisensteinRing) : ℂ)) →
      z = c1 + scalar m lam * yy := by
    intro z yy h
    ext i
    rw [AddMonoidAlgebra.coeff_add, Finsupp.add_apply, GeneralNorm.QuaternaryQ.scalar_mul_coeff]
    simp only [Subalgebra.coe_add, Subalgebra.coe_mul, Subalgebra.coe_sub, Subalgebra.coe_one,
      coe_omegaE, lam, h]
  have hx1 : cw = c1 + scalar m lam * y1 := by
    apply hdiff
    intro i
    simp only [cw, c1, y1, ω, ofCoeffs_apply, Subalgebra.coe_add, Subalgebra.coe_mul,
      Subalgebra.coe_sub, Subalgebra.coe_pow, Subalgebra.coe_one, coe_omegaE]
    linear_combination (omega - 2 + (omega - 1) * (W 2 i : ℂ)) * omega_sq
  have hx2 : cw2 = c1 + scalar m lam * y2 := by
    apply hdiff
    intro i
    simp only [cw2, c1, y2, ω, ofCoeffs_apply, Subalgebra.coe_add, Subalgebra.coe_mul,
      Subalgebra.coe_sub, Subalgebra.coe_pow, Subalgebra.coe_one, coe_omegaE]
    linear_combination (omega - 2 + (omega - 1) * (W 1 i : ℂ)) * omega_sq
  have h2c : ((2 : EisensteinRing) : ℂ) = 2 := map_ofNat EisensteinRing.val 2
  have hsum : ∑ j ∈ (Finset.univ : Finset (Fin 2)), ![y1, y2] j =
      scalar m (2 * ω ^ 2) * Jelem m + scalar m lam * WW := by
    rw [Fin.sum_univ_two]
    ext i
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
    rw [AddMonoidAlgebra.coeff_add, Finsupp.add_apply, AddMonoidAlgebra.coeff_add,
      Finsupp.add_apply, GeneralNorm.QuaternaryQ.scalar_mul_coeff,
      GeneralNorm.QuaternaryQ.scalar_mul_coeff]
    simp only [y1, y2, WW, ofCoeffs_apply, BinaryCoefficients.JOver_apply, mul_one]
    simp only [ω, lam, Subalgebra.coe_add, Subalgebra.coe_mul, Subalgebra.coe_sub,
      Subalgebra.coe_pow, Subalgebra.coe_one, coe_omegaE, h2c]
    linear_combination ((omega - 1) * ((W 1 i : ℂ) + (W 2 i : ℂ))) * omega_sq
  apply master_lemma m hm h3 (Finset.univ : Finset (Fin 2)) c1 h1 ![cw, cw2]
    (fun j _ => by fin_cases j <;> assumption) (fun _ => lam)
    (fun _ _ h => by
      have h' := congrArg (fun z : EisensteinRing => (z : ℂ)) h
      simp only [lam, Subalgebra.coe_sub, Subalgebra.coe_one, coe_omegaE,
        ZeroMemClass.coe_zero] at h'
      exact omega_ne_one (sub_eq_zero.mp h').symm)
    (fun _ _ => ⟨1, (mul_one _).symm⟩) ![y1, y2]
    (fun j _ => by fin_cases j <;> assumption) (2 * ω ^ 2)
    ⟨2, 2 * ω ^ 2, by
      apply Subtype.ext
      simp only [ω, lam, Subalgebra.coe_add, Subalgebra.coe_mul, Subalgebra.coe_sub,
        Subalgebra.coe_pow, Subalgebra.coe_one, coe_omegaE, SubringClass.coe_intCast, h2c]
      push_cast
      linear_combination (2 * (omega - 1)) * omega_sq, by decide⟩
    WW hsum

end CirculantHadamard.EisensteinNorm.TernaryT3

end

end OAI
