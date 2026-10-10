import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremQ

namespace OAI

/-!
# Theorem T1: no perfect binary array / Menon–Hadamard difference set in `C₂ × C₂ × C_N`

A function `f : ZMod 2 × ZMod 2 × ZMod N → {±1}` whose nontrivial periodic
autocorrelations all vanish is a perfect binary array, equivalently a
Menon–Hadamard difference set in `C₂ × C₂ × C_N`. For odd `N > 1` none exists.
(Counting forces `N = u²`; the case `C₂ × C₂ × C_{u²}` is the non-trivial case of T1.)

Proof. Write `h = H₀₀ + X H₁₀ + Y H₀₁ + XY H₁₁` and `c_{αβ} = h(α, β) / 2`
for `α, β ∈ {±1}`; these are integral of norm `N`. With `x₀ = c₁₁` the three
differences are `c_{αβ} - c₁₁ = -2 b_j` with `∑ b_j ≡ J (mod 2)`, so the Master
Lemma (three indices, `t = -2 = (1+i)(i-1)`, `e = 1`) applies.
-/

noncomputable section

namespace CirculantHadamard.GeneralNorm.MenonT1

open scoped BigOperators
open CyclotomicRings CyclicRing QuaternaryQ

theorem sum_zmod_two {M : Type*} [AddCommMonoid M] (F : ZMod 2 → M) :
    ∑ a, F a = F 0 + F 1 :=
  Fin.sum_univ_two F

theorem zmod_two_one_add_one : (1 : ZMod 2) + 1 = 0 := by decide

variable {N : ℕ} [NeZero N]

/-- Split autocorrelation. -/
def A (H : ZMod 2 → ZMod 2 → ZMod N → ℤ) (s t : ZMod 2) (k : ZMod N) : ℤ :=
  ∑ a, ∑ b, ∑ i, H a b i * H (a + s) (b + t) (i - k)

/-- Character-evaluated block sum. -/
def P (H : ZMod 2 → ZMod 2 → ZMod N → ℤ) (α β : ℤ) (i : ZMod N) : ℤ :=
  H 0 0 i + α * H 1 0 i + β * H 0 1 i + α * β * H 1 1 i

theorem autocorr_P (H : ZMod 2 → ZMod 2 → ZMod N → ℤ) (α β : ℤ)
    (hα : α = 1 ∨ α = -1) (hβ : β = 1 ∨ β = -1) (k : ZMod N) :
    ∑ i, P H α β i * P H α β (i - k) =
      A H 0 0 k + α * A H 1 0 k + β * A H 0 1 k + α * β * A H 1 1 k := by
  simp only [A, P, sum_zmod_two, add_zero, zero_add, zmod_two_one_add_one, Finset.mul_sum,
    ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i _
  rcases hα with rfl | rfl <;> rcases hβ with rfl | rfl <;> ring

theorem A_eq (f : ZMod 2 × ZMod 2 × ZMod N → ℤ) (s t : ZMod 2) (k : ZMod N) :
    A (fun a b i => f (a, b, i)) s t k = ∑ g, f g * f (g + (s, t, -k)) := by
  rw [Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro a _
  rw [Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro b _
  apply Finset.sum_congr rfl
  intro i _
  simp only [Prod.mk_add_mk, sub_eq_add_neg]

/-- If `2 X = P` coefficientwise (integers) and `P` has autocorrelation
`4N δ`, then `X X^* = N`. -/
theorem norm_of_doubled (Xz Pz : ZMod N → ℤ) (hP : ∀ i, 2 * Xz i = Pz i)
    (hauto : ∀ k : ZMod N, ∑ i, Pz i * Pz (i - k) = if k = 0 then 4 * (N : ℤ) else 0) :
    (ofCoeffs fun i => ((Xz i : ℤ) : GaussianRing)) *
        ringStar (ofCoeffs fun i => ((Xz i : ℤ) : GaussianRing)) =
      scalar N ((N : ℕ) : GaussianRing) := by
  ext k
  rw [mul_ringStar_apply, scalar_apply]
  simp only [ofCoeffs_apply, star_intCast]
  have hz : ∑ i, Xz i * Xz (i - k) = if k = 0 then (N : ℤ) else 0 := by
    have h4 : 4 * ∑ i, Xz i * Xz (i - k) = ∑ i, Pz i * Pz (i - k) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      rw [← hP i, ← hP (i - k)]
      ring
    rw [hauto k] at h4
    split_ifs at h4 ⊢ <;> omega
  have hc : (∑ i, ((Xz i : ℤ) : GaussianRing) * ((Xz (i - k) : ℤ) : GaussianRing)) =
      ((∑ i, Xz i * Xz (i - k) : ℤ) : GaussianRing) := by
    push_cast
    rfl
  rw [hc, hz]
  split_ifs <;> simp

/-- **Theorem T1 (general odd N).** There is no `±1`-valued function on
`ZMod 2 × ZMod 2 × ZMod N` (`N > 1` odd) all of whose nontrivial periodic
autocorrelations vanish. -/
theorem not_exists_perfect_binary_array (N : ℕ) [NeZero N] (hN : 1 < N) (hodd : Odd N) :
    ¬ ∃ f : ZMod 2 × ZMod 2 × ZMod N → ℤ,
      (∀ g, f g = 1 ∨ f g = -1) ∧
      ∀ τ : ZMod 2 × ZMod 2 × ZMod N, τ ≠ 0 → ∑ g, f g * f (g + τ) = 0 := by
  classical
  rintro ⟨f, hf, hauto⟩
  let H : ZMod 2 → ZMod 2 → ZMod N → ℤ := fun a b i => f (a, b, i)
  -- values of the split autocorrelations
  have hA0 : ∀ k, A H 0 0 k = if k = 0 then 4 * (N : ℤ) else 0 := by
    intro k
    rw [A_eq]
    by_cases hk : k = 0
    · subst hk
      simp only [neg_zero, ↓reduceIte]
      have h1 : ∀ g : ZMod 2 × ZMod 2 × ZMod N, f g * f (g + (0, 0, 0)) = 1 := by
        intro g
        have : g + ((0 : ZMod 2), (0 : ZMod 2), (0 : ZMod N)) = g := by
          rw [Prod.mk_zero_zero, Prod.mk_zero_zero, add_zero]
        rw [this]
        rcases hf g with h | h <;> rw [h] <;> norm_num
      rw [Finset.sum_congr rfl (fun g _ => h1 g), Finset.sum_const, Finset.card_univ]
      simp only [Fintype.card_prod, ZMod.card, nsmul_eq_mul, mul_one]
      push_cast
      ring
    · simp only [hk, ↓reduceIte]
      apply hauto
      intro h
      apply hk
      have := congrArg (fun g : ZMod 2 × ZMod 2 × ZMod N => g.2.2) h
      simpa using this
  have hAne : ∀ s t k, ¬ (s = 0 ∧ t = 0) → A H s t k = 0 := by
    intro s t k hst
    rw [A_eq]
    apply hauto
    intro h
    apply hst
    have h1 := congrArg (fun g : ZMod 2 × ZMod 2 × ZMod N => g.1) h
    have h2 := congrArg (fun g : ZMod 2 × ZMod 2 × ZMod N => g.2.1) h
    exact ⟨h1, h2⟩
  have hPauto : ∀ α β : ℤ, (α = 1 ∨ α = -1) → (β = 1 ∨ β = -1) → ∀ k,
      ∑ i, P H α β i * P H α β (i - k) = if k = 0 then 4 * (N : ℤ) else 0 := by
    intro α β hα hβ k
    rw [autocorr_P H α β hα hβ k, hA0, hAne 1 0 k (by decide), hAne 0 1 k (by decide),
      hAne 1 1 k (by decide)]
    ring
  -- f = 1 + 2 w
  let w : ZMod 2 → ZMod 2 → ZMod N → ℤ := fun a b i => (H a b i - 1) / 2
  have hw : ∀ a b i, H a b i = 1 + 2 * w a b i := by
    intro a b i
    rcases hf (a, b, i) with h | h <;> simp only [w, H, h] <;> norm_num
  -- the four integral character values and three differences
  let c11 : ZMod N → ℤ := fun i => 2 + w 0 0 i + w 1 0 i + w 0 1 i + w 1 1 i
  let c1m : ZMod N → ℤ := fun i => w 0 0 i + w 1 0 i - w 0 1 i - w 1 1 i
  let cm1 : ZMod N → ℤ := fun i => w 0 0 i - w 1 0 i + w 0 1 i - w 1 1 i
  let cmm : ZMod N → ℤ := fun i => w 0 0 i - w 1 0 i - w 0 1 i + w 1 1 i
  let b1 : ZMod N → ℤ := fun i => 1 + w 0 1 i + w 1 1 i
  let b2 : ZMod N → ℤ := fun i => 1 + w 1 0 i + w 1 1 i
  let b3 : ZMod N → ℤ := fun i => 1 + w 1 0 i + w 0 1 i
  let E : (ZMod N → ℤ) → Elem GaussianRing N := fun z => ofCoeffs fun i => ((z i : ℤ) : GaussianRing)
  have hnorm : ∀ (z : ZMod N → ℤ) (α β : ℤ), (α = 1 ∨ α = -1) → (β = 1 ∨ β = -1) →
      (∀ i, 2 * z i = P H α β i) → E z * ringStar (E z) = scalar N ((N : ℕ) : GaussianRing) :=
    fun z α β hα hβ hz => norm_of_doubled z (P H α β) hz (hPauto α β hα hβ)
  have h11 := hnorm c11 1 1 (Or.inl rfl) (Or.inl rfl) (fun i => by
    simp only [c11, P, hw]; ring)
  have h1m := hnorm c1m 1 (-1) (Or.inl rfl) (Or.inr rfl) (fun i => by
    simp only [c1m, P, hw]; ring)
  have hm1 := hnorm cm1 (-1) 1 (Or.inr rfl) (Or.inl rfl) (fun i => by
    simp only [cm1, P, hw]; ring)
  have hmm := hnorm cmm (-1) (-1) (Or.inr rfl) (Or.inr rfl) (fun i => by
    simp only [cmm, P, hw]; ring)
  have hdiff : ∀ (z b : ZMod N → ℤ), (∀ i, z i = c11 i - 2 * b i) →
      E z = E c11 + scalar N ((-2 : ℤ) : GaussianRing) * E b := by
    intro z b hzb
    ext i
    rw [AddMonoidAlgebra.coeff_add, Finsupp.add_apply, scalar_mul_coeff]
    simp only [E, ofCoeffs_apply, hzb]
    push_cast
    ring
  let x : Fin 3 → Elem GaussianRing N := ![E c1m, E cm1, E cmm]
  let y : Fin 3 → Elem GaussianRing N := ![E b1, E b2, E b3]
  let W : Elem GaussianRing N := ofCoeffs fun i =>
    (1 - gaussianI) * ((1 + w 1 0 i + w 0 1 i + w 1 1 i : ℤ) : GaussianRing)
  have hsum : ∑ j ∈ (Finset.univ : Finset (Fin 3)), y j =
      scalar N (1 : GaussianRing) * Jelem N + scalar N (1 + gaussianI) * W := by
    ext i
    rw [Fin.sum_univ_three]
    simp only [y, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.tail_cons, Matrix.head_cons]
    rw [AddMonoidAlgebra.coeff_add, Finsupp.add_apply, AddMonoidAlgebra.coeff_add,
      Finsupp.add_apply, AddMonoidAlgebra.coeff_add, Finsupp.add_apply,
      scalar_mul_coeff, scalar_mul_coeff]
    simp only [E, W, ofCoeffs_apply, BinaryCoefficients.JOver_apply, one_mul, b1, b2, b3]
    simp only [Subalgebra.coe_add, Subalgebra.coe_mul, Subalgebra.coe_sub,
      Subalgebra.coe_one, coe_gaussianI, SubringClass.coe_intCast]
    push_cast
    linear_combination (1 + (w 1 0 i : ℂ) + (w 0 1 i : ℂ) + (w 1 1 i : ℂ)) * Complex.I_sq
  apply master_lemma N hN hodd (Finset.univ : Finset (Fin 3)) (E c11) h11 x
    (fun j _ => by fin_cases j <;> assumption)
    (fun _ => ((-2 : ℤ) : GaussianRing))
    (fun _ _ h => by
      have h' := congrArg (fun z : GaussianRing => (z : ℂ)) h
      simp only [SubringClass.coe_intCast, ZeroMemClass.coe_zero] at h'
      norm_num at h')
    (fun _ _ => ⟨gaussianI - 1, by
      apply Subtype.ext
      simp only [SubringClass.coe_intCast, Subalgebra.coe_mul,
        Subalgebra.coe_add, Subalgebra.coe_sub, Subalgebra.coe_one, coe_gaussianI]
      push_cast
      linear_combination (-1 : ℂ) * Complex.I_sq⟩)
    y
    (fun j _ => by
      fin_cases j
      · exact hdiff c1m b1 (fun i => by simp only [c1m, c11, b1]; ring)
      · exact hdiff cm1 b2 (fun i => by simp only [cm1, c11, b2]; ring)
      · exact hdiff cmm b3 (fun i => by simp only [cmm, c11, b3]; ring))
    1 ⟨0, by simp⟩ W hsum

/-- **Theorem T1, non-trivial case.** For odd `u > 1` there is no Menon–Hadamard
difference set (perfect binary array) in `C₂ × C₂ × C_{u²}`. -/
theorem not_exists_menon_C2_C2_square (u : ℕ) [NeZero u] (hu : 1 < u) (hodd : Odd u) :
    ¬ ∃ f : ZMod 2 × ZMod 2 × ZMod (u ^ 2) → ℤ,
      (∀ g, f g = 1 ∨ f g = -1) ∧
      ∀ τ : ZMod 2 × ZMod 2 × ZMod (u ^ 2), τ ≠ 0 → ∑ g, f g * f (g + τ) = 0 := by
  exact not_exists_perfect_binary_array (u ^ 2)
    (Nat.one_lt_pow (by norm_num) hu) hodd.pow

end CirculantHadamard.GeneralNorm.MenonT1

end

end OAI
