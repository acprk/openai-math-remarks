import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.AlternatingTorsion
import OAI.LinearAlgebra.CirculantHadamard.RankOne.Fourier

/-!
# Theorem F (rank-one decomposition), Gaussian two-prime squarefree case

Let `p ≠ q` be odd primes and `x ∈ ℤ[i][C_{pq}]` with `x x^* = pq`. Then

  `ε(x) · x(g) = pr_p(x)(g mod p) · pr_q(x)(g mod q)`   for every `g ∈ ZMod (p*q)`,

where `ε` is the augmentation and `pr_p(x) ∈ ℤ[i][C_p]` is the push-forward
(marginal) along `ZMod (p*q) → ZMod p`.

Proof outline (all steps formalized here, no new hypotheses):
1. `alternatingProduct_eq_one_of_two_primes` (Δ = 1) for the fixed primitive
   characters `ρ_p`, `ρ_q`, rewritten as `x̂(1) x̂(χψ) = x̂(χ) x̂(ψ)`.
2. Every other nontrivial character pair `(χ_p^k, χ_q^l)` is reached by applying
   step 1 to the relabelled element `a ↦ x(w⁻¹ a)` for the unit `w` with CRT
   coordinates `(k, l)`; relabelling by a unit preserves `x x^* = N`.
   (This replaces Galois transport.) Pairs with a trivial component are tautologies.
3. Fourier inversion on `C_p × C_q` (`RankOne.rankOne_of_transform`).
4. The identity in `ℂ` descends to `ℤ[i]` since `ℤ[i] → ℂ` is injective.
-/

namespace OAI

noncomputable section

namespace CirculantHadamard.RankOne

open Finset
open CyclicRing CyclotomicRings PrimeComponentsN PrimeCharacterEvaluationsN

/-! ### Relabelling the cyclic group by a unit -/

section Twist

variable {R : Type*} [Semiring R] {n : ℕ} [NeZero n]

/-- Relabel by the unit `w`: `(twist w x)(a) = x(w⁻¹ a)`. -/
def twist (w : (ZMod n)ˣ) (x : Elem R n) : Elem R n :=
  ofCoeffs fun a => x.coeff (((w⁻¹ : (ZMod n)ˣ) : ZMod n) * a)

@[simp] theorem twist_apply (w : (ZMod n)ˣ) (x : Elem R n) (a : ZMod n) :
    (twist w x).coeff a = x.coeff (((w⁻¹ : (ZMod n)ˣ) : ZMod n) * a) := rfl

theorem twist_mul (w : (ZMod n)ˣ) (x y : Elem R n) :
    twist w (x * y) = twist w x * twist w y := by
  ext a
  rw [twist_apply, mul_apply, mul_apply]
  symm
  apply Fintype.sum_equiv (Units.mulLeft w⁻¹)
  intro j
  simp only [twist_apply, Units.mulLeft_apply, mul_sub]

theorem twist_ringStar [StarRing R] (w : (ZMod n)ˣ) (x : Elem R n) :
    twist w (ringStar x) = ringStar (twist w x) := by
  ext a
  simp only [twist_apply, ringStar_apply, mul_neg]

theorem twist_scalar (w : (ZMod n)ˣ) (r : R) :
    twist w (scalar n r) = scalar n r := by
  ext a
  simp only [twist_apply, scalar_apply, Units.mul_right_eq_zero]

theorem twist_norm [StarRing R] (w : (ZMod n)ˣ) (x : Elem R n) (r : R)
    (hx : x * ringStar x = scalar n r) :
    twist w x * ringStar (twist w x) = scalar n r := by
  rw [← twist_ringStar, ← twist_mul, hx, twist_scalar]

end Twist

/-! ### Arithmetic of `N = p * q` -/

section TwoPrimes

variable {N : ℕ} (p q : PrimeIndex N) (hpq : p ≠ q) (hN : N = p.val * q.val)
include hpq hN

omit hpq in
theorem primeIndex_eq_or (r : PrimeIndex N) : r = p ∨ r = q := by
  have hr := prime_of_mem r.property
  have hdvd : r.val ∣ p.val * q.val := hN ▸ dvd_of_mem r.property
  rcases (Nat.Prime.dvd_mul hr).mp hdvd with h | h
  · exact Or.inl (Subtype.ext ((Nat.prime_dvd_prime_iff_eq hr (prime_of_mem p.property)).mp h))
  · exact Or.inr (Subtype.ext ((Nat.prime_dvd_prime_iff_eq hr (prime_of_mem q.property)).mp h))

theorem order_eq_left : order N p.val = p.val := by
  have hp := prime_of_mem p.property
  have hq := prime_of_mem q.property
  have hne : p.val ≠ q.val := fun h => hpq (Subtype.ext h)
  have hexp : exponent N p.val = 1 := by
    have h : (p.val * q.val).factorization p.val = 1 := by
      rw [Nat.factorization_mul hp.ne_zero hq.ne_zero, Finsupp.add_apply,
        hp.factorization_self, Nat.factorization_eq_zero_of_not_dvd]
      intro h
      exact hne ((Nat.prime_dvd_prime_iff_eq hp hq).mp h)
    rw [← hN] at h
    exact h
  rw [order, hexp, pow_one]

omit hpq in
theorem univ_eq_pair : (Finset.univ : Finset (PrimeIndex N)) = {p, q} := by
  ext r
  simp only [Finset.mem_univ, Finset.mem_insert, Finset.mem_singleton, true_iff]
  exact primeIndex_eq_or p q hN r

theorem coordinate_bijective (hu : 0 < N) :
    Function.Bijective (fun g : ZMod N => (coordinate N hu p g, coordinate N hu q g)) := by
  have : NeZero N := ⟨hu.ne'⟩
  have hinj : Function.Injective
      (fun g : ZMod N => (coordinate N hu p g, coordinate N hu q g)) := by
    intro a b hab
    simp only [Prod.mk.injEq] at hab
    apply (crt N hu).injective
    funext r
    rcases primeIndex_eq_or p q hN r with rfl | rfl
    · exact hab.1
    · exact hab.2
  refine (Fintype.bijective_iff_injective_and_card _).mpr ⟨hinj, ?_⟩
  rw [Fintype.card_prod, ZMod.card, ZMod.card, ZMod.card, order_eq_left p q hpq hN,
    order_eq_left q p (Ne.symm hpq) (hN.trans (Nat.mul_comm _ _))]
  exact hN

end TwoPrimes

theorem coordinate_mul (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) (a b : ZMod N) :
    coordinate N hu p (a * b) = coordinate N hu p a * coordinate N hu p b := by
  change crtRing N hu (a * b) p = crtRing N hu a p * crtRing N hu b p
  rw [map_mul]
  rfl

theorem coordinate_one (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    coordinate N hu p 1 = 1 := by
  change crtRing N hu 1 p = 1
  rw [map_one]
  rfl

/-! ### Characters of `ZMod m` with `m = P` -/

theorem addChar_eq_rho_pow {m P : ℕ} [NeZero m] (hm : m = P)
    (χ : AddChar (ZMod m) ℂ) :
    ∃ k : ZMod m, ∀ c : ZMod m, χ c = rho P ^ (k * c).val := by
  subst hm
  have hP : m ≠ 0 := NeZero.ne m
  have hroot := rho_isPrimitiveRoot m hP
  have h1 : χ 1 ^ m = 1 := by
    rw [← AddChar.map_nsmul_eq_pow, nsmul_one, ZMod.natCast_self,
      AddChar.map_zero_eq_one]
  obtain ⟨i, hi, hiχ⟩ := hroot.eq_pow_of_pow_eq_one h1
  refine ⟨(i : ZMod m), fun c => ?_⟩
  have hc : χ c = χ 1 ^ c.val := by
    rw [← AddChar.map_nsmul_eq_pow, nsmul_one, ZMod.natCast_zmod_val]
  rw [hc, ← hiχ, ← pow_mul, ZMod.val_mul, ZMod.val_natCast, Nat.mod_eq_of_lt hi,
    ← pow_eq_pow_mod _ (rho_pow m hP)]

theorem addChar_eq_one_of_zero {m P : ℕ} [NeZero m] (χ : AddChar (ZMod m) ℂ)
    (hχ : ∀ c : ZMod m, χ c = rho P ^ ((0 : ZMod m) * c).val) : χ = 1 := by
  ext c
  simp [hχ c]

/-! ### Δ = 1 as a product identity -/

/-- For `N = p q`, `Δ = 1` gives `x̂(1) x̂(ρ_p ρ_q) = x̂(ρ_p) x̂(ρ_q)` in `B N`. -/
theorem subsetEvaluation_pair_identity (N : ℕ) [NeZero N] (hu : 1 < N) (huodd : Odd N)
    (x : Elem GaussianRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : GaussianRing))
    (p q : PrimeIndex N) (hpq : p ≠ q) (hN : N = p.val * q.val) :
    subsetEvaluation N (lt_trans Nat.zero_lt_one hu) ∅ x *
        subsetEvaluation N (lt_trans Nat.zero_lt_one hu) {p, q} x =
      subsetEvaluation N (lt_trans Nat.zero_lt_one hu) {p} x *
        subsetEvaluation N (lt_trans Nat.zero_lt_one hu) {q} x := by
  classical
  set hu0 : 0 < N := lt_trans Nat.zero_lt_one hu
  let K := FractionRing (B N)
  let v : Finset (PrimeIndex N) → K := fun S =>
    algebraMap (B N) K (subsetEvaluation N hu0 S x)
  have hv : ∀ S, v S ≠ 0 := by
    intro S hS
    apply subsetEvaluation_ne_zero N hu0 S x hx
    apply IsFractionRing.injective (B N) K
    simpa only [map_zero] using hS
  have hΔ := GeneralNorm.alternatingProduct_eq_one_of_two_primes N hu huodd x hx p q hpq
  change alternatingProduct (Finset.univ : Finset (PrimeIndex N)) v = 1 at hΔ
  rw [alternatingProduct_pair _ v p (Finset.mem_univ p)] at hΔ
  have herase : (Finset.univ : Finset (PrimeIndex N)).erase p = {q} := by
    rw [univ_eq_pair p q hN, Finset.erase_insert (by simpa using hpq)]
  have hne : (∅ : Finset (PrimeIndex N)) ≠ {q} := (Finset.singleton_ne_empty q).symm
  have hpow : ({q} : Finset (PrimeIndex N)).powerset = {∅, {q}} := by
    ext S
    simp [Finset.subset_singleton_iff]
  rw [herase, hpow, Finset.prod_pair hne, subsetSign_empty,
    zpow_one] at hΔ
  have hsq : subsetSign ({q} : Finset (PrimeIndex N)) = -1 := by
    simp [subsetSign]
  rw [hsq, zpow_neg_one, insert_empty_eq] at hΔ
  -- hΔ : v ∅ / v {p} * (v {q} / v (insert p {q}))⁻¹ = 1
  have hkey : v ∅ * v {p, q} = v {p} * v {q} := by
    have h := hΔ
    rw [inv_div, div_mul_div_comm, div_eq_one_iff_eq (mul_ne_zero (hv _) (hv _))] at h
    exact h
  apply IsFractionRing.injective (B N) K
  simpa only [map_mul] using hkey

/-! ### Main theorem -/

/-- **Theorem F, coordinate form.** For `N = p q` with `p ≠ q` odd primes and
`x ∈ ℤ[i][C_N]` with `x x^* = N`, the augmentation-weighted element is the product
of its two CRT marginals. -/
theorem theoremF_coordinates (N : ℕ) [NeZero N] (hu : 1 < N) (huodd : Odd N)
    (x : Elem GaussianRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : GaussianRing))
    (p q : PrimeIndex N) (hpq : p ≠ q) (hN : N = p.val * q.val) (g : ZMod N) :
    augmentation N x * x.coeff g =
      (∑ h with coordinate N (lt_trans Nat.zero_lt_one hu) p h =
          coordinate N (lt_trans Nat.zero_lt_one hu) p g, x.coeff h) *
      (∑ h with coordinate N (lt_trans Nat.zero_lt_one hu) q h =
          coordinate N (lt_trans Nat.zero_lt_one hu) q g, x.coeff h) := by
  classical
  set hu0 : 0 < N := lt_trans Nat.zero_lt_one hu
  have hbij := coordinate_bijective p q hpq hN hu0
  have hordp := order_eq_left p q hpq hN
  have hordq := order_eq_left q p (Ne.symm hpq) (hN.trans (Nat.mul_comm _ _))
  set f : ZMod N → ℂ := fun b => ((x.coeff b : GaussianRing) : ℂ) with hfdef
  -- Every subset evaluation of a relabelled element, as a sum over the original group.
  have heval : ∀ (w : (ZMod N)ˣ) (S : Finset (PrimeIndex N)),
      ((subsetEvaluation N hu0 S (twist w x) : B N) : ℂ) =
        ∑ b, f b * ∏ r ∈ S, rho r.val ^ (coordinate N hu0 r ((w : ZMod N) * b)).val := by
    intro w S
    rw [subsetEvaluation, evaluate_apply, AddSubmonoidClass.coe_finsetSum]
    symm
    apply Fintype.sum_equiv (Units.mulLeft w)
    intro b
    simp only [Units.mulLeft_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, twist_apply,
      subsetCharacter, rootsCharacter_apply, Units.inv_mul_cancel_left, MulMemClass.coe_mul,
      SubmonoidClass.coe_finsetProd, SubmonoidClass.coe_pow, coe_rhoB, hfdef]
    rw [coe_gaussianToB]
  have hF : ∀ (χ : AddChar (Component N p) ℂ) (ψ : AddChar (Component N q) ℂ),
      transform (coordinate N hu0 p) (coordinate N hu0 q) f χ ψ *
          transform (coordinate N hu0 p) (coordinate N hu0 q) f 1 1 =
        transform (coordinate N hu0 p) (coordinate N hu0 q) f χ 1 *
          transform (coordinate N hu0 p) (coordinate N hu0 q) f 1 ψ := by
    intro χ ψ
    obtain ⟨k, hk⟩ := addChar_eq_rho_pow hordp χ
    obtain ⟨l, hl⟩ := addChar_eq_rho_pow hordq ψ
    by_cases hk0 : k = 0
    · subst hk0
      rw [addChar_eq_one_of_zero χ hk]
      ring
    by_cases hl0 : l = 0
    · subst hl0
      rw [addChar_eq_one_of_zero ψ hl]
    -- the unit of `ZMod N` with CRT coordinates `(k, l)`
    have hcopk : Nat.Coprime k.val (order N p.val) := by
      have key : ∀ n : ℕ, 0 < n → n < order N p.val → Nat.Coprime n (order N p.val) := by
        intro n h1 h2
        rw [hordp] at h2 ⊢
        exact (Nat.coprime_of_lt_prime h1.ne' h2 (prime_of_mem p.property)).symm
      exact key _ (ZMod.val_pos.mpr hk0) (ZMod.val_lt k)
    have hcopl : Nat.Coprime l.val (order N q.val) := by
      have key : ∀ n : ℕ, 0 < n → n < order N q.val → Nat.Coprime n (order N q.val) := by
        intro n h1 h2
        rw [hordq] at h2 ⊢
        exact (Nat.coprime_of_lt_prime h1.ne' h2 (prime_of_mem q.property)).symm
      exact key _ (ZMod.val_pos.mpr hl0) (ZMod.val_lt l)
    set ku := ZMod.unitOfCoprime k.val hcopk
    set lu := ZMod.unitOfCoprime l.val hcopl
    have hku : (ku : Component N p) = k := by
      simp [ku, ZMod.coe_unitOfCoprime]
    have hlu : (lu : Component N q) = l := by
      simp [lu, ZMod.coe_unitOfCoprime]
    obtain ⟨w, hw⟩ := hbij.2 (k, l)
    obtain ⟨w', hw'⟩ := hbij.2 (((ku⁻¹ : (Component N p)ˣ) : Component N p),
      ((lu⁻¹ : (Component N q)ˣ) : Component N q))
    simp only [Prod.mk.injEq] at hw hw'
    have hww' : w * w' = 1 := by
      apply hbij.1
      simp only [Prod.mk.injEq, coordinate_mul, coordinate_one, hw.1, hw.2, hw'.1, hw'.2]
      rw [← hku, ← hlu]
      exact ⟨ku.mul_inv, lu.mul_inv⟩
    let wu : (ZMod N)ˣ := ⟨w, w', hww', by rw [mul_comm]; exact hww'⟩
    have hy := twist_norm wu x _ hx
    have hid := subsetEvaluation_pair_identity N hu huodd (twist wu x) hy p q hpq hN
    have hidC := congrArg (fun z : B N => (z : ℂ)) hid
    simp only [MulMemClass.coe_mul] at hidC
    rw [heval, heval, heval, heval] at hidC
    simp only [Finset.prod_empty, Finset.prod_singleton, Finset.prod_pair hpq, coordinate_mul,
      wu, hw.1, hw.2, mul_one] at hidC
    simp only [transform, hk, hl, AddChar.one_apply, mul_one, mul_assoc] at hidC ⊢
    linear_combination hidC
  have hC := rankOne_of_transform (coordinate N hu0 p) (coordinate N hu0 q) hbij f hF g
  apply Subtype.ext
  rw [augmentation_apply]
  push_cast
  simpa only [hfdef] using hC

/-! ### Headline form with explicit primes and push-forwards -/

/-- Push-forward (marginal) of `x ∈ R[C_N]` along the reduction `ZMod N → ZMod m`. -/
def pushforward {R : Type*} [Semiring R] {N : ℕ} [NeZero N] (m : ℕ) [NeZero m] (x : Elem R N) :
    Elem R m :=
  ofCoeffs fun a => ∑ h with (ZMod.cast h : ZMod m) = a, x.coeff h

theorem pushforward_apply {R : Type*} [Semiring R] {N : ℕ} [NeZero N] (m : ℕ) [NeZero m]
    (x : Elem R N) (a : ZMod m) :
    (pushforward m x).coeff a = ∑ h with (ZMod.cast h : ZMod m) = a, x.coeff h := rfl

private theorem cast_eq_iff_of_eq {N m P : ℕ} (hm : m = P) (h g : ZMod N) :
    ((ZMod.cast h : ZMod m) = ZMod.cast g) ↔ ((ZMod.cast h : ZMod P) = ZMod.cast g) := by
  subst hm
  rfl

/-- **Theorem F** (Gaussian, squarefree two-prime case). For distinct odd primes
`p ≠ q` and `x ∈ ℤ[i][C_{pq}]` with `x x^* = pq`,
`ε(x) · x = pr_p(x) ⊗ pr_q(x)` coefficientwise:
`ε(x) · x(g) = pr_p(x)(g mod p) · pr_q(x)(g mod q)`. -/
theorem theoremF_gaussian {p q : ℕ} (hp : p.Prime) (hq : q.Prime) (hp2 : p ≠ 2)
    (hq2 : q ≠ 2) (hpq : p ≠ q) [NeZero p] [NeZero q] [NeZero (p * q)]
    (x : Elem GaussianRing (p * q))
    (hx : x * ringStar x = scalar (p * q) (((p * q : ℕ)) : GaussianRing))
    (g : ZMod (p * q)) :
    augmentation (p * q) x * x.coeff g =
      (pushforward p x).coeff (ZMod.cast g) * (pushforward q x).coeff (ZMod.cast g) := by
  classical
  have hN0 : p * q ≠ 0 := Nat.mul_ne_zero hp.ne_zero hq.ne_zero
  let P : PrimeIndex (p * q) := ⟨p, mem_primes.mpr ⟨hp, dvd_mul_right p q, hN0⟩⟩
  let Q : PrimeIndex (p * q) := ⟨q, mem_primes.mpr ⟨hq, dvd_mul_left q p, hN0⟩⟩
  have hPQ : P ≠ Q := fun h => hpq (congrArg Subtype.val h)
  have hu : 1 < p * q := Nat.one_lt_mul_iff.mpr ⟨hp.pos, hq.pos, Or.inl hp.one_lt⟩
  have hodd : Odd (p * q) := (hp.odd_of_ne_two hp2).mul (hq.odd_of_ne_two hq2)
  have hmain := theoremF_coordinates (p * q) hu hodd x hx P Q hPQ rfl g
  have hordP := order_eq_left P Q hPQ rfl
  have hordQ := order_eq_left Q P hPQ.symm (Nat.mul_comm _ _)
  have hcoord : ∀ (R : PrimeIndex (p * q)) (h : ZMod (p * q)),
      coordinate (p * q) (lt_trans Nat.zero_lt_one hu) R h =
        (ZMod.cast h : ZMod (order (p * q) R.val)) := fun R h =>
    crt_apply (p * q) _ h R
  rw [pushforward_apply, pushforward_apply, hmain]
  congr 1
  · apply Finset.sum_congr _ (fun _ _ => rfl)
    ext h
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, hcoord]
    exact cast_eq_iff_of_eq hordP h g
  · apply Finset.sum_congr _ (fun _ _ => rfl)
    ext h
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, hcoord]
    exact cast_eq_iff_of_eq hordQ h g

end CirculantHadamard.RankOne

end

end OAI
