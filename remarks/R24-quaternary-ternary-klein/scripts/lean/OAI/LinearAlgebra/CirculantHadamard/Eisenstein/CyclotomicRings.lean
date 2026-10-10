import OAI.LinearAlgebra.CirculantHadamard.CyclotomicRings

namespace OAI

/-!
# Eisenstein coefficient rings `ℤ[ω, ρ_p : p ∣ u]`

Generalization of `OAI.LinearAlgebra.CirculantHadamard.CyclotomicRings`,
where the extra root `i` (order 4) is replaced by `ω = ρ_3` (order 3), and the
hypothesis "`u` odd" by "`3 ∤ u`". The prime roots `ρ_p` and the conjugation lemmas
for roots of unity are imported unchanged.
-/

noncomputable section

namespace CirculantHadamard.CyclotomicRingsW

open scoped BigOperators
open CirculantHadamard.CyclotomicRings (rho rho_isPrimitiveRoot rho_pow star_rho rho_mul_star
  star_eq_pow_pred_of_pow_eq_one rho_zero)

/-- The primitive cube root of unity `ω = exp(2πi/3)`. -/
def omega : ℂ := rho 3

theorem omega_isPrimitiveRoot : IsPrimitiveRoot omega 3 :=
  rho_isPrimitiveRoot 3 (by norm_num)

theorem omega_pow_three : omega ^ 3 = 1 := omega_isPrimitiveRoot.pow_eq_one

theorem omega_sq : omega ^ 2 = -1 - omega := by
  have h := omega_isPrimitiveRoot.geom_sum_eq_zero (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, pow_one, zero_add] at h
  linear_combination h

theorem star_omega : star omega = omega ^ 2 := star_rho 3 (by norm_num)

def EisensteinRing : Subalgebra ℤ ℂ := Algebra.adjoin ℤ ({omega} : Set ℂ)

def omegaE : EisensteinRing := ⟨omega, Algebra.subset_adjoin (by simp)⟩

@[simp] theorem coe_omegaE : (omegaE : ℂ) = omega := rfl

theorem omegaE_sq : omegaE ^ 2 = -1 - omegaE :=
  Subtype.ext omega_sq

/-- A finite list of prime directions, keeping `ω` in the coefficient ring. -/
def rootRing (s : Finset ℕ) : Subalgebra ℤ ℂ :=
  Algebra.adjoin ℤ (({omega} : Set ℂ) ∪ rho '' (s : Set ℕ))

/-- The literal global ring in the alternating products. -/
abbrev B (u : ℕ) : Subalgebra ℤ ℂ := rootRing u.primeFactors

/-- The literal coefficient ring with the `p`-direction omitted. -/
abbrev awayRing (u p : ℕ) : Subalgebra ℤ ℂ := rootRing (u.primeFactors.erase p)

theorem omega_mem_rootRing (s : Finset ℕ) : omega ∈ rootRing s :=
  Algebra.subset_adjoin (Or.inl (by simp))

theorem rho_mem_rootRing (s : Finset ℕ) (p : ℕ) (hp : p ∈ s) : rho p ∈ rootRing s :=
  Algebra.subset_adjoin (Or.inr ⟨p, hp, rfl⟩)

theorem omega_mem_B (u : ℕ) : omega ∈ B u := omega_mem_rootRing _

theorem rho_mem_B (u p : ℕ) (hp : p ∈ u.primeFactors) : rho p ∈ B u :=
  rho_mem_rootRing _ _ hp

theorem omega_mem_awayRing (u p : ℕ) : omega ∈ awayRing u p := omega_mem_rootRing _

theorem rho_mem_awayRing (u p q : ℕ) (hq : q ∈ u.primeFactors) (hne : q ≠ p) :
    rho q ∈ awayRing u p :=
  rho_mem_rootRing _ _ (Finset.mem_erase.mpr ⟨hne, hq⟩)

theorem rootRing_mono {s t : Finset ℕ} (hst : s ⊆ t) : rootRing s ≤ rootRing t := by
  apply Algebra.adjoin_mono
  exact Set.union_subset_union_right _ (Set.image_mono hst)

theorem EisensteinRing_le_rootRing (s : Finset ℕ) : EisensteinRing ≤ rootRing s :=
  Algebra.adjoin_mono Set.subset_union_left

theorem awayRing_le_B (u p : ℕ) : awayRing u p ≤ B u :=
  rootRing_mono (Finset.erase_subset _ _)

def awayInclusion (u p : ℕ) : awayRing u p →ₐ[ℤ] B u :=
  Subalgebra.inclusion (awayRing_le_B u p)

def eisensteinToB (u : ℕ) : EisensteinRing →ₐ[ℤ] B u :=
  Subalgebra.inclusion (EisensteinRing_le_rootRing _)

def eisensteinToAway (u p : ℕ) : EisensteinRing →ₐ[ℤ] awayRing u p :=
  Subalgebra.inclusion (EisensteinRing_le_rootRing _)

@[simp] theorem coe_awayInclusion (u p : ℕ) (x : awayRing u p) :
    (awayInclusion u p x).val = x.val := by
  exact Subalgebra.coe_inclusion (R := ℤ) (A := ℂ)
    (S := awayRing u p) (T := B u) (awayRing_le_B u p) x

@[simp] theorem coe_eisensteinToB (u : ℕ) (x : EisensteinRing) :
    (eisensteinToB u x).val = x.val := by
  exact Subalgebra.coe_inclusion (R := ℤ) (A := ℂ)
    (S := EisensteinRing) (T := B u) (EisensteinRing_le_rootRing _) x

@[simp] theorem coe_eisensteinToAway (u p : ℕ) (x : EisensteinRing) :
    (eisensteinToAway u p x).val = x.val := by
  exact Subalgebra.coe_inclusion (R := ℤ) (A := ℂ)
    (S := EisensteinRing) (T := awayRing u p) (EisensteinRing_le_rootRing _) x

theorem eisenstein_mem (s : Finset ℕ) (a b : ℤ) :
    (a : ℂ) + (b : ℂ) * omega ∈ rootRing s :=
  (rootRing s).add_mem ((rootRing s).algebraMap_mem a)
    ((rootRing s).mul_mem ((rootRing s).algebraMap_mem b) (omega_mem_rootRing s))

private theorem star_mem_adjoin {S : Set ℂ}
    (hS : ∀ z ∈ S, star z ∈ Algebra.adjoin ℤ S) {z : ℂ}
    (hz : z ∈ Algebra.adjoin ℤ S) : star z ∈ Algebra.adjoin ℤ S := by
  refine Algebra.adjoin_induction (p := fun z _ => star z ∈ Algebra.adjoin ℤ S)
    hS ?_ ?_ ?_ hz
  · intro a
    simp
  · intro x y _ _ hx hy
    simpa only [star_add] using (Algebra.adjoin ℤ S).add_mem hx hy
  · intro x y _ _ hx hy
    simpa only [star_mul] using (Algebra.adjoin ℤ S).mul_mem hy hx

theorem star_mem_EisensteinRing {z : ℂ} (hz : z ∈ EisensteinRing) :
    star z ∈ EisensteinRing := by
  apply star_mem_adjoin (S := ({omega} : Set ℂ)) ?_ hz
  rintro z (rfl : z = omega)
  rw [star_omega]
  exact EisensteinRing.pow_mem (Algebra.subset_adjoin (by simp)) 2

theorem star_mem_rootRing (s : Finset ℕ) {z : ℂ} (hz : z ∈ rootRing s) :
    star z ∈ rootRing s := by
  apply star_mem_adjoin ?_ hz
  rintro z (hI | ⟨p, hp, rfl⟩)
  · rcases hI with rfl
    rw [star_omega]
    exact (rootRing s).pow_mem (omega_mem_rootRing s) 2
  · by_cases hp0 : p = 0
    · simpa only [hp0, rho_zero, star_one, rootRing] using (rootRing s).one_mem
    · rw [star_rho p hp0]
      exact (rootRing s).pow_mem (rho_mem_rootRing s p hp) _

instance EisensteinRingStarRing : StarRing EisensteinRing where
  star x := ⟨star (x : ℂ), star_mem_EisensteinRing x.property⟩
  star_involutive x := Subtype.ext (star_star (x : ℂ))
  star_mul x y := Subtype.ext (star_mul (x : ℂ) (y : ℂ))
  star_add x y := Subtype.ext (star_add (x : ℂ) (y : ℂ))

instance rootRingStarRing (s : Finset ℕ) : StarRing (rootRing s) where
  star x := ⟨star (x : ℂ), star_mem_rootRing s x.property⟩
  star_involutive x := Subtype.ext (star_star (x : ℂ))
  star_mul x y := Subtype.ext (star_mul (x : ℂ) (y : ℂ))
  star_add x y := Subtype.ext (star_add (x : ℂ) (y : ℂ))

@[simp] theorem coe_star_EisensteinRing (x : EisensteinRing) :
    ((star x : EisensteinRing) : ℂ) = star (x : ℂ) := rfl

@[simp] theorem coe_star_rootRing (s : Finset ℕ) (x : rootRing s) :
    ((star x : rootRing s) : ℂ) = star (x : ℂ) := rfl

@[simp] theorem eisensteinToB_star (u : ℕ) (x : EisensteinRing) :
    eisensteinToB u (star x) = star (eisensteinToB u x) := by
  apply Subtype.ext
  exact (coe_eisensteinToB u (star x)).trans
    ((coe_star_EisensteinRing x).trans
      ((congrArg (star : ℂ → ℂ) (coe_eisensteinToB u x).symm).trans
        (coe_star_rootRing u.primeFactors (eisensteinToB u x)).symm))

@[simp] theorem eisensteinToAway_star (u p : ℕ) (x : EisensteinRing) :
    eisensteinToAway u p (star x) = star (eisensteinToAway u p x) := by
  apply Subtype.ext
  exact (coe_eisensteinToAway u p (star x)).trans
    ((coe_star_EisensteinRing x).trans
      ((congrArg (star : ℂ → ℂ) (coe_eisensteinToAway u p x).symm).trans
        (coe_star_rootRing (u.primeFactors.erase p) (eisensteinToAway u p x)).symm))

private theorem rho_finite_order (p : ℕ) : ∃ n : ℕ, n ≠ 0 ∧ rho p ^ n = 1 := by
  by_cases hp : p = 0
  · exact ⟨1, by decide, by simp [hp]⟩
  · exact ⟨p, hp, rho_pow p hp⟩

private theorem map_star_of_pow_eq_one (S : Subalgebra ℤ ℂ) [StarRing S]
    (hstar : ∀ x : S, ((star x : S) : ℂ) = star (x : ℂ))
    (σ : S →+* ℂ) (x : S) (n : ℕ) (hn : n ≠ 0) (hx : (x : ℂ) ^ n = 1) :
    σ (star x) = star (σ x) := by
  have hxS : x ^ n = 1 := Subtype.ext hx
  have hsource : star x = x ^ (n - 1) := by
    apply Subtype.ext
    simpa only [hstar, Subalgebra.coe_pow] using star_eq_pow_pred_of_pow_eq_one hn hx
  have htarget : σ x ^ n = 1 := by rw [← map_pow, hxS, map_one]
  rw [hsource, map_pow, star_eq_pow_pred_of_pow_eq_one hn htarget]

private theorem map_star_adjoin (S : Set ℂ) [StarRing (Algebra.adjoin ℤ S)]
    (hstar : ∀ x : Algebra.adjoin ℤ S, ((star x : Algebra.adjoin ℤ S) : ℂ) = star (x : ℂ))
    (hS : ∀ z ∈ S, ∃ n : ℕ, n ≠ 0 ∧ z ^ n = 1)
    (σ : Algebra.adjoin ℤ S →+* ℂ) (x : Algebra.adjoin ℤ S) :
    σ (star x) = star (σ x) := by
  rcases x with ⟨x, hx⟩
  refine Algebra.adjoin_induction
    (p := fun z hz => σ (star (⟨z, hz⟩ : Algebra.adjoin ℤ S)) =
      star (σ (⟨z, hz⟩ : Algebra.adjoin ℤ S))) ?_ ?_ ?_ ?_ hx
  · intro z hz
    obtain ⟨n, hn, hzn⟩ := hS z hz
    exact map_star_of_pow_eq_one _ hstar σ ⟨z, Algebra.subset_adjoin hz⟩ n hn hzn
  · intro a
    change σ (star (a : Algebra.adjoin ℤ S)) = star (σ (a : Algebra.adjoin ℤ S))
    simp
  · intro z w hz hw hzp hwp
    change σ (star ((⟨z, hz⟩ : Algebra.adjoin ℤ S) + ⟨w, hw⟩)) =
      star (σ ((⟨z, hz⟩ : Algebra.adjoin ℤ S) + ⟨w, hw⟩))
    simpa only [star_add, map_add] using congrArg₂ (· + ·) hzp hwp
  · intro z w hz hw hzp hwp
    change σ (star ((⟨z, hz⟩ : Algebra.adjoin ℤ S) * ⟨w, hw⟩)) =
      star (σ ((⟨z, hz⟩ : Algebra.adjoin ℤ S) * ⟨w, hw⟩))
    simpa only [star_mul, map_mul] using congrArg₂ (· * ·) hwp hzp

/-- Every complex-valued ring map respects the actual conjugation. -/
theorem map_star_rootRing (s : Finset ℕ) (σ : rootRing s →+* ℂ) (x : rootRing s) :
    σ (star x) = star (σ x) := by
  apply @map_star_adjoin (({omega} : Set ℂ) ∪ rho '' (s : Set ℕ))
    (rootRingStarRing s) (coe_star_rootRing s) ?_ σ x
  rintro z (hI | ⟨p, _, rfl⟩)
  · rcases hI with rfl
    exact ⟨3, by decide, omega_pow_three⟩
  · exact rho_finite_order p

theorem map_star_B (u : ℕ) (σ : B u →+* ℂ) (x : B u) :
    σ (star x) = star (σ x) := map_star_rootRing _ σ x

theorem map_star_awayRing (u p : ℕ) (σ : awayRing u p →+* ℂ) (x : awayRing u p) :
    σ (star x) = star (σ x) := map_star_rootRing _ σ x

theorem map_star_EisensteinRing (σ : EisensteinRing →+* ℂ) (x : EisensteinRing) :
    σ (star x) = star (σ x) := by
  apply @map_star_adjoin ({omega} : Set ℂ) EisensteinRingStarRing
    coe_star_EisensteinRing ?_ σ x
  rintro z (rfl : z = omega)
  exact ⟨3, by decide, omega_pow_three⟩

private theorem integral_of_finite_order {z : ℂ} (h : ∃ n : ℕ, n ≠ 0 ∧ z ^ n = 1) :
    IsIntegral ℤ z := by
  obtain ⟨n, hn, hz⟩ := h
  apply IsIntegral.of_pow (Nat.pos_of_ne_zero hn)
  rw [hz]
  exact isIntegral_one

instance rootRing_finite (s : Finset ℕ) : Module.Finite ℤ (rootRing s) := by
  apply Algebra.finite_adjoin_of_finite_of_isIntegral
    ((Set.finite_singleton omega).union (s.finite_toSet.image rho))
  rintro z (hI | ⟨p, _, rfl⟩)
  · rcases hI with rfl
    exact integral_of_finite_order ⟨3, by decide, omega_pow_three⟩
  · exact integral_of_finite_order (rho_finite_order p)

instance EisensteinRing_finite : Module.Finite ℤ EisensteinRing :=
  Algebra.finite_adjoin_simple_of_isIntegral
    (integral_of_finite_order ⟨3, by decide, omega_pow_three⟩)

private theorem primitive_mul {z w : ℂ} {m n : ℕ}
    (hz : IsPrimitiveRoot z m) (hw : IsPrimitiveRoot w n) (hcop : m.Coprime n) :
    IsPrimitiveRoot (z * w) (m * n) := by
  have horder := (Commute.all z w).orderOf_mul_eq_mul_orderOf_of_coprime
    (by simpa only [← hz.eq_orderOf, ← hw.eq_orderOf] using hcop)
  rw [← hz.eq_orderOf, ← hw.eq_orderOf] at horder
  rw [← horder]
  exact IsPrimitiveRoot.orderOf _

theorem prime_roots_prod_primitive (s : Finset ℕ) (hs : ∀ p ∈ s, Nat.Prime p) :
    IsPrimitiveRoot (∏ p ∈ s, rho p) (∏ p ∈ s, p) := by
  classical
  revert hs
  induction s using Finset.induction_on with
  | empty =>
      intro _
      simp
  | @insert p s hp ih =>
      intro hs
      have hpp := hs p (Finset.mem_insert_self p s)
      have hsp : ∀ q ∈ s, Nat.Prime q := fun q hq => hs q (Finset.mem_insert_of_mem hq)
      rw [Finset.prod_insert hp, Finset.prod_insert hp]
      apply primitive_mul (rho_isPrimitiveRoot p hpp.ne_zero) (ih hsp)
      apply Nat.coprime_prod_right_iff.mpr
      intro q hq
      apply (Nat.coprime_primes hpp (hsp q hq)).mpr
      intro heq
      exact hp (heq.symm ▸ hq)

/-- The single generator of a finite prime-root ring. -/
def etaOf (s : Finset ℕ) : ℂ := omega * ∏ q ∈ s, rho q

def orderOfSet (s : Finset ℕ) : ℕ := 3 * ∏ q ∈ s, q

theorem etaOf_isPrimitiveRoot (s : Finset ℕ)
    (hs : ∀ p ∈ s, Nat.Prime p) (hodd : ∀ p ∈ s, ¬ 3 ∣ p) :
    IsPrimitiveRoot (etaOf s) (orderOfSet s) := by
  apply primitive_mul omega_isPrimitiveRoot (prime_roots_prod_primitive s hs)
  apply Nat.coprime_prod_right_iff.mpr
  intro p hp
  exact (Nat.Prime.coprime_iff_not_dvd Nat.prime_three).mpr (hodd p hp)

theorem orderOfSet_pos (s : Finset ℕ) (hs : ∀ p ∈ s, Nat.Prime p) :
    0 < orderOfSet s :=
  Nat.mul_pos (by decide) (Finset.prod_pos (fun p hp => (hs p hp).pos))

theorem rootRing_eq_singleton (s : Finset ℕ)
    (hs : ∀ p ∈ s, Nat.Prime p) (hodd : ∀ p ∈ s, ¬ 3 ∣ p) :
    rootRing s = Algebra.adjoin ℤ ({etaOf s} : Set ℂ) := by
  classical
  have hη := etaOf_isPrimitiveRoot s hs hodd
  let : NeZero (orderOfSet s) := ⟨(orderOfSet_pos s hs).ne'⟩
  have hηmem : etaOf s ∈ rootRing s :=
    (rootRing s).mul_mem (omega_mem_rootRing s)
      ((rootRing s).prod_mem (fun q hq => rho_mem_rootRing s q hq))
  have hpower {z : ℂ} (hz : z ^ orderOfSet s = 1) :
      z ∈ Algebra.adjoin ℤ ({etaOf s} : Set ℂ) := by
    obtain ⟨j, _, hj⟩ := hη.eq_pow_of_pow_eq_one hz
    rw [← hj]
    exact (Algebra.adjoin ℤ ({etaOf s} : Set ℂ)).pow_mem
      (Algebra.subset_adjoin (by simp)) j
  apply le_antisymm
  · apply Algebra.adjoin_le
    rintro z (hI | ⟨q, hq, rfl⟩)
    · rcases hI with rfl
      apply hpower
      exact (omega_isPrimitiveRoot.pow_eq_one_iff_dvd _).mpr (dvd_mul_right 3 _)
    · apply hpower
      apply ((rho_isPrimitiveRoot q (hs q hq).ne_zero).pow_eq_one_iff_dvd _).mpr
      exact (Finset.dvd_prod_of_mem (fun q : ℕ => q) hq).trans (dvd_mul_left _ 3)
  · apply Algebra.adjoin_le
    simpa using hηmem

def eta (u p : ℕ) : ℂ := etaOf (u.primeFactors.erase p)

def awayOrder (u p : ℕ) : ℕ := orderOfSet (u.primeFactors.erase p)

theorem eta_isPrimitiveRoot (u p : ℕ) (hu : ¬ 3 ∣ u) :
    IsPrimitiveRoot (eta u p) (awayOrder u p) := by
  apply etaOf_isPrimitiveRoot
  · intro q hq
    exact Nat.prime_of_mem_primeFactors (Finset.mem_of_mem_erase hq)
  · intro q hq
    exact fun h3 => hu (h3.trans (Nat.dvd_of_mem_primeFactors (Finset.mem_of_mem_erase hq)))

theorem awayOrder_pos (u p : ℕ) : 0 < awayOrder u p :=
  orderOfSet_pos _ (fun _q hq => Nat.prime_of_mem_primeFactors (Finset.mem_of_mem_erase hq))

theorem awayRing_eq_singleton (u p : ℕ) (hu : ¬ 3 ∣ u) :
    awayRing u p = Algebra.adjoin ℤ ({eta u p} : Set ℂ) := by
  apply rootRing_eq_singleton
  · intro q hq
    exact Nat.prime_of_mem_primeFactors (Finset.mem_of_mem_erase hq)
  · intro q hq
    exact fun h3 => hu (h3.trans (Nat.dvd_of_mem_primeFactors (Finset.mem_of_mem_erase hq)))

theorem prime_not_dvd_awayOrder (u p : ℕ) (hu : ¬ 3 ∣ u) (hp : p ∈ u.primeFactors) :
    ¬p ∣ awayOrder u p := by
  have hpp := Nat.prime_of_mem_primeFactors hp
  apply hpp.coprime_iff_not_dvd.mp
  change p.Coprime (3 * ∏ q ∈ u.primeFactors.erase p, q)
  apply Nat.coprime_mul_iff_right.mpr
  constructor
  · have h3 : ¬ 3 ∣ p := fun h3 => hu (h3.trans (Nat.dvd_of_mem_primeFactors hp))
    exact ((Nat.Prime.coprime_iff_not_dvd Nat.prime_three).mpr h3).symm
  · apply Nat.coprime_prod_right_iff.mpr
    intro q hq
    apply (Nat.coprime_primes hpp
      (Nat.prime_of_mem_primeFactors (Finset.mem_of_mem_erase hq))).mpr
    exact (Finset.ne_of_mem_erase hq).symm

theorem B_eq_adjoin_away (u p : ℕ) (hp : p ∈ u.primeFactors) :
    B u = (Algebra.adjoin (awayRing u p) ({rho p} : Set ℂ)).restrictScalars ℤ := by
  change Algebra.adjoin ℤ (({omega} : Set ℂ) ∪ rho '' (u.primeFactors : Set ℕ)) =
    (Algebra.adjoin
      (Algebra.adjoin ℤ (({omega} : Set ℂ) ∪ rho '' ((u.primeFactors.erase p : Finset ℕ) : Set ℕ)))
      ({rho p} : Set ℂ)).restrictScalars ℤ
  rw [← Algebra.adjoin_union_eq_adjoin_adjoin]
  congr 1
  ext z
  constructor
  · rintro (hI | ⟨q, hq, rfl⟩)
    · exact Or.inl (Or.inl hI)
    · by_cases heq : q = p
      · exact Or.inr (by simp [heq])
      · exact Or.inl (Or.inr ⟨q, Finset.mem_erase.mpr ⟨heq, hq⟩, rfl⟩)
  · rintro ((hI | ⟨q, hq, rfl⟩) | heq)
    · exact Or.inl hI
    · exact Or.inr ⟨q, Finset.mem_of_mem_erase hq, rfl⟩
    · exact Or.inr ⟨p, hp, (Set.mem_singleton_iff.mp heq).symm⟩

end CirculantHadamard.CyclotomicRingsW

end

end OAI
