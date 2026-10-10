import OAI.LinearAlgebra.CirculantHadamard.CyclicEvaluation
import OAI.LinearAlgebra.CirculantHadamard.CyclicPolynomial
import OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.PrimeComponents
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.CyclotomicRings
import Mathlib.Algebra.Star.BigOperators
import Mathlib.Tactic.Linarith

namespace OAI

universe uR uS uT uA

/-!
# Prime-component character evaluations over ℤ[ω] (Eisenstein coefficients)

Generalization of
`OAI.LinearAlgebra.CirculantHadamard.PrimeCharacterEvaluations`, which is stated
only for `ZMod (u ^ 2)` and norm `u ^ 2`. Here the group is `ZMod N` and the norm
is `N`; the complex absolute value of every character value is `Real.sqrt N`.

The partial-evaluation homomorphism retains one cyclic coordinate and evaluates
the others. Its multiplication and subsequent evaluation laws follow from its
construction, rather than assumptions attached to the factors of a norm equation.
-/

noncomputable section

namespace CirculantHadamard.PrimeCharacterEvaluationsW

open scoped BigOperators
open CyclicRing

section WeightedProjection

variable {R : Type uR} {S : Type uS} {T : Type uT} [Semiring R] [CommSemiring S] [CommSemiring T]
variable {n m : ℕ}

/-- A character with values in the actual residual cyclic group algebra. -/
def weightedCharacter (f : ZMod n →+ ZMod m)
    (χ : Multiplicative (ZMod n) →* S) :
    Multiplicative (ZMod n) →* Elem S m where
  toFun g := AddMonoidAlgebra.single (f g.toAdd) (χ g)
  map_one' := by
    change AddMonoidAlgebra.single (f 0) (χ 1) = 1
    simp only [map_zero, map_one, AddMonoidAlgebra.one_def]
  map_mul' g h := by
    change AddMonoidAlgebra.single (f (g.toAdd + h.toAdd)) (χ (g * h)) = _
    simp only [map_add, map_mul, AddMonoidAlgebra.single_mul_single]

@[simp] theorem weightedCharacter_apply (f : ZMod n →+ ZMod m)
    (χ : Multiplicative (ZMod n) →* S) (g : ZMod n) :
    weightedCharacter f χ (Multiplicative.ofAdd g) =
      AddMonoidAlgebra.single (f g) (χ (Multiplicative.ofAdd g)) := rfl

/-- Coefficients are mapped by `φ`, the coordinate `f` is kept as a variable,
and all other character values become coefficients. -/
def weightedProjection (φ : R →+* S) (f : ZMod n →+ ZMod m)
    (χ : Multiplicative (ZMod n) →* S) : Elem R n →+* Elem S m :=
  evaluate (AddMonoidAlgebra.singleZeroRingHom.comp φ) (weightedCharacter f χ)

@[simp] theorem weightedProjection_single (φ : R →+* S)
    (f : ZMod n →+ ZMod m) (χ : Multiplicative (ZMod n) →* S)
    (g : ZMod n) (r : R) :
    weightedProjection φ f χ (AddMonoidAlgebra.single g r) =
      AddMonoidAlgebra.single (f g) (φ r * χ (Multiplicative.ofAdd g)) := by
  simp only [weightedProjection, evaluate_single, RingHom.comp_apply,
    AddMonoidAlgebra.singleZeroRingHom_apply, AddMonoidAlgebra.singleAddHom,
    weightedCharacter_apply,
    AddMonoidAlgebra.single_mul_single, zero_add]

@[simp] theorem weightedProjection_scalar (φ : R →+* S)
    (f : ZMod n →+ ZMod m) (χ : Multiplicative (ZMod n) →* S) (r : R) :
    weightedProjection φ f χ (scalar n r) = scalar m (φ r) := by
  rw [scalar, weightedProjection_single, map_zero]
  change AddMonoidAlgebra.single 0 (φ r * χ 1) = _
  simp only [map_one, mul_one, scalar]

private theorem toMultiplicative_apply (f : ZMod n →+ ZMod m) (a : ZMod n) :
    f.toMultiplicative (Multiplicative.ofAdd a) = Multiplicative.ofAdd (f a) := rfl

/-- Subsequent evaluation is exactly the product of the evaluated coefficients
and the character of the retained coordinate. -/
theorem evaluate_weightedProjection (φ : R →+* S) (ψ : S →+* T)
    (f : ZMod n →+ ZMod m) (χ : Multiplicative (ZMod n) →* S)
    (θ : Multiplicative (ZMod m) →* T) :
    (evaluate ψ θ).comp (weightedProjection φ f χ) =
      evaluate (ψ.comp φ)
        ((ψ.toMonoidHom.comp χ) * (θ.comp f.toMultiplicative)) := by
  apply AddMonoidAlgebra.ringHom_ext
  · intro r
    simp only [RingHom.comp_apply, weightedProjection_single, evaluate_single,
      map_zero, map_mul, MonoidHom.mul_apply, MonoidHom.comp_apply,
      toMultiplicative_apply, RingHom.toMonoidHom_eq_coe]
    simp
  · intro g
    simp only [RingHom.comp_apply, weightedProjection_single, evaluate_single,
      map_one, one_mul, MonoidHom.mul_apply, MonoidHom.comp_apply,
      toMultiplicative_apply, RingHom.toMonoidHom_eq_coe]
    rfl

/-- Norm preservation is a consequence of the constructed ring homomorphism.
In particular it does not require the residual group algebra to be a domain. -/
theorem weightedProjection_product (φ : R →+* S) (f : ZMod n →+ ZMod m)
    (χ : Multiplicative (ZMod n) →* S) (x y : Elem R n) (r : R)
    (hxy : x * y = scalar n r) :
    weightedProjection φ f χ x * weightedProjection φ f χ y = scalar m (φ r) := by
  rw [← map_mul, hxy, weightedProjection_scalar]

end WeightedProjection

open PrimeComponentsN CyclotomicRingsW
open CyclotomicRings (rho rho_isPrimitiveRoot rho_pow rho_mul_star)

/-- The actual CRT coordinate, without a chosen-decomposition hypothesis. -/
def coordinate (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    ZMod N →+ Component N p where
  toFun a := crt N hu a p
  map_zero' := by simp
  map_add' a b := by simp

@[simp] theorem coordinate_generator (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    coordinate N hu p (primaryGenerator N hu p) = 1 :=
  crt_primaryGenerator_same N hu p

@[simp] theorem coordinate_generator_ne (N : ℕ) (hu : 0 < N)
    (p q : PrimeIndex N) (h : q ≠ p) :
    coordinate N hu q (primaryGenerator N hu p) = 0 :=
  crt_primaryGenerator_ne N hu p q h

section RootCharacters
variable {A : Type uA} [CommMonoid A]

theorem root_order_pow (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (r : A) (hr : r ^ p.val = 1) : r ^ order N p.val = 1 := by
  obtain ⟨t, ht⟩ := prime_dvd_order hu p.property
  rw [ht, pow_mul, hr, one_pow]

def componentCharacter (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (r : A) (hr : r ^ p.val = 1) : Multiplicative (ZMod N) →* A :=
  (CyclicPolynomial.powerCharacter (order N p.val) r (root_order_pow N hu p r hr)).comp
    (coordinate N hu p).toMultiplicative

@[simp] theorem componentCharacter_apply (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (r : A) (hr : r ^ p.val = 1) (a : ZMod N) :
    componentCharacter N hu p r hr (Multiplicative.ofAdd a) =
      r ^ (coordinate N hu p a).val := rfl

@[simp] theorem componentCharacter_generator (N : ℕ) (hu : 0 < N)
    (p : PrimeIndex N) (r : A) (hr : r ^ p.val = 1) :
    componentCharacter N hu p r hr (Multiplicative.ofAdd (primaryGenerator N hu p)) = r := by
  change CyclicPolynomial.powerCharacter (order N p.val) r (root_order_pow N hu p r hr) (Multiplicative.ofAdd
    (coordinate N hu p (primaryGenerator N hu p))) = r
  rw [coordinate_generator, CyclicPolynomial.powerCharacter_generator]

@[simp] theorem componentCharacter_generator_ne (N : ℕ) (hu : 0 < N)
    (p q : PrimeIndex N) (h : q ≠ p) (r : A) (hr : r ^ q.val = 1) :
    componentCharacter N hu q r hr (Multiplicative.ofAdd (primaryGenerator N hu p)) = 1 := by
  rw [componentCharacter_apply, coordinate_generator_ne N hu p q h]
  simp

/-- A finite product of the actual primary-coordinate characters. -/
def rootsCharacter (N : ℕ) (hu : 0 < N) (r : PrimeIndex N → A)
    (hr : ∀ p, r p ^ p.val = 1) (S : Finset (PrimeIndex N)) :
    Multiplicative (ZMod N) →* A :=
  ∏ p ∈ S, componentCharacter N hu p (r p) (hr p)

@[simp] theorem rootsCharacter_apply (N : ℕ) (hu : 0 < N) (r : PrimeIndex N → A)
    (hr : ∀ p, r p ^ p.val = 1) (S : Finset (PrimeIndex N)) (a : ZMod N) :
    rootsCharacter N hu r hr S (Multiplicative.ofAdd a) =
      ∏ p ∈ S, r p ^ (coordinate N hu p a).val := by
  simp [rootsCharacter]

@[simp] theorem rootsCharacter_empty (N : ℕ) (hu : 0 < N) (r : PrimeIndex N → A)
    (hr : ∀ p, r p ^ p.val = 1) : rootsCharacter N hu r hr ∅ = 1 := by
  simp [rootsCharacter]

theorem rootsCharacter_generator (N : ℕ) (hu : 0 < N) (r : PrimeIndex N → A)
    (hr : ∀ p, r p ^ p.val = 1) (S : Finset (PrimeIndex N)) (p : PrimeIndex N) :
    rootsCharacter N hu r hr S (Multiplicative.ofAdd (primaryGenerator N hu p)) =
      if p ∈ S then r p else 1 := by
  classical
  simp only [rootsCharacter, MonoidHom.finsetProd_apply]
  have heq : ∀ q ∈ S, componentCharacter N hu q (r q) (hr q)
      (Multiplicative.ofAdd (primaryGenerator N hu p)) = if q = p then r p else 1 := by
    intro q _
    by_cases h : q = p
    · subst q
      simp only [componentCharacter_generator, ite_true]
    · simp [h, componentCharacter_generator_ne N hu p q h]
  rw [Finset.prod_congr rfl heq]
  simp

end RootCharacters

/-- The chosen primitive prime root in the literal global coefficient ring. -/
def rhoB (N : ℕ) (p : PrimeIndex N) : B N :=
  ⟨rho p.val, rho_mem_B N p.val p.property⟩

@[simp] theorem coe_rhoB (N : ℕ) (p : PrimeIndex N) :
    (rhoB N p : ℂ) = rho p.val := rfl

theorem rhoB_pow (N : ℕ) (p : PrimeIndex N) : rhoB N p ^ p.val = 1 :=
  Subtype.ext (rho_pow p.val (prime_of_mem p.property).ne_zero)

theorem rhoB_ne_one (N : ℕ) (p : PrimeIndex N) : rhoB N p ≠ 1 := by
  intro h
  exact (rho_isPrimitiveRoot p.val (prime_of_mem p.property).ne_zero).ne_one
    (prime_of_mem p.property).one_lt (congrArg Subtype.val h)

def subsetCharacter (N : ℕ) (hu : 0 < N) (S : Finset (PrimeIndex N)) :
    Multiplicative (ZMod N) →* B N :=
  rootsCharacter N hu (rhoB N) (rhoB_pow N) S

@[simp] theorem subsetCharacter_empty (N : ℕ) (hu : 0 < N) :
    subsetCharacter N hu ∅ = 1 := rootsCharacter_empty _ _ _ _

theorem subsetCharacter_generator (N : ℕ) (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (p : PrimeIndex N) :
    subsetCharacter N hu S (Multiplicative.ofAdd (primaryGenerator N hu p)) =
      if p ∈ S then rhoB N p else 1 := rootsCharacter_generator _ _ _ _ _ _

theorem subsetCharacter_nonempty (N : ℕ) (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (hS : S.Nonempty) :
    ∃ a : ZMod N, subsetCharacter N hu S (Multiplicative.ofAdd a) ≠ 1 := by
  obtain ⟨p, hp⟩ := hS
  refine ⟨primaryGenerator N hu p, ?_⟩
  rw [subsetCharacter_generator, ite_eq_left hp]
  exact rhoB_ne_one N p

def subsetEvaluation (N : ℕ) (hu : 0 < N) (S : Finset (PrimeIndex N)) :
    Elem EisensteinRing N →+* B N :=
  evaluate (eisensteinToB N).toRingHom (subsetCharacter N hu S)

/-- The omitted root is replaced by one; all remaining roots lie in the actual away ring. -/
def rhoAway (N : ℕ) (p q : PrimeIndex N) : awayRing N p.val := by
  classical
  exact if h : q = p then 1 else
    ⟨rho q.val, rho_mem_awayRing N p.val q.val q.property
      (fun he => h (Subtype.ext he))⟩

theorem rhoAway_pow (N : ℕ) (p q : PrimeIndex N) : rhoAway N p q ^ q.val = 1 := by
  classical
  by_cases h : q = p
  · simp [rhoAway, h]
  · apply Subtype.ext
    simpa [rhoAway, h] using rho_pow q.val (prime_of_mem q.property).ne_zero

def awayCharacter (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) : Multiplicative (ZMod N) →* awayRing N p.val :=
  rootsCharacter N hu (rhoAway N p) (rhoAway_pow N p) (S.erase p)

/-- Leave precisely the full p-primary cyclic coordinate unevaluated. -/
def partialEvaluation (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) :
    Elem EisensteinRing N →+* Elem (awayRing N p.val) (order N p.val) :=
  weightedProjection (eisensteinToAway N p.val).toRingHom (coordinate N hu p)
    (awayCharacter N hu p S)

theorem partialEvaluation_product (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) (x y : Elem EisensteinRing N)
    (hxy : x * y = scalar N ((N : ℕ) : EisensteinRing)) :
    partialEvaluation N hu p S x * partialEvaluation N hu p S y =
      scalar (order N p.val) ((N : ℕ) : awayRing N p.val) := by
  unfold partialEvaluation
  simpa only [map_natCast] using weightedProjection_product
    (eisensteinToAway N p.val).toRingHom (coordinate N hu p) (awayCharacter N hu p S)
    x y ((N : ℕ) : EisensteinRing) hxy

theorem awayCharacter_map (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) :
    (awayInclusion N p.val).toMonoidHom.comp (awayCharacter N hu p S) =
      subsetCharacter N hu (S.erase p) := by
  classical
  apply MonoidHom.ext
  intro a
  change awayInclusion N p.val
    (rootsCharacter N hu (rhoAway N p) (rhoAway_pow N p) (S.erase p)
      (Multiplicative.ofAdd a.toAdd)) =
    rootsCharacter N hu (rhoB N) (rhoB_pow N) (S.erase p)
      (Multiplicative.ofAdd a.toAdd)
  simp only [rootsCharacter_apply, map_prod, map_pow]
  apply Finset.prod_congr rfl
  intro q hq
  have hqp := (Finset.mem_erase.mp hq).1
  congr 1
  apply Subtype.ext
  simp [rhoAway, hqp, rhoB]

theorem subsetCharacter_insert (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) :
    subsetCharacter N hu (insert p S) =
      subsetCharacter N hu (S.erase p) * componentCharacter N hu p (rhoB N p) (rhoB_pow N p) := by
  classical
  have hi : insert p (S.erase p) = insert p S := by
    ext q
    by_cases h : q = p <;> simp [h]
  calc
    subsetCharacter N hu (insert p S) = subsetCharacter N hu (insert p (S.erase p)) :=
      congrArg (subsetCharacter N hu) hi.symm
    _ = _ := by
      simp only [subsetCharacter, rootsCharacter, Finset.prod_insert (Finset.notMem_erase p S),
        mul_comm]

/-- Evaluate the remaining p-primary variable at one. -/
def partialAtOne (N : ℕ) (p : PrimeIndex N) :
    Elem (awayRing N p.val) (order N p.val) →+* B N :=
  evaluate (awayInclusion N p.val).toRingHom 1

/-- Evaluate the remaining full p-primary variable at the selected pth root. -/
def partialAtRoot (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    Elem (awayRing N p.val) (order N p.val) →+* B N :=
  evaluate (awayInclusion N p.val).toRingHom
    (CyclicPolynomial.powerCharacter (order N p.val) (rhoB N p)
      (root_order_pow N hu p (rhoB N p) (rhoB_pow N p)))

private theorem eisenstein_away_comp (N p : ℕ) :
    (awayInclusion N p).toRingHom.comp (eisensteinToAway N p).toRingHom =
      (eisensteinToB N).toRingHom := by
  apply RingHom.ext
  intro r
  apply Subtype.ext
  exact (coe_awayInclusion N p (eisensteinToAway N p r)).trans
    ((coe_eisensteinToAway N p r).trans (coe_eisensteinToB N r).symm)

theorem partialEvaluation_at_one (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) :
    (partialAtOne N p).comp (partialEvaluation N hu p S) =
      subsetEvaluation N hu (S.erase p) := by
  unfold partialAtOne partialEvaluation
  rw [evaluate_weightedProjection]
  rw [eisenstein_away_comp]
  have ht : (1 : Multiplicative (Component N p) →* B N).comp
      (coordinate N hu p).toMultiplicative = 1 := by ext; rfl
  rw [ht, mul_one]
  exact congrArg (evaluate (eisensteinToB N).toRingHom) (awayCharacter_map N hu p S)

theorem partialEvaluation_at_root (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) :
    (partialAtRoot N hu p).comp (partialEvaluation N hu p S) =
      subsetEvaluation N hu (insert p S) := by
  unfold partialAtRoot partialEvaluation
  rw [evaluate_weightedProjection]
  rw [eisenstein_away_comp]
  change evaluate (eisensteinToB N).toRingHom
    (((awayInclusion N p.val).toMonoidHom.comp (awayCharacter N hu p S)) *
      componentCharacter N hu p (rhoB N p) (rhoB_pow N p)) = _
  rw [awayCharacter_map, ← subsetCharacter_insert]
  rfl

@[simp] theorem partialEvaluation_at_one_apply (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) (hp : p ∉ S) (x : Elem EisensteinRing N) :
    partialAtOne N p (partialEvaluation N hu p S x) = subsetEvaluation N hu S x := by
  have h := congrArg (fun f : Elem EisensteinRing N →+* B N => f x)
    (partialEvaluation_at_one N hu p S)
  simpa only [RingHom.comp_apply, Finset.erase_eq_of_notMem hp] using h

@[simp] theorem partialEvaluation_at_root_apply (N : ℕ) (hu : 0 < N) (p : PrimeIndex N)
    (S : Finset (PrimeIndex N)) (x : Elem EisensteinRing N) :
    partialAtRoot N hu p (partialEvaluation N hu p S x) =
      subsetEvaluation N hu (insert p S) x :=
  congrArg (fun f : Elem EisensteinRing N →+* B N => f x)
    (partialEvaluation_at_root N hu p S)

theorem rhoB_mul_star (N : ℕ) (p : PrimeIndex N) : rhoB N p * star (rhoB N p) = 1 :=
  Subtype.ext (rho_mul_star p.val)

theorem subsetCharacter_mul_star (N : ℕ) (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (a : ZMod N) :
    subsetCharacter N hu S (Multiplicative.ofAdd a) *
      star (subsetCharacter N hu S (Multiplicative.ofAdd a)) = 1 := by
  unfold subsetCharacter
  rw [rootsCharacter_apply, star_prod, ← Finset.prod_mul_distrib]
  apply Finset.prod_eq_one
  intro p _
  rw [star_pow, ← mul_pow, rhoB_mul_star, one_pow]

theorem subsetCharacter_neg (N : ℕ) (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (a : ZMod N) :
    subsetCharacter N hu S (Multiplicative.ofAdd (-a)) =
      star (subsetCharacter N hu S (Multiplicative.ofAdd a)) := by
  let χ := subsetCharacter N hu S
  have hi : χ (Multiplicative.ofAdd (-a)) * χ (Multiplicative.ofAdd a) = 1 := by
    rw [← map_mul]
    change χ (Multiplicative.ofAdd (-a + a)) = 1
    simp
  have hs := subsetCharacter_mul_star N hu S a
  change χ (Multiplicative.ofAdd a) * star (χ (Multiplicative.ofAdd a)) = 1 at hs
  calc
    χ (Multiplicative.ofAdd (-a)) = χ (Multiplicative.ofAdd (-a)) * 1 := (mul_one _).symm
    _ = χ (Multiplicative.ofAdd (-a)) *
        (χ (Multiplicative.ofAdd a) * star (χ (Multiplicative.ofAdd a))) := by rw [hs]
    _ = star (χ (Multiplicative.ofAdd a)) := by rw [← mul_assoc, hi, one_mul]

theorem subsetEvaluation_star (N : ℕ) [NeZero N] (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (x : Elem EisensteinRing N) :
    subsetEvaluation N hu S (ringStar x) = star (subsetEvaluation N hu S x) := by
  exact evaluate_ringStar (eisensteinToB N).toRingHom (subsetCharacter N hu S)
    (eisensteinToB_star N) (subsetCharacter_neg N hu S) x

theorem subsetEvaluation_norm (N : ℕ) [NeZero N] (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing)) :
    subsetEvaluation N hu S x * star (subsetEvaluation N hu S x) = ((N : ℕ) : B N) := by
  simpa only [subsetEvaluation, map_natCast] using evaluate_norm_eq
    (eisensteinToB N).toRingHom (subsetCharacter N hu S)
    (eisensteinToB_star N) (subsetCharacter_neg N hu S) x ((N : ℕ) : EisensteinRing) hx

theorem subsetEvaluation_ne_zero (N : ℕ) [NeZero N] (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing)) :
    subsetEvaluation N hu S x ≠ 0 := by
  intro h
  have hn := subsetEvaluation_norm N hu S x hx
  rw [h, zero_mul] at hn
  exact (Nat.cast_ne_zero.mpr hu.ne') hn.symm

/-- Every complex embedding has the required norm. Compatibility with conjugation
is derived from the actual adjoined roots in `map_star_B`. -/
theorem subsetEvaluation_map_norm (N : ℕ) [NeZero N] (hu : 0 < N)
    (S : Finset (PrimeIndex N)) (x : Elem EisensteinRing N)
    (hx : x * ringStar x = scalar N ((N : ℕ) : EisensteinRing))
    (σ : B N →+* ℂ) : ‖σ (subsetEvaluation N hu S x)‖ = Real.sqrt (N : ℝ) := by
  have hn := congrArg (fun z : B N => ‖σ z‖) (subsetEvaluation_norm N hu S x hx)
  simp only [map_mul, map_star_B, norm_mul, norm_star, map_natCast] at hn
  rw [RCLike.norm_natCast] at hn
  rw [← hn, Real.sqrt_mul_self (norm_nonneg _)]

end CirculantHadamard.PrimeCharacterEvaluationsW

end

end OAI
