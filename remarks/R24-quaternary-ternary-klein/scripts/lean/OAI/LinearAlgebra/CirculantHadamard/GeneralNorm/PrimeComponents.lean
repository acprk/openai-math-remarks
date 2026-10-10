import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Data.ZMod.QuotientRing
import Mathlib.Algebra.Ring.Parity
import Mathlib.Algebra.Group.Pi.Lemmas

namespace OAI

/-!
# Primary components of the cyclic group of an arbitrary order N

Generalization of `OAI.LinearAlgebra.CirculantHadamard.PrimeComponents`,
which is stated only for the square order `u ^ 2`. This copy replaces `u ^ 2`
by an arbitrary positive `N`; the proofs are otherwise unchanged.

The prime set is the literal finite set of prime factors of `N`. The exponent
of each component is `N.factorization p`.
Factorization and coprimality supply the hypotheses of the finite Chinese
remainder theorem. Its coordinates are the canonical reduction maps.

The generator of each primary component is the inverse image of its unit
coordinate vector. Its additive order is proved to be the full prime power,
and the complementary integer factor of N is proved prime to p.
-/

noncomputable section

namespace CirculantHadamard.PrimeComponentsN

open scoped BigOperators

def primes (N : ℕ) : Finset ℕ := N.primeFactors

abbrev PrimeIndex (N : ℕ) := {p : ℕ // p ∈ primes N}

def exponent (N p : ℕ) : ℕ := N.factorization p

def order (N p : ℕ) : ℕ := p ^ exponent N p

abbrev Component (N : ℕ) (p : PrimeIndex N) := ZMod (order N p.val)

abbrev PrimaryProduct (N : ℕ) := (p : PrimeIndex N) → Component N p

theorem mem_primes {N p : ℕ} : p ∈ primes N ↔ p.Prime ∧ p ∣ N ∧ N ≠ 0 :=
  Nat.mem_primeFactors

theorem prime_of_mem {N p : ℕ} (hp : p ∈ primes N) : p.Prime :=
  (mem_primes.mp hp).1

theorem dvd_of_mem {N p : ℕ} (hp : p ∈ primes N) : p ∣ N :=
  (mem_primes.mp hp).2.1

theorem primes_nonempty {N : ℕ} (hu : 1 < N) : (primes N).Nonempty :=
  Nat.nonempty_primeFactors.mpr hu

theorem prime_ne_two_of_odd {N p : ℕ} (hu : Odd N) (hp : p ∈ primes N) : p ≠ 2 := by
  intro h
  apply hu.not_two_dvd_nat
  simpa only [h] using dvd_of_mem hp

theorem odd_of_mem {N p : ℕ} (hu : Odd N) (hp : p ∈ primes N) : Odd p :=
  (prime_of_mem hp).odd_of_ne_two (prime_ne_two_of_odd hu hp)

theorem exponent_pos {N p : ℕ} (hu : 0 < N) (hp : p ∈ primes N) :
    0 < exponent N p := by
  exact (prime_of_mem hp).factorization_pos_of_dvd hu.ne' (dvd_of_mem hp)

theorem order_pos {N p : ℕ} (hp : p ∈ primes N) : 0 < order N p :=
  pow_pos (prime_of_mem hp).pos _

instance componentOrder_neZero (N : ℕ) (p : PrimeIndex N) : NeZero (order N p.val) :=
  ⟨(order_pos p.property).ne'⟩

theorem prime_dvd_order {N p : ℕ} (hu : 0 < N) (hp : p ∈ primes N) :
    p ∣ order N p := dvd_pow_self p (exponent_pos hu hp).ne'

theorem order_dvd_square (N p : ℕ) : order N p ∣ N :=
  Nat.ordProj_dvd N p

/-- This is the full primary factor: the next power does not divide N. -/
theorem next_power_not_dvd {N p : ℕ} (hu : 0 < N) (hp : p ∈ primes N) :
    ¬p ^ (exponent N p + 1) ∣ N :=
  Nat.pow_succ_factorization_not_dvd hu.ne' (prime_of_mem hp)

def cofactor (N p : ℕ) : ℕ := N / order N p

theorem order_mul_cofactor (N p : ℕ) : order N p * cofactor N p = N := by
  exact Nat.mul_div_cancel' (order_dvd_square N p)

theorem cofactor_not_dvd {N p : ℕ} (hu : 0 < N) (hp : p ∈ primes N) :
    ¬p ∣ cofactor N p :=
  Nat.not_dvd_ordCompl (prime_of_mem hp) hu.ne'

theorem cofactor_coprime {N p : ℕ} (hu : 0 < N) (hp : p ∈ primes N) :
    Nat.Coprime p (cofactor N p) :=
  Nat.coprime_ordCompl (prime_of_mem hp) hu.ne'

theorem prod_order_finset (N : ℕ) (hu : 0 < N) :
    (∏ p ∈ primes N, order N p) = N := by
  have h := Nat.prod_factorization_pow_eq_self hu.ne'
  simpa only [Nat.prod_factorization_eq_prod_primeFactors,
    primes, order, exponent] using h

theorem prod_order (N : ℕ) (hu : 0 < N) :
    (∏ p : PrimeIndex N, order N p.val) = N := by
  change (∏ p : (primes N), order N p.val) = N
  rw [Finset.prod_coe_sort]
  exact prod_order_finset N hu

theorem pairwise_coprime (N : ℕ) :
    Pairwise (fun p q : PrimeIndex N => Nat.Coprime (order N p.val) (order N q.val)) := by
  intro p q hpq
  exact Nat.coprime_pow_primes _ _ (prime_of_mem p.property) (prime_of_mem q.property)
    (fun h => hpq (Subtype.ext h))

/-- The actual finite Chinese remainder equivalence, with no chosen-map premise. -/
def crtRing (N : ℕ) (hu : 0 < N) : ZMod N ≃+* PrimaryProduct N :=
  (ZMod.ringEquivCongr (prod_order N hu).symm).trans
    (ZMod.prodEquivPi (fun p : PrimeIndex N => order N p.val) (pairwise_coprime N))

def crt (N : ℕ) (hu : 0 < N) : ZMod N ≃+ PrimaryProduct N :=
  (crtRing N hu).toAddEquiv

@[simp] theorem crt_intCast (N : ℕ) (hu : 0 < N) (z : ℤ) (p : PrimeIndex N) :
    crt N hu (z : ZMod N) p = (z : Component N p) := by
  change (crtRing N hu (z : ZMod N)) p = _
  exact congrArg (fun y : PrimaryProduct N => y p) (map_intCast (crtRing N hu) z)

@[simp] theorem crt_natCast (N : ℕ) (hu : 0 < N) (z : ℕ) (p : PrimeIndex N) :
    crt N hu (z : ZMod N) p = (z : Component N p) := by
  change (crtRing N hu (z : ZMod N)) p = _
  exact congrArg (fun y : PrimaryProduct N => y p) (map_natCast (crtRing N hu) z)

/-- Every CRT coordinate is literal reduction from ZModN to the full p-power. -/
theorem crt_apply (N : ℕ) (hu : 0 < N) (x : ZMod N) (p : PrimeIndex N) :
    crt N hu x p = (ZMod.cast x : Component N p) := by
  obtain ⟨z, rfl⟩ := ZMod.intCast_surjective x
  rw [crt_intCast, ZMod.cast_intCast (order_dvd_square N p.val)]

def primaryGenerator (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) : ZMod N :=
  (crt N hu).symm (Pi.single p (1 : Component N p))

@[simp] theorem crt_primaryGenerator (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    crt N hu (primaryGenerator N hu p) = Pi.single p (1 : Component N p) :=
  (crt N hu).apply_symm_apply _

@[simp] theorem crt_primaryGenerator_same (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    crt N hu (primaryGenerator N hu p) p = 1 := by
  rw [crt_primaryGenerator, Pi.single_eq_same]

@[simp] theorem crt_primaryGenerator_ne (N : ℕ) (hu : 0 < N)
    (p q : PrimeIndex N) (hqp : q ≠ p) :
    crt N hu (primaryGenerator N hu p) q = 0 := by
  rw [crt_primaryGenerator, Pi.single_eq_of_ne hqp]

/-- The chosen primary generator has exactly the full primary order. -/
theorem primaryGenerator_addOrderOf (N : ℕ) (hu : 0 < N) (p : PrimeIndex N) :
    addOrderOf (primaryGenerator N hu p) = order N p.val := by
  rw [primaryGenerator, AddEquiv.addOrderOf_eq]
  have h := addOrderOf_injective (AddMonoidHom.single (Component N) p)
    (Pi.single_injective p) (1 : Component N p)
  exact h.trans (ZMod.addOrderOf_one (order N p.val))

end CirculantHadamard.PrimeComponentsN

end

end OAI
