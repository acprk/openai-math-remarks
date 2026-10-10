import OAI.LinearAlgebra.CirculantHadamard.AlternatingField
import OAI.LinearAlgebra.CirculantHadamard.AlternatingLocalization
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.CyclotomicTorsion
import OAI.LinearAlgebra.CirculantHadamard.Eisenstein.AlternatingLocalization
import OAI.LinearAlgebra.CirculantHadamard.ResidueTorsion

namespace OAI

/-!
# Prop″ (ℤ[ω] case): the alternating character product, norm N with 3 ∤ N

Generalization of `OAI.LinearAlgebra.CirculantHadamard.AlternatingTorsion`,
which is stated only for `x : Elem EisensteinRing (u ^ 2)` with `x * x^* = u ^ 2`.
Here `x : Elem EisensteinRing N` (i.e. `ℤ[ω][ZMod N]`), `x * x^* = N`, `N > 1`, `3 ∤ N`
(N may be even).
Conclusion: the alternating product `Δ(x) = ∏_S χ_S(x)^{(-1)^{|S|}}` is a root of
unity whose order is a power of `p` for every prime `p ∣ N`.

The index is the full set of prime divisors of N. The local comparison and
all-maximal-ideal argument are applied to the actual partial evaluations of
the given constant-norm group-ring element. No edge-comparison, integrality,
complex-norm, or torsion conclusion is assumed.
-/

noncomputable section

namespace CirculantHadamard.EisensteinNorm

open CyclotomicRingsW PrimeComponentsN PrimeCharacterEvaluationsW
open AlternatingLocalization AlternatingLocalizationW

/-- The single fixed field element used for every prime-direction comparison. -/
abbrev alternatingValue (N : ℕ) (hu : 0 < N)
    (x : CyclicRing.Elem EisensteinRing N) : FractionRing (B N) :=
  alternatingProduct (Finset.univ : Finset (PrimeIndex N))
    (fun S => algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu S x))

/-- Denominator removal and all complex homomorphisms are discharged for the
actual character values before invoking Kronecker's criterion. -/
theorem exists_integral_alternating_root_of_unity (N : ℕ) [NeZero N]
    (hu : 1 < N) (huodd : ¬ 3 ∣ N) (x : CyclicRing.Elem EisensteinRing N)
    (hx : x * CyclicRing.ringStar x =
      CyclicRing.scalar N ((N : ℕ) : EisensteinRing)) :
    ∃ z : B N,
      algebraMap (B N) (FractionRing (B N)) z =
        alternatingValue N (lt_trans Nat.zero_lt_one hu) x ∧
      ∃ n : ℕ, 0 < n ∧ z ^ n = 1 := by
  classical
  let hu0 : 0 < N := lt_trans Nat.zero_lt_one hu
  obtain ⟨z, hz⟩ := exists_integral_alternatingProduct N hu0 huodd x hx
  refine ⟨z, hz, BW_root_of_unity N z ?_⟩
  intro σ
  have ht : (Finset.univ : Finset (PrimeIndex N)).Nonempty := by
    obtain ⟨p, hp⟩ := primes_nonempty hu
    exact ⟨⟨p, hp⟩, Finset.mem_univ _⟩
  have hne : ∀ S ⊆ (Finset.univ : Finset (PrimeIndex N)),
      algebraMap (B N) (FractionRing (B N)) (subsetEvaluation N hu0 S x) ≠ 0 := by
    intro S _ hS
    apply subsetEvaluation_ne_zero N hu0 S x hx
    apply IsFractionRing.injective (B N) (FractionRing (B N))
    simpa only [map_zero] using hS
  exact norm_hom_integral_alternatingProduct
    (algebraMap (B N) (FractionRing (B N)))
    (IsFractionRing.injective (B N) (FractionRing (B N)))
    Finset.univ ht (fun S => subsetEvaluation N hu0 S x) z hne hz σ (Real.sqrt (N : ℝ))
    (Real.sqrt_pos.mpr (by exact_mod_cast hu0)) (fun S _ => subsetEvaluation_map_norm N hu0 S x hx σ)

/-- The same alternating product is killed by a power of every actual prime
dividing N. The maximal ideal and characteristic-p residue field are built
from the concrete integral coefficient ring. -/
theorem alternatingProduct_prime_power_torsion (N : ℕ) [NeZero N]
    (hu : 1 < N) (huodd : ¬ 3 ∣ N) (x : CyclicRing.Elem EisensteinRing N)
    (hx : x * CyclicRing.ringStar x =
      CyclicRing.scalar N ((N : ℕ) : EisensteinRing))
    (p : PrimeIndex N) :
    ∃ a : ℕ, (alternatingValue N (lt_trans Nat.zero_lt_one hu) x) ^ (p.val ^ a) = 1 := by
  classical
  let hu0 : 0 < N := lt_trans Nat.zero_lt_one hu
  obtain ⟨z, hz, hfinite⟩ := exists_integral_alternating_root_of_unity N hu huodd x hx
  have hp : p.val.Prime := prime_of_mem p.property
  obtain ⟨m, hm, _, hpm⟩ := exists_maximal_ideal_over_nat_prime_of_charZero (B N) hp
  let M : MaximalSpectrum (B N) := ⟨m, hm⟩
  let : Algebra (B N) (AtMaximal (B N) M) := BW_local_algebra N M
  obtain ⟨v, hv, hvres⟩ :=
    exists_local_unit_alternatingProduct N hu0 huodd x hx p M hpm
  have hlocal : algebraMap (B N) (AtMaximal (B N) M) z =
      (v : AtMaximal (B N) M) :=
    local_representative_eq (FractionRing (B N)) M z v (hz.trans hv.symm)
  let : CharP (IsLocalRing.ResidueField (AtMaximal (B N) M)) p.val :=
    residue_charP M p.val hp hpm
  let f : B N →+* IsLocalRing.ResidueField (AtMaximal (B N) M) :=
    (IsLocalRing.residue (AtMaximal (B N) M)).comp
      (algebraMap (B N) (AtMaximal (B N) M))
  have hres : f z = 1 := by
    change IsLocalRing.residue (AtMaximal (B N) M)
      (algebraMap (B N) (AtMaximal (B N) M) z) = 1
    rw [hlocal]
    exact hvres
  obtain ⟨a, ha⟩ := exists_prime_power_of_residue_one f z hp hres hfinite
  refine ⟨a, ?_⟩
  have h := congrArg (algebraMap (B N) (FractionRing (B N))) ha
  simpa only [map_pow, map_one, hz] using h

/-- The alternating product has a positive torsion exponent prime to `3`. -/
theorem alternatingProduct_coprime_three_torsion (N : ℕ) [NeZero N]
    (hu : 1 < N) (huodd : ¬ 3 ∣ N) (x : CyclicRing.Elem EisensteinRing N)
    (hx : x * CyclicRing.ringStar x =
      CyclicRing.scalar N ((N : ℕ) : EisensteinRing)) :
    ∃ m : ℕ, 0 < m ∧ m.Coprime 3 ∧
      (alternatingValue N (lt_trans Nat.zero_lt_one hu) x) ^ m = 1 := by
  classical
  obtain ⟨p, hp⟩ := primes_nonempty hu
  obtain ⟨a, ha⟩ := alternatingProduct_prime_power_torsion N hu huodd x hx ⟨p, hp⟩
  refine ⟨p ^ a, pow_pos (prime_of_mem hp).pos _, ?_, ha⟩
  apply Nat.Coprime.pow_left
  apply (Nat.coprime_primes (prime_of_mem hp) Nat.prime_three).mpr
  intro h
  exact huodd (h ▸ dvd_of_mem hp)

/-- With two distinct prime directions the single fixed product is one. -/
theorem alternatingProduct_eq_one_of_two_primes (N : ℕ) [NeZero N]
    (hu : 1 < N) (huodd : ¬ 3 ∣ N) (x : CyclicRing.Elem EisensteinRing N)
    (hx : x * CyclicRing.ringStar x =
      CyclicRing.scalar N ((N : ℕ) : EisensteinRing))
    (p q : PrimeIndex N) (hpq : p ≠ q) :
    alternatingValue N (lt_trans Nat.zero_lt_one hu) x = 1 := by
  obtain ⟨a, ha⟩ := alternatingProduct_prime_power_torsion N hu huodd x hx p
  obtain ⟨b, hb⟩ := alternatingProduct_prime_power_torsion N hu huodd x hx q
  have hcop : Nat.Coprime (p.val ^ a) (q.val ^ b) :=
    Nat.coprime_pow_primes a b (prime_of_mem p.property) (prime_of_mem q.property)
      (fun h => hpq (Subtype.ext h))
  have hd := Nat.dvd_gcd (orderOf_dvd_of_pow_eq_one ha) (orderOf_dvd_of_pow_eq_one hb)
  have hd1 : orderOf (alternatingValue N (lt_trans Nat.zero_lt_one hu) x) ∣ 1 := by
    simpa only [hcop.gcd_eq_one] using hd
  exact orderOf_eq_one_iff.mp (Nat.dvd_one.mp hd1)

end CirculantHadamard.EisensteinNorm

end

end OAI
