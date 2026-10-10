# VERIFICATION — family 238 constants C_u and ε_n

All runs were made on 2026-10-09 with CPython 3.10.12, NumPy 2.2.6, mpmath 1.3.0 and SageMath 9.5 (route A
only). Each script was single-threaded, except `exact_tables.py`, which uses 3 processes.
`scripts/run_all.sh` reproduced all six expected outputs byte for byte; see §2.

## 1. Who checked what

| Item | Primary (route A) | Independent (route B) | Adversarial review | Status |
|---|---|---|---|---|
| Exact values Z_u(m*) for u = 1, 1/2, 3/10, 1/5, 3/20 | `exact_tables.py` (hook formula, m ≤ 55) + MPFI sums in `cu_bound.py` | `independent/indep_Cu.py`: branching-rule f, m ≤ 50, f^{−p/q} enclosed by exact integer q-th roots at 2⁻⁸⁰; Z₁(6) = 473/144 as an exact Fraction | see §3 | identical 10-digit enclosures |
| Peeling inequality | ballot-word Lemma A; unit test `check_lemma.py` (4507 partitions, 0 violations) | Lemma B: exact first-row hook identity + Chebyshev/AM–GM; unit test inside `indep_Cu.py` (exact rationals, 7313 shapes with m ≤ 24, 0 violations; identity exact) | see §3 | two different lemmas |
| Recursive region m > M0 | 3-index DP to M_E, then 1-D recursion to M_1 (floats inflated by 1+2⁻⁴⁰; MPFI weights) | 3-index DP with Lemma B weights to M_E, then 1-D recursion to M_1 (floats inflated by 1+2⁻³⁰; 2188 weights re-checked in `mpmath.iv`, 0 violations) | see §3 | all majorants are below Z_u(m*) |
| Tail m > M_1 | Lemma "tail" (Catalan series), T(M_1) in MPFI | Lemma T (Robbins bound, three ranges), T(M_1) in `mpmath.iv` | see §3 | both close |
| ε_n exponent | `eps_n_exponent.py` (mpmath.iv) | `independent/indep_eps.py` (integer comparisons 2^{495q−p} vs 31^{99q}) | see §3 | e* ∈ (3.03456, 3.03457) |
| C_{1/10} ≤ 1267.27 | route A (M_E = 450, M_1 = 2500), from the internal notes | **not reproduced**: route B's 1-D step blows up at the M_E it can afford | — | **not claimed** |
| C_{1/10} ≥ 562.0737 | route A exact region | route B exact region | see §3 | confirmed |

Provenance. Route A is the earlier internal computation (paperH notes, 2026-10-07). Before this build its
tail had only one group of checkers. The build session wrote route B from scratch to supply a second check
by a different method, as required for inclusion. Route B imports nothing from route A.

## 2. Commands and key outputs

```
sage -python exact_tables.py 55        -> done 55 2868521 distinct (r,c,f) classes
python3 check_lemma.py                 -> checked 4507 violations 0 equality cases 163
sage -python cu_bound.py <u> 55 ME M1  (route A)
  u=1    m*=6   Z in [3.2847222222, 3.2847222223]  max m>55: 2.039271 (m=56)  T(400)  <= 2.054831
  u=1/2  m*=8   Z in [5.9061465497, 5.9061465498]  max m>55: 2.428254 (m=56)  T(400)  <= 2.280570
  u=3/10 m*=14  Z in [13.0732169915, 13.0732169916] max m>55: 3.871009 (m=56) T(400)  <= 3.083811
  u=1/5  m*=22  Z in [34.5870288005, 34.5870288006] max m>55: 11.561959 (m=56) T(400) <= 10.322035
  u=3/20 m*=30  Z in [89.3331732895, 89.3331732896] max m>55: 69.593086 (m=56) T(1000)<= 10.822273
python3 independent/indep_Cu.py        (route B)
  u=1    Z(6)  in [3.2847222222, 3.2847222223]   DP max 2.043553 (m=51), 1-D max 2.017246, T(600)  <= 2.003361
  u=1/2  Z(8)  in [5.9061465497, 5.9061465498]   DP max 2.462676, 1-D max 2.336627, T(600)  <= 2.093084
  u=3/10 Z(14) in [13.0732169915, 13.0732169916] DP max 4.053970, 1-D max 3.857438, T(600)  <= 2.646706
  u=1/5  Z(22) in [34.5870288005, 34.5870288006] DP max 11.753292, 1-D max 8.059078, T(1000) <= 4.725881
  u=3/20 Z(30) in [89.3331732895, 89.3331732896] DP max 52.276652, 1-D max 46.623290, T(1000) <= 26.745931
  u=1/10 Z(50) in [562.0737331731, 562.0737331732]
  ALL VERIFIED
python3 eps_n_exponent.py   -> 198*log2(32/31) = 9.069130543; e* = 3.034565272; e^{-n/256} term subdominant for d >= 15
python3 independent/indep_eps.py -> e* > 3.03456: True; e* > 3.03457: False
```
Both routes now print the same outward-rounded 10-digit enclosures. (Before the adversarial review,
`cu_bound.py` rounded both endpoints to nearest, so it printed e.g. [5.9061465498, 5.9061465498]. The
internal MPFI interval was correct, but the printed lower endpoint was above the true value 5.906146549796….
The print statement now rounds outward. This is a formatting change only.)

`./run_all.sh` (2026-10-09, final full run after the review corrections; wall time 10 min 29 s):
```
OK   exact_tables
OK   check_lemma
OK   cu_bound
OK   eps_n_exponent
OK   independent/indep_Cu
OK   independent/indep_eps
```

## 3. Adversarial review

A separate AI session was told to refute the claim. It was not allowed to edit the remark.
**Verdict: CONFIRMED WITH CORRECTIONS.** It found nothing critical or major, and all the minor points are fixed.

What it checked:
- **Proofs, line by line.** Lemma A (ballot peeling, including the case k > L). Lemma B: the hook
  identity, AM–GM, the direction of Chebyshev's inequality, and the monotonicity of t·log(1+k/t). The DP
  recursions: conjugation fill, orientation factor 2, constraint ranges, ⌈√m⌉ = `isqrt(m−1)+1`. Both tail
  lemmas: range coverage, Robbins' bound in the MacWilliams–Sloane form, concavity of
  −(1−t)log(1−t), and all four monotonicity steps. The induction logic.
- **Code against proof.** The slicing `E[kp,1:bmax+1,1:cmax]` ↔ b = 1..bmax, c−1 = 1..cmax−1;
  `CS[k,min(L,k),min(L−1,k)]` = R(k,L). The tail hypotheses (3K ≤ M, u(K+1) ≥ 1+1/M, ρ < 1,
  M ≥ (1+u/2)/(u c_τ)) are all asserted in code. The integer-root enclosures are correct. Route A's sums have
  fewer than 2¹² terms, so its 1+2⁻⁴⁰ inflation is rigorous.
- **Third independent computation.** A Frobenius-formula enumeration for m ≤ 62 at 40 digits, with
  Σf² = m! checked. It reproduces all five 10-digit values and Z_{1/10}(50) = 562.07373317315…. Within
  m ≤ 62 the arg-max is unique and equals m* in every case. Z₁(6) = 473/144 was also checked by hand.
- **Majorants are ≥ the true values.** For m = 51..62, route B's majorants are at least the exact Z_u(m) for
  every u. The minimum relative margin ranges from 4.6·10⁻⁹ (u = 1, where Lemma B is almost exact on hooks)
  to 1.9·10⁻². Route A's margin for u = 3/10, m = 56..62, is at least 6.2·10⁻⁴.
- **Brute-force lemma test** over all 28,598 shapes with m ≤ 30, in exact rationals: Lemmas A and B (all
  three inequalities) have 0 violations.
- **ε_n.** It confirmed e* = 3.0345652717…. The threshold d = 15 is sharp: at d = 14 the e^{−n/256} term
  still dominates. The paper's O(n^{−5/2}) is correct but loose, not an error.

Corrections it requested, all applied:
1. The note claimed "τ = 1/10, K = 50". In fact the script picks the best K ≤ 50: K = 17, 44, 50, 50, 50. Fixed.
2. The note claimed route B's float margin exceeded the error "by >10³". The measured worst case was about
   230×, dominated by an fsum-accumulated ln n! table. The tables are now computed at 40 digits, and the
   output is unchanged byte for byte. The note now states the IEEE/libm assumption explicitly and says that
   route B's recursive region is not a formal interval computation.
3. Route A printed endpoints rounded to nearest. Fixed (see §2).
4. "Θ(n^{−3.0346})" is not literally true. Replaced by Θ(n^{−e*}), e* = 3.03456….
5. The description of route A's tail term in note §2 was inaccurate. Corrected.
6. T(M₁) for route A at u = 1/5 is now 10.322 (was rounded to 10.33).
8. The note now says that route B's independence is at the level of code, lemmas and arithmetic, and that
   the overall architecture is shared.

Two further points: (7) the expected outputs were still being generated when the review began, and they are
now frozen; (9) Gamburd 2006 could not be checked offline.

## 4. What is numeric and what is exact

- Exact: the integer values f_λ (both routes); Z₁(6) = 473/144 (Fraction); the exact-region enclosures in
  route B (integer roots). The ε_n exponent bounds in `indep_eps.py` are integer comparisons.
- Interval: route A uses MPFI for its exact sums, weights and tail; route B uses `mpmath.iv` for its tail and
  for the spot-check of weights.
- IEEE double with explicit outward inflation: both recursive regions. Route A inflates by 1+2⁻⁴⁰ per stage.
  Route B inflates by 1+2⁻³⁰ (ln ≈ 9.3·10⁻¹⁰) per weight and per stage. Its ln n! and H_n tables are
  computed at 40 digits, so the log-weight error is ≲ 5·10⁻¹³, a margin of more than 10³.
  This assumes standard IEEE-754 behaviour and glibc/NumPy `exp`/`log`/`log1p` errors of a few ulp. To guard
  against that assumption, route B re-checks a seeded random sample of 2188 weights in `mpmath.iv`; none
  were too small.
- Proof-level inputs that are not computed: Lemma A/B, the tail lemmas, Chebyshev's sum inequality, AM–GM,
  Robbins' bound binom(m,L) ≥ (8L(m−L)/m)^{−1/2} e^{mH(L/m)} (MacWilliams–Sloane Ch. 10, Lemma 7).

## 5. GAPs

- C_{1/10}: upper bound not independently confirmed, so not claimed.
- Growth of C_u as u → 0: only the elementary lower bound is stated.
- No Lean.
