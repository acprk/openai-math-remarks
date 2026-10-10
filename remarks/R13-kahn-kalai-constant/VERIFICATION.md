# R13 verification record (2026-10-09)

## 1. Primary computation (topt-lean session; `scripts/primary/kk176_verify.py`)
| Item | Command | Result |
|---|---|---|
| Main | `kk176_verify.py 6.15 1.4 0.4 6` | B = 0.402469, K = 7.20698, λ = 17.6514, ar1 = 1.68017, ar2 = 3.40960 → C1 = 100 183.8, C2 = 203 304.5 |
| Positive control | `kk176_verify.py 50 10 1 16 2/3` | B = 1.234·10⁻⁴ (the paper claims < 1/4) |

Arithmetic: mpmath intervals plus exact `Fraction` binomial tails. For s > 150 it uses ⌈s/p⌉ (Kaas–Buhrman). For C2 it uses a crude tail 2/p + 1/16 for n ≥ 2¹⁶; this crude tail is what fixes ar2.

## 2. Independent re-implementation (adversarial reviewer; `scripts/independent/kk176_indep.py`)
This was written from the manuscript and shares no code with the primary script. The reviewer had read the primary script beforehand to learn its parameter conventions. It uses a different library (Sage RealBallField/Arb, 256 bits), exact integer binomial tails, and a tighter but still rigorous prefix-sum bound S(m) = Σᵢ⌊m/gⁱ⌋.
| Item | Result |
|---|---|
| B (our S(m)) | 0.3959071 |
| B (primary's relaxation ⌊m/g⌋·g/(g−1), recomputed independently) | 0.40246859389831449…, agrees with primary to all printed digits |
| K = ∏_{v≥0}(1−e^{−0.4·2^v})^{−1} | 7.2069771 (agrees) |
| λ = 7⁷/6⁶ check | 2·C(7r,r) < λ^r exact for r ≤ 3000; Robbins factor 2√(7/(12π)) < 0.87 |
| N(s), s = 1..8 (target 1/2) | 1, 3, 5, 6, 8, 10, 11, 13 (agrees with primary for s ≤ 6); N(s) ≤ ⌈s/p⌉ checked for every s ≤ 200 |
| sup N(s)/s | 5/3 (s = 3); tail for s > 200 is ≤ 1/p + 1/201 = 1.660 |
| C1 | **99 378.34** (primary 100 183.8: both are valid upper bounds) |
| C2 (sup over s ≤ 400 at the least n with each s, tail 2/p + 1/200) | **197 984.8** (primary 203 304.5) |
| Target 2/3 (exact for s ≤ 200, Chernoff–KL in ball arithmetic for s ≤ 5000, Hoeffding with ε = 0.014 beyond) | C1 ≤ 119 254.1, C2 ≤ 212 768.1 |
| Positive control at the paper's parameters (50, 10, 1, 16) | B = 1.2341·10⁻⁴ (identical to primary), K = 1.8643 < 4 |

The note states C1 = 10⁵ and C2 = 2·10⁵. Both are above the independent certificates. The primary's C1 = 100 184 is slightly above 10⁵, so the headline value relies on the tighter independent bookkeeping.

## 3. Adversarial review of the proof (by hand)
Verdicts on the four questions in the brief:

**(a) Are the listed constraints really all the places where constants are fixed? CONFIRMED.**
We re-read `resampling.tex` (148 lines), `trees.tex` (279), `graphs.tex` (127) and `introduction.tex` (189) line by line. Every number in the proof traces to one of these sites:
1. local lemma margin e^{−50}, the Z-cut e^{−10m}, and the rounding to e^{−9m} / e^{11m} (`resampling.tex` L63, L68, L117–147): parameters L and α;
2. pruning ε = e^{−m}, the potential recurrence e^{−8m}/e^{12m}, and the geometric bound Σ e^{−7m} < 1/4 (`trees.tex` L98–168): parameters β and g;
3. the spread-inflation product < 4 (L244–252): parameter β;
4. attempts 4s, Markov, 2/3, and r* ≤ 4sρ (L257–277): target and N(s);
5. the hierarchy factor 17, 128 = σ/q, and the bound (17e/128)^r (`graphs.tex` L9, L60–72): parameters g and λ;
6. the final step 1+log₂h ≤ 3log₂n (`graphs.tex` L123–125): replaced by the exact sup over n.

We found no hidden "sufficiently large/small" requirement. The only side conditions are m ≥ 1 (capacities are positive integers), 0 < ρ < 1 (if ρ ≥ 1 the bound is trivial because of the min{1, ·}), distinct capacities (implied by m_{j+1} ≥ g·m_j with g ≥ 2), ⌊gx⌋ ≥ g⌊x⌋ (true for integer g), and n ≥ |V(H)|. The tail rates r1 = 0.6396 and r2 = 0.2215 are positive, as both scripts assert. One remark: the primary's "θ = 1/2 fixed" is not a constraint. It is the manuscript's halving rule, which fixes s = ⌊log₂ℓ⌋+1, and we left it unchanged.

**(b) Is it legitimate to change the success probability from 2/3 to 1/2? CONFIRMED.**
The definition (`introduction.tex` L17–19) reads: "pc(n,H) = inf{p∈[0,1] : P(N(G(n,p),H)>0) ≥ 1/2}". Lean `criticalThreshold` (`SecondKahnKalai.lean`) reads: `sInf {p | p ∈ Icc 0 1 ∧ (1/2) ≤ containmentProbability n p H}`. Probability ≥ 1/2 at density r therefore gives p_c ≤ r. The 2/3 in Theorem 3.1 is slack. For s > 150 (primary) or s > 200 (ours), the bound uses the Kaas–Buhrman theorem (every median of Bin(N,p) is ≥ ⌊Np⌋), which is valid for target 1/2. We also give a 2/3 variant with its own rigorous tail (§2) as a cross-check.

**(c) Are the λ inequality and the convergence of K correct? CONFIRMED.**
- λ: Robbins gives C(n,k) ≤ √(n/(2πk(n−k)))·nⁿ/(kᵏ(n−k)^{n−k}), because r_n − r_k − r_{n−k} < 1/(12n) − 1/(12k+1) − 1/(12(n−k)+1) < 0. Hence 2C(7r,r) ≤ 2√(7/(12πr))λ^r < 0.87λ^r. The reduction "r ≥ h/7 ⇒ C(h,r) ≤ C(7r,r)" is valid because C(·,r) is monotone. Exact check for r ≤ 3000.
- K: the dyadic replacement is new compared with the manuscript, so we checked it separately. The capacities met by one level are ⌊ℓ/2^u⌋ for u < t, all ≥ 1. Because ⌊y⌋ ≥ 2^v⌊y/2^v⌋, the v-th smallest is ≥ 2^v, and (1−e^{−βm})^{−1} decreases in m. The infinite product converges: the code truncates once the next term is below 10⁻⁴⁰ and bounds the remainder by 4e^{−β2^v}.
- Stochastic domination by Bin(N,p): induction on N, conditioning on the first attempt. Monotonicity of the binomial tail in s is used.

**(d) Independent recomputation of C1 and C2 with a positive control: CONFIRMED** (see §2).

**Overall R13 verdict: CONFIRMED.**

## 4. Duplicate search (2026-10-09)
- GitHub: `gh search repos` for "kahn kalai", "second kahn kalai" and "openai math" (about 35 community repos updated 10/06–10/09, among them aster2024/openai-math-followups, Rogerhu12/openai-math-extensions, mahrud/openai-math-comments, fortnow/AI-Papers, piiiico/openai-math-counterexamples, Kohze/…review). Tree listings of the six most relevant repos contain no Kahn–Kalai, 176 or constant files. Issue search found nothing.
- arXiv/web: Tran 2609.20546 (v2, 2026-09-24) proves q_f(F_H) ≤ C q_H log(2e(H)) with an unspecified C. Dubroff–Kahn–Park 2508.14269 say they "make no effort to optimize constant factors", and 2510.02005 covers cliques, cycles and trees with universal L. G. Kalai's blog post (2026-10-07) and Fortnow's TCS selection page mention family 176 but give no constant.
- Conclusion: we found no published explicit constant other than 2048e^50.

## GAPs
- No human review. No Lean (by policy, this is a re-parametrisation).
- We did not independently re-prove the manuscript's correctness; the premise is stated.
- The grid search over (L, α, β, g) was coarse, so the constant is not claimed to be optimal for the method.
