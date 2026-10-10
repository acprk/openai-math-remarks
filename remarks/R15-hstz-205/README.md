# R15: the staircase tensor square in the alternating group (HSTZ at n = 28, 36, 45, 66, 78)

**One sentence.** OpenAI proved Saxl's conjecture for S_n. The motivating conjecture of Heide–Saxl–Tiep–Zalesski (HSTZ) for **A_n** is still open. We show three things. First, by a known parity fact (Bessenrodt–Bowman 2023, Cor. 6.2), the staircase split character can only be a witness when m ≢ 2 (mod 4). Second, S_n-universality of a self-conjugate label does **not** pass to A_n; the first non-hook instance is (4,4,4,3) ⊢ 15. Third, a direct modular computation of all Kronecker coefficients g(ρ_m, ρ_m, λ) **proves HSTZ for A_28, A_36, A_45, A_66 and A_78**. Previously published checks stopped at n ≤ 25.

## Precise statement (Before / After)

| | Before | After (this remark) |
|---|---|---|
| S_n | Saxl's conjecture and the tensor-square conjecture for S_n, n ∉ {2,4,9} (OpenAI, family 205) | unchanged |
| A_n, HSTZ conjecture | open. Computer checks: O'Brien n ≤ 17 (HSTZ 2013, Remark 1.3(3)); Bessenrodt n ≤ 25 (Alg. Comb. 2018, §5) | additionally **n = 28, 36, 45, 66, 78** (m = 7, 8, 9, 11, 12), with witness a constituent χ₊ of χ^{ρ_m}\|_{A_n} |
| Staircase as A_n witness | χ₊ is real iff m ≢ 2 (mod 4) (Bessenrodt–Bowman, arXiv:2202.03066, Thm 6.1/Cor 6.2) | recorded consequence: impossible for m ≡ 2 (mod 4), since χ₊² then misses the trivial character |
| Does S_n-universality of λ = λ′ give A_n-universality of χ^λ₊? | not addressed | **No.** For λ = ((q+s)^q, q^s), χ₊ is never universal. (3,1,1) ⊢ 5 and (4,4,4,3) ⊢ 15 are S_n-universal; for (4,4,4,3), χ₊² misses the 14-dimensional standard character of A_15 |

Theorem (after). For m ∈ {3,4,5,7,8,9,11,12}, every irreducible character of A_{m(m+1)/2} is a constituent of χ₊², where χ₊ is either irreducible constituent of the restricted staircase character. The cases n = 6, 10, 15 were already known.

## Attribution

Saxl's conjecture and the S_n tensor-square theorem are **OpenAI's work** (family 205). The preprints are `preprints/A-Cyclic-Polytabloid-Proof-of-Saxls-Conjecture-September-24-2026/` and `preprints/Universal-Tensor-Squares-for-Symmetric-Groups-September-24-2026/` in https://github.com/openai/math at commit `fd4aeeb2e` (Lean: `lean/ComparatorChallenges/Saxl.lean`). The square formula (Proposition 1 of the note) is classical in substance: Bessenrodt–Behns 2004; Pak–Panova–Vallejo, Main Lemma. We include it because we need the exact multiplicities. The realness criterion, the m ≡ 2 (mod 4) fact and g(λ,λ,(n−1,1)) = #distinct parts − 1 are in Bessenrodt–Bowman (arXiv:2202.03066). This remark contributes:
* the A_n bookkeeping (Prop. 1(3) and the criterion Cor. 2);
* the family of non-universal split characters (Prop. 4);
* the computer proof for five new n;
* an analysis of how far the multiplicity bounds hidden in OpenAI's cyclic module W_m reach (§4 of the note).

**Theorem 5 does not depend on OpenAI's proof.** The direct computation also certifies g > 0 wherever a = 0. Section 4 records a second, conditional route through 205's Propositions 3.2/4.2/5.3.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are not human peer review.

## Proof outline

1. **Square formula.** Let λ = λ′, ε = (−1)^{(n−d)/2}, and a_ν = χ^ν(ĥ(λ)). The multiplicity of Res χ^ν (ν ≠ ν′) in χ₊² is (g(λ,λ,ν) + εa_ν)/2. For ν = ν′ ≠ λ it is (g + εa_ν)/4. For χ_± it is (g(λ,λ,λ) + 1 ± 2ε)/4. Hence, for ε = +1, χ₊ is universal iff g ≠ −a whenever a ≤ 0, and g(λ,λ,λ) ≠ 1.
2. **ε for staircases.** (n−d)/2 ≡ d(m−1) (mod 2), so ε = −1 iff m ≡ 2 (mod 4).
3. **Counterexample family.** g(λ,λ,(n−1,1)) = (#distinct parts) − 1 = 1. Moreover ĥ has no part 1, so a = −1, and the multiplicity is 0.
4. **Certificate.** Since ρ_m is a 2-core, Σ_λ g s_λ = Σ_{μ odd} χ^ρ(μ)² z_μ⁻¹ p_μ. `kron_stair.cpp` evaluates this mod p = 2³¹−1 in the Schur basis, for all p(N_m) partitions (12,132,164 at m = 12). It uses Murnaghan–Nakayama on β-sets, a head/tail split of μ, and Horner on the tail trie. It also computes the exact a(λ). No λ with g ≡ −a (mod p) and a ≤ 0 occurs, and g(ρ,ρ,ρ) ≢ 1.

## Reproduce

From `scripts/` (g++ ≥ 9 with OpenMP, GAP 4, SageMath 9.5, Python 3):

```
./run_all.sh            # everything; m=12 needs ~24 min (3 threads) and ~11 GB RAM
# or individually:
g++ -O3 -march=native -fopenmp -o kron_stair kron_stair.cpp
./kron_stair 9 2147483647 7 3 out_m9.txt           # -> expected_output/kron_stair_m9.txt (seconds)
./kron_stair 12 2147483647 9 3 out_m12.txt         # -> expected_output/kron_stair_m12.txt
gap -q -o 8g an_square_check.g                     # -> expected_output/an_square_check.out
python3 crt_exact.py 9 ks_m9_P*.txt.full           # -> expected_output/crt_exact_m9.txt (needs 4 prime runs)
sage w205.sage hooks 16                            # -> expected_output/w205_hooks_m1-16.txt
```

Key lines to check in `expected_output/kron_stair_m{7,8,9,11,12}.txt`: `eps=1`, `... g == |a| (mod P): 0`, `a==0 and g == 0 (mod P): 0`, `g(rho,rho,rho) mod P ≠ 1`, and both `check ... : OK`.

## What is NOT claimed

* HSTZ for any non-triangular n, for m ≡ 2 (mod 4) (n = 21, 55, …), or for any infinite family. Section 4 of the note explains why we see no way to carry OpenAI's cyclic-module bounds to all m (heuristic obstructions, not an impossibility proof). On the hook (m, 1^{N−m}) they give multiplicity exactly 1 for every m ≥ 5; the true value is, e.g., 84465 at m = 8. Nothing in the construction controls the margin g − |a|.
* For n = 66 and 78 the certificate is a single computation modulo p = 2³¹ − 1, with internal checksums: Σ g·dim = dim(ρ)², Σ g·a = 1, and threshold-independence at m = 10, 11. It is not an exact-integer computation. For n ≤ 45 we computed the exact integers (four primes + CRT, exact identity Σ g·dim = dim(ρ)² over ℤ).
* No Lean formalization (per project rules this is a computational remark).
* Proposition 1 is not new in substance; only the bookkeeping is spelled out.
* Novelty: to our knowledge (search of GitHub, including the openai/math commentary aggregators cloned locally — eliottcassidy2000, VibeMathed, satyajitghana, supsup, mahrud, piiiico, Rogerhu12, aster2024, zjukop3, SafeVerifyAgent, blisspixel, jkbennitt — GitHub code/repo search, and arXiv up to 2026-10-09), no one has stated HSTZ for these n or the (4,4,4,3)-type failure. arXiv:2609.17319 (15 Sep 2026) still calls the A_n case open.

## Premises

* Correctness of the C++ certificate `kron_stair.cpp`. It is cross-validated against GAP for m ≤ 6, by exact identities for m ≤ 9, by mod-p identities for m ≤ 12, and against an earlier independent hook computation.
* Standard character theory of A_n (the split-character values on the classes of type ĥ(λ)).
* Theorem 5 does **not** assume OpenAI's 205 propositions. The conditional route of §4 does, together with the correctness of the LR code (Sage lrcalc).
