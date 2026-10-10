# R19 — Jacobsthal's function in F_q[x]: K_q(1) = q and K_q(2) = 2q − 1

## Premise and confidence (read first)

- **What this is.** OpenAI's family 021 bounds the *integer* Jacobsthal function. Here we record the first two values of the analogous quantity in the polynomial ring F_q[x], for every prime power q. Before → After: *not stated anywhere we could find* → **K_q(1) = q, K_q(2) = 2q − 1**, with the extremal configurations characterised for q ≥ 3, plus a short table of exact small values.
- **This is a direct analogue, and it is elementary.** Nothing from OpenAI's integer proof is used. The m = 1 case is immediate, and m = 2 is a two-line incidence count of lines in F_q². We record it as a remark because the function-field version appears not to have been written down. We make no claim about the integer problem or about the growth of K_q(m).
- **Confidence: high.** Two independent exact computations agree on every tabulated value: CP-SAT, and a separate branch-and-bound with its own field arithmetic. A SageMath check factors the explicit construction. An adversarial AI session found no error (see `VERIFICATION.md`).

## Statement

Let V_m = {g ∈ F_q[x] : deg g < m}. Define K_q(m) as the least ω(M) (number of distinct monic irreducible factors) such that, for some f, every element of the "interval" f + V_m shares a factor with M. Equivalently, the function-field Jacobsthal function h_q(k) = max{j(M) : ω(M) ≤ k} satisfies h_q(k) ≤ m ⇔ K_q(m) > k.

| | Before | After |
|---|---|---|
| K_q(1) | — | q |
| K_q(2) | — | 2q − 1 (every prime power q) |
| h_q(k) | — | 1 for k ≤ q−1, 2 for q ≤ k ≤ 2q−2, ≥ 3 at k = 2q−1 |
| equality at m = 2, q ≥ 3 | — | exactly: all q linear primes with concurrent residue lines (x^q − x divides f after a shift within the coset), plus q − 1 primes of degree ≥ 2, each dividing a different one of the q − 1 uncovered elements, and no other prime factor |

Exact small values (two independent programs; q = 4 by branch-and-bound and Sage only):

| q \ m | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 2 | 2 | 3 | 4 | 6 | 9 |
| 3 | 3 | 5 | 8 | | |
| 4 | 4 | 7 | | | |
| 5 | 5 | 9 | | | |
| 7 | 7 | 13 | | | |

K_2(4) = 6, K_2(5) = 9 and K_3(3) = 8 exceed m(q−1)+1, so the pattern of m ≤ 2 does not continue.

## Attribution

The integer problem and its solution are OpenAI's work: family 021 of
[github.com/openai/math](https://github.com/openai/math), commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/A-quadratic-bound-for-Jacobsthals-function-September-25-2026/`
(h(k) ≪ k²/(log log 3k)²). This remark studies only the F_q[x] analogue in the smallest cases. It does
not use and does not modify any line of that preprint.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are
not human peer review.

## Proof sketch

Full proofs are in `note.pdf`.
1. **Covering model.** A prime of degree ≥ m divides at most one element of f + V_m, because a difference of two elements is nonzero of degree < m. Conversely, by the CRT on f, any one prescribed element can be made divisible by a fresh such prime. A prime P of degree < m covers the residue class g ≡ −f (mod P), and these residues are independent and arbitrary (CRT again). So K_q(m) is a finite set-cover minimum.
2. **m = 2.** Identify g = c_0 + c_1x with (c_0, c_1) ∈ F_q². The prime x − a covers the line c_0 + a·c_1 = −f(a), and different a give non-parallel lines. ℓ lines cover at most q + (ℓ−1)(q−1) points, so the cost is ≥ q² − 1 − ℓ(q−2) ≥ 2q − 1. The bound is attained by f = x^q − x, with all q linear primes and one irreducible factor of each x^q − x + c, c ≠ 0.
3. **Equality (q ≥ 3).** All q lines are needed, and each new line must meet the union of the previous ones in exactly one point, which forces concurrency.

## Reproduce

```bash
cd scripts && ./run_all.sh     # python3 + ortools + python-sat (primary), SageMath (check); single thread; < 1 min
```
It prints `OK` for each of three outputs.
- `primary/run_table.py`: CP-SAT model (from the internal review's `fqx_Kqm.py`, single worker), q ∈ {2,3,5,7}.
- `independent/indep_Kqm.py`: exhaustive branch-and-bound with its own GF(q) (q ∈ {2,3,4,5,7}), its own polynomial reduction and its own irreducibility test. No solver.
- `independent/construction_check.sage`: factors x^q − x + c and checks that all q² elements of f + V_2 are non-coprime to M, for q ∈ {2,3,4,5,7,8,9}. It also enumerates all residue vectors for q ∈ {3,4,5} and confirms that the optimum 2q − 1 is attained exactly by affine (concurrent) ones.

## Scope and limitations (what is NOT claimed)

- Nothing about the integer Jacobsthal function, and nothing about OpenAI's proof.
- No formula for m ≥ 3 and no asymptotics in m. The integer analogue of the latter would be a function-field version of Iwaniec or of OpenAI's bound, which we do not attempt.
- The table beyond m = 2 is computational: exact, and agreed by two programs, but only for the listed (q, m).
- No Lean.
- **Novelty.** To our knowledge (search of GitHub and arXiv up to 2026-10-09: web searches for "Jacobsthal function" with "F_q[x]", "F_q[t]", "polynomial ring", "function field", and "covering" variants; GitHub repository searches for "jacobsthal", "jacobsthal finite field", "jacobsthal polynomial"), the function-field Jacobsthal function has not been studied. Nearby but different:
  - Kalmynin–Konyagin, arXiv:2302.00459, study shifts of integer polynomial sequences;
  - arXiv:2310.02495 studies powerfree polynomials in short intervals over F_q.
  - Covering systems of F_q[x] (Azlin 2011; Li–Wang–Wang–Yi 2024) concern distinct-moduli coverings of the whole ring, which is a different problem.
  - Jamison (1977) and Brouwer–Schrijver (1978) on affine blocking sets are close in spirit to the m = 2 count. We cite them in the note.
  - The community repository [jdb19937/jacobsthal](https://github.com/jdb19937/jacobsthal) (HEAD a91ffe3, 2026-10-08) improves the *integer* bound of family 021 by a factor (log k)^{2/25} and does not consider F_q[x]. We checked its full TeX source.
