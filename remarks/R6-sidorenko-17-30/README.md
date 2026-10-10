# R6 — Sidorenko counterexample: 35 vertices / 66 edges → 17 / 30 (conditional)

> **Premise and confidence (read first).**
> - The result is **conditional**. It assumes that Lemma 3.1 (rank-layer limit) of the original paper, together
>   with its proof in Sections 4–7, holds for the 7-point complex below by the same argument after the
>   complex-dependent constants are re-instantiated. We did not re-derive Sections 4–7 for this complex.
> - The original 35-vertex theorem is formalised in Lean, but that formalisation is hard-coded to the
>   13-point / 22-face complex (`Fin 13`, `Fin 22`, `Fin 33`). **Our instance carries over only at the
>   level of the TeX argument and has no machine certificate. Its reliability is lower than the original's.**
> - Checked exactly and unconditionally: conditions (1)–(3) for the new complex, its colour-preserving
>   homomorphisms, the coefficients C(f), and the generalised pinning lemma (a short proof, given below).

**One-line result.** Drop condition (4) (discrete colour refinement) from OpenAI's construction and replace it by
"every colour-preserving homomorphism of H is an automorphism". Then the point–face incidence graph H₇ of a
7-point, 10-face triangulated sphere, with **17 vertices and 30 edges**, is a counterexample to Sidorenko's
conjecture, conditionally on the premise above. The original graph has 35 vertices and 66 edges.

## Statement (Before / After)

| | Before (OpenAI, family 161) | After (this remark) |
|---|---|---|
| Complex | 13 points, 22 faces, conditions (1)–(4) (`complex.tex` L73–89) | 7 points, 10 faces, (1)–(3) plus "colour-preserving homs = Aut" |
| H | 35 vertices, 66 edges | **17 vertices, 30 edges** |
| Status | unconditional; Lean-formalised | conditional on premise (P); no Lean |

**The complex K₇.** Points 0–6. Faces f0–f9 = 012, 016, 023, 034, 045, 056, 125, 156, 234, 245.

**Edge list of H₇** (point–face incidences):
```
f0: p0 p1 p2    f1: p0 p1 p6    f2: p0 p2 p3    f3: p0 p3 p4    f4: p0 p4 p5
f5: p0 p5 p6    f6: p1 p2 p5    f7: p1 p5 p6    f8: p2 p3 p4    f9: p2 p4 p5
```
Point degrees are (6,4,5,3,4,5,3) for p0…p6; every face has degree 3.

**Bit certificate for (2)–(3).** Bits are read left to right.

| face | f0 | f1 | f2 | f3 | f4 | f5 | f6 | f7 | f8 | f9 |
|---|---|---|---|---|---|---|---|---|---|---|
| bits | 000 | 001 | 001 | 011 | 111 | 011 | 100 | 110 | 100 | 110 |

Exposure orders (each later face shares a pair with an earlier one and adds exactly one new point; the
new points are in brackets):

| bit = 0 | bit = 1 |
|---|---|
| b1: 0,1,2,3,5 [6,3,4,5] | b1: 4,9,6,7,8 [2,1,6,3] |
| b2: 0,1,2,6,8 [6,3,5,4] | b2: 3,4,5,7,9 [5,6,1,2] |
| b3: 0,6,7,9,8 [5,6,4,3] | b3: 1,5,4,3,2 [5,4,3,2] |

Separations: a_e ≥ 1/3 for every pair. a₁₆ = a₃₄ = 1 and a₂₃ = a₅₆ = 2/3; every other a_e is 1/3.

**Why (4) fails, and why it does not matter here.** The stable colouring has point classes {0}, {1,4}, {2,5}, {3,6}
and face classes {012,045}, {016,034}, {023,056}, {125,245}, {156,234}. Aut(H₇) = {id, φ}, where φ = (1 4)(2 5)(3 6)
on points. Of the 2⁶·2¹⁰ = 65 536 colour-preserving maps, exactly two are homomorphisms: id and φ.

## Proof sketch: generalised pinning (note §3)

1. **Ties (Lemma 3.1).** The scores of the original Lemma 3.4 still satisfy Q(f) ≤ Q(id). Equality holds exactly
   for maps preserving the *stable* colouring. Discreteness was only used to say that this set is {id}.
2. **Off homomorphisms (Lemma 3.2).** If f is not a homomorphism, then C(f) = 0, because a non-incident type pair
   has S = ∅.
3. **Automorphism invariance (Lemma 3.3).** For φ ∈ Aut(H), C_σ(φ) = C_{σ∘φ}(id) = (∏σ)(θ/2)^{|𝒞|} = C(id). The support
   forcing of the original Lemma 3.3 holds for every sign vector, and ∏σ is invariant under permutations.
4. Hence the maximal-score part equals C(id) + C(φ) = −2(θ/2)³⁰ < 0. The original steps apply from here on:
   ρ large, dilution with exponents ε⁷ and ε¹⁰, the crossed product, and the finite host (Prop. 3.5, Lemma 3.6).
   All of them use Lemma 3.1, hence premise (P).

Lines changed relative to the original: `complex.tex` Table 1 and Prop. 2.1(4); `kernel.tex` L383–450 (the
pinning lemma, generalised) and L454–510 (its application: one maximiser → two). The constants 13/22/66 become
7/10/30 throughout. In the tail section, 618 → 242, 22C* → 10C*, and 1575 → 407 (note §4).

## Controls: why not smaller

| complex | outcome |
|---|---|
| V=4 (tetrahedron), V=5 (bipyramid), V=6 #2 (octahedron), V=7 #4 | H is already Sidorenko: a reflection gluing reduces it to a half with a part of size ≤ 4 or a dominating vertex. No construction can work. |
| V=6 #1 | Pinning fails. 8 folding homomorphisms have C = +θ²²/2²⁴, which dominates the automorphisms' −θ²⁴/2²⁴. They pair up with vec(f) + vec(f′) = 2 vec(id), so by AM–GM no priors can isolate the identity. Whether this H is Sidorenko is open. |
| V=7 #5 | Generalised pinning is inconclusive: 270 maximising homomorphisms versus \|Aut\| = 6. |
| **V=7 #3 = K₇** | works (this remark); V=8 #6 and V=9 #10 also work |

## H₇ is not in a known Sidorenko class

| class | why H₇ is not covered |
|---|---|
| trees, even cycles, complete bipartite, one part ≤ 4 (Sidorenko) | parts of sizes 7 and 10; contains cycles; minimum degree 3 |
| vertex complete to the other part (Conlon–Fox–Sudakov) | maximum degrees 6 < 10 and 3 < 7 |
| weakly norming (Hatami: needs regular parts) | point degrees 3,3,4,4,5,5,6 |
| Coregliano degree-count criterion | 10 < C(7,3) = 35; the reverse orientation also fails |
| Conlon–Lee divisibility / \|B\| ≥ C(\|A\|, r) | fails (checked by the review session) |
| blow-ups (Conlon–Lee) | no twin vertices |
| Cartesian products | 17 is prime |
| subdivisions / theta graphs / even-path substitution | minimum degree 3 |
| symmetric gluing (reflection) | Aut = ℤ/2; the involution fixes only p0, and H₇ − p0 is connected |
| strongly tree-decomposable (CKLL), reflection trees (Li–Szegedy), tree-arrangeable (KLL) | exhaustive first-step leaf test: 56 candidate leaf sets, none passes |
| Szegedy's "thick" class | **not decided** |

In addition, 69 numerical step-kernel searches all converged to the constant kernel (see the recorded logs).

## Relation to rohanarun/sidorenko-template-refinements

The repository [rohanarun/sidorenko-template-refinements](https://github.com/rohanarun/sidorenko-template-refinements)
(created 2026-10-09; commit 4e1b9b6 at our check) **keeps condition (4)** and proves that, inside that template,
the 13-point / 35-vertex complex is minimal. For V ≤ 8 there is no discrete complex; for V = 9–12 they give exact
cut obstructions after a degree-≤6 filter. They also classify V = 13. Their identity Σ_{e∋v} a_e = 2 forces point
degree ≤ 6. K₇ meets this with equality at p0.

We do not claim that minimality result. Our remark **relaxes (4)**, which is exactly what takes us below their
template minimum. Our own exact LP dual certificates for V ≤ 10 (appendix; scripts `appendix_lp/`) independently
recompute their conclusion in that range:
- 155 discrete complexes, all infeasible, each with an exact `Fraction` certificate;
- after restricting to degree ≤ 6, the counts are 0, 1, 4 at V = 8, 9, 10, matching theirs.

An earlier primary run (`primary/`, plantri 5.5) covered V ≤ 12, i.e. 7881 complexes. Its V = 11, 12 part has not
been re-run independently.

## Reproduce

```bash
cd scripts && PY=python3 ./run_all.sh     # ~1-3 min, single thread; prints OK for 12 outputs
```
Requirements: Python 3.10 with networkx, numpy and scipy (HiGHS); a C compiler; network access for plantri.

| script | role |
|---|---|
| `certificate/check_v7.py` | **Main certificate.** Written from the TeX independently of the construction code. Checks (1)–(3) with the bits above, the stable colouring, Aut, the pinning maximisers over all homomorphisms (32 033 point maps extend), and the exact C = −(θ/2)³⁰ for each maximiser. |
| `construction/verify3.py`, `pinexact.py` | Primary checks. Enumerates the 65 536 colour-preserving maps (2 are homomorphisms) and runs exact max-score elimination over all part-preserving maps and over all homomorphisms. |
| `construction/survey.py`, `fold1.py` | Controls: Sidorenko-by-reflection cases; the V=6 folding obstruction. |
| `classes/invariants.py`, `stdcheck.py` | Invariants for the class table; the strong tree-decomposability first-step test. |
| `appendix_lp/r6_verify.py` | Appendix: V ≤ 10 enumeration, exact dual certificates, and the paper's 13-point complex as positive control. |
| `numerics/wsearch*.py` | Heuristic kernel search (PyTorch, hours). **Not re-run**; logs in `expected_output/*_RECORDED.log`. |

## Third-party software

plantri (G. Brinkmann, B. D. McKay et al.), Apache License 2.0, https://users.cecs.anu.edu.au/~bdm/plantri/. It is
**not redistributed**. `scripts/get_plantri.sh` downloads, sha256-checks and builds it into the ignored `scripts/_build/`:
- 5.8: `e78a944116fec9f2c9f5e484206276cc2b0043bae803e9815f4b2683614629b8`
- 5.5: `911cdf5bcca7294eb80f8f79fefc148183f7ba81da15b3aa4d6d2401a3bc7ded`

## Attribution

The original construction and proof are OpenAI's work: family 161,
[`preprints/A-counterexample-to-Sidorenkos-conjecture-September-23-2026`](https://github.com/openai/math/tree/main/preprints/A-counterexample-to-Sidorenkos-conjecture-September-23-2026),
[github.com/openai/math](https://github.com/openai/math) at commit `08fcd618213aa43b0fcaff4526a8a48d66d02e55`.
This remark only replaces the finite complex and generalises the pinning lemma. The new mathematical content is
small: Lemmas 3.1–3.3 of the note, essentially "C is invariant under automorphisms". We regard it as a
parameter-level variation of the original method that crosses a template boundary.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under [AUTHOR NAME]'s direction;
independent adversarial reviews were run by separate AI sessions and are not human peer review.

## What is NOT claimed

- Not an unconditional theorem: premise (P) is assumed, and there is no Lean.
- Not minimality of 17/30 among all Sidorenko counterexamples. It is only the smallest instance in the relaxed
  template for which *this* argument works: V=6 #1 and V=7 #5 are undecided.
- Not the template minimality with (4) kept, which is rohanarun's result; we only cross-check it for V ≤ 10.
- No explicit host graph G and no numerical violation: the host comes from the asymptotic rank-layer construction,
  as in the original.
- Novelty: to our knowledge (search of GitHub and arXiv up to 2026-10-09: GitHub repository search "sidorenko" /
  "sidorenko conjecture", code search "Sidorenko counterexample", the rohanarun and piiiico/openai-math-counterexamples
  repositories, arXiv API query "Sidorenko AND counterexample"), no Sidorenko counterexample with fewer than 35 vertices has been
  published.

## Open: relaxing (4)

This remark is one way to relax (4): colour-preserving homomorphisms must all be automorphisms. Weaker
replacements are not settled, for example handling folding homomorphisms (V=6 #1) by different priors or
activation laws. Nor is it settled whether the rank-layer premise (P) can be machine-checked for general complexes.
