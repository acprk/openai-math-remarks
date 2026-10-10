# VERIFICATION — R6 (Sidorenko, conditional 17/30)

All runs: 2026-10-09, CPython 3.10.12, networkx 3.4.2, scipy 1.15.3, one core.
`scripts/run_all.sh` reproduced all 14 expected outputs (`OK` ×14, about 95 s wall time).

## 1. Who checked what

| Item | Primary (construction session) | Independent re-implementation | Adversarial review (separate session) | Status |
|---|---|---|---|---|
| (1)–(3) for K₇, bit certificate, exposure orders | `construction/verify3.py` | `certificate/check_v7.py`, written from the TeX; the exposure orders were recomputed once more in this build | re-derived the bit table | exact, holds |
| Stable colouring; Aut(H₇) = ℤ/2 | `cplx.py` | `check_v7.py`; `classes/invariants.py` (networkx isomorphisms) | yes | exact, holds |
| Colour-preserving maps: 65 536, of which 2 are homomorphisms | `verify3.py` (brute force) | `pinexact.py` (max-plus elimination: 65 536 maximisers among all maps, 2 among homomorphisms); `check_v7.py` (32 033 hom point-maps; 2 maximisers) | yes | exact, holds |
| C(id) = C(φ) = −(θ/2)³⁰ | `ccoef.py` | `check_v7.py` (exact `Fraction` polynomial in θ) | proved C(φ) = C(id) in general (∏σ is invariant) | exact, holds |
| Generalised pinning (ties); C = 0 off homomorphisms; automorphism invariance | construction notes | written out in note §3 by this build | accepted | proof written; short |
| Location outside known Sidorenko classes | `survey.py` (reflection test) | `invariants.py`, `stdcheck.py` (56 leaf sets, none passes) | full table incl. Conlon–Lee and Coregliano | not covered by any class checked; **Szegedy thick class undecided** |
| Controls: V=4,5, V=6 #2, V=7 #4 Sidorenko; V=6 #1 folding; V=7 #5 inconclusive | `survey.py`, `fold1.py`, `pinexact.py` | — | V=6 and V=5 controls re-run | consistent |
| Kernel search (69 starts) → constant kernel | — | `numerics/wsearch*.py` (logs recorded, not re-run in this build) | yes | heuristic only |
| Generalised pinning, quantitative margin | `pinexact.py` | `certificate/pin_margin.py` (LP + integer check over all 7⁷ point maps): δ* = 3/22 | — | exact check of the rounded scores, holds |
| Premise (P): Sections 4–7 transfer to K₇ | — | `premise_audit/rank_layer_audit.py` recomputes every complex-dependent constant; calibrated on the paper's complex | line-by-line audit (PREMISE-AUDIT.md): 21 places, none uses structure beyond (1)–(3) and constants (618→242, 1575→407, 1641→437, 44C*→20C*) | **assumed; audited, not proved or machine-checked** |
| Appendix: with (4) kept, V ≤ 10 infeasible | `primary/` (plantri 5.5, V ≤ 12, 7881 complexes) | `appendix_lp/r6_verify.py` (plantri 5.8, dart face-tracing, Hamiltonian-cycle splits cross-checked by brute force for V ≤ 10, `Fraction` certificates) | — | V ≤ 10: 155/155 certified; V = 11, 12 primary only |
| Positive control: the paper's 13-point complex | `primary/sidsearch` | `r6_verify.py`: (1) holds, (4) discrete, 1756 splits, paper bits admissible, t* = 1/3 feasible | — | passes |

## 2. Key outputs

```
check_v7.py:  (1) |E|=15, Euler 2, connected; (2) six exposure orders; (3) 3a_e: 1 (x11), 2 (x2), 3 (x2)
              stable classes points [[3,6],[1,4],[2,5],[0]]; |Aut| = 2
              point maps admitting a hom: 32033; maximising homs: 2 (id, phi)
              C(id) = C(phi) = -(theta/2)^30; sum over maximisers = -2 (theta/2)^30
pinexact.py 3 (homs only): #maximisers = 2 = |Aut|      pinexact.py 3 (all maps): 65536 maximisers
pinexact.py 1 (V=6 #1): 16 maximisers vs |Aut| = 4;  4 (V=7 #4): 48520 vs 20;  5 (V=7 #5): 270 vs 6
fold1.py:     8 fold maps, 8 midpoint pairs, sigma-dependence violations 0
stdcheck.py:  candidate leaf sets 56, passing the copy test 0
r6_verify.py: discrete counts V=8,9,10 = 2,16,137; exact certificates 155/155; worst bounds 1/5, 2/7, 2/7
```

## 3. Exact versus numeric

Everything in §2 except the kernel search uses integer or `Fraction` arithmetic. The LP is solved in floating
point (HiGHS) only to *find* dual vectors; every certificate is then checked in exact arithmetic against every
Hamiltonian cycle.

## 4. Duplicate search (2026-10-09)

- rohanarun/sidorenko-template-refinements (commit 4e1b9b6): keeps (4); proves 35 is minimal within that
  template. This is complementary to our result; we cross-check V ≤ 10.
- piiiico/openai-math-counterexamples: row 161 only re-checks the original 35/66 complex.
- GitHub repository and code searches for "sidorenko" and "Sidorenko counterexample", and an arXiv API query for
  "Sidorenko AND counterexample", found no smaller counterexample.

## 5. Open gaps

- **GAP-1 (main):** premise (P). Lemma 3.1 and Sections 4–7 have not been re-derived or formalised for K₇. The
  official Lean is hard-coded to the 13/22/33 complex. A line-by-line AI audit (PREMISE-AUDIT.md) found no
  dependence on the 13-point complex beyond conditions (1)–(3) and constants; this narrows the gap but does not
  close it.
- **GAP-2:** whether H₇ lies in Szegedy's thick-graph class is not decided. A proof that H₇ is Sidorenko (by any route) would refute (P) or the argument here.
- **GAP-3:** V = 11, 12 of the appendix (970 + 6756 complexes) rely on the primary computation only.
- **GAP-4:** V=6 #1 and V=7 #5 are undecided. So 17/30 is the smallest size for which *this argument* works,
  not a proven minimum.
- **GAP-5:** the numerical kernel search was not re-run in this build.
