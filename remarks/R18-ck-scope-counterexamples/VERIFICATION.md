# VERIFICATION — R18 (scope of the family-119 bounds)

All runs: 2026-10-09, CPython 3.10.12, mpmath 1.3.0, sympy 1.14.0, NumPy 2.2.6, single core.
`scripts/run_all.sh` reproduces all three outputs byte-for-byte (`OK` ×3).

## 1. Who checked what

| Item | Primary | Independent re-implementation | Adversarial review | Status |
|---|---|---|---|---|
| Extended Cor. 1.2 / Thm 1.1 fail for small λ when w ≠ ½ (resampling noise) | internal review; `primary/qary_refined_boolean.py` scan, max gap 4.6e-3 | `independent/indep_scope.py`: interval certificates at q=4, q=3, p=0.3; exact sympy series | separate AI session: sympy series incl. λ³ term −m²λ³/(6w(1−w)); direct 4×4 joint-law sum for q=4 matches the closed form to 40 digits; windows reproduced | holds |
| Extended Hellinger bound fails (0.005464 > 0.005013) | internal review (CKB-ROUND1, Prop. 4); `primary/hellinger_counter.py` | interval certificate [0.00546448802612, …] vs [0.00501256289338, …] | reproduced: 0.00546449 > 0.00501256; window root 0.3227258 | holds |
| Plain CK curve not violated by these examples | — | interval check in each witness | grid w = 0.01..0.99 × 199 values of λ: never violated (equality only at w = ½) | observation only |
| Uniform binary cube cannot produce a witness | — | AND control under bit flips satisfies both bounds | same control computed by the reviewer | consistent with 119 |

## 2. Key outputs

```
[q=4 uniform, |A|=1] w=0.25 lam=0.1: I in [0.00481700696933, ...], psi(lam r_m) <= 0.00444575410525 -> (A) violated
   Hellinger loss in [0.00546448802612, ...] vs 1-sqrt(1-lam^2) <= 0.00501256289338 -> (B) violated
Taylor: I = 1/2 lam^2 ; psi(lam r) = r^2/2 lam^2 ; Hellinger loss = 1/(2 sqrt(1-m^2)) lam^2 ; 1-sqrt(1-lam^2) = 1/2 lam^2
Windows: w=0.25: (A) 0<lam<0.43647647, (B) 0<lam<0.32272579
Control AND, rho=0.1: I = 0.0033451 <= 0.00444575 ; Hellinger 0.00386403 <= 0.00501256
```

## 3. Adversarial review (separate AI session; not human peer review)

Verdict: **CONFIRMED-WITH-CORRECTIONS**.
- (medium, framing; applied) The noise model decides the outcome, not the product space alone.
  - The counterexample stands under resampling: the Bonami–Beckner p-biased T_ρ, and the q-ary symmetric channel parametrised by its correlation λ = 1 − εq/(q−1).
  - It gets stronger under the parametrisation ρ = 1 − 2ε.
  - It disappears under ρ = Pr[Y = X], which is degenerate at ρ = 0.
  - Under independent bit flips on p-biased bits: an exhaustive check of all Boolean functions for n = 2, 3, at 19 values of p and 25 values of ρ, found no violation (max LHS − RHS = 4.4e-16).
  - README and note now say "under the resampling noise operator".
- (low; applied) The mechanism should be stated plainly: at n = 1 all of f's variance sits in the λ-eigenspace, whereas the level-1 inequality prevents this on the uniform cube. The AND control was added.
- (scope; no overreach found) Both preprints fix the uniform measure: Sharp-binary main.tex l. 53 and l. 113; Hellinger Thm 1.1 "with respect to the uniform probability measure", and "Throughout the paper the input measure is uniform".
- (literature) No prior statement was found in the papers listed in the README. The reviewer notes the computation is folklore-level.

## 4. Duplicate search (2026-10-09)

- **openai/math**: family 119's two preprints (above) and its Lean docs (`lean/docs/119.md`, "uniform discrete cube").
- **GitHub**: repository searches for "courtade kumar", "courtade-kumar" and "hellinger contraction boolean". We cloned and grepped dpwoodru/general-courtade-kumar-lean, hessammahdavifar/most-informative-bit and deepinfra/courtade-kumar. All three work on uniform inputs. The deepinfra notes say the work is not "an extension to biased inputs/asymmetric noise".
- **arXiv** (by the reviewer, via web): 2609.26444, 2609.24184, 2609.24931, 1506.03167, 2601.09679, 2410.10147, 1807.11289.

## 5. Open gaps

- **GAP-1:** whether the plain Courtade–Kumar inequality I ≤ ψ(λ) holds on general finite product spaces with resampling noise is open. The internal review has numerical evidence with no violation, which is not part of this remark.
- **GAP-2:** for bit-flip noise on p-biased cubes, only n ≤ 3 has been checked, and nothing is claimed.
- **GAP-3:** no Lean.
