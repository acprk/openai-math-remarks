# VERIFICATION — R17 (latent gadget rank threshold, family 102)

All runs: 2026-10-09, CPython 3.10.12, NumPy 2.2.6, mpmath 1.3.0, each pinned to one core.
`scripts/run_all.sh` reproduces both outputs byte-for-byte (`OK` ×2).

## 1. Who checked what

| Item | Primary | Independent re-implementation | Adversarial review | Status |
|---|---|---|---|---|
| Lower bound r_* ≥ 1/(4p) (Prop. 2.1 of the note) | internal review session (C-complexity, "Theorem N1"); `primary/gadget_bounds.py` Part 1: floating-point Hadamard computation on 30 random equivariant maps (n=7, l=3): r=1 identity to 1.1e-16, r=2 ratio ≤ 0.777 | `independent/indep_gadget.py` Part A: exact rationals, direct summation (no Hadamard), 12 random maps, r = 1, 2, 3: r=1 identity exact, max ratio 0.7806 | separate AI session re-derived every step (Fourier support, Parseval, p' step, union bound, quantifier order over l): **correct** | holds |
| Upper bound r_* ≤ ⌈1.5p^-3⌉+1 (Thm 3.1 of the note) | this build session (proof); internal review had an earlier O(p^-6) version | Part B: exact rational rank chains confirm the supermartingale property and the extinction bound; Part D: worst-case chains at the note's parameters (extinction 0.22–0.23 against the bound 0.65–0.67) | same session, two rounds: **CONFIRMED-WITH-CORRECTIONS** (see §3); all corrections applied | holds |
| "Before" thresholds (paper main text / effective recipe) | `primary/gadget_bounds.py` Part 3 (effective: 2^647.8 at p=1e-8) | Part C (main 2^449.8, effective 2^647.8) | reviewer recomputed both independently: 449.8 and 647.8 | agree |
| After threshold 2^80.3 at p = 1e-8 | — (formula ⌈1.5p^-3⌉+1) | Part C | reviewer: log2(1.5e24) = 80.31 | agree |
| Route limit 2 ln(8/7) p^-3 | — | Part D | the reviewer corrected our first value, ln(4/3), to 2 ln(8/7) by accounting for the leaf bound | corrected |

## 2. Commands and key outputs

```
python3 scripts/independent/indep_gadget.py
  Part A: E1 = p'/2 exactly in all 12 trials; max E_r/(r p'/2), r = 2,3: 0.7806
  Part B: e.g. p=3/10, r0=40: Pr[extinct] = 3.35e-1 <= e^(3 theta s)/(r0+1) = 9.06e-1; E[W_j] nonincreasing
  Part C (p = 1e-8): A (paper main) 2^449.8 | B (paper effective) 2^647.8 | C (harmonic, no enlargement) 2^223.3
                     | D (this note) 2^80.3 | L (lower) 2^24.6 | M (route limit) 2^77.8
  Part D: p = 0.5..0.2 at r0 = ceil(1.5 p^-3) - 1: exact extinction 0.22-0.23, supermartingale bound 0.65-0.67
python3 scripts/primary/gadget_bounds.py      (internal review script, earlier O(p^-6) claim)
  Part 1: |E1 - p'/2| = 1.11e-16; r=2 ratio 0.7771
  Part 3: p=1e-08 paper log2 r* ~ 647.8 ; "new" 166.7 (superseded O(p^-6) analysis) ; lower log2 r* >= 24.6
```

## 3. Adversarial review (separate AI session; not human peer review)

**Round 1** (the earlier claim r_* = Θ(p^-6) and the lower bound): CONFIRMED-WITH-CORRECTIONS.
- (medium) The parameters were circular: ε_0 and r_0 depended on θ, which depends on d. Fixed: ε_0 and r_0 now depend on p only.
- (low) "Before" must be described as the threshold certified by the preprint's union bound. Done.
- (low) r_* = b+1 already suffices in the paper's enlargement. This is superseded, since the enlargement step is no longer used.
- (low) The primary script printed log2(3/(16p)) as the lower bound. Fixed to log2(1/(4p)).
- The reviewer's own checks:
  - exhaustive checks over F_4 and F_8 of the single-block change probability 1−θ and of the kernel criterion (Lemma 3.4);
  - a Monte Carlo run of an adversarial rank chain that includes nongeneric wipe-outs; it stayed under the bound in all settings;
  - an independent recomputation of the thresholds.

**Duplicate search, then round 2.** Our search found `lgboim/unique-games-short-proof` (2026-10-09). It already removes the enlargement step for ℓ ∈ 3ℤ, which gives r_* = r_0 ≈ p^-8. With the supermartingale, the Θ(p^-6) claim was therefore dominated, and it was strengthened to r_* ≤ ⌈1.5p^-3⌉+1. Round 2 verdict on the revised claim: **CONFIRMED-WITH-CORRECTIONS**.
- (medium) Arithmetic: (2/3)e^{0.09} = 0.7295, not < 0.7273. The conclusion still held, with a margin of only 0.004. Fixed by requiring θ ≤ 0.01: (2/3)e^{0.03} + 1/60 < 0.704.
- (low) "Every d ≥ d_0" must be justified through the explicit union bound of Lemma 3.6, because the true nongeneric fraction need not be monotone in d. Fixed in the note.
- (low) Route limit: ln(4/3) should be 2 ln(8/7) ≈ 0.267, and the result should be labelled a limit of the method only. Fixed.
- The reviewer verified:
  - descent from an arbitrary (non-full, non-generic) domain S, since a lift may have any domain and Lemma 3.8 uses neither property;
  - treating the step-0 loss as a bad event;
  - fresh orientation for j ≥ 1;
  - the bad-event sum;
  - the supermartingale;
  - the padding (equivariance, change probability, dim(G_V|_B) ≥ r_0, g(a) = g_V(a_V));
  - the number 80.31.
- Duplicate check by the reviewer: lgboim states its contract only for ℓ ∈ 3ℤ and uses the harmonic analysis. The supermartingale, the p^-3 bound and the padding to every ℓ are not there.

## 4. Duplicate search (2026-10-09)

- **openai/math** (HEAD 08fcd61): Section 3 of the preprint gives no explicit r_*. The companion preprints of family 102 (Max-Cut, Vertex Cover, Min-UnCut, DFVS) do not use this gadget.
- **GitHub**:
  - repository searches: "unique games", "unique games theorem openai", "latent alphabet gadget";
  - `lgboim/unique-games-short-proof` read in full in the relevant section (acknowledged).
- **Web/arXiv**: a search for explicit bounds on the latent-gadget threshold found only news coverage and the VibeMathed summary (no explicit bound). The arXiv API was unreachable from this machine during the session, so the arXiv check went through web search only.

## 5. Open gaps

- **GAP-1:** the order of r_*(p) is open between p^-1 and p^-3.
- **GAP-2:** we rely on the preprint's Lemmas 3.2–3.8 and Prop. 3.9. They are not formalised upstream at this level; the upstream Lean covers the final reduction.
- **GAP-3:** no Lean for this remark.
- **GAP-4:** the primary script predates the final claim. The final threshold (2^80.3) comes from a closed-form formula. It was evaluated by the independent script and by the reviewer, and the proof is in the note.
