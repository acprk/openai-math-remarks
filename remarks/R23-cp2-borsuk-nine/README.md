# R23 — Family 156: the complex projector set of CP² has Borsuk number 9 (no complex-projective Borsuk counterexample below dimension 9)

## Premise and confidence (read first)

**This is a negative result, and the mathematics in it is elementary.** OpenAI's family-156 preprint
proves that the projector image of RP³ in ℝ⁹ cannot be covered by ten sets of smaller diameter. A natural
follow-up is to complexify: does the projector image of **CP²**, which lives in ℝ⁸, fail to be
coverable by nine sets? It does not. Nine explicit caps suffice. The lower bound 9 is the preprint's
own odd-map argument transferred verbatim to Hermitian matrices (a parallel transfer, not new). The upper
bound is a covering-radius computation. Its only non-routine part is that we certify it **globally on all of CP²**
with two independent rigorous computations. We record it so that nobody spends time on this route. CP² is the only complex projective projector set in dimension < 9 that is not a round sphere: CPⁿ for n ≥ 3 lives in ℝ^{n²+2n} ⊇ ℝ¹⁵, and CP¹, HP¹, OP¹ are round spheres. So the route "complexify to go below dimension 9" is closed. **Caveat on the analogy:** CP² is the complex analogue of RP² (k = 3), not of the RP³ witness (k = 4). For k = 3 the preprint does not exclude the equality case either, so b = 9 here is the expected outcome. The complex analogue of the witness itself, CP³ ⊂ ℝ¹⁵, is untouched and irrelevant for the dimension record.
It is not a research contribution in its own right.

**One line.** Before: the 156 lower-bound argument applied to CP² gives b ≥ 9 = d + 1 in ℝ⁸. Whether nine
sets actually suffice was not settled (an internal exploration had conjectured they do not). After: **b(X_CP²) = 9.** Nine caps {|⟨u,p_i⟩|² ≥ 21/40} cover the set, each of diameter
≤ √(399/200) < √2. So CP² gives no Borsuk counterexample in dimension 8.

- **Correctness:** high. The cover inequality is certified on the whole parameter domain by (a)
  outward-rounded interval arithmetic and (b) an independent exact-integer branch and bound. A negative control
  at T = 0.53 fails with an explicit witness, as it must.
- **Novelty:** low. To our knowledge (search of GitHub and arXiv up to 2026-10-09) the Borsuk number of
  the CP² projector set has not been written down. It is, however, a direct question about a classical
  object, and the answer is unsurprising.
- **History of this line.** An earlier internal exploration (code name K-A) conjectured the
  *opposite*, namely that CP² might be an 8-dimensional counterexample, and spent effort on
  Clifford-equivariant SAT searches. That conjecture was judged **wrong** in our third-round audit,
  once a 9-cap cover was found. The cover was first checked only by a floating-point branch and bound with a
  Lipschitz bound. This remark replaces that check with rigorous global certificates.

## Statement (Before / After)

X = { uu* : u ∈ ℂ³, |u| = 1 } ⊂ { A ∈ Herm₃(ℂ) : tr A = 1 } ≅ ℝ⁸, with the Frobenius metric.
Then ‖uu* − vv*‖² = 2 − 2|⟨u,v⟩|², so diam X = √2, attained exactly at orthogonal pairs.

| | Before | After (this remark) |
|---|---|---|
| Lower bound | 156's argument (real case) gives m ≥ k(k+1)/2 = d+1. Its complex form gives b(X) ≥ 9 (not written in the preprint) | written out: b(X) ≥ 9 (note, Prop. 2) |
| Upper bound | none rigorous; float-only B&B in an internal audit | **b(X) ≤ 9**: explicit 9 caps, threshold T = 21/40 (also certified at 0.526), each cap of diameter ≤ √(399/200) ≈ 1.41244 |
| Consequence | "CP² as an 8-dim Borsuk counterexample" open as a route | closed: b(X) = 9 = d+1 |

## Attribution

The original construction and proof are OpenAI's work: family 156 of
[github.com/openai/math](https://github.com/openai/math), commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`,
`preprints/A-nine-dimensional-counterexample-to-Borsuks-covering-assertion-September-23-2026/`.
This remark only (i) transfers the preprint's strict-gap reduction (`build/sections/geometry.tex`,
`lem:cover-to-map`, l. 92–123), its odd extension (`build/sections/extension.tex`,
`lem:positive-extension` / `lem:signed-extension`, l. 24–185) and the first sentence of its
`thm:label-bound` (`build/sections/topology.tex`, l. 28–39) from real symmetric to complex Hermitian
matrices. The changes are: U(3)-invariant measure on S(ℂ³), ∫uu* dμ = I/3, dim_ℝ Herm₃ = 9, and continuity
suffices. It (ii) adds a certified nine-cap cover. It takes the correctness of the preprint's real-case
arguments as given only as a template; the complex transfer is re-proved in the note.

Author: [AUTHOR NAME]. Developed with AI assistance (Anthropic Claude, Opus 5.5) under
[AUTHOR NAME]'s direction; independent adversarial reviews were run by separate AI sessions and are
not human peer review.

## Proof sketch

1. **Lower bound** (note §2). Cover by m sets of diameter < √2 ⇒ continuous f : CP² → Δ^{m−1} with
   disjoint supports on orthogonal lines ⇒ E(Q) = 3∫H(Q^{1/2}u)dμ, F(A) = E(A₊) − E(A₋) is odd with
   ‖F(A)‖₁ = tr|A| ⇒ odd map S⁸ → S^{m−1} ⇒ m ≥ 9 (Borsuk–Ulam).
2. **Cap lemma** (note Lemma 4). Unit x, y with |⟨x,p⟩|², |⟨y,p⟩|² ≥ T > 1/2 satisfy |⟨x,y⟩| ≥ 2T − 1. (Write
   x = ap + x′ and use |a||b| − |x′||y′|.) With T = 21/40 every cap has diameter² ≤ 2 − 2/400.
3. **Cover** (note Prop. 5). For all x ≠ 0, max_i (|⟨x,p_i⟩|² − T|p_i|²|x|²) > 0. Charts x_j = 1, |x_k| ≤ 1
   ⊂ [−1,1]⁴ cover CP². Each h_i is a real quadratic, bounded on boxes rigorously.

## Reproduce

```bash
cd scripts && PY=python3 ./run_all.sh     # python3 + numpy (+ scipy for numerics/); single thread; ~1.5 min
```
It checks the SHA-256 of `P9.txt`, then prints `OK` for eight outputs:
- `cert_interval.py P9.txt 21 40` and `263 500`: **primary certificate.** numpy interval arithmetic with
  `nextafter` outward rounding after every IEEE operation, 16-way branch and bound on three charts.
  For 21/40: 1 897 683 accepted leaves, minimum half-width 2⁻¹².
- `control_53_100`: the same code at T = 0.53 must fail. It stops with the witness
  x = (1, 0.67578125+0.30859375i, −0.26171875+0.08984375i), whose normalized max overlap is 0.52950.
- `control_drop9`: the same code with the ninth point deleted (`controls_P8_drop9.txt`) must fail, since by the
  lower bound no eight caps can work. It fails at the first level.
- `independent/cert_exact.py P9.txt 21 40` and `263 500`: **independent certificate**, standard library only.
  It uses exact integers, an exact second-order Taylor bound, and bisection along the widest coordinate. For 21/40:
  16 876 leaves, depth 36.
- `numerics/minimum_search.py`, `numerics/cap_lemma_check.py`: illustrations only. The numerical minimum of the
  normalized max overlap is 0.526055, and a random test of the cap lemma is included.
- `provenance/`: the heuristic search that found the points, and the conversion to `P9.txt`. Neither is run and neither is part of the proof.

## Scope and limitations (what is NOT claimed)

- **No new method.** The lower bound is the preprint's argument with ℝ replaced by ℂ. The upper bound is an elementary
  cap cover with a computer check.
- **Only this set.** Nothing is claimed about Borsuk's problem in dimensions ≤ 8 in general. We show only that the
  CP² projector set is not a counterexample.
- **A symmetric configuration does not suffice.** For the Hesse SIC, e₁ has overlap 0 or 1/2 with every SIC line, so its max
  overlap is exactly 1/2 (borderline, diameter √2). The certified points are a perturbed configuration found by computer.
- **The optimal threshold is not determined.** The best 9-point covering threshold is ≥ 0.526. Whether a symmetric
  (SIC/MUB/Clifford-equivariant) nine-cover exists is not settled here. An earlier exploratory SAT computation
  suggested that equivariant covers have margin ≤ 0.005; it is not used here.
- **Floating-point model.** The primary certificate assumes IEEE-754 binary64 round-to-nearest (numpy's default).
  The independent certificate makes no floating-point assumption.
- **No Lean.** Nothing is formalized (per project policy, negative/elementary results are not formalized).

## Prior work searched

GitHub and arXiv up to 2026-10-09: "Borsuk" with "complex projective", "Hermitian projector", "Fubini–Study";
the openai/math repository and its family-156 discussion; Kalai's survey (arXiv:1505.04952, §2.3, the real
elliptic-space version via x ↦ x⊗x) and his blog post of 2026-10-07 on the OpenAI release. We found no
statement of b(X_CP²). The 63-dimensional record (arXiv:2608.12561) and the Bondarenko/Jenrich
constructions over PG(2,16) with a Hermitian form are finite configurations and are unrelated to this set.
