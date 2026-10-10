# R4 verification record: Smith–Toda V(4) exists if and only if p ≥ 11 (conditional on the p = 1009 manuscript)

> **Precondition: conditional on the correctness of the original p = 1009 manuscript, which we did not re-verify.**

Claim: every finite (May E₁ page / inequality) statement in OpenAI's *A finite Smith–Toda complex V(4) at the prime 1009* holds at every prime p ≥ 11. At p ≥ 19 the inequalities hold as written in symbolic form. At p = 11, 13, 17 three inequalities fail, and the statements they were used for are checked by exact enumeration. Hence, given that proof, V(4) exists at every prime p ≥ 11. Layer 4 adds the unconditional nonexistence at p = 2, 3, 5, 7, so that V(4) exists ⇔ p ≥ 11.

Environment: Python 3.10.12 (`/usr/bin/python3`, numpy 2.2.6, sympy 1.14.0). For `independent/` we used a Python 3.12 venv with z3-solver 4.x (version string 5.1.0) and numpy. All scripts are single-threaded and all arithmetic is exact. Re-run on 2026-10-09 inside this repository; outputs are in `scripts/expected_output/`.

## Layer 1: main computation (rescan session; re-run by the coordinating agent and again here)

* `main/may308.py p` uses a direct recursive enumeration of E₁ monomials and computes d₁ ranks over 𝔽_p. It runs these checks:
  * (A) Lemma `low:algebraic-vanishing`: every lower obstruction bidegree is empty;
  * (B) the v₄-target lists;
  * (C) Table `may:n5-table` and d₁-injectivity;
  * (D) the Frobenius quotient-group table;
  * (E) Lemma `top:small-stems`;
  * (F) the page bound.

  Results:

| p | result | runtime |
|---|---|---|
| 7 | `FAILS ['A']`: product obstructions e = (0,2,0,2), k = 3 (`y A^16` at (33,119)) and e = (0,2,1,2), k = 2 (`yb A^9`) | 1 s |
| 11 | `ALL MAY-PAGE / FINITE CLAIMS HOLD` (v₄ targets non-empty at k = 1, 2, 10, 11) | 2 s |
| 13 | ALL HOLD | 1 s |
| 17 | ALL HOLD | 2 s |
| 19, 23, 29, 31, 37, 41, 43, 47 | ALL HOLD | — |

* `main/may308_dp.py 7 11 13 17` cross-checks with knapsack-DP counts. At p = 11, 13, 17 there are 0 non-empty lower bidegrees and the v₄-target sizes are {1:2, 2:2, p−1:3, p:2}, equal to the paper's. At p = 7 the lower bidegrees are non-empty.
* `main/ineq308.py` encodes the 33 p-size inequalities of the manuscript, with line numbers. Only 02 L109–111 and 02 L154–155 (from 17) and 03 L161–163 (from 19) fail at p = 11. The literal threshold is 19.
* `main/spot1009.py` is the positive control. At p = 1009 the v₄-target lists at k = 1, 2, 1008, 1009 reproduce the paper's monomials verbatim.

## Layer 2: independent re-implementation (separate session, `independent/`)

* `mayenum.py` enumerates exactly with z3 (with blocking clauses), plus a DP count.
* `sweep.py p` checks every bidegree of check (A), all v₄ targets up to the cutoff, and the §5 small stems, using an int64 DP mod 2⁶¹−1 to locate non-empty bidegrees and z3 to list them. It also checks (q+1, u₅+1) (DP count 5 = z3 count 5), A^N-divisibility of the k = p target, and w < u_{n+1} for the products. At p = 11, 13, 17 the only non-empty bidegrees are the four v₄ targets and the two small-stem bidegrees the paper expects. At p = 7 it additionally finds the two product obstructions above. That matches layer 1.
* `d1check.py p` checks:
  * d₁ rank 4 of 4 on the m ≤ 5 span at (q+1, u₅+1);
  * acb is not a cycle and acybz is a cycle;
  * the (6, W+1) d₁-cohomology is 0;
  * the Frobenius table matches;
  * S < qp+1 (p = 11: 205 < 221).

  All of these pass at p = 7, 11, 13, 17.
* `t1009.py` is the positive control at p = 1009. The k = 1, 2, p−1, p lists agree with the paper and with `spot1009.py`. The k = 3 and k = p+1 lists are empty.

## Layer 3: adversarial review (independent session, 2026-10-09). Verdict: **CONFIRMED (conditional)**

* The reviewer read sections 01–05, 02a and the appendix line by line. The appendix has no p-dependence. The reviewer found no implicit argument of the form "p large ⇒ some group vanishes"; sparseness and convergence use only explicit inequalities.
* The three failing inequalities agree with the original report. The reviewer added one point: 03 L161–163 still fails at p = 17 (it becomes an equality), so p = 11, 13 and 17 all need the enumeration, while p ≥ 19 works verbatim.
* Consistency signal: the p = 7 failure sits exactly at the product lemma, e = (0,2,0,2), k = 3. That matches Toda's result that V(3) is a ring spectrum for p ≥ 11 (Ravenel, green book 2nd ed., §5.6, p. 184: "According to Toda [1] it is a ring spectrum for p ≥ 11"; the theorem number "Toda 1971, Thm 4.4" given by the reviewer was not checked, because Toda's paper was not retrieved) and Nave 2010.
* Literature: there is no conflict. Nave at p = 11 rules out only V(6), and shows V(5) is not a ring spectrum. HMM rules out V(9). For n ≥ 4, V(n) has not been constructed at any prime (Ravenel 5.6.13; Ausoni–Rognes, arXiv 2204.05890).
* This session also grepped every number ≥ 100 in the manuscript (74 hits). Each is a citation or a numerical instantiation listed in README table (b).


## Layer 4: sharp classification (literature check, 2026-10-09). Verdict: **"⇒" holds unconditionally; "⇐" is conditional**

Sources retrieved and read. Local copies were made in a scratch directory and are not included in this release.
* Nave, *The Smith-Toda complex V((p+1)/2) does not exist*, Ann. of Math. 171 (2010) 491–509, journal PDF from annals.math.princeton.edu.
* Nave, *On the nonexistence of Smith-Toda complexes*, arXiv:math/9810178v1 (1998).
* Ravenel, *Complex cobordism and stable homotopy groups of spheres*, 2nd ed., author's PDF.
* Salch, arXiv:1607.01108.

**Definition.** Nave 2010, Definition 1.1 (p. 492): "A Smith-Toda complex is a finite spectrum V(k) satisfying BP_*V(k) ≅ BP_*/(p, v₁, …, v_k) as a BP_*-module, hence as a BP_*BP-comodule." The 1998 version (intro) says the same and attributes the problem to Smith [Smi70] and Toda [Tod71]. Ravenel 2nd ed. p. 285 also uses "V(3) [satisfying BP_*(V(3)) = BP_*/I₄]". **The definition does not presuppose V(k−1).**

**Reduction V(k) ⇒ V(k−1) (odd p).** Nave 2010, Theorem 1.2 (p. 493): "Let p > 2 and X be a finite spectrum. Then BP_*X ≅ BP_*/(p, v₁, …, v_k) as a BP_*-module if and only if H_*(X; F_p) ≅ E(τ₀, τ₁, …, τ_k) as an A_*-comodule." The proof is in Appendix A. Directly after it: "Such spectra need not be unique. However, it follows from Theorem 1.2 that if V(k) exists, it arises as a cofiber Σ^{|v_k|}V(k−1)′ → V(k−1) → V(k), where f induces multiplication by v_k in BP-homology. Note that f is not necessarily an honest self-map, i.e., the spectra V(k−1)′, and V(k−1) need not be equivalent." The 1998 version, §3, gives the same with "This point of view allows one to produce V(k−1) as a skeleton of V(k)". We re-checked the degree count: monomials without τ_k have degree ≤ Σ_{i<k}(2p^i−1) < 2p^k−1 = |τ_k|. So the sub-|τ_k| skeleton of a minimal cell structure has homology E(τ₀..τ_{k−1}). **Hence, at odd p, nonexistence of V(k−1) implies nonexistence of V(k). No alternative definition evades this.** At p = 2 Theorem 1.2 is not stated, so we do not use this reduction there.

| p | claim | source (exact) | status |
|---|---|---|---|
| 2 | no V(4) | our Lemma 5 (note.pdf §3). It extends Ravenel 2nd ed. §1.3, p. 19: "The nonexistence of V(1) for p = 2 follows easily from the structure of the Steenrod algebra; if it existed its cohomology would contradict the Adem relation Sq²Sq² = Sq¹Sq²Sq¹." The lemma uses only BP_iX for i ≤ 3, so it applies to every k ≥ 1. | proof written and self-checked; **not human-reviewed** |
| 3 | no V(2), hence no V(4) | Nave 2010 p. 493: "The first negative results were obtained by Toda, who showed that V(1) cannot be constructed when p = 2, and likewise, for V(2) when p = 3." Ravenel §1.3, p. 19: Toda's argument uses "the nonassociativity of the mod (3) Moore spectrum"; second proof via β₄ (5.5.1, pp. 180–181). | published, via secondary sources; **GAP: Toda 1971 (Topology 10, 53–65) not retrieved, theorem number unknown** |
| 5 | no V(3), hence no V(4) | Ravenel 2nd ed. **Theorem 7.6.1** (p. 284): "For p = 5, β₁¹⁷ ≠ 0 and there are Adams–Novikov differentials d₃₃(γ₃) = β₁¹⁸. Consequently the Smith–Toda complex V(3) does not exist, and V(2) is not a ring spectrum." 1st-ed. numbering: [Rav86, 7.5.1], cited by Nave 2010 p. 493: "Later, Ravenel [Rav86, 7.5.1] showed that V(3) does not exist when p = 5." | published |
| 7 | no V(4) | Nave 2010 (Annals) **Theorem 1.3** (p. 493): "If p ≥ 7, V((p+1)/2) does not exist. If V((p−1)/2) exists, it is not a ring spectrum." Abstract: "for p a prime greater than 5". **arXiv 1998 Theorem 1.1 is weaker**: "If p > 5, V((p+3)/2) does not exist. If V((p+1)/2) exists, it is not a ring spectrum." At p = 7 the 1998 version gives only V(5), and says V(4) is not a ring spectrum. The Annals version (received 2005-11-16) is required. | published |
| ≥ 11 | V(4) exists | Layers 1–3 + OpenAI 308 | **conditional** |

Note: Nave 2010 also states that Theorem 1.3 is "the best result possible using our techniques" (p. 494). It cites HMM (unpublished) for the nonexistence of V(p−2), which at p ≥ 11 only concerns V(9) and higher. No published result excludes V(4) at p ≥ 11. Salch (arXiv 1607.01108, final section) writes: "I have been told that M. Mahowald conjectured that V(3) is the last Smith-Toda complex to exist, that is, V(n) does not exist for n > 3 and for all primes". The conditional "⇐" contradicts this reported conjecture.

**Dedup (2026-10-09).** The following searches found no statement of "V(4) exists ⇔ p ≥ 11" or "V(4) exists for p ≥ 11":
* `gh search code` for "Smith-Toda V(4)", "Smith–Toda", '"V(4)" "p ≥ 11"', '"V(4)" "p >= 11"', 'Smith-Toda "every prime p"';
* `gh search repos` for "Smith-Toda";
* the arXiv API, for abstracts containing "Smith-Toda"/"Toda-Smith", newest first;
* web and MathOverflow searches.

The hits were only the two openai/math READMEs, an overlap index (mrconter1/vibemathed), and generic conjecture lists (The-Last-Math-Competition 00000004934: "for sufficiently large primes all exist", no bound). The OpenAI p = 1009 introduction says only that Nave's obstructions "leave room for V(4) at a sufficiently large prime" (00-introduction.tex L37–40).

## GAPs

1. **The correctness of the p = 1009 manuscript itself is not verified.** That covers §3 adjusted-filtration injectivity, §4 Frobenius left-hand side, §5 extension of powers, and the appendix pairing. The value of the remark rests entirely on it.
2. "All p-dependent steps are listed" rests on human-style reading by one adversarial reviewer plus grep. There was no second reviewer for this specific point.
3. Nonexistence at p = 3 rests on Toda 1971, cited through Nave 2010 and Ravenel; the original was not retrieved and its theorem numbers are unchecked. Nonexistence at p = 5 rests on Ravenel's Theorem 7.6.1 computation, which is published but was not re-derived here.
4. Policy note: by the project's internal criterion this is a parameter generalisation of the original argument. It is reported under the T-OPT (constant/threshold optimisation) heading.
5. The p = 2 lemma (note.pdf Lemma 5) is our own extension of Ravenel's one-line V(1) argument. It is elementary and was self-checked, but no human has reviewed it.
