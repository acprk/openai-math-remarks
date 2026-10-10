# R4 verification record

> **Precondition: conditional on the correctness of the original p = 1009 manuscript, which we did not re-verify.**

Claim: every finite (May E₁ page / inequality) statement in OpenAI's *A finite Smith–Toda complex V(4) at the prime 1009* holds at every prime p ≥ 11. At p ≥ 19 the inequalities hold as written in symbolic form. At p = 11, 13, 17 three inequalities fail, and the statements they were used for are checked by exact enumeration. Hence, given that proof, V(4) exists at every prime p ≥ 11.

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
* Consistency signal: the p = 7 failure sits exactly at the product lemma, e = (0,2,0,2), k = 3. That matches Toda 1971, Theorem 4.4 (a unital multiplication on V(3) needs p ≥ 11), and Nave 2010.
* Literature: there is no conflict. Nave at p = 11 rules out only V(6), and shows V(5) is not a ring spectrum. HMM rules out V(9). For n ≥ 4, V(n) has not been constructed at any prime (Ravenel 5.6.13; Ausoni–Rognes, arXiv 2204.05890).
* This session also grepped every number ≥ 100 in the manuscript (74 hits). Each is a citation or a numerical instantiation listed in README table (b).

## GAPs

1. **The correctness of the p = 1009 manuscript itself is not verified.** That covers §3 adjusted-filtration injectivity, §4 Frobenius left-hand side, §5 extension of powers, and the appendix pairing. The value of the remark rests entirely on it.
2. "All p-dependent steps are listed" rests on human-style reading by one adversarial reviewer plus grep. There was no second reviewer for this specific point.
3. p = 3 and 5 are not addressed.
4. Policy note: by the project's internal criterion this is a parameter generalisation of the original argument. It is reported under the T-OPT (constant/threshold optimisation) heading.
