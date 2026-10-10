# VERIFICATION — E1 (documentation, packaging and reproducibility errata)

Upstream state: github.com/openai/math `main` = `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`
(re-fetched 2026-10-09; the October 7 update, merge of PR #1). All runs: 2026-10-09, CPython 3.10.12,
NumPy 2.2.6, mpmath 1.3.0, sympy 1.14.0, one thread.
`scripts/run_all.sh` (full mode) reproduced all 10 expected outputs: `OK` ×10, `ALL OK` (764 s CPU, 924 s wall on a CPU-throttled host). `SKIP_LONG=1` mode: 8 × `OK` in under a minute.

## 1. Who checked what

Three layers, as for the other remarks: (a) primary check, (b) independent re-implementation written
separately, (c) adversarial review by two separate AI sessions told to refute each item against the
current upstream `main` (not against our earlier notes).

| Item | (a) Primary | (b) Independent | (c) Adversarial review | Result |
|---|---|---|---|---|
| P1 266 `verification/` absent | `check_upstream.sh` (`git ls-tree`, basename counts) | `indep_check_paths.py` (`git cat-file -e` at the commit and all ancestors) | ls-tree + `git log --all` + deleted-file log; also confirmed the binary64 / a-priori error wording (`rounding.tex` l.22–45; no interval library or directed rounding in the TeX) | **confirmed**; correction adopted: `full-run-01/…` are archive-internal paths |
| P2 095 `verify_certificate.py` absent | same | same | confirmed; same-named files exist only under *Snaky-in-21-Maker-moves*; minors are F₇ blocks of size 28–66 | **confirmed**; mitigation (paper says the lemma does not rely on the program) stated in README |
| P3 318 `support/` absent | same + `318_ce_ranks.py` | `indep_318.py` (cochain formula, different basis, sympy rank) | confirmed; line numbers SUPPORT-COVERAGE.md l.7/11/13/15 supplied by reviewer | **confirmed**; ranks (0,3,0,0), (0,8,20,35) agree |
| P4 222 certificate file absent | same + extraction + SHA-256 + full run | hash recomputed by reviewer separately (also: without the final newline the hash differs, so the published digest is of the newline-terminated file) | confirmed | **confirmed**; certificate reproduces |
| P5 238 `support/network_bounds.py` absent | same + `238_network_bounds_recheck.py` | `indep_238.py` (from the definitions of h_j and H; generic knapsack; own tail bound; interval logs) | confirmed; noted that the five inequalities are hand-checkable | **confirmed**; all five inequalities hold |
| P6 227 `numerics-*.py` absent | same | same | confirmed, lowest severity (paper says not needed) | **confirmed** |
| P7 266-B listing path | same + run of shipped `verification/run.py` | `indep_check_paths.py` | confirmed; the checksum in paper.tex l.259–262 matches the shipped file; **correction adopted**: the fix must be `../verification/verify.py` relative to `build/` | **confirmed** |
| D1 107 docs paper list | `check_upstream.sh` | reviewer's `git grep` over the 2.258 manuscript (no 9/4 proof there) | confirmed; "README has no Lean link" **dropped** (no manuscript README links Lean) | **confirmed (narrowed)** |
| D2 026 docs scope | same | reviewer read `largeGapIndices` (C·log p_n < p_{n+1} − p_n) | confirmed: understatement | **confirmed** |
| D3 142 CONTENTS | same | reviewer checked 029's theorem numbering (Thm 1.2 = `thm:hecke-zero-free`) | "mostly refuted": CONTENTS l.3436 (abstract row) already states the dependency; only the summary line l.3427 lacks it | **kept as an optional transparency suggestion, with that fact stated** |
| W1 103 wording | same | reviewer checked q_fr (l.68–70, exact frame probability) and that the additive reading would be vacuous | "partly": a one-word slip, no gap | **kept as wording** |
| W2 140 mixture constant | `140_mixture_constant.py` (binary64 grid, D' ≤ 40, 9 880 triples, 0 violations) | reviewer's own mpmath check (40 digits, D' ≤ 15, tails to 10⁻²⁹): 0 violations; max of sup/2^{(b−1)/2} = 1/√2 | each proof step re-derived (beta-product identity, monotonicity of B, cases j = 1 and b = 2); downstream use l.253–256 and c_* unaffected; Lean `MixtureEstimate` states the same loose 2^{finrank E} | **confirmed (optional, no change needed)** |
| W3 099 docs | `check_upstream.sh` | reviewer read all four comparators (`EditDistance`, `FiniteCircle`, `BinaryEditLower`, `TreeEdit`) and the Tree paper abstract | **refuted as a contradiction**: l.19 refers to the binary comparators and reuses the Tree paper's abstract | **downgraded to an optional wording note** |

Not included after review: 238 ε_n = O(n^{−5/2}) (a constant slack, left to a constants remark);
families 205 and 269 (handled in other remarks).

## 2. Key outputs

```
check_upstream.sh       openai/math commit: fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
                        every cited missing file: "files named … in tree: 0"
                        (verify_certificate.py: 3, all under Snaky-in-21-Maker-moves)
indep_check_paths.py    all cited paths present@commit=False present@any-ancestor=False,
                        except verification/verify.py (266-B, present)
222_extract.txt         sha256 of listing = sha256 in paper = 7bd5dbd918d9…a85dc7f, match: True
222_certificate_run.txt 24 checkpoints, first "0 410000", last "283000 304000"; exit 0
266B_verification_run   rank table: 8871 909 897 / 46041 3429 2767 / -613302797399911/6480 … (= Appendix A)
238 (both scripts)      H(1/8)<.791402, H(13/50)<1.098069, H(1)<1.958069,
                        bracket(7/25)<.904054, bracket(1)<.515374: all OK
318 (both scripts)      sl2 (0,3,0,0), sl3 (0,8,20,35); H^0..H^3 = 1,0,0,1
140                     9880 triples, 0 violations of log2 C <= (b-1)/2
```

## 3. Exact versus numeric

- P3, P5: exact (`Fraction`, sympy). The independent 238 script uses `mpmath.iv` interval arithmetic for
  the logarithms. Both give upper bounds below the claimed constants.
- P4: the certificate is the paper's own binary64 code; we only re-ran it.
- W2: binary64 grid evaluation, illustration only; the claim rests on the paper's own proof, re-derived
  by the reviewer.
- Everything else is a file-existence or text check.

## 4. Open gaps

- **GAP-1.** For P1 (N(6) = 3), P2 (095 rank certificate) and P6 (227) we confirmed only that the cited
  artifacts are missing. We did not re-implement those computations.
- **GAP-2.** Line numbers are for `fd4aeeb2e`. If upstream changes, `check_upstream.sh` will show the
  differences in its output (that is the intended failure mode).
- **GAP-3.** The duplicate search covered GitHub only (openai/math PRs, mahrud/openai-math-comments
  threads, repository search). GitHub code search hit a rate limit for some file names.
