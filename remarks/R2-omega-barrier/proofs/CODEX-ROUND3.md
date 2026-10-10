# Codex round 3: certified interval lasers and a zero/rest barrier

Read against `BRIEF-codex-R3.md`, both earlier `CODEX-ROUND*.md`,
`../allfields/BARRIER.md`, and Sections 2–5 of `../allfields/paper94.tex`.
The field and character framework are those of the preceding rounds: complex
tensors and nonnegative real-valued characters.

**Result.** No tested laser inequality refutes the balanced assignment. I prove
an **all-size outer-partition barrier for a legal zero/rest laser partition
with rank-feasible child values**, even if its truncated-convolution block
is given its character upper bound. Consecutive-interval calculations are certified for every `n <= 30` in
the recursive family specified below, with an additional exhaustive search over
the two-interval cut positions at the root for `n <= 12`. The ordinary
support/quantum functionals do not close the unrestricted multi-matrix route:
their matrix normalization is different, and their values on convolution can
be computed exactly without producing the required bound.

**I have not proved a barrier for arbitrary laser extractions, arbitrary
partitions of higher powers, or arbitrary changes of basis.** The remaining
route in the brief must not be marked closed. Only `CODEX-ROUND3.md` was written
in the project; the verification program is included here.

## 1. Three quantities that must not be identified

Put

\[
 \phi(a)=\prod_{j=1}^{a-1}\frac{3j+1}{3j},\qquad
 F(a,b)=\phi(\min(a,b))\left(\max(a,b)+\frac{\min(a,b)-1}{2}\right),
 \qquad \tau=\frac34.
\]

For an actual character whose average dot exponent is \(\tau\), define

\[
 G_\lambda(T)=\left(\prod_{\pi\in S_3}\lambda(T^\pi)\right)^{1/6}.
\]

Thus \(G_\lambda(C(a,b))=P(a,b)^\tau\), and
\(G_\lambda(M(e,h,l))=(ehl)^\tau\).

I use \(V^{\rm sym}_\tau(T)\) for the operational matrix-extraction value with
all leg permutations balanced. Concretely, take powers of
\(\bigotimes_{\pi\in S_3}T^\pi\), extract full sums of matrix products, sum
their \((ehl)^\tau\) values, and normalize by the number of original factors.
This is the usual symmetrized value convention with parameter
\(\rho=3\tau=9/4\). It permits permuted source factors; it must not be
silently described as a calculation using only the original orientation of
\(C(n,n)^{\otimes k}\).

There are two different kinds of lower certificate below:

* **Operational:** dot products, scalars, and recursively certified matrix
  extractions. These lower-bound \(V^{\rm sym}_\tau\) and hence
  \(G_\lambda\).
* **Profile-assisted:** also use the paper's proved lower bounds on
  \(G_\lambda(C(a,b))\). These lower-bound actual character values. They
  are not, merely for that reason, constructive lower bounds on
  \(V^{\rm sym}_\tau\).

Neither kind is an evaluation of the unknown unrestricted value. In particular,
a number in the table below is not a claim that \(V^{\rm sym}_\tau\) equals it.

### 1.1 Why the profile seed is valid, including unequal input lengths

The paper's diagonal recurrence proves the stronger bound
\(P(a,a)\geq(3a-1)\phi(a)/2\). For fixed \(a\), write
\(g(h)=P(a,h)\) and \(\Delta_h=g(h+1)-g(h)\). Concavity and tripling imply

\[
 2g(h)\leq g(3h+a-1)-g(h)\leq(2h+a-1)\Delta_h.
\]

Consequently

\[
 g(h+1)\geq g(h)\frac{2h+a+1}{2h+a-1}.
\]

For \(b\geq a\), telescoping from \(h=a\) proves

\[
 P(a,b)\geq P(a,a)\frac{2b+a-1}{3a-1}\geq F(a,b).
 \tag{1}
\]

Symmetry handles the other ordering. Thus **\(F(a,b)^\tau\) is a proved
lower bound on the symmetrized character value**, conditional on the paper's
character/profile lemmas. It is not a proved lower bound on each individual
\(\lambda_\pi(C(a,b))\), and is not an operational matrix-extraction seed.
The program keeps these distinctions.

## 2. The entropy optimization and its certificates

For a tight block support \(S\), a distribution \(p\) on \(S\), and certified
block values \(v_s\), the expression used is

\[
 \mathcal L_S(p;v)=
 \exp\left(\frac{H(p_X)+H(p_Y)+H(p_Z)}3
       +\sum_s p_s\log v_s-\Gamma_S(p)\right),
 \tag{2}
\]

where
\(\Gamma_S(p)=\max_{q:\,q_X=p_X,q_Y=p_Y,q_Z=p_Z}H(q)-H(p)\).
This is the tight-support value bound in
[Le Gall, *Powers of Tensors and Fast Matrix Multiplication*, Theorem 4.1](https://arxiv.org/pdf/1401.7714).

For character certificates the same block-extraction proof applies before
converting the blocks into matrices: use the same frequencies in the six
permuted sources. Their marginal entropies average to the coefficient in (2),
and the block products give \(\prod_sG_\lambda(T_s)^{p_s}\). This explains
why symmetrized profile seeds can be used in this version without assuming
that each unsymmetrized character is balanced.

**Compatibility is checked, not omitted.** Every finite optimization in the
program has an injective marginal map. Exact rational row reduction verifies
that its incidence matrix has full column rank. Its compatible distribution is
therefore unique and \(\Gamma_S(p)=0\). The singleton-convolution seed below
has a separate maximum-entropy proof; its generally nonunique compatibility
class is not passed through that rank test.

Numerical optimization only proposes a probability vector. Round it to positive
rational probabilities and evaluate (2) using directed interval arithmetic.
For an independent upper certificate, take any positive rational reference
marginals \(r_X,r_Y,r_Z\). Relative entropy gives

\[
 \max_p\mathcal L_S(p;v)
 \leq\max_{s\in S}
 \frac{v_s}{[r_X(s_X)r_Y(s_Y)r_Z(s_Z)]^{1/3}}.
 \tag{3}
\]

The program uses lower block endpoints in (2) and upper endpoints in (3).
Thus both recursive lower bounds and upper bounds on the **specified
optimization recipe** are certified. Solver termination, floating tolerances,
and a small gradient are not used as proofs of optimality. `mpmath.iv` uses
40 decimal digits; its binary endpoints are converted exactly to rationals.

## 3. Consecutive-interval partitions and recursive windows

Represent a window by

\[
 W(a,b,c;s)=\sum_{\substack{0\leq i<a,\ 0\leq j<b,\ 0\leq k<c\\i+j+k=s}}
 x_i y_j z_k.
 \tag{4}
\]

Reversing the output identifies \(C(n,n)\) with
\(W(n,n,2n-1;2n-2)\). A consecutive interval on each leg produces another
window, with its lengths and degree shifted exactly. Unused leading/trailing
coordinates are trimmed. Leg permutations and simultaneous reversal of all
three index sets give a canonical state. A state with one dimension 1 is an
oriented dot product after trimming; its symmetrized value is exactly
\(m^{3/4}\).

At each nonterminal window, split each index interval at either its floor
half or its ceiling half. Let \(b_X,b_Y,b_Z\in\{0,1\}\) be the sector labels.
For each

\[
 (w_X,w_Y,w_Z)\in\{-2,-1,1,2\}^3,
\]

retain all supported blocks minimizing
\(w_Xb_X+w_Yb_Y+w_Zb_Z\). Subtracting that minimum makes every supported
coefficient's weight nonnegative. The retained support is tight because each
nonconstant two-sector label has a nonzero weight. Repeated retained sets are
deduplicated. This explicitly enforces monomial-degeneration legality; it does
not delete inconvenient blocks independently.

The child values are obtained recursively. The operational recursion also
includes:

* the singleton laser seed
  \([ab(a+b-1)]^{1/3}\) for a full \(C(a,b)\), proved in Section 7;
* the explicit truncated-convolution to big-CW degeneration from round 2,
  with its three scalar and three dot-product blocks, optimizing (2).

The character recursion additionally includes (1) at full-convolution children.
An upper endpoint for this recursion bounds the values delivered by these
rules, **not the actual character or unrestricted tensor value**.

For the root \(C(n,n)\), I record the best new half-interval calculation
before adding its own profile seed. Otherwise the maximum would simply report
\(F(n,n)^{3/4}\), concealing how much weaker the added calculations were.

### 3.1 Additional enumeration of cut positions

For `2 <= n <= 12`, every root cut

\[
 1\leq c_X,c_Y<n,\qquad 1\leq c_Z<2n-1
\]

was examined, with the same 64 weight vectors and recursive half-interval child
rules. This is **8,712 cut triples and 29,915 distinct retained root cases**.
It is exhaustive over these cut positions and weight choices, not over all
partitions, weights, or higher-power laser analyses.

| n | root cut triples | retained cases | certified new bound interval | lower endpoint for \(F(n,n)^{3/4}\) |
|---:|---:|---:|---:|---:|
| 2 | 2 | 15 | [1.9999, 2.0001] | 2.4669 |
| 3 | 16 | 81 | [3.4669, 3.4670] | 3.9396 |
| 4 | 54 | 235 | [4.9396, 4.9397] | 5.4138 |
| 5 | 128 | 509 | [6.4138, 6.4139] | 6.8885 |
| 6 | 250 | 943 | [7.8885, 7.8886] | 8.3635 |
| 7 | 432 | 1569 | [9.3635, 9.3636] | 9.8386 |
| 8 | 686 | 2427 | [10.8386, 10.8387] | 11.3138 |
| 9 | 1024 | 3549 | [12.3138, 12.3139] | 12.7891 |
| 10 | 1458 | 4975 | [13.7891, 13.7892] | 14.2645 |
| 11 | 2000 | 6737 | [15.2645, 15.2646] | 15.7398 |
| 12 | 2662 | 8875 | [16.7398, 16.7400] | 17.2152 |

These intervals contain \(1+F(n-1,n-1)^{3/4}\). The source of that
particularly simple lower bound is an actual full direct sum:

\[
 C(n,n)\rightsquigarrow 1\oplus C(n-1,n-1).
 \tag{5}
\]

Give \(x_0,y_0\) weight 1, the other inputs weight 0, \(z_0\) weight
\(-2\), and all other outputs weight 0. The cross terms disappear, the
scalar survives, and the remaining convolution uses outputs 2 through
\(2n-2\). Hence

\[
 P(n,n)^\tau\geq1+P(n-1,n-1)^\tau.
 \tag{6}
\]

This family survives the assignment at **every** size. Indeed,
\(\phi(n)^3/n\) is nondecreasing, since

\[
 n(3n+1)^3-27n^3(n+1)=9n^2+n>0.
\]

Consequently \(f_n/n\) is nondecreasing for
\(f_n=F(n,n)^{3/4}\), because

\[
 \frac{f_n}{n}=
 \left(\frac{\phi(n)^3}{n}\right)^{1/4}
 \left(\frac{3-1/n}{2}\right)^{3/4}.
\]

Since \(f_1=1\), \(f_n-f_{n-1}\geq1\). This also handles repeated scalar
peeling. The numerical intervals do not constitute a proof that (5) is the
optimal root partition at every size.

## 4. The zero/rest split: what can legally survive

The informal split into \(x=0\), \(y=0\), and the remaining
\(C(n-1,n-1)\) is not already a tight block decomposition. The dots use
outputs \(1,\ldots,n-1\), while the remaining convolution uses
\(2,\ldots,2n-2\). For \(n\geq3\) they overlap. Assigning separate output
labels to these overlapping occurrences would duplicate variables.

There are at least two legal choices: discard the cross dots to obtain (5),
or preserve the dots and truncate the convolution as follows. Put \(q=n-1\).
Partition

\[
 X,Y:\ \{0\}\mid[1,n-1],\qquad
 Z:\ \{0\}\mid[1,n-1]\mid[n,2n-2].
\]

Give these input sectors weights \((0,1)\) and output sectors weights
\((0,-1,-2)\). The possible supported sectors and their weights are

| sector | block | weight |
|---|---|---:|
| 000 | scalar | 0 |
| 011 | \(B_X(q)\) | 0 |
| 101 | \(B_Y(q)\) | 0 |
| 111 | middle part of the remaining convolution, empty when \(n=2\) | 1 |
| 112 | upper window of the remaining convolution | 0 |

After reversing indices, block 112 is

\[
 t_q=\sum_{i+j<q}x_i y_j z_{i+j},
\]

the **truncated** convolution, not full \(C(q,q)\). The surviving support is

\[
 S_0=\{000,011,101,112\}.
 \tag{7}
\]

It is tight, with labels \(b_X+b_Y-b_Z=0\), and its marginal map has full
column rank. Thus (2) applies with block values

\[
 (1,q^{3/4},q^{3/4},v_q).
 \tag{8}
\]

The calculation uses recursive operational lower bounds for \(v_q\). For
an upper bound on the **outer optimization with rank-feasible child values**,
one may instead use \(v_q=q\). The character upper bound is legitimate:
the cyclic group tensor of
order \(q\), whose rank over \(\mathbb C\) is \(q\), degenerates to
\(t_q\) by weights \(i,j,-k\) with \(k=i+j\pmod q\). Supported weights
are 0 or \(q\). Hence the border rank, and every character value, of
\(t_q\) is at most \(q\); its flattenings show equality of the border rank.

## 5. All-size theorem for the corrected zero/rest laser

**Theorem.** For every \(n\geq2\), the optimized entropy bound (2) on (7),
even with the upper block values \((1,q^{3/4},q^{3/4},q)\), is strictly
less than \(F(n,n)^{3/4}\). This includes any improvement of the recursive
block value compatible with its border-rank bound.

This is an outer-partition theorem with the hypothesis \(v_q\leq q\),
not an unconditional upper bound on \(V^{\rm sym}_{3/4}(t_q)\). An operational
certificate \(V^{\rm sym}_{3/4}(t_q)>q\) would already prove
\(\omega<9/4\) using the truncated tensor itself. That inner possibility is
not ruled out here. Likewise, a character lower bound exceeding its rank
upper bound would show that no character at the assumed exponent exists;
one must not discard such a contradiction by capping the lower bound.

Here is a one-variable dual proof. Choose identical input reference marginals
\((u,1-u)\). Choosing output reference masses proportional to

\[
 u^{-2},\quad\frac{q^{9/4}}{u(1-u)},\quad\frac{q^3}{(1-u)^2}
\]

makes every bound in (3) equal. Therefore the optimized value is at most

\[
 U(q)=\inf_{0<u<1}
 \left[u^{-2}+\frac{q^{9/4}}{u(1-u)}
                    +\frac{q^3}{(1-u)^2}\right]^{1/3}.
 \tag{9}
\]

In fact equality holds in this dual problem. One can symmetrize the input
reference marginals and minimize over the output marginals in the entropy
dual. Only the upper bound is needed for the barrier proof.

For \(q=1\), take \(u=1/2\), giving \(U(1)^3\leq12\). For
\(2\leq q\leq15\), take \(u=1/4\), giving

\[
 U(q)^3\leq R_q:=16+\frac{16}{3}q^{9/4}+\frac{16}{9}q^3.
 \tag{10}
\]

The program certifies \(12^4<F(2,2)^9\) and
\(R_q^4<F(q+1,q+1)^9\) for that finite range with rational interval
endpoints. There are no unresolved comparisons.

For \(q\geq16\), \(q^{3/4}\geq8\) and \(16\leq q^3/256\), so

\[
 R_q\leq\left(\frac{22}{9}+\frac1{256}\right)q^3,
 \qquad
 \left(\frac{22}{9}+\frac1{256}\right)^4<\left(\frac{27}{8}\right)^3.
\]

Meanwhile the exact induction bound \(\phi(q+1)^3\geq q+1\) gives

\[
 F(q+1,q+1)^3
 \geq(q+1)\left(\frac{3q+2}{2}\right)^3
 >\frac{27}{8}q^4.
\]

Cubing the last inequality and comparing with the fourth power of (10)
proves the strict comparison for every remaining \(q\). This completes
the proof without fitting or extrapolating numerical data.

The optimized bound has rate

\[
 U(q)=q+\frac{2\sqrt2}{3}q^{5/8}+o(q^{5/8}).
 \tag{11}
\]

For example, substitute \(u=q^{-3/8}/\sqrt2\) for the matching upper
expansion. For the lower expansion of the expression in (9), use
\((1-u)^{-2}\geq1+2u\), \([u(1-u)]^{-1}\geq u^{-1}\), and
\(2q^3u+q^{9/4}/u\geq2\sqrt2q^{21/8}\). In contrast,

\[
 F(n,n)^{3/4}\sim
 \left(\frac{3}{2\Gamma(4/3)}\right)^{3/4}n,
\]

whose coefficient is about 1.47549, while the rank is \(2n-1\).
The all-size proof above, not these decimal/asymptotic comparisons, certifies
survival of the assignment.

## 6. Certified results through n = 30 and their implications

In this table:

* \(B_{\rm op}\) is a rounded-down **operational lower bound** on
  \(V^{\rm sym}_{3/4}(C(n,n))\): the maximum of the singleton seed,
  recursive interval extraction, and zero/rest extraction.
* \(I_{\rm char}\) is a rounded-up upper certificate for the best **new
  half-interval recipe**, allowing profile-assisted child values.
* \(Z_{\rm op}\) is the rounded-down zero/rest operational lower bound.
* \(U_0\) is a rounded-up entropy-dual upper certificate for zero/rest
  even with \(v_q=q\).
* The displayed \(f_n\) is a lower endpoint: its exact value lies between
  this decimal and this decimal plus \(10^{-4}\), inclusive at \(n=1\).

The upper columns are bounds on the named recipes, not on unrestricted
\(V^{\rm sym}_{3/4}\). All strict comparisons were performed on rational
endpoints before formatting.

| n | \(B_{\rm op}\geq\) | \(I_{\rm char}\leq\) | \(Z_{\rm op}\geq\) | \(U_0\leq\) | \(f_n\) lower | rank |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1 |
| 2 | 2.2894 | 2.0001 | 2.2894 | 2.2895 | 2.4669 | 3 |
| 3 | 3.5568 | 3.0582 | 3.5483 | 3.6309 | 3.9396 | 5 |
| 4 | 4.8318 | 3.9668 | 4.8318 | 4.9664 | 5.4138 | 7 |
| 5 | 6.0822 | 5.1500 | 6.0678 | 6.2866 | 6.8885 | 9 |
| 6 | 7.3434 | 6.1138 | 7.3359 | 7.5906 | 8.3635 | 11 |
| 7 | 8.6042 | 7.2934 | 8.5721 | 8.8798 | 9.8386 | 13 |
| 8 | 9.8648 | 8.2201 | 9.7737 | 10.1557 | 11.3138 | 15 |
| 9 | 11.1252 | 9.4014 | 10.9415 | 11.4201 | 12.7891 | 17 |
| 10 | 12.3856 | 10.4002 | 12.0778 | 12.6742 | 14.2645 | 19 |
| 11 | 13.6458 | 11.5767 | 13.1854 | 13.9193 | 15.7398 | 21 |
| 12 | 14.9060 | 12.5578 | 14.2668 | 15.1565 | 17.2152 | 23 |
| 13 | 16.1662 | 13.7325 | 15.3245 | 16.3863 | 18.6906 | 25 |
| 14 | 17.4263 | 14.6870 | 16.3604 | 17.6097 | 20.1660 | 27 |
| 15 | 18.6864 | 15.8619 | 17.4388 | 18.8271 | 21.6415 | 29 |
| 16 | 19.9465 | 16.7878 | 18.3744 | 20.0392 | 23.1169 | 31 |
| 17 | 21.2065 | 17.9643 | 19.5716 | 21.2464 | 24.5923 | 33 |
| 18 | 22.4666 | 18.8620 | 20.4054 | 22.4490 | 26.0678 | 35 |
| 19 | 23.7266 | 20.0409 | 21.6444 | 23.6473 | 27.5432 | 37 |
| 20 | 24.9866 | 20.9120 | 22.4533 | 24.8418 | 29.0187 | 39 |
| 21 | 26.2466 | 22.0939 | 23.6629 | 26.0326 | 30.4941 | 41 |
| 22 | 27.5066 | 22.9403 | 24.4496 | 27.2200 | 31.9696 | 43 |
| 23 | 28.7666 | 24.1254 | 25.6321 | 28.4042 | 33.4451 | 45 |
| 24 | 30.0266 | 24.9491 | 26.3992 | 29.5855 | 34.9205 | 47 |
| 25 | 31.2866 | 26.1375 | 27.5569 | 30.7637 | 36.3960 | 49 |
| 26 | 32.5465 | 26.9404 | 28.3064 | 31.9393 | 37.8715 | 51 |
| 27 | 33.8065 | 28.1322 | 29.4412 | 33.1124 | 39.3470 | 53 |
| 28 | 35.0665 | 28.9960 | 30.2169 | 34.2830 | 40.8224 | 55 |
| 29 | 36.3265 | 30.1889 | 31.4103 | 35.4514 | 42.2979 | 57 |
| 30 | 37.5864 | 30.8770 | 32.0873 | 36.6175 | 43.7734 | 59 |

The default run checked 237 recursion states, 8,012 retained recursive cases,
and 4,914 entropy optimizations, with zero nonunique compatibility cases.
The exhaustive-root run, started with an empty cache, checked 365 child states
and made 28,674 entropy optimizations. It also found no violation.
An additional check of the profile-assisted truncated recipes for
\(2\leq q\leq29\) certified every upper endpoint strictly below \(q\);
no inner rank contradiction was hidden by the use of \(v_q\leq q\).
This additional check runs after the default experiment's statistics are printed.

### 6.1 What each possible comparison would mean

For a valid new lower bound \(L\) at \(\tau=3/4\):

1. **\(L>F(n,n)^{3/4}\)** refutes the particular balanced assignment.
   By itself it does not prove a growth exponent greater than \(4/3\), or
   rule out another feasible profile with that exponent.
2. **An operational \(L>2n-1\)** gives \(\omega<9/4\) through the
   asymptotic sum inequality. The strict margin can be retained at a slightly
   smaller exponent in a finite extraction certificate.
3. **\(L\leq F(n,n)^{3/4}\)** supplies no refutation. A weak lower bound
   below this threshold does not upper-bound the true value.

Here the new recipe upper certificates themselves are below the assignment,
so this is stronger than merely observing weak feasible lower values. It is
still restricted to the enumerated recipes and the all-size families proved
above. The condition triggering Task 2 did not occur. I did not modify or run
an allegedly strengthened profile LP: no violated inequality was found to add.
The new all-size family (6) also has an explicit surviving witness.

## 7. Attempted general barrier: exact support/quantum values and the gap

There is an explicit useful calculation here. For \(C(a,b)\), assign supported
coefficient \((i,j,i+j)\) probability

\[
 p_{ij}=\frac{\binom{a-1}{i}\binom{b-1}{j}}
 {(a+b-1)\binom{a+b-2}{i+j}}.
 \tag{12}
\]

All three marginals are uniform, of sizes \(a,b,a+b-1\). One proof is
to choose \(u\) uniformly on \([0,1]\), then conditionally sample
independent \(\operatorname{Bin}(a-1,u)\) and
\(\operatorname{Bin}(b-1,u)\) variables. Integrating the binomial probabilities
gives (12); their sum is conditionally \(\operatorname{Bin}(a+b-2,u)\).
The verification code also checks the three rational marginals for all
\(1\leq a,b\leq30\), 900 exact cases.

Moreover, \(\log p_{ij}\) is a sum of a function of \(i\), a function of
\(j\), and a function of \(i+j\). It therefore maximizes joint entropy
among distributions with those marginals: the relative entropy from any
compatible distribution to \(p\) equals \(H(p)-H(q)\). Thus
\(\Gamma(p)=0\). Applying the singleton tight-support laser proves

\[
 V^{\rm sym}_{3/4}(C(a,b))\geq[ab(a+b-1)]^{1/3}.
 \tag{13}
\]

This is a scalar-extraction certificate from the symmetrized source, not an
assertion that a single oriented \(C(a,b)\) has subrank above its smallest
flattening. Round 2, inequality (5), already proves that this seed is at most
\(F(a,b)^{3/4}\) for every size.

The same distribution computes the ordinary quantum functionals. Rescale the
three bases by the factors whose product is \(\sqrt{p_{ij}}\). All factors
are nonzero. The resulting unit tensor has diagonal, uniform reduced density
matrices: two fixed indices determine the third, so the off-diagonal entries
vanish. Every marginal entropy attains its dimension bound. Hence, writing
\(\mathcal Q_\theta\) to avoid confusion with the profile \(F\),

\[
 \boxed{\mathcal Q_\theta(C(a,b))
 =a^{\theta_X}b^{\theta_Y}(a+b-1)^{\theta_Z}.}
 \tag{14}
\]

This uses the orbit-entropy characterization of the quantum functionals in
[Christandl–Vrana–Zuiddam, *Universal Points in the Asymptotic Spectrum of Tensors*](https://ir.cwi.nl/pub/27873/27873.pdf).
The standard inequality \(\mathcal Q_\theta\leq\zeta^\theta\), together
with (12) and the coordinate dimension upper bound, gives the same value for
the upper support functional. Thus even minimizing over bases cannot improve
these particular functionals on convolution.

However, for a matrix block,

\[
 \mathcal Q_\theta(M(e,h,l))
 =(eh)^{\theta_X}(hl)^{\theta_Y}(le)^{\theta_Z}.
 \tag{15}
\]

Its total degree in \(e,h,l\) is **2**, since \(\sum\theta_L=1\).
At uniform weights it is \((ehl)^{2/3}\), whereas the target in this brief
is \((ehl)^{3/4}\). In particular, for \(T_m\) these values are
\(m^2\) and \(m^{9/4}\). An upper bound on
\(\sum_s(e_sh_sl_s)^{2/3}\) does not bound the required sum with exponent
\(3/4\) by the same number. Raising the functional to \(9/8\) changes
its value on a scalar direct sum \(\Delta_r\) to \(r^{9/8}\), destroying
the additive normalization needed for that argument.

There is a concrete illustration in this round's certified data:
\(V^{\rm sym}_{3/4}(C(4,4))\geq4.8318\), whereas its uniform quantum
functional is \(112^{1/3}<4.8203\). Treating that functional as an upper
bound on the \(3/4\)-value would already contradict this small example.

One can combine this resource bound with round 2's single-block monomial
bound, but the elementary combination is insufficient. In the original
convolution coordinates, put \(R=2n-1\) and
\(D=[n^2R]^{1/3}\). For a monomial extraction from \(k\) copies of
\(C(n,n)\), each matrix volume is at most \(R^k\), by the abelian-group
argument of round 2. The uniform quantum functional gives

\[
 \sum_s v_s^{3/4}
 \leq R^{k/12}\sum_s v_s^{2/3}
 \leq[R^{1/12}D]^k.
 \tag{16}
\]

Already at \(n=2\), this upper bound exceeds the assignment: its twelfth
power is \(2^8 3^5=62208\), whereas the assignment's twelfth power is
\((10/3)^9<62208\). As \(n\) grows, (16) has order \(n^{13/12}\),
also too large. This is a valid multi-block bound within the monomial scope,
not the desired barrier. See also the distinction between monomial and
unrestricted comparisons in
[Christandl–Vrana–Zuiddam, *Barriers for Fast Matrix Multiplication from Irreversibility*](https://theoryofcomputing.org/articles/v017a002/v017a002.pdf).

The conclusion is specific: the zero/rest outer optimization with
rank-feasible child values, and the scalar-peeling inequalities, survive at
all sizes; the stated interval searches are certified negative;
the suggested standard functional route does not furnish the missing
\(3/4\)-value upper bound. Arbitrary higher-power repartitioning, merging,
other decompositions, and non-monomial basis changes remain unclosed.

## 8. Reproduction

The following single program contains both experiments and the exact checks.
Save the block to a temporary file, for example `<local path>`, and
run with the system Python. The default run prints the seven numeric columns
of Section 6 and verifies the all-size theorem's finite head and rational tail.
`--all-cuts` runs the additional root-cut enumeration with a fresh cache.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 <local path>
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 <local path> --all-cuts
```

SciPy supplies proposals only; all comparison certificates use rational
arithmetic and interval endpoints. A harmless SLSQP bounds-clipping warning
can occur while generating proposals. No acceptance criterion depends on it.

```python
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import log, comb
import sys
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.iv.dps=40
iv=mp.iv
DEN=10**9
stats=dict(states=0,partitions=0,optimizations=0,nonunique=0)
records={}
@lru_cache(None)
def phi(a):return Q(1) if a==1 else phi(a-1)*Q(3*a-2,3*a-3)
def F(a,b):
 a,b=sorted((a,b));return phi(a)*(b+Q(a-1,2))
def I(q):
 q=Q(q);return iv.mpf(q.numerator)/q.denominator
def dyadic(t):
 sign,man,e,bc=t
 assert bc>=0
 return (-1 if sign else 1)*Q(man)*(Q(2)**e)
def enclose(v):return tuple(dyadic(t) for t in v._mpi_)
@lru_cache(None)
def canonical(a,b,c,s):
 ds=(a,b,c);lo=[max(0,s-(sum(ds)-ds[i]-2)) for i in range(3)]
 hi=[min(d-1,s) for d in ds]
 if any(l>h for l,h in zip(lo,hi)):return None
 ds=tuple(sorted(h-l+1 for l,h in zip(lo,hi)));s-=sum(lo)
 return (*ds,min(s,sum(ds)-3-s))
def blocks(T,cuts):
 a,b,c,s=T;D=(a,b,c)
 out={}
 for bits in product(range(2),repeat=3):
  lows=[cuts[j] if bits[j] else 0 for j in range(3)]
  lens=[D[j]-cuts[j] if bits[j] else cuts[j] for j in range(3)]
  B=canonical(*lens,s-sum(lows))
  if B:out[bits]=B
 return out
@lru_cache(None)
def faces(S):
 ans=set()
 for w in product((-2,-1,1,2),repeat=3):
  scores=[sum(x*y for x,y in zip(w,s)) for s in S];m=min(scores)
  A=tuple(s for s,v in zip(S,scores) if v==m)
  ans.add(A)
 return tuple(sorted(ans))
def exact_rank(A):
 A=[[Q(x) for x in row] for row in A];r=0
 for j in range(len(A[0])):
  piv=next((i for i in range(r,len(A)) if A[i][j]),None)
  if piv is None:continue
  A[r],A[piv]=A[piv],A[r];d=A[r][j];A[r]=[x/d for x in A[r]]
  for i in range(r+1,len(A)):
   d=A[i][j];A[i]=[x-d*y for x,y in zip(A[i],A[r])]
  r+=1
 return r
def entropy_solver(S,values):
 stats['optimizations']+=1
 A=np.array([[int(s[j]==k) for s in S] for j in range(3) for k in sorted(set(s[j] for s in S))],dtype=float)
 # Check injectivity of the marginal map by exact rational elimination.
 if exact_rank(A.astype(int).tolist())!=len(S):
  stats['nonunique']+=1;raise ValueError(S)
 l=np.array([log(float((a+b)/2)) for a,b in values])
 def fun(p):
  q=A@p;return sum(x*log(x) for x in q if x>0)/3-l@p
 def jac(p):return A.T@(np.log(np.maximum(A@p,1e-150))+1)/3-l
 k=len(S)
 if k==1:p=np.ones(1)
 else:
  r=minimize(fun,np.full(k,1/k),jac=jac,bounds=[(1e-12,1)]*k,constraints={'type':'eq','fun':lambda p:p.sum()-1,'jac':lambda p:np.ones(k)},method='SLSQP',options={'ftol':1e-12,'maxiter':150})
  p=r.x
 counts=[max(1,round(float(v)*DEN)) for v in p];total=sum(counts);P=[Q(v,total) for v in counts]
 M=[[sum(P[i] for i,s in enumerate(S) if s[j]==k) for k in sorted(set(s[j] for s in S))] for j in range(3)]
 H=-sum(I(q)*iv.ln(I(q)) for row in M for q in row)/3
 LB=iv.exp(H+sum(I(p)*iv.ln(I(v[0])) for p,v in zip(P,values)))
 ub=[]
 for i,s in enumerate(S):
  prodq=Q(1)
  for j in range(3):prodq*=M[j][sorted(set(v[j] for v in S)).index(s[j])]
  ub.append(I(values[i][1])/iv.exp(iv.ln(I(prodq))/3))
 lower=enclose(LB)[0];upper=max(enclose(v)[1] for v in ub)
 return lower,upper
@lru_cache(None)
def solve(T,mode='character'):
 stats['states']+=1
 a,b,c,s=T
 if a==1:return enclose(iv.exp(iv.ln(I(b**3))/4))
 seeds=[(Q(1),Q(1))]
 if c==a+b-1:seeds.append(enclose(iv.exp(iv.ln(I(a*b*c))/3)))
 if a==b==c and s==a-1 and a>=3:
  CW=((2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0))
  dot=enclose(iv.exp(iv.ln(I((a-2)**3))/4))
  seeds.append(entropy_solver(CW,[(Q(1),Q(1))]*3+[dot]*3))
 if c==a+b-1 and mode=='character':seeds.append(enclose(iv.exp(iv.ln(I(F(a,b)**3))/4)))
 bestlo=max(v[0] for v in seeds);besthi=max(v[1] for v in seeds)
 # Consecutive interval halves; all floor/ceil choices on all three legs.
 cutchoices=[sorted(set((d//2,(d+1)//2))) for d in (a,b,c)]
 seen=set();newlo=Q(1);newhi=Q(1)
 for cuts in product(*cutchoices):
  bs=blocks(T,cuts)
  for S in faces(tuple(bs)):
   kids=tuple(bs[x] for x in S)
   key=(S,kids)
   if key in seen:continue
   seen.add(key);stats['partitions']+=1
   if len(S)==1:L,U=solve(kids[0],mode)
   else:L,U=entropy_solver(S,[solve(k,mode) for k in kids])
   bestlo=max(bestlo,L);besthi=max(besthi,U);newlo=max(newlo,L);newhi=max(newhi,U)
 records[(T,mode)]=(newlo,newhi)
 return bestlo,besthi
def floor4(q):return f'{(q.numerator*10000//q.denominator)/10000:.4f}'
def ceil4(q):return floor4(Q(-((-q.numerator*10000)//q.denominator),10000))
def all_cuts():
 for n in range(2,13):
  T=(n,n,2*n-1,2*n-2);seen=set();L=U=Q(1);ncuts=0
  for cuts in product(range(1,n),range(1,n),range(1,2*n-1)):
   ncuts+=1;bs=blocks(T,cuts)
   for S in faces(tuple(bs)):
    kids=tuple(bs[x] for x in S);key=(S,kids)
    if key in seen:continue
    seen.add(key)
    if len(S)==1:l,u=solve(kids[0])
    else:l,u=entropy_solver(S,[solve(k) for k in kids])
    L=max(L,l);U=max(U,u)
  f=enclose(iv.exp(iv.ln(I(F(n,n)**3))/4))[0]
  assert U<f
  print(n,ncuts,len(seen),floor4(L),ceil4(U),floor4(f),flush=True)
 print(stats)

if __name__=='__main__':
 if '--all-cuts' in sys.argv:
  all_cuts();sys.exit(0)
 Z=((0,0,0),(0,1,1),(1,0,1),(1,1,2))
 for n in range(1,31):
  T=(n,n,2*n-1,2*n-2);C=solve(T);O=solve(T,'operational')
  fl,fu=enclose(iv.exp(iv.ln(I(F(n,n)**3))/4))
  new=records.get((T,'character'),(Q(1),Q(1)))
  if n==1:zlo=zhi=zupper=Q(1)
  else:
   q=n-1;D=enclose(iv.exp(iv.ln(I(q**3))/4))
   tr=solve((q,q,q,q-1),'operational')
   zlo,zhi=entropy_solver(Z,[(Q(1),Q(1)),D,D,tr])
   _,zupper=entropy_solver(Z,[(Q(1),Q(1)),D,D,(Q(q),Q(q))])
   assert new[1]<fl and zupper<fl
  print(n,floor4(max(O[0],zlo)),ceil4(new[1]),floor4(zlo),ceil4(zupper),floor4(fl),2*n-1,flush=True)
 print(stats)
 # All-size zero/rest proof: finite head and rational tail.
 for q in range(1,16):
  if q==1:R=Q(12)
  else:R=enclose(I(16)+I(Q(16,3))*iv.exp(iv.ln(I(q**9))/4)+I(Q(16,9)*q**3))[1]
  assert R**4<F(q+1,q+1)**9
 assert (Q(22,9)+Q(1,256))**4 < Q(27,8)**3
 print('zero/rest finite-head and all-size tail: PASS')

 for a,b in product(range(1,31),repeat=2):
  p={(i,j):Q(comb(a-1,i)*comb(b-1,j),
               (a+b-1)*comb(a+b-2,i+j))
     for i in range(a) for j in range(b)}
  assert all(sum(p[i,j] for j in range(b))==Q(1,a) for i in range(a))
  assert all(sum(p[i,j] for i in range(a))==Q(1,b) for j in range(b))
  assert all(sum(p[i,j] for i in range(a) for j in range(b) if i+j==k)
             ==Q(1,a+b-1) for k in range(a+b-1))
 print('Uniform marginals: 900 exact checks PASS')
 assert 2**8*3**5 > F(2,2)**9
 print('Quantum/group combined bound does not exclude F at n=2: PASS')
 for q in range(2,30):
  assert solve((q,q,q,q-1),'character')[1]<q
 print('Truncated character recipes q=2..29: all upper endpoints below q PASS')
```

Observed terminal certificates:

```text
Default: states=237, partitions=8012, optimizations=4914, nonunique=0
All cuts: states=365, partitions=11539, optimizations=28674, nonunique=0
zero/rest finite-head and all-size tail: PASS
Uniform marginals: 900 exact checks PASS
Quantum/group combined bound does not exclude F at n=2: PASS
Truncated character recipes q=2..29: all upper endpoints below q PASS
```
