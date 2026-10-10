# Codex round 1: polynomial-multiplication profile inequalities

Read against `BRIEF-codex.md` and `../allfields/paper94.tex`, especially Sections 2–5. All computations below use exact rational arithmetic; the one family involving fourth roots uses certified rational enclosing intervals. The executable verification code is included at the end of this file. No other project file was changed.

**Result:** I found valid additional degenerations, including a two-input-sector construction with shared output and a genuine bivariate degeneration to a full direct sum. None of the valid inequalities tested excludes P*. I can prove a barrier for a specified enlarged constraint system, including nonsymmetrized entropy optimization. I cannot prove that the entire polynomial-multiplication approach has a 9/4 barrier.

The distinction matters: excluding this particular P* would not, by itself, prove a better exponent. Another profile with growth exponent 4/3 might survive. A finite violation also need not produce a uniform exponent gap.

## 1. Conventions and how the inequalities are obtained

Write

\[
 F(a,b)=P^*(a,b)=\phi(\min(a,b))\left(\max(a,b)+\frac{\min(a,b)-1}{2}\right),
 \qquad \phi(a)=\prod_{j=1}^{a-1}\frac{3j+1}{3j}.
\]

All profile tests use t=3/4. For example,

\[
 F(2,2)=10/3,\quad F(2,3)=14/3,\quad
 F(3,3)=56/9,\quad F(4,4)=770/81.
\]

Here and below a sum with no summands means zero. A convolution with a nonpositive input length is omitted, not treated as a positive tensor.

The paper's separation proposition requires **already matched sectors on two legs**, with no terms between unmatched sectors (`paper94.tex`, lines 277–289). It cannot be applied directly to a grid of arbitrary input blocks. Its Fourier weights use the same existing sector label on both of those legs (lines 303–332).

If a degeneration leaves blocks S_i with shared leg L, its character consequence is

\[
 \lambda(\text{source})\geq e^{p_L H(q)}\prod_i\lambda(S_i)^{q_i}
 \tag{E}
\]

for every probability vector q. For p_L>0, optimization gives

\[
 \lambda(\text{source})\geq
 \left(\sum_i\lambda(S_i)^{1/p_L}\right)^{p_L}.
 \tag{E-opt}
\]

For p_L=0 the corresponding bound is the maximum of the block values. Multiplying (E) for the six leg permutations with a **common q**, then optimizing q, gives a sum of profile values. This argument works with shared Z just as with shared X, since the sum of the six permuted p_Z values is also 6t.

For a **full direct sum**, additivity instead gives an entropy coefficient 1 before symmetrization. The resulting profile inequality has sums of P^t, not sums of P. This stronger inequality is tested separately below.

## 2. Candidate A: transfer length between the two inputs

This directly realizes the requested two-input-sector/shared-output construction, using an auxiliary B_Z(m).

### Statement

Set r=m−1 and

\[
 I=\{i:0\leq i<\min(a,m),\quad 0\leq r-i<\min(b,m)\}.
\]

There is a degeneration

\[
 C(a,b)\otimes B_Z(m)\rightsquigarrow
 \sum_{i\in I}^{\text{shared }Z}
 C(a+m-1-2i,\ b-m+1+2i).
 \tag{A-tensor}
\]

The X and Y sectors are disjoint; all outputs have length a+b−1 and use the same Z space. Consequently

\[
 m^{p_Z}\lambda(C(a,b))\geq
 e^{p_ZH(q)}\prod_{i\in I}
 \lambda(C(a+m-1-2i,b-m+1+2i))^{q_i},
 \tag{A-character}
\]

and

\[
 \boxed{mP(a,b)\geq\sum_{i\in I}P(a+m-1-2i,b-m+1+2i).}
 \tag{A}
\]

In particular, for a,b≥2,

\[
 \boxed{2P(a,b)\geq P(a+1,b-1)+P(a-1,b+1).}
 \tag{A2}
\]

This is concavity along constant a+b, a different direction from the paper's separate concavity.

### Explicit construction for m=2

Let V_e be binary forms of degree e in u,v, and use s,w for a second two-dimensional space. Interpret B_Z(2), after an invertible change of basis, as the alternating pairing on span(s,w). The bilinear map is

\[
 (f_s s+f_w w,\ g_s s+g_w w)\longmapsto f_sg_w-f_wg_s.
\]

On each input use the quotient/kernel decomposition under (s,w)=(u,v), with D=uw−vs, as in the paper's determinant sequence. A kernel element is Dh=−vh\,s+uh\,w. Thus:

* kernel × kernel is zero;
* quotient × kernel is multiplication of the quotient's diagonal image by h;
* kernel × quotient is the same multiplication with a minus sign.

Give each input quotient weight +1 and each input kernel weight −1; give Z weight zero. The quotient × quotient part has weight 2 and disappears. The apparent weight −2 part is identically zero. The two weight-zero blocks are exactly C(a+1,b−1) and C(a−1,b+1), up to a removable sign.

### General m: explicit filtration and weights

Put W_{e,f}=V_e(u,v)\otimes V_f(s,w). Filter each input by powers of D:

\[
 F_iW_{a-1,r}=D^iW_{a-1-i,r-i},\qquad
 F_i/F_{i+1}\simeq V_{a+r-1-2i}.
\]

The second formula's degree is a+m−2−2i. Quotient lifts can be chosen explicitly: for 0≤ℓ≤e+f−2i, choose α=min(ℓ,e−i), β=ℓ−α and lift the diagonal monomial by

\[
 D^i u^{e-i-\alpha}v^\alpha s^{f-i-\beta}w^\beta.
 \tag{lifts}
\]

Pair the W_r factors apolarly, with normalization

\[
 \langle s^{r-h}w^h,s^{r-k}w^k\rangle
 =\begin{cases}(-1)^h/\binom rh,&h+k=r,\\0,&h+k\ne r.\end{cases}
\]

This is a nondegenerate pairing, so the auxiliary tensor is B_Z(m). At any fixed nonzero (u,v), two forms divisible by D^i and D^j pair to zero if i+j>r. For i+j=r the induced quotient pairing is a nonzero scalar times multiplication of the diagonal images. The scalar in these lifts is (−1)^i/binom(r,i), and can be absorbed into X_i.

Assign weights −i to X_i, −j to Y_j, and r to every Z coordinate. All nonzero terms have weight r−i−j≥0. Exactly i+j=r survives, giving (A-tensor). The included coefficient check verifies the zero and leading blocks for 1≤a,b,m≤4: **780 block-basis-pair checks**.

### P* test and an all-size proof of survival

At a=b=m=2, (A) is 20/3≥6, with slack 2/3. Exact enumeration for 1≤a,b≤50 and 1≤m≤30 gives **75,000 cases, zero violations**.

More strongly, F satisfies (A) for all positive a,b,m. Here is a proof, not an extrapolation of that enumeration.

For s≥2, let f_s(x)=F(x,s−x) for integer 1≤x≤s−1, and set f_s(0)=f_s(s)=0. This sequence is concave. For 2≤a≤b, the slack in (A2), divided by φ(a), is

\[
 \begin{cases}
 2/(3a-2),&b=a\text{ or }a+1,\\
 2(2a+b-1)/(3a(3a-2)),&b\geq a+2.
 \end{cases}
\]

At the endpoint a=1,b≥2 the slack is 2(b+1)/3; s=2 is immediate. Hence the piecewise-linear extension of f_s is concave on [0,s].

If m<s, the first coordinates A_i=a+m−1−2i of the surviving blocks, supplemented by max(m−a,0) copies of 0 and max(m−b,0) copies of s, form a multiset of **m numbers whose sum is ma**. Jensen's inequality gives (A), since both endpoint values vanish. If m≥s there are no surviving blocks. Thus this whole family cannot exclude F.

## 3. Candidate B: higher-dimensional determinant filtration on one input

The paper's B_X(2) construction extends to every m:

\[
 C(a,b)\otimes B_X(m)\rightsquigarrow
 \sum_{i=0}^{\min(b,m)-1,\ \text{shared }X}
 C(a,b+m-1-2i).
\]

Use the D-adapted lifts above on the second input W_{b−1,m−1} and output W_{a+b−2,m−1}. Multiplication by a form in V_{a−1} preserves every power of D. Give second-input grade i weight −i and output-dual grade j weight +j. The surviving diagonal blocks are the displayed convolutions. Therefore

\[
 \boxed{mP(a,b)\geq\sum_{i=0}^{\min(b,m)-1}P(a,b+m-1-2i).}
 \tag{B}
\]

The character version has m^{p_X} on the left and p_XH(q) on the right, exactly as in (A-character).

Exact test: **75,000 cases, zero violations**, on the same box as (A). For a=b=m=2, both sides are 20/3, so this family includes tight points.

There is also an all-size reason. For fixed a, extend F(a,0)=0. The resulting sequence is concave, including at zero: 2F(a,1)≥F(a,2). The min(b,m) positive lengths in (B) sum to bm. Adding m−min(b,m) zeros gives m lengths with mean b. Jensen proves (B). Thus merely replacing B_X(2) by a larger dot product provides no escape.

## 4. Candidate C: bivariate multiplication, then an anti-diagonal grade

This is a degeneration of a tensor product, rather than the carry-free Kronecker restriction already in the brief.

### Statement

Put r_X=min(a,c), r_Y=min(b,d), A_i=a+c−1−2i, and B_j=b+d−1−2j. For any 0≤k≤r_X+r_Y−2,

\[
 \boxed{P(a,b)P(c,d)\geq
 \sum_{\substack{0\leq i<r_X,\ 0\leq j<r_Y\\i+j=k}}
 P(A_i,B_j).}
 \tag{C}
\]

Before symmetrization, the inequality is

\[
 \lambda(C(a,b))\lambda(C(c,d))\geq
 e^{p_ZH(q)}\prod_{i+j=k}\lambda(C(A_i,B_j))^{q_{ij}}.
 \tag{C-character}
\]

### Construction

View the product as multiplication

\[
 W_{a-1,c-1}\times W_{b-1,d-1}
 \longrightarrow W_{a+b-2,c+d-2}.
\]

Use the same D-adapted bases on all three spaces. Products of input filtration grades i,j lie in output filtration grade at least i+j. Weights (−i,−j,+k) therefore have nonnegative total weight, and retain exactly the associated-graded multiplication with k=i+j. On the quotient, diagonal substitution identifies the block with ordinary multiplication

\[
 V_{A_i-1}\times V_{B_j-1}\longrightarrow V_{A_i+B_j-2}.
\]

Now keep only output grade k. Each X_i is matched to exactly one Y_{k−i}; these sectors share that output space. Apply separation with shared Z. This succeeds precisely because the preliminary filtration has created a genuine grading without carries.

The included exact basis-matrix check verifies the associated-graded coefficients for all 1≤a,b,c,d≤3, including absence of negative-weight terms.

### P* result

For 1≤a,b,c,d≤9 and every allowed k: **39,609 inequalities, zero violations**. At a=b=c=d=2, the left side is 100/9; the three k values give right sides 56/9, 6, and 1.

I have a proof of the tensor inequality, but **not an all-size proof that F satisfies (C)**. The finite check is not such a proof. This family remains a legitimate target for a larger or symbolic search.

## 5. Candidate D: extract a full direct sum from the bivariate grading

Do not discard output grades after the first degeneration in Section 4. On its support k=i+j, give the three grades the further weights

\[
 2i^2,\quad 2j^2,\quad -k^2.
\]

Their sum is (i−j)^2. Thus only i=j survives, and k=2i identifies the same sector on all three legs. With R=min(a,b,c,d),

\[
 C(a,b)\otimes C(c,d)\rightsquigarrow
 \bigoplus_{i=0}^{R-1} C(a+c-1-2i,b+d-1-2i).
 \tag{D-tensor}
\]

The two stages can be combined into one integer weight by multiplying the first-stage weights by an integer larger than the absolute values of all second-stage total weights. This also checks compatibility with the paper's stated degeneration convention.

Additivity gives the nonsymmetrized inequality

\[
 \lambda(C(a,b))\lambda(C(c,d))\geq
 \sum_{i=0}^{R-1}\lambda(C(a+c-1-2i,b+d-1-2i)).
\]

Applying generalized Hölder to the six sums gives the stronger profile statement

\[
 \boxed{[P(a,b)P(c,d)]^t\geq
 \sum_{i=0}^{R-1}P(a+c-1-2i,b+d-1-2i)^t.}
 \tag{D}
\]

For F and t=3/4, exact rational interval certification on 1≤a,b,c,d≤10 gives **10,000 cases, zero violations and zero unresolved comparisons**. Single-summand cases are compared directly as rationals. Multiple-summand cases use rational lower/upper bounds for fourth roots of exact rational cubes; no floating-point tolerance is used.

For example, at a=b=c=d=2, the comparison is

\[
 (100/9)^{3/4}\geq (56/9)^{3/4}+1.
\]

The left side is greater than 6 and the right side less than 5: these claims follow exactly from (100/9)^3>6^4 and (56/9)^3<4^4. As with (C), I do not have an all-size survival proof for F.

## 6. Fourier filtering on the original input blocks: cost accounting

Partition the inputs of C(Ma,Mb) into M consecutive blocks of lengths a and b. Label them h,j∈{0,…,M−1}. A desired anti-diagonal consists of h+j=M−1.

Here is a valid Fourier restriction with an explicit auxiliary cost. Put L=2M−1 and tensor the source with B_Z(L). If ρ labels the dot-product index, substitute

\[
 x_{ha+i,\rho}\mapsto\zeta^{\rho h}X_{h,i},\qquad
 y_{jb+\ell,\rho}\mapsto L^{-1}\zeta^{\rho(j-M+1)}Y_{j,\ell}.
\]

Leave Z unchanged. Summing over ρ gives the exact indicator h+j=M−1, since the possible mismatch has absolute value at most M−1. Each surviving block is a translated C(a,b), and the blocks have disjoint X,Y sectors. Their Z supports may overlap, which is allowed for the shared leg.

Separation consequently gives

\[
 L^{p_Z}\lambda(C(Ma,Mb))\geq M^{p_Z}\lambda(C(a,b)),
 \qquad
 \boxed{LP(Ma,Mb)\geq MP(a,b).}
 \tag{F}
\]

This survives F for all parameters simply by monotonicity and L≥M. The code also checks 1≤M,a,b≤10: **1,000 exact cases**. For M=a=b=2 it reads 770/27≥20/3.

Using L full scalar copies instead gives the weaker factor L^{1/t} on the profile side. Neither copy cost may be silently deleted. This construction is not a proof that every more elaborate Fourier-assisted approach must incur this cost.

### A precise obstruction to direct coordinate weights

For C(4,4), try to retain both full anti-diagonal 2×2 input blocks. Write W_ij=α_i+β_j+γ_{i+j}. All eight desired terms must have weight zero. But the exact linear identity

\[
 W_{11}+W_{22}=W_{21}+W_{02}+W_{13}-W_{03}=0
\]

forces the sum of two undesired weights to be zero. They cannot both be strictly positive. This rules out that coordinate-weight construction, not arbitrary changes of basis or assisted degenerations. It also avoids the overbroad claim that every anti-diagonal is impossible: singleton blocks can be isolated by output restriction.

## 7. Asymmetry: valid additional bounds, and a balanced feasible lift

The outer-product tensor B_Y(a)⊗B_X(b) has terms x_i y_j z_{ij}. The substitution z_{ij}↦z_{i+j} restricts it to C(a,b). Thus

\[
 \lambda(C(a,b))\leq a^{p_Y}b^{p_X},\qquad
 \boxed{P(a,b)\leq ab.}
 \tag{O}
\]

Setting one input to its first coordinate gives B_X(b) or B_Y(a); retaining output index min(a,b)−1 gives B_Z(min(a,b)). Consequently

\[
 \lambda(C(a,b))\geq
 \max\{b^{p_X},a^{p_Y},\min(a,b)^{p_Z}\}.
 \tag{L}
\]

These are useful nonsymmetric constraints, but do not exclude the following formal assignment:

\[
 p_X=p_Y=p_Z=t=3/4,\qquad
 L_\pi(a,b)=F(a,b)^{3/4}\quad\text{for every }\pi.
 \tag{balanced lift}
\]

Its fourth powers are exactly rational: L_π(a,b)^4=F(a,b)^3. It has the correct dot-product boundary values and the required symmetrized profile. The rank upper bound becomes F(a,b)^3≤(a+b−1)^4; (O) becomes F(a,b)≤ab; (L) becomes F(a,b)≥max(a,b). All pass the exact grid checks below, and all hold globally. For (O), φ(a)≤2a/3 for a≥2, by induction from a=2; for a≤b this implies F(a,b)≤ab. The a=1 case is equality. The lower bound follows immediately from φ≥1.

In this lift, optimizing (E) **separately for every permutation** reduces exactly to

\[
 S^{3/4}\geq\left(\sum_i F_i\right)^{3/4}
 \quad\Longleftrightarrow\quad S\geq\sum_iF_i,
\]

where S is the corresponding source profile, including auxiliary factors. Thus unequal optimizing q's cannot help against this assignment. In the paper's tripling construction all three block values are equal in the lift, so the optimal q is uniform. The same observation applies to the new shared-leg constructions. For full direct sums it reduces to (D), which was tested with the actual 3/4 power.

This is **a feasible assignment to these constraints, not a constructed tensor character**. In particular, its consistency with every other tensor relation has not been shown. It is not permissible to argue that actual characters have different permuted values and hence discard the balanced case without a theorem forcing that difference at the relevant t.

## 8. A restricted barrier theorem, and the remaining gap

**Restricted barrier.** Consider the constraint system consisting of the paper's profile axioms and their nonsymmetrized entropy forms, all the B_X(m) filtrations (B), all the B_Z(m) transfer filtrations (A), the outer-product and dot-product bounds (O),(L), and the explicitly costed Fourier inequalities (F). Include all permutations, all probability vectors, and consequences obtained by multiplying and combining these numerical inequalities. This system admits the balanced lift at t=3/4. It therefore cannot force t<3/4.

Proof: the brief's original P* witness handles the original axioms. Sections 2–3 prove all-size survival for both additional filtration families, Section 6 handles Fourier filtering, and Section 7 handles the nonsymmetrized constraints and independent entropy optimization. A feasible assignment remains feasible under logical consequences of the constraints. This is a barrier for a specified inequality system; it does not assert realizability in the tensor semiring.

For completeness, the original witness's relevant all-size checks can also be made directly. For a≤b, away from the diagonal its first-coordinate concavity slack is

\[
 2F(a,b)-F(a-1,b)-F(a+1,b)
 =\phi(a)\frac{2(b-a)}{3a(3a-2)}\quad(a\geq2,b\geq a+1).
\]

The cross-diagonal slack is zero and the long-coordinate tail is linear. For tripling with h<a, the required comparison is

\[
 \phi(a)(a+2h-1)\geq\phi(h)(2a+h-1).
\]

The ratio on the left after division by 2a+h−1 is nondecreasing in a≥h: the cross-multiplied difference for its successive ratio is 2(a−h)(a−h+1)≥0. For h≥a, tripling is equality. In fact **h=a−1 is also equality**, an additional tight case worth preserving in searches.

For the rank upper bound, at fixed a the ratio F(a,b)^3/(a+b−1)^4 decreases for b≥a, so it suffices to check the diagonal. The diagonal ratio decreases from 1 at a=1. Its successive-ratio comparison has numerator

\[
 (3a(3a-1))^3(2a+1)^4-((3a+1)(3a+2))^3(2a-1)^4,
\]

whose expansion in x=a−1 has coefficients, from degree 8 down to 0,

\[
 7776,\ 69552,\ 266544,\ 570348,\ 743614,\ 603680,\
 297618,\ 81476,\ 9496.
\]

They are all positive. This proves the required upper bound without floating-point powers.

I do **not** extend the barrier theorem to (C),(D) for all sizes, to all possible multivariate degenerations, or to arbitrary polynomial-multiplication tensor products. The finite searches found no violation of (C),(D), but global feasibility remains unproved. Conversely, no valid inequality in this round refutes the 9/4 witness. A next search can focus on richer extraction patterns from the associated-graded bivariate tensor, with every full-direct-sum gain and auxiliary cost retained.

### Rejected tempting improvement: remove Kronecker carries for free

The claimed degeneration C(4,4)⇝C(2,2)⊗C(2,2) would imply

\[
 770/81=F(4,4)\geq F(2,2)^2=900/81,
\]

so it really would exclude F, with deficit 130/81. However, its output flattening ranks are 7 and 9. Flattening rank cannot increase under degeneration. The degeneration is impossible. The code checks the exact deficit; the rank obstruction is independent of P*.

## 9. Reproducible exact verification

`barrier_check.py` uses floating-point powers and tolerances; its numerical output alone is not an exact certificate. `barrier_extra.py` uses Fractions for its additional tests. The independent script below uses Fractions throughout, SymPy exact matrices for the bivariate basis check, and integer fourth-root bounds for (D).

**Observed execution:** the embedded script completed with exit code 0: **239,009 exact inequality comparisons**, **780 apolar coefficient checks**, and **1,296 bivariate basis-product checks**, all passing. There were no unresolved fourth-root comparisons. These counts describe the finite verification, independently of the all-size proofs explicitly given above.

Run the embedded script without creating another file:

```sh
/usr/bin/python3 - <<'PY'
from pathlib import Path
s = Path('CODEX-ROUND1.md').read_text()
code = s.split('```python\n', 1)[1].split('\n```', 1)[0]
exec(compile(code, 'CODEX-ROUND1.md:verification', 'exec'))
PY
```

```python
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb
from sympy import Matrix, Rational, integer_nthroot, symbols, expand

@lru_cache(None)
def phi(a):
    return Q(1) if a == 1 else phi(a-1) * Q(3*a-2, 3*a-3)

@lru_cache(None)
def F(a, b):
    a, b = sorted((a, b))
    assert a >= 1
    return phi(a) * (b + Q(a-1, 2))

counts = {}
def ge(name, lhs, rhs):
    counts[name] = counts.get(name, 0) + 1
    assert lhs >= rhs, (name, lhs, rhs)

for a, b in product(range(1, 81), repeat=2):
    assert F(a, b) == F(b, a)
    ge('rank upper (cubed)', (a+b-1)**4, F(a, b)**3)
    ge('outer product upper', a*b, F(a, b))
    ge('dot product lower', F(a, b), max(a, b))
    ge('monotonicity', F(a, b+1), F(a, b))
    if b >= 2:
        ge('separate concavity', 2*F(a, b), F(a, b-1)+F(a, b+1))
    ge('tripling', F(a, 3*b+a-1), 3*F(a, b))
    if b >= a-1:
        assert F(a, 3*b+a-1) == 3*F(a, b)
for a in range(1, 81):
    assert F(1, a) == a
    ge('concavity at zero', 2*F(a, 1), F(a, 2))

for a, b in product(range(1, 51), repeat=2):
    for m in range(1, 31):
        transfer = [F(a+m-1-2*i, b-m+1+2*i)
                    for i in range(min(a, m))
                    if 0 <= m-1-i < min(b, m)]
        higher = [F(a, b+m-1-2*i) for i in range(min(b, m))]
        ge('A transfer', m*F(a, b), sum(transfer, Q(0)))
        ge('B higher dot', m*F(a, b), sum(higher, Q(0)))
        # Exact endpoint/mean identity used in the all-size Jensen proof.
        if m < a+b:
            coords = [a+m-1-2*i for i in range(min(a, m))
                      if 0 <= m-1-i < min(b, m)]
            coords += [0]*max(m-a, 0) + [a+b]*max(m-b, 0)
            assert len(coords) == m and sum(coords) == m*a

for a, b, c, d in product(range(1, 10), repeat=4):
    for k in range(min(a, c)+min(b, d)-1):
        terms = [F(a+c-1-2*i, b+d-1-2*(k-i))
                 for i in range(min(a, c)) if 0 <= k-i < min(b, d)]
        ge('C bivariate shared output', F(a, b)*F(c, d), sum(terms, Q(0)))

SCALE = 10**24
@lru_cache(None)
def power_interval(q):
    # Bound q^(3/4) using only integer arithmetic and rational comparisons.
    cube = q**3
    n = integer_nthroot(cube.numerator*SCALE**4 // cube.denominator, 4)[0]
    lo = Q(int(n), SCALE)
    if lo**4 == cube:
        return lo, lo
    hi = Q(int(n)+1, SCALE)
    assert lo**4 < cube < hi**4
    return lo, hi

for a, b, c, d in product(range(1, 11), repeat=4):
    lhs = F(a, b)*F(c, d)
    terms = [F(a+c-1-2*i, b+d-1-2*i) for i in range(min(a, b, c, d))]
    if len(terms) == 1:
        ge('D full direct sum', lhs, terms[0])
    else:
        # Failure here would be reported, not accepted by a tolerance.
        ge('D full direct sum', power_interval(lhs)[0],
           sum((power_interval(x)[1] for x in terms), Q(0)))

for m, a, b in product(range(1, 11), repeat=3):
    ge('F Fourier with cost', (2*m-1)*F(m*a, m*b), m*F(a, b))

# Exact coefficients of adapted D-filtration bases.
# Keys (i,j) denote u^(r-i)v^i s^(s-j)w^j.
def lift(r, s, k, ell):
    i = min(ell, r-k)
    j = ell-i
    return {(i+h, j+k-h): Q((-1)**h * comb(k, h)) for h in range(k+1)}

@lru_cache(None)
def basis(r, s):
    return [(k, ell, lift(r, s, k, ell))
            for k in range(min(r, s)+1)
            for ell in range(r+s-2*k+1)]

apolar_checks = 0
for a, b, m in product(range(1, 5), repeat=3):
    r = m-1
    for i, ell, x in basis(a-1, r):
        for j, n, y in basis(b-1, r):
            if i+j < r:
                continue
            out = {}
            for (v, w), cx in x.items():
                for (vv, ww), cy in y.items():
                    if w+ww == r:
                        out[v+vv] = out.get(v+vv, Q(0)) + cx*cy*Q((-1)**w, comb(r, w))
            out = {k: v for k, v in out.items() if v}
            expected = {} if i+j > r else {ell+n: Q((-1)**i, comb(r, i))}
            assert out == expected, (a, b, m, i, j, out, expected)
            apolar_checks += 1

@lru_cache(None)
def inverse_basis(r, s):
    bs = basis(r, s)
    size = (r+1)*(s+1)
    mat = Matrix.zeros(size, size)
    for col, (_, _, p) in enumerate(bs):
        for (i, j), coeff in p.items():
            mat[i*(s+1)+j, col] = Rational(coeff.numerator, coeff.denominator)
    return mat.inv()

bivariate_checks = 0
for a, b, c, d in product(range(1, 4), repeat=4):
    rr, ss = a+b-2, c+d-2
    out_basis = basis(rr, ss)
    inv = inverse_basis(rr, ss)
    for i, ell, x in basis(a-1, c-1):
        for j, n, y in basis(b-1, d-1):
            raw = Matrix.zeros((rr+1)*(ss+1), 1)
            for (v, w), cx in x.items():
                for (vv, ww), cy in y.items():
                    coeff = cx*cy
                    raw[(v+vv)*(ss+1)+w+ww] += Rational(coeff.numerator, coeff.denominator)
            transformed = inv*raw
            for coeff, (k, deg, _) in zip(transformed, out_basis):
                if k < i+j:
                    assert coeff == 0
                elif k == i+j:
                    assert coeff == int(deg == ell+n)
            bivariate_checks += 1

# The coordinate-weight obstruction is an identity in all leg weights.
coeffs = [0]*15  # 4 X weights, 4 Y weights, 7 Z weights.
for i, j, sign in [(1,1,1), (2,2,1), (2,1,-1), (0,2,-1), (1,3,-1), (0,3,1)]:
    coeffs[i] += sign
    coeffs[4+j] += sign
    coeffs[8+i+j] += sign
assert coeffs == [0]*15

assert F(2,2) == Q(10,3) and F(2,3) == Q(14,3)
assert F(3,3) == Q(56,9) and F(4,4) == Q(770,81)
assert 2*F(2,2) - F(3,1) - F(1,3) == Q(2,3)
assert F(2,2)**2 - F(4,4) == Q(130,81)
assert (F(2,2)**2)**3 > 6**4 and F(3,3)**3 < 4**4
assert 7 < 9  # Output flattenings of the rejected carry-free degeneration.

# Symbolic certificates used in the all-size survival proof.
x = symbols('x')
aa = x+1
upper_numerator = (3*aa*(3*aa-1))**3*(2*aa+1)**4 - ((3*aa+1)*(3*aa+2))**3*(2*aa-1)**4
positive_coefficients = [9496, 81476, 297618, 603680, 743614, 570348, 266544, 69552, 7776]
assert expand(upper_numerator - sum(c*x**i for i, c in enumerate(positive_coefficients))) == 0
aa, hh = symbols('aa hh')
tripling_numerator = (3*aa+1)*(aa+2*hh)*(2*aa+hh-1) - 3*aa*(aa+2*hh-1)*(2*aa+hh+1)
assert expand(tripling_numerator - 2*(aa-hh)*(aa-hh+1)) == 0

assert counts['A transfer'] == counts['B higher dot'] == 75000
assert counts['C bivariate shared output'] == 39609
assert counts['D full direct sum'] == 10000
assert counts['F Fourier with cost'] == 1000
assert sum(counts.values()) == 239009
assert apolar_checks == 780 and bivariate_checks == 1296
for name, count in counts.items():
    print(f'{name}: {count} exact comparisons, PASS')
print(f'apolar leading/zero coefficients: {apolar_checks}, PASS')
print(f'bivariate adapted-basis products: {bivariate_checks}, PASS')
print('All displayed witnesses and the coordinate-weight identity: PASS')
print('No floating-point comparisons; no unresolved root intervals.')
```

The external literature check was limited to the spectral framework: [Christandl–Vrana–Zuiddam, *Universal points in the asymptotic spectrum of tensors*](https://arxiv.org/abs/1709.07851). The constructions and P* calculations above are derived from the supplied local paper, with self-contained proofs; they do not rely on a claimed literature novelty result.
