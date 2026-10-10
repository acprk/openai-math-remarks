# Codex round 2: all-size survival and further extraction barriers

Read against `BRIEF-codex-R2.md`, `CODEX-ROUND1.md`, `BRIEF-codex.md`, and Sections 2–5 of `../allfields/paper94.tex`.

**Results.** Both (C) and (D) survive for **all positive integer sizes**, with proofs below. There are also all-size barriers for (i) extracting one matrix-multiplication tensor by monomial degenerations in the standard convolution coordinates, with arbitrary tensor powers and dot-product auxiliaries, and (ii) degenerating a square matrix-multiplication tensor into any full direct sum of polynomial multiplications, even with mixed leg orientations and arbitrary changes of basis. An explicit degeneration to the small Coppersmith–Winograd tensor gives a standard laser inequality that also survives at every size.

I found **no valid inequality with a power-law violation of the proposed assignment**. This does not prove a barrier for arbitrary laser extractions into multiple matrix products, arbitrary basis changes in such extractions, or arbitrary mixed `T ⊗ C` constructions. The precise enlarged barrier theorem is in Section 7. It is a feasible assignment to a specified system, not a character on the actual tensor semiring.

All finite computations use rational/integer arithmetic. Section 8 contains the executable verification code. Only this project file was written.

## 1. Notation and useful exact bounds

Put

\[
 F(a,b)=P^*(a,b)=\phi(u)\left(v+\frac{u-1}{2}\right),
 \quad u=\min(a,b),\quad v=\max(a,b),\qquad
 \phi(n)=\prod_{j=1}^{n-1}\frac{3j+1}{3j}.
\]

The assignment being tested is

\[
 t=p_X=p_Y=p_Z=\frac34,\qquad
 L(C(a,b)^\pi)=F(a,b)^{3/4},\quad
 L(B_L(m))=m^{3/4},\quad L(T_m)=m^{9/4}.
 \tag{1}
\]

For rectangular multiplication, use

\[
 M(e,h,l)=\sum_{i=1}^{e}\sum_{j=1}^{h}\sum_{k=1}^{l}
 x_{ij}y_{jk}z_{ki}
 \cong B_Y(e)\otimes B_Z(h)\otimes B_X(l).
\]

Thus the actual character value is

\[
 \lambda(M(e,h,l))=e^{p_Y}h^{p_Z}l^{p_X},
 \qquad L(M(e,h,l))=(ehl)^{3/4}.
 \tag{2}
\]

The orientation of these exponents matters before specializing to (1).

### 1.1 Two elementary induction bounds

For every integer \(n\geq1\),

\[
 \boxed{n\leq\phi(n)^3\leq\frac{3n-1}{2}.}
 \tag{3}
\]

Both bounds are equalities at 1. The induction steps are certified by

\[
 n(3n+1)^3-27n^3(n+1)=9n^2+n>0,
\]

\[
 27n^3(3n+2)-(3n+1)^3(3n-1)=6n+1>0.
\]

In particular, \(F(n,n)\geq n^{4/3}\), while
\(F(n,n)^{3/4}\leq(3n-1)/2\). These suffice for the power-law comparisons here; no numerical Gamma-function constant is needed.

We also have

\[
 \boxed{F(a,b)\geq a+b-1.}
 \tag{4}
\]

For \(a=1\leq b\), this is equality. For \(2\leq a\leq b\), use \(\phi(a)\geq4/3\): the resulting difference from \(a+b-1\) is \((b-a+1)/3>0\).

### 1.2 A sandwich by flattening dimensions

Let \(d_X=a,d_Y=b,d_Z=a+b-1\) be the flattening ranks of \(C(a,b)\). Then

\[
 \boxed{(d_Xd_Yd_Z)^{1/3}\leq F(a,b)^{3/4}
 \leq\frac{3(a+b)-2}{4}
 =\frac38(d_X+d_Y+d_Z)-\frac18.}
 \tag{5}
\]

Here is a proof that uses only (3) and arithmetic–geometric mean inequalities. Assume \(a\leq b\) and write \(v=b+(a-1)/2\). For the lower bound,

\[
 F(a,b)^9\geq a^3v^9
 \geq a^4b^4(a+b-1)^4,
\]

because \(v\geq a\) and \(v^2\geq b(a+b-1)\). For the upper bound,

\[
 F(a,b)^{3/4}\leq
 \left(\frac{3a-1}{2}v^3\right)^{1/4}
 \leq\frac{(3a-1)/2+3v}{4}
 =\frac{3(a+b)-2}{4}.
\]

Since \((3(a+b)-2)/4\leq a+b-1\), this also gives a shorter all-size proof of the rank constraint \(F(a,b)^{3/4}\leq a+b-1\).

### 1.3 A stretching lemma

If \(c,d\geq r\geq1\) and \(u\geq0\) is an integer, then

\[
 \frac{F(c,d+u)}{F(c,d)}\leq1+\frac{2u}{3r-1}.
 \tag{6}
\]

To see this, use the separate concavity proved in round 1 and bound the first forward increment. If \(d\geq c\), the relative increment is
\(1/(d+(c-1)/2)\leq2/(3r-1)\). If \(d<c\), the increment is

\[
 F(c,d+1)-F(c,d)=\phi(d)\frac{c+2d}{3d}.
\]

Its ratio to \(F(c,d)\) is at most \(2/(3d-1)\), since the cross-multiplied difference is

\[
 6d\left(c+\frac{d-1}{2}\right)-(c+2d)(3d-1)
 =(3d+1)(c-d)\geq0.
\]

Concavity bounds every subsequent increment by this first one.

Consequently, if \(r\leq b,c,d\),

\[
 \boxed{F(r,b)F(c,d)\geq F(r,r)F(c,b+d-r).}
 \tag{7}
\]

Indeed, \(F(r,b)/F(r,r)=1+2(b-r)/(3r-1)\), and (6) applies with \(u=b-r\). This is the useful four-parameter reduction.

## 2. Task 1: all-size proof for (C)

Recall (C): with \(A_i=a+c-1-2i\) and \(B_j=b+d-1-2j\),

\[
 F(a,b)F(c,d)\geq
 \sum_{\substack{0\leq i<\min(a,c)\\0\leq j<\min(b,d)\\i+j=k}}
 F(A_i,B_j).
 \tag{C}
\]

Exchanging the two factors and/or simultaneously exchanging their inputs preserves this inequality, so assume \(a=r=\min(a,b,c,d)\). Put \(y=b+d-r\geq r\).

Each allowed \(i\) occurs at most once, and \(B_j\leq b+d-1\). Monotonicity, padding by the omitted positive terms, and separate concavity give

\[
 \begin{aligned}
 \text{RHS(C)}
 &\leq\sum_{i=0}^{r-1}F(c+r-1-2i,b+d-1)\\
 &\leq rF(c,b+d-1).
 \end{aligned}
 \tag{8}
\]

The mean of the first coordinates is exactly \(c\). Extend \(F(c,0)=0\), as in round 1. Concavity then implies

\[
 F(c,y+r-1)\leq\frac{y+r-1}{y}F(c,y)
 \leq\frac{2r-1}{r}F(c,y).
\]

Finally, (7) and (4) on the diagonal imply

\[
 F(r,b)F(c,d)
 \geq F(r,r)F(c,y)
 \geq(2r-1)F(c,y)
 \geq\text{RHS(C)}.
\]

This proves (C) for every size and every allowed \(k\), including \(r=1\). It does not depend on a large-size approximation.

## 3. Task 1: all-size proof for (D)

Recall (D), at \(t=3/4\):

\[
 [F(a,b)F(c,d)]^{3/4}\geq
 \sum_{i=0}^{r-1}F(a+c-1-2i,b+d-1-2i)^{3/4},
 \quad r=\min(a,b,c,d).
 \tag{D}
\]

Again assume \(a=r\), and set \(y=b+d-r\). Write

\[
 \delta_i=r-1-2i,\qquad
 S=\sum_{i=0}^{r-1}F(c+\delta_i,y+\delta_i),\qquad G=F(c,y).
\]

We will prove the stronger rational inequality

\[
 \boxed{[F(r,b)F(c,d)]^3\geq rS^3.}
 \tag{9}
\]

Hölder then gives
\(\sum_i F_i^{3/4}\leq r^{1/4}S^{3/4}\), which proves (D).

### 3.1 A rational bound on the diagonal sum

The increments of \(\phi\) decrease, because

\[
 \phi(n+1)-\phi(n)=\frac{\phi(n)}{3n},\qquad
 \frac{\phi(n+1)/(3(n+1))}{\phi(n)/(3n)}
 =\frac{3n+1}{3n+3}<1.
\]

Thus its forward slope is a supporting slope at the integer \(c\):

\[
 \phi(c+\delta)\leq\phi(c)\left(1+\frac{\delta}{3c}\right)
 \quad(c+\delta\geq1).
 \tag{10}
\]

By symmetry assume \(c\leq y\); both are at least \(r\). Put \(v=y+(c-1)/2\). The ordering stays the same after adding \(\delta_i\), so (10) gives

\[
 F(c+\delta_i,y+\delta_i)
 \leq\phi(c)\left(1+\frac{\delta_i}{3c}\right)
                  \left(v+\frac32\delta_i\right).
\]

The second factor is positive even at the smallest index. Since

\[
 \sum_i\delta_i=0,\qquad
 \sum_i\delta_i^2=\frac{r(r^2-1)}3,
\]

we obtain

\[
 \begin{aligned}
 S&\leq rG\left(1+\frac{r^2-1}{6cv}\right)\\
  &\leq rG\left(1+\frac{r^2-1}{3r(3r-1)}\right)
  \leq\frac98rG.
 \end{aligned}
 \tag{11}
\]

The last step is the exact inequality
\(8(r^2-1)\leq3r(3r-1)\), equivalent to
\(r^2-3r+8>0\).

### 3.2 Finish

For \(r\geq2\), (3) gives

\[
 F(r,r)=\phi(r)\frac{3r-1}{2}
 \geq\frac54r^{4/3}.
\]

Combine this with (7) and (11):

\[
 [F(r,b)F(c,d)]^3
 \geq\frac{125}{64}r^4G^3
 >\frac{729}{512}r^4G^3
 \geq rS^3.
\]

The strict rational comparison is \(1000>729\). For \(r=1\), there is one summand and (9) follows directly from (7). This closes the all-size gap for (D).

## 4. Task 2(a): tightness, cyclic groups, and a concrete laser extraction

### 4.1 Tightness and scalar extraction

The weights \(i,j,-k\) show that \(C(a,b)\) is tight. Tightness is a property of its support; it does not itself specify the sizes or multiplicities of matrix blocks in a laser extraction.

In fact the ordinary subrank of \(C(a,b)\) is exactly \(\min(a,b)\). For the lower bound, evaluate both input polynomials at that many distinct points and retain those evaluations of the product. The two evaluation maps have full row rank, and suitable input bases identify the result with that many independent scalar products. The opposite bound is the smaller input flattening rank. Thus scalar extraction from powers of this single tensor cannot beat the minimum-leg dimension, and (1) survives because

\[
 F(a,b)^{3/4}\geq F(r,r)^{3/4}\geq r,
 \qquad r=\min(a,b).
\]

There is a more general scalar-extraction barrier. Let \(S\) be any tensor product of permuted \(C\)'s, oriented \(B\)'s, and rectangular matrix products, with arbitrary sizes. For each generator, its assigned value is at least the geometric mean of its three flattening ranks: this is (5) for \(C\), \(m^{3/4}\geq m^{2/3}\) for \(B\), and \((ehl)^{3/4}\geq(ehl)^{2/3}\) for a matrix product. Flattening ranks multiply. Consequently, for **any** restriction or degeneration \(S\rightsquigarrow\Delta_N\),

\[
 N\leq\min_L\operatorname{rank}_L(S)
 \leq\left(\prod_L\operatorname{rank}_L(S)\right)^{1/3}
 \leq L(S).
 \tag{12}
\]

This includes arbitrary changes of basis and arbitrary tensor powers. It is a barrier for scalar targets, not for matrix targets valued above their flattening ranks.

### 4.2 The wrap-around relation

Let \(G_n\) denote multiplication in \(\mathbb C[\mathbb Z_n]\). Reducing output exponents modulo \(n\) is a restriction

\[
 C(n,n)\geq G_n.
\]

The complex Fourier transform identifies \(G_n\) with the diagonal tensor \(\Delta_n\), so **every** character has \(\lambda(G_n)=n\). The resulting inequality is

\[
 P(n,n)^t\geq n.
 \tag{13}
\]

It survives (1) by (3). The group-algebra value is not \(n^t\). Conversely, the no-wrap coordinate embedding is
\(C(a,b)\leq G_{a+b-1}\), and supplies the familiar rank upper bound. These are different group sizes and different restriction directions.

Fourier diagonalization is a genuine linear change of basis. One cannot subsequently claim a barrier for arbitrary extractions from \(G_n\) just from a barrier about its group-coordinate support. General extractions from \(\Delta_n\) include general bilinear algorithms.

### 4.3 An explicit degeneration to Coppersmith–Winograd

For \(n\geq3\), put \(q=n-2\). First restrict the output of \(C(n,n)\) to degrees \(0,\ldots,n-1\), giving the truncated multiplication tensor

\[
 t_n=\sum_{i+j<n}x_i y_j z_{i+j}.
\]

Define \(w(0)=0\), \(w(i)=1\) for \(1\leq i\leq n-2\), and \(w(n-1)=2\). Give \(x_i,y_i\) weight \(w(i)\), and \(z_k\) weight \(-w(k)\). All supported weights are nonnegative. The leading tensor is exactly

\[
 x_0y_0z_0+
 \sum_{i=1}^{n-1}(x_0y_i+x_i y_0)z_i+
 \sum_{i=1}^{n-2}x_i y_{n-1-i}z_{n-1}.
 \tag{14}
\]

This is the big \(CW_q\): the middle multiplication pairing is the nondegenerate symmetric reversal form. An exact congruence over \(\mathbb C\), accompanied by the dual output change, puts it in the usual identity-pairing basis. Reversing the output grading gives the usual symmetric trilinear convention. Setting \(x_{n-1}=y_{n-1}=z_0=0\) leaves the small tensor

\[
 cw_q=\sum_{i=1}^q(x_0y_i z_i+x_i y_0z_i+x_i y_i z_0)
\]

up to that change of basis. Thus

\[
 C(n,n)\rightsquigarrow CW_{n-2}\geq cw_{n-2}.
 \tag{15}
\]

The displayed construction distinguishes full polynomial multiplication, truncated multiplication, and cyclic convolution. Their dimensions and character values cannot be interchanged.

### 4.4 The standard small-CW laser inequality, with its cost

The three small-CW blocks are \(B_X(q),B_Y(q),B_Z(q)\). At frequency \((1/3,1/3,1/3)\), the usual tight-support laser extraction from powers gives, for \(k\) divisible by 3, a full direct sum of matrix blocks \(T_{q^{k/3}}\), with the number of blocks having exponential rate

\[
 \alpha=e^{H(1/3,2/3)}=\left(\frac{27}{4}\right)^{1/3}.
\]

More precisely, the count is at least \(\alpha^{k-o(k)}\); this is an asymptotic extraction statement, not a claim that \(\alpha^k\) is an integer finite multiplicity. Applying additivity, multiplicativity, and degeneration monotonicity before taking \(k\)-th roots gives

\[
 \lambda(cw_q)\geq\alpha q^t.
\]

This is the standard small-CW calculation; see the displayed small-CW formula in [Landsberg, *Secant varieties and the complexity of matrix multiplication*, Section 4](https://rendiconti.dmi.units.it/volumi/54/011.pdf). Its block-count argument applies to characters, rather than only to the rank upper bound usually used to state the algorithmic consequence.

The small tensor is invariant under leg permutations. Hence (15), applied to all six characters, yields the valid new profile inequality

\[
 \boxed{P(n,n)^t\geq\left(\frac{27}{4}\right)^{1/3}(n-2)^t.}
 \tag{16}
\]

At \(t=3/4\), this survives **for every \(n\geq3\)**. Indeed \(\alpha<2\), while weighted AM–GM gives

\[
 n^4=((n-2)+2)^4\geq\frac{512}{27}(n-2)^3
 >16(n-2)^3.
\]

Together with \(F(n,n)^{3/4}\geq n\), this proves
\(F(n,n)^{3/4}>\alpha(n-2)^{3/4}\). An entirely rational certificate for each parameter is

\[
 F(n,n)^9\geq\left(\frac{27}{4}\right)^4(n-2)^9.
\]

For the proposed scaling, the left side of (16) has order \(n\), while the right side has order \(n^{3/4}\). Thus this route is unbalanced in the direction that makes it weaker. Repeating this extraction or tensoring both sides with auxiliary \(B\)'s cannot reverse its already-proved all-size inequality.

This calculation does not exhaust higher-power repartitioning of \(CW_q\), \(t_n\), or \(C(n,n)\).

## 5. A structural barrier for one matrix block in convolution coordinates

Here “monomial” means deletion, rescaling, and injective relabeling of coordinate variables, followed by integer-weight degeneration. It excludes general linear combinations, variable identifications, and changes to a different basis.

### 5.1 Abelian-group lemma

**Lemma.** If the standard multiplication tensor of a finite abelian group \(G\) monomially degenerates to one rectangular matrix product \(M(e,h,l)\), then

\[
 \boxed{ehl\leq |G|.}
 \tag{17}
\]

Use the symmetric support convention \(x+y+z=0\) in \(G\). Write the selected labels of the matrix variables as \(X_{ij},Y_{jk},Z_{ki}\). The identities
\(X_{ij}+Y_{jk}+Z_{ki}=0\) for all \(i,j,k\) imply that there are group elements \(U_i,V_j,W_k\) such that

\[
 X_{ij}=U_i+V_j,\quad Y_{jk}=-V_j+W_k,\quad Z_{ki}=-U_i-W_k.
 \tag{18}
\]

For example, fixing a reference index in the identities shows that each difference of two rows of \(X\) is constant along the row; the other two formulas then follow. The same argument applies to the integer variable weights, since every desired matrix term has weight zero. Denote the corresponding weight potentials by \(\widehat U_i,\widehat V_j,\widehat W_k\).

Suppose

\[
 U_i+V_j+W_k=U_{i'}+V_{j'}+W_{k'}.
\]

Then both coordinate triples

\[
 (X_{ij},Y_{j'k},Z_{k'i'}),\qquad
 (X_{i'j'},Y_{jk'},Z_{ki})
\]

are supported in the group tensor. Their weights are opposite numbers, namely
\(\pm(\widehat U_i+\widehat V_j+\widehat W_k-\widehat U_{i'}-\widehat V_{j'}-\widehat W_{k'})\).
Both must be nonnegative; consequently both are zero and must belong to the surviving matrix tensor. Matching its three matrix indices forces \(i=i',j=j',k=k'\). The map \((i,j,k)\mapsto U_i+V_j+W_k\) is therefore injective, proving (17).

For restrictions this is the familiar abelian triple-product obstruction. The opposite-weight argument is why the stated version also covers monomial degenerations. The group-theoretic context is described by [Cohn–Kleinberg–Szegedy–Umans, *Group-theoretic Algorithms for Matrix Multiplication*](https://people.csail.mit.edu/virgi/6.890/papers/groupmatrix.pdf); the proof above is the exact argument used here.

### 5.2 Application to arbitrary products with auxiliaries

Let

\[
 S=\bigotimes_{s=1}^k C(a_s,b_s)^{\pi_s}
 \otimes\bigotimes_v B_{L_v}(m_v),
 \qquad
 N=\prod_s(a_s+b_s-1)\prod_v m_v.
\]

Each convolution is a coordinate restriction of the cyclic group tensor of order \(a_s+b_s-1\). Each oriented dot product is a coordinate restriction of a cyclic group tensor of order \(m_v\). Their tensor product embeds in the tensor of the direct-product abelian group of order \(N\). Permutations of the three legs are harmless in the symmetric group-support convention. Matrix-product auxiliaries are included by their factorization into three dot products.

Therefore, if \(S\) monomially degenerates to **one** \(M(e,h,l)\) in these coordinates, (17) and (4) give

\[
 L(M(e,h,l))=(ehl)^{3/4}
 \leq N^{3/4}
 \leq\left(\prod_sF(a_s,b_s)\prod_v m_v\right)^{3/4}
 =L(S).
 \tag{19}
\]

This proves survival for all sizes, all powers \(k\), rectangular targets, all leg permutations, and all such auxiliaries. The actual character inequality being tested is

\[
 \prod_s\lambda(C(a_s,b_s)^{\pi_s})\prod_v m_v^{p_{L_v}}
 \geq e^{p_Y}h^{p_Z}l^{p_X}.
\]

The theorem is stronger than a degree-counting observation: the full constants satisfy the inequality.

**Scope limitation:** the collision proof is internal to one matrix block. It does not prove \(\sum_s e_sh_sl_s\leq |G|\) for multiple blocks. Do not use that unproved extension. Full direct sums and subsequent combinations are precisely a remaining laser-method possibility. Likewise, Fourier basis changes take the problem outside (19).

## 6. Task 2(b): upper-bound relations from matrix multiplication

### 6.1 Arbitrary full sums of polynomial multiplications: an all-size barrier

**Theorem.** Every restriction or degeneration

\[
 T_m\rightsquigarrow\bigoplus_{s=1}^N C(a_s,b_s)^{\pi_s}
 \tag{20}
\]

satisfies the proposed character assignment:

\[
 \boxed{m^{9/4}\geq\sum_s F(a_s,b_s)^{3/4}.}
 \tag{21}
\]

This statement permits arbitrary linear substitutions and arbitrary mixtures of leg permutations in the summands.

For each of the three legs, flattening monotonicity and direct-sum additivity give
\(\sum_s d_{s,L}\leq m^2\). Apply (5), summed over the blocks:

\[
 \sum_sF(a_s,b_s)^{3/4}
 \leq\frac38\sum_{s,L}d_{s,L}-\frac N8
 \leq\frac98m^2-\frac N8.
 \tag{22}
\]

For \(m\geq2\), \(m^{1/4}>9/8\), as witnessed exactly by
\(2\cdot8^4>9^4\). This proves (21). For \(m=1\), flattening ranks permit at most one nonzero scalar block, so (21) also holds. Zero output is immediate.

The nonsymmetrized character relation for (20) is

\[
 m^{3t}\geq\sum_s\lambda(C(a_s,b_s)^{\pi_s}).
\]

Under (1) it is exactly (21); under general asymmetry it must retain the individual character values. The theorem shows that optimizing which such degeneration to construct cannot exclude the balanced assignment.

If all summands have the same output leg, there is the simpler stronger bound
\(\sum_sF(a_s,b_s)^{3/4}\leq\sum_s(a_s+b_s-1)\leq m^2\).

In scaling terms, all these target sums are at most \((9/8)m^2\); the proposed source value grows as \(m^{9/4}\). This entire upper-bound route is unbalanced in a harmless direction.

### 6.2 What remains possible with a `T ⊗ C` source

For an arbitrary degeneration

\[
 T_m\otimes C(a,b)\rightsquigarrow\bigoplus_s C(A_s,B_s)^{\pi_s},
\]

the same calculation gives

\[
 \sum_sF(A_s,B_s)^{3/4}
 \leq\frac38m^2(2a+2b-1).
 \tag{23}
\]

Thus **all** such relations survive whenever

\[
 \boxed{mF(a,b)^3\geq
 \left[\frac38(2a+2b-1)\right]^4.}
 \tag{24}
\]

This is a sufficient condition, not a necessary one. For \(a=b=n\), it holds for **every \(m\geq2\)**. For \(n\geq4\), (3) gives

\[
 2F(n,n)^3\geq2n\left(\frac{3n-1}{2}\right)^3
 \geq\frac{1331}{256}n^4
 >\frac{1296}{256}n^4
 =\left(\frac32n\right)^4
 >\left[\frac38(4n-1)\right]^4.
\]

For \(n=1,2,3\), substitute \(F(n,n)=1,10/3,56/9\); all three comparisons also hold exactly. Increasing \(m\) preserves (24). Hence a balanced polynomial factor together with any nontrivial square matrix factor cannot contradict (1) through targets of the form in (23).

This proof does not settle the pure-convolution case \(m=1\), highly unequal \(a,b\) outside (24), or targets that are products of multiple \(C\)'s or contain matrix factors. These are explicit remaining cases, rather than covered by the theorem through an appeal to scaling.

### 6.3 Simple mixed constructions with all costs included

There are several valid relations, but they add no constraint at (1):

* \(M(a,1,b)\geq C(a,b)\), by sending the outer-product output \(z_{ki}\) to \(z_{i+k}\). The character inequality is \(\lambda(C(a,b))\leq a^{p_Y}b^{p_X}\), already covered by round 1.
* Tensoring gives \(M(ma,m,mb)\geq T_m\otimes C(a,b)\); canceling the \(T_m\) value returns exactly that same outer-product bound.
* Interpolation gives \((a+b-1)T_m\geq T_m\otimes C(a,b)\); canceling returns the original rank bound.
* Restricting the matrix indices to its diagonal gives \(T_m\otimes C(a,b)\geq\bigoplus_{i=1}^m C(a,b)\). The value comparison is \(m^{3t}\geq m\), which (1) satisfies.

These relations are actual restrictions. None supplies an uncanceled matrix exponent capable of contradicting the assignment. The broader result (20)–(22) also rules out improving the last type merely by finding more full convolution blocks inside \(T_m\).

## 7. Task 3: the enlarged barrier, and why it is not a global character

**Enlarged barrier theorem.** The balanced assignment (1) satisfies simultaneously:

1. The numerical system in the restricted barrier theorem of round 1: the paper's axioms, all nonsymmetrized shared-leg entropy forms and their permutations, all higher-dot and transfer filtrations (A),(B), the outer-product/dot-product bounds, and the explicitly costed Fourier family.
2. Every size and every parameter of (C) and (D), including their nonsymmetrized character versions. At (1), the independently optimized shared-leg entropy form of (C) is exactly (C) to the power \(3/4\); the full direct-sum form of (D) is exactly (D).
3. Wrap-around scalar extraction (13), and the standard small-CW laser inequalities (16).
4. Every scalar-target relation (12), from any tensor product of the listed generators.
5. Every single-matrix-target monomial relation (19) in the specified convolution/group coordinates, with arbitrary powers and dot-product or matrix auxiliaries.
6. Every matrix-to-full-convolution-sum relation (20), and the mixed-source relations covered by (24).

It continues to satisfy every consequence obtained by addition, multiplication, tensoring these verified relations, and combining their numerical inequalities. In particular, optimizing probabilities separately on different leg permutations or mixing these constraints in an LP cannot force \(t<3/4\).

The proof is the round-one theorem together with Sections 1–6 above. Unlike round 1, (C),(D) are now part of the all-size statement.

One may formalize the additive/multiplicative part as a positive homomorphism on the **presented semiring of formal generator expressions**, ordered by the verified relations and their additive/multiplicative closure. Its generator values are (1), and rectangular matrix factors obey (2). This is a precise formal consistency statement. It does not establish that evaluation is well-defined modulo every actual tensor isomorphism, mutual restriction, or degeneration.

There is a substantial reason not to assert the latter casually. An actual character with \(\lambda(T_m)=m^{9/4}\) would imply

\[
 \operatorname{rank}(T_m)\geq m^{9/4}\quad\text{for every }m,
 \qquad \omega\geq\frac94.
\]

Together with the supplied paper's claimed upper bound, it would establish equality. Even a character on the actual ordered subsemiring generated by the dots and scalars, with all their real restriction relations retained, would already have this consequence, because that subsemiring contains every \(T_m\). The current feasibility proofs do not establish this matrix-multiplication lower bound.

The remaining directions are therefore specific:

* Full direct sums of multiple matrix products extracted from powers of \(C\), especially after useful non-monomial changes of basis. The single-block theorem cannot be summed across blocks.
* Higher-power repartitioning that is not merely the standard three-sector small-CW extraction.
* Mixed `T ⊗ C` relations outside (24), or with more complicated targets than full sums of individual convolutions.
* Further relations involving asymmetry that actually rule out the balanced lift, rather than assuming that a character must be asymmetric.

Finally, the two asymptotic parameters must be kept distinct. A comparison can be degree-balanced as input length \(n\) grows while a fixed-size constant violation, if one exists, amplifies exponentially with tensor-power length \(k\). The results here prove satisfaction with constants for the families included in the barrier. They do not dismiss untested families merely because their dimensional exponents match.

For the distinction between unrestricted and monomial tensor comparisons, and between a barrier for a method and a conclusion about \(\omega\), see [Christandl–Vrana–Zuiddam, *Barriers for Fast Matrix Multiplication from Irreversibility*, Sections 2–3](https://theoryofcomputing.org/articles/v017a002/v017a002.pdf). No numerical barrier from that paper is being used as a substitute for the explicit proofs here.

## 8. Reproducible exact verification

These checks supplement the all-size proofs; they are not their justification. The script needs only Python's standard library and can be run through `/usr/bin/python3` without creating a project-side script or cache.

The (D) check tests the stronger rational certificate (9), so it requires no fourth-root approximations. The CW check verifies every original supported coefficient's weight and the complete surviving support, before any basis normalization.

```python
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product

@lru_cache(None)
def phi(a):
    assert a >= 1
    return Q(1) if a == 1 else phi(a-1)*Q(3*a-2, 3*a-3)

@lru_cache(None)
def F(a, b):
    a, b = sorted((a, b))
    return phi(a)*(b+Q(a-1, 2))

counts = dict(union=0, C=0, D_sufficient=0,
              phi_tangent=0, cw_uniform=0, sandwich=0, cw_support=0,
              mixed_square=0)

for a, b, c, d in product(range(1, 17), repeat=4):
    lhs = F(a, b)*F(c, d)
    for k in range(min(a, c)+min(b, d)-1):
        rhs = sum(F(a+c-1-2*i, b+d-1-2*(k-i))
                  for i in range(min(a, c))
                  if 0 <= k-i < min(b, d))
        assert lhs >= rhs, ('C', a, b, c, d, k)
        counts['C'] += 1
    r = min(a, b, c, d)
    s = sum(F(a+c-1-2*i, b+d-1-2*i) for i in range(r))
    assert lhs**3 >= r*s**3, ('D certificate', a, b, c, d)
    counts['D_sufficient'] += 1
    if a == r:
        assert lhs >= F(a, a)*F(c, b+d-a)
        counts['union'] += 1

for a in range(1, 33):
    assert phi(a)**3 >= a
    assert F(a, a) >= 2*a-1
    for z in range(1, 64):
        assert phi(z) <= phi(a)*(1+Q(z-a, 3*a))
        counts['phi_tangent'] += 1

for a, b in product(range(1, 101), repeat=2):
    r = min(a, b)
    assert r <= phi(r)**3 <= Q(3*r-1, 2)
    assert F(a, b)**3 <= Q(3*(a+b)-2, 4)**4
    assert F(a, b)**9 >= (a*b*(a+b-1))**4
    counts['sandwich'] += 1

for n in range(3, 65):
    assert F(n, n)**9 >= Q(27, 4)**4*(n-2)**9
    counts['cw_uniform'] += 1
    def w(i):
        return 0 if i == 0 else 2 if i == n-1 else 1
    original = {(i, j, i+j) for i in range(n) for j in range(n)
                if i+j < n}
    assert all(w(i)+w(j)-w(k) >= 0 for i, j, k in original)
    leading = {(i, j, k) for i, j, k in original
               if w(i)+w(j)-w(k) == 0}
    expected = ({(0, j, j) for j in range(n)}
                | {(i, 0, i) for i in range(n)}
                | {(i, n-1-i, n-1) for i in range(1, n-1)})
    assert leading == expected and len(leading) == 3*n-3
    small = {(i, j, k) for i, j, k in leading
             if i < n-1 and j < n-1 and k != 0}
    assert len(small) == 3*(n-2)
    counts['cw_support'] += 1

for n in range(1, 101):
    assert 2*F(n, n)**3 >= (Q(3, 8)*(4*n-1))**4
    counts['mixed_square'] += 1

# Polynomial identities: coefficient-vector operations, not sampling.
def mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def sub(a, b):
    return [(a[i] if i < len(a) else 0)
            - (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]

cube = mul(mul([1, 3], [1, 3]), [1, 3])
assert sub(mul([0, 1], cube), mul([0, 0, 0, 27], [1, 1])) == [0, 1, 9, 0, 0]
assert sub(mul([0, 0, 0, 27], [2, 3]), mul(cube, [-1, 3])) == [1, 6, 0, 0, 0]
assert 1000 > 729
assert 2*8**4 > 9**4
assert Q(512, 27) > 16
assert Q(27, 4) < 8
assert 1331 > 1296

assert counts == dict(union=18496, C=700416, D_sufficient=65536,
                     phi_tangent=2016, cw_uniform=62,
                     sandwich=10000, cw_support=62, mixed_square=100)
for name, count in counts.items():
    print(f'{name}: {count} exact checks, PASS')
print('Polynomial identities and rational constants: PASS')
print('No floating-point arithmetic or unresolved comparisons.')
```

Observed output:

```text
union: 18496 exact checks, PASS
C: 700416 exact checks, PASS
D_sufficient: 65536 exact checks, PASS
phi_tangent: 2016 exact checks, PASS
cw_uniform: 62 exact checks, PASS
sandwich: 10000 exact checks, PASS
cw_support: 62 exact checks, PASS
mixed_square: 100 exact checks, PASS
Polynomial identities and rational constants: PASS
No floating-point arithmetic or unresolved comparisons.
```
