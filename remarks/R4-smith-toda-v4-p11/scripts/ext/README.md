# Actual Ext over the reduced-power algebra (layer 5)

`ext308.rs` is an example program for the `ext` crate of
[spectralsequences/sseq](https://github.com/spectralsequences/sseq). It computes a minimal resolution
of F_p over the reduced-power algebra P (Milnor basis, profile `q_part = 0`, i.e. all Bocksteins removed).
P is dual to F_p[t_1, t_2, …] with the additive coproduct, so Ext_P(F_p, F_p) is the additive cobar
cohomology of the manuscript. Lemma `fw:additive-comparison` identifies it with the R/I_4 cobar
cohomology at every v_4-target bidegree.

Usage: `ext308 p s_max w_max s1,w1 s2,w2 …`. Here the weight is w = t/(2p−2). The program prints dim Ext^{s,w}
for each requested bidegree, after six sanity values (h_{10}, h_{11}, b_{10}, b_{11} have dimension 1; h_{10}^2 = 0).

## Build

Use sseq commit `bf938a74d3344f050fb0a56615e67f6f20dc0da1`. The next commit, `b496e8ee`, bit-packs
Milnor basis elements and caps the degree at 2045, which is too small for p ≥ 7.

```
git clone https://github.com/spectralsequences/sseq && cd sseq && git checkout bf938a74d3344f050fb0a56615e67f6f20dc0da1
cp <this dir>/ext308.rs ext/examples/
cd ext && cargo +nightly build --release --features concurrent --example ext308
```

## Runs (8 threads)

```
ext308 5 42 162 9,157 17,158 25,159 33,160 41,161               # ~10 s       -> ../expected_output/ext308_p5.txt
ext308 7 86 408 13,401 25,402 37,403 49,404 61,405 73,406 85,407 # 4662 s, 9.3 GB -> ../expected_output/ext308_p7.txt
```
