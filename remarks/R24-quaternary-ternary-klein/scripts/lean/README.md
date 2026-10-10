# Lean files for this remark

**Base.** [openai/math](https://github.com/openai/math) at commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, directory `lean/`. The toolchain is `leanprover/lean4:v4.34.1`
(`lean-toolchain`, copied here), and dependencies are as pinned in the upstream `lake-manifest.json`
(Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612`). No upstream file is modified, and Mathlib has no
local changes.

**Files added** (17 files, 4085 lines; checksums in `SHA256SUMS`):

| Directory | Files | Content |
|---|---|---|
| `OAI/LinearAlgebra/CirculantHadamard/GeneralNorm/` | 7 | Prop′ (ℤ[i], odd N, norm N), Master Lemma over ℤ[i], Theorems Q and T1 |
| `OAI/LinearAlgebra/CirculantHadamard/Eisenstein/` | 8 | ℤ[ω] copies of the cyclotomic and localization layers, Prop″, Master Lemma over ℤ[ω], Theorem T3 |
| `OAI/LinearAlgebra/CirculantHadamard/RankOne/` | 2 | Fourier inversion on A × B; Theorem F (ℤ[i], C_{pq}) |

They import only Mathlib and 15 unmodified upstream modules of `OAI/LinearAlgebra/CirculantHadamard/`
(AlternatingField, AlternatingLocalization, BinaryCharacterSum, BinaryLocalRing,
BinaryProductTransport, BinaryResidue, CyclicEvaluation, CyclicPolynomial, CyclotomicLocalBase,
CyclotomicRings, CyclotomicTorsion, FirstOrderResidue, Kronecker, LocalizationMaps, ResidueTorsion).

**Check.**
```bash
git clone https://github.com/openai/math && cd math && git checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
cp -r <this-dir>/OAI lean/ && cd lean && lake exe cache get   # Mathlib cache
lake build OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremQ \
           OAI.LinearAlgebra.CirculantHadamard.GeneralNorm.TheoremT1 \
           OAI.LinearAlgebra.CirculantHadamard.Eisenstein.TheoremT3 \
           OAI.LinearAlgebra.CirculantHadamard.RankOne.TheoremF
lake env lean <this-dir>/Checks.lean
```
Expected output: `expected_output/axioms.txt`. All nine `#print axioms` lines list only
`[propext, Classical.choice, Quot.sound]`.

**Difference from the built copy.** Two comment lines of `GeneralNorm/TheoremT1.lean` (l. 11 and the
docstring at l. 210) were reworded after the build, changing "new case" to "non-trivial case". No
code changed. All other files are byte-identical to the copies that were built.
