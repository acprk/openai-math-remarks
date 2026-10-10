"""Independent cross-check of the 'referenced but not shipped' items (does not use
check_upstream.sh).  For each path that a manuscript cites, resolve it relative to
the manuscript directory (and, for build-relative listings, to build/) and ask git
whether the object exists at COMMIT (`git cat-file -e COMMIT:path`).  Also checks
every historical commit reachable from COMMIT, to show the file was never shipped.
Usage: python3 indep_check_paths.py OPENAI_MATH COMMIT        (python3, git; seconds)
"""
import subprocess, sys
repo, commit = sys.argv[1], sys.argv[2]
P = "preprints/"
cited = {
 "266-A": (P + "The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026",
           ["verification/code/seed.cpp", "verification/code/covout.cpp", "verification/code/ub.cpp",
            "verification/code/patchx.cpp", "verification/code/atlas.cpp", "verification/code/stageB.cpp",
            "verification/code/stageC.cpp", "verification/code/extraction.json",
            "verification/code/proof_run.py", "verification/computation/README.md",
            "verification/REPRODUCE.md"]),
 "266-B": (P + "Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026",
           ["certificate/verify.py", "build/certificate/verify.py", "verification/verify.py"]),
 "095":   (P + "Hyperbolicity-Cones-Without-Semidefinite-Lifts-October-5-2026",
           ["verify_certificate.py", "build/verify_certificate.py", "build/source/verify_certificate.py"]),
 "318":   (P + "A-rational-obstruction-to-strong-chromatic-splitting-at-height-three-September-25-2026",
           ["support/r-ce/ce_certificate.py", "support/r-ce/ce_certificate.json",
            "support/r-ce/ce_certificate.txt", "support/historical/r-ce", "BUILDING.md",
            "SUPPORT-COVERAGE.json"]),
 "222":   (P + "Microscopic-jamming-in-the-negative-spherical-perceptron-September-24-2026",
           ["certificate/original_certificate.py", "certificate"]),
 "238":   (P + "Routing-densities-and-representation-contraction-for-Thorp-sweeps-September-26-2026",
           ["support/network_bounds.py", "support"]),
 "227":   (P + "Cutoff-throughout-the-high-temperature-Sherrington-Kirkpatrick-phase-September-24-2026",
           ["numerics-interval.py", "numerics-cross-errors.py",
            "build/numerics-interval.py", "build/appendices/numerics-interval.py"]),
}
def exists(rev, path):
    return subprocess.run(["git", "-C", repo, "cat-file", "-e", f"{rev}:{path}"],
                          capture_output=True).returncode == 0
revs = subprocess.run(["git", "-C", repo, "rev-list", commit], check=True,
                      capture_output=True, text=True).stdout.split()
print("commit", revs[0], "history length", len(revs))
for key, (d, files) in cited.items():
    assert exists(commit, d), d
    for f in files:
        now = exists(commit, f"{d}/{f}")
        ever = any(exists(r, f"{d}/{f}") for r in revs)
        print(f"{key:6s} {f:42s} present@commit={now}  present@any-ancestor={ever}")
