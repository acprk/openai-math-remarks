"""Extract the Appendix A listing of family 222 ('Microscopic jamming in the
negative spherical perceptron', build/main.tex) from a clone of openai/math,
check its SHA-256 against the digest printed in the paper (main.tex l.9496-9497),
and write it to ./222_original_certificate.py.

Usage: python3 extract_222_certificate.py OPENAI_MATH COMMIT
The listing is the text strictly between the unique \\begin{lstlisting} and
\\end{lstlisting} of build/main.tex, taken byte for byte (with final newline).
"""
import hashlib, re, subprocess, sys

repo, commit = sys.argv[1], sys.argv[2]
path = ("preprints/Microscopic-jamming-in-the-negative-spherical-perceptron-"
        "September-24-2026/build/main.tex")
tex = subprocess.run(["git", "-C", repo, "show", f"{commit}:{path}"],
                     check=True, capture_output=True).stdout.decode("utf-8")
lines = tex.split("\n")
b = [i for i, l in enumerate(lines) if l.strip() == r"\begin{lstlisting}"]
e = [i for i, l in enumerate(lines) if l.strip() == r"\end{lstlisting}"]
assert len(b) == 1 and len(e) == 1, (b, e)
code = "\n".join(lines[b[0] + 1:e[0]]) + "\n"
digest = hashlib.sha256(code.encode("utf-8")).hexdigest()
# digest printed in the paper, split by \allowbreak
m = re.search(r"([0-9a-f]{33})\\allowbreak\s*\n\s*([0-9a-f]{31})\.", tex)
published = m.group(1) + m.group(2)
print(f"listing lines {b[0]+2}-{e[0]} of main.tex ({len(code.splitlines())} lines)")
print("sha256 of listing :", digest)
print("sha256 in paper   :", published)
print("match:", digest == published)
open("222_original_certificate.py", "w").write(code)
sys.exit(0 if digest == published else 1)
