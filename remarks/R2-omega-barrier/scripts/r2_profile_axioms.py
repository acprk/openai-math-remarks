#!/usr/bin/env python3
"""
R2 core check: the least profile P* of OpenAI's 9/4 argument.

Environment: Python >= 3.8, standard library only (fractions) for sections 1-5;
section 6 uses mpmath (interval arithmetic) for the optional Wendel check.
Run:  /usr/bin/python3 r2_profile_axioms.py [N]        (default N = 200)
Expected output: expected_output/r2_profile_axioms.txt  (runtime ~1 min, 1 thread)

All comparisons are exact rationals (fractions.Fraction); no floating point
is used for any PASS/FAIL decision.

Notation (OpenAI preprint, build/paper.tex at openai/math adc7f1241b42):
  phi(1) = 1, phi(m) = prod_{j=1}^{m-1} (1 + 1/(3j))
  P*(a,b) = phi(m) * (M + (m-1)/2),   m = min(a,b), M = max(a,b)
Axioms actually used by the paper (Lemma 5.1 hypotheses + final step):
  (Pos) P > 0                         paper.tex l.611
  (S)   P(a,b) = P(b,a)               l.427, l.611
  (B)   P(1,b) = b                    l.427, l.611
  (C)   2P(a,b) >= P(a,b+1)+P(a,b-1)  (a>=1, b>=2)   l.444-447, l.612
  (T)   P(a,3h+a-1) >= 3P(a,h)        (a,h>=1)       l.528-530, l.613
  (U_t) P(a,b) <= (a+b-1)^(1/t)       l.428 ; only the diagonal is used, l.657
At t = 3/4, (U) is checked as P*(a,b)^3 <= (a+b-1)^4 (exact integers/rationals).
"""
import sys
from fractions import Fraction as Fr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
BIG = 4 * N + 10

# phi via the recurrence phi(m+1) = phi(m) * (3m+1)/(3m)   (independent of other scripts)
PHI = [None, Fr(1)]
for m in range(1, BIG + 1):
    PHI.append(PHI[m] * Fr(3 * m + 1, 3 * m))

def Pstar(a, b):
    m, M = (a, b) if a <= b else (b, a)
    return PHI[m] * (M + Fr(m - 1, 2))

counts = {}
def ok(name, cond):
    counts.setdefault(name, [0, 0])
    counts[name][0] += 1
    if not cond:
        counts[name][1] += 1

# ---- 1. the five axioms, exhaustive on a finite box -----------------------
for a in range(1, N + 1):
    for b in range(1, N + 1):
        p = Pstar(a, b)
        ok("(Pos) P*>0", p > 0)
        ok("(S) symmetry", p == Pstar(b, a))
        ok("(U) P*^3 <= (a+b-1)^4", p ** 3 <= (a + b - 1) ** 4)
        ok("(U) equality only at (1,1)", (p ** 3 == (a + b - 1) ** 4) == (a == b == 1))
    for b in range(2, 4 * N):
        ok("(C) concavity in b", 2 * Pstar(a, b) >= Pstar(a, b + 1) + Pstar(a, b - 1))
        ok("monotone in b (not an axiom)", Pstar(a, b) >= Pstar(a, b - 1))
    for h in range(1, N + 1):
        ok("(T) tripling", Pstar(a, 3 * h + a - 1) >= 3 * Pstar(a, h))
        ok("(T) equality iff h >= a-1", (Pstar(a, 3 * h + a - 1) == 3 * Pstar(a, h)) == (h >= a - 1))
for b in range(1, BIG):
    ok("(B) P*(1,b)=b", Pstar(1, b) == b)

# ---- 2. tightness of OpenAI's Lemma 5.1 chain on P* -----------------------
# H_a := 2 P(a,a)/(3a-1); Delta_{a,h} := P(a,h+1)-P(a,h)
for a in range(1, N + 1):
    D = Pstar(a, a)
    H = 2 * D / (3 * a - 1)
    ok("H_a = phi(a) (eq. normalized-growth tight)", H == PHI[a])
    ok("Delta_{a,a} = H_a", Pstar(a, a + 1) - D == H)
    # 2D_a <= P(a,4a-1)-D_a <= (3a-1) Delta_{a,a}: both equalities
    ok("2D_a = P(a,4a-1)-D_a", 2 * D == Pstar(a, 4 * a - 1) - D)
    ok("P(a,4a-1)-D_a = (3a-1)Delta_{a,a}", Pstar(a, 4 * a - 1) - D == (3 * a - 1) * (Pstar(a, a + 1) - D))
    if a >= 2:
        # (diagonal-rec): D_a - D_{a-1} = Delta_{a-1,a-1} + Delta_{a,a-1} >= H_{a-1}+H_a, equality
        ok("eq. diagonal-rec equality", D - Pstar(a - 1, a - 1) == PHI[a - 1] + PHI[a])
    # the two final simplifications are strict for a >= 2 (constant loss only)
    if a >= 2:
        ok("strict: phi(a)^3 > a", PHI[a] ** 3 > a)

# ---- 3. the elementary bounds used in the proof of (U) -------------------
#   m <= phi(m)^3 <= 2m-1   (induction; no Gamma function needed)
for m in range(1, BIG):
    ok("phi(m)^3 <= 2m-1", PHI[m] ** 3 <= 2 * m - 1)
    ok("phi(m)^3 >= m", PHI[m] ** 3 >= m)
    # induction-step certificates (polynomial identities checked pointwise):
    #   (2m+1)(3m)^3 - (2m-1)(3m+1)^3 = 9m^2 + 7m + 1      (upper bound step)
    #   m(3m+1)^3 - (m+1)(3m)^3       = 9m^2 + m           (lower bound step)
    ok("upper step identity = 9m^2+7m+1",
       (2 * m + 1) * (3 * m) ** 3 - (2 * m - 1) * (3 * m + 1) ** 3 == 9 * m * m + 7 * m + 1)
    ok("lower step identity = 9m^2+m",
       m * (3 * m + 1) ** 3 - (m + 1) * (3 * m) ** 3 == 9 * m * m + m)

# ---- 4. least element: the closed forms in the minimality proof ----------
# Row recursion forced by (C)+(T): P(a,h+1) >= P(a,h)*(2h+a+1)/(2h+a-1);
# on P* (b >= a) this is an equality; telescoping from P(a,a) gives P*(a,b).
for a in range(1, N + 1):
    for h in range(a, 3 * N):
        ok("row recursion tight on P* (h>=a)",
           Pstar(a, h + 1) == Pstar(a, h) * Fr(2 * h + a + 1, 2 * h + a - 1))
# (T) slack formula for h < a (from the note): L-R telescopes with increments
#   phi(a)(2+h/a) - 3phi(h) >= 0
for h in range(1, N):
    for a in range(h, N):
        inc = PHI[a] * (2 + Fr(h, a)) - 3 * PHI[h]
        ok("(T) h<a increment >= 0", inc >= 0)

# ---- 5. summary -----------------------------------------------------------
bad = 0
for name, (n, f) in counts.items():
    print(f"{name:55s} checks={n:8d} failures={f}")
    bad += f
print("D_1..D_6 =", [str(Pstar(a, a)) for a in range(1, 7)])

# ---- 6. optional: Wendel-type bound phi(m) <= m^(1/3)/Gamma(4/3) (interval) --
try:
    import mpmath
    iv = mpmath.iv
    iv.dps = 50
    g = iv.gamma(iv.mpf(4) / 3)
    v = 1 / g ** 3                  # interval containing Gamma(4/3)^(-3)
    sgn, man, ex, _ = v._mpi_[1]     # certified upper endpoint
    c3 = Fr(man) * Fr(2) ** ex       # exact rational >= Gamma(4/3)^-3
    assert sgn == 0
    w_bad = sum(1 for m in range(1, BIG) if PHI[m] ** 3 > c3 * m)
    print(f"Wendel check phi(m)^3 <= m/Gamma(4/3)^3 (m<{BIG}), failures={w_bad}; "
          f"Gamma(4/3)^-3 <= {float(c3):.12f}")
except ImportError:
    print("mpmath not available: Wendel check skipped (not needed by the proof)")

print("ALL PASS" if bad == 0 else f"FAILURES: {bad}")
sys.exit(1 if bad else 0)
