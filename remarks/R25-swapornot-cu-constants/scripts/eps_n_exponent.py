#!/usr/bin/env python3
"""
Family 238, main paper (Optimal-order mixing of the Thorp shuffle), Prop. 4.3 / eq. (h-8), main.tex l.381-389:
    eps_n = (n^{3/2}/2) * sqrt((31/32)^{198 d} + 32 e^{-n/256}),   n = 2^d,   stated "= O(n^{-5/2})".
Since (31/32)^{198d} = n^{-gamma} with gamma = 198 log2(32/31), the first term gives exactly n^{-(gamma/2 - 3/2)}.
This script certifies, in mpmath interval arithmetic:
  * gamma/2 - 3/2 = 99 log2(32/31) - 3/2 in [3.03456, 3.03457] (3.035 to three decimals);
  * for d >= 15 the e^{-n/256} term is below the first one, so  (1/2) n^{-e*} <= eps_n <= (sqrt2/2) n^{-e*};
  * the table eps_n * n^{5/2} and eps_n * n^{e*} for d = 15..40.
Run: python3 eps_n_exponent.py   (mpmath >= 1.3; < 1 s)
"""
from mpmath import iv
iv.dps = 40
g = 198 * iv.log(iv.mpf(32) / 31) / iv.log(2)
e = g / 2 - iv.mpf(3) / 2
print(f"198*log2(32/31) in [{float(g.a):.9f}, {float(g.b):.9f}]  (> 8: {g.a > 8})")
print(f"exponent e* = 99*log2(32/31) - 3/2 in [{float(e.a):.9f}, {float(e.b):.9f}]")
assert e.a > iv.mpf('3.03456').b and e.b < iv.mpf('3.03457').a
print("paper: O(n^{-5/2}); actual order of the stated eps_n: n^{-3.03456...} (Theta, not only O)")
ok = []
for d in range(1, 200):
    n = iv.mpf(2) ** d
    t1 = (iv.mpf(31) / 32) ** (198 * d); t2 = 32 * iv.exp(-n / 256)
    ok.append(t2.b <= t1.a)
dmin = next(d for d in range(1, 200) if all(ok[d - 1:]))
print(f"e^(-n/256) term <= (31/32)^(198d) term for every d >= {dmin} (checked up to d=199; beyond, n/256 grows faster than d)")
print(" d   eps_n*n^(5/2)        eps_n*n^(e*)")
for d in [15, 16, 20, 24, 32, 40]:
    n = iv.mpf(2) ** d
    eps = n ** iv.mpf(1.5) / 2 * iv.sqrt((iv.mpf(31) / 32) ** (198 * d) + 32 * iv.exp(-n / 256))
    r1 = eps * n ** iv.mpf(2.5); r2 = eps * n ** e
    print(f"{d:2d}   {float(r1.b):.6e}     [{float(r2.a):.6f}, {float(r2.b):.6f}]")
