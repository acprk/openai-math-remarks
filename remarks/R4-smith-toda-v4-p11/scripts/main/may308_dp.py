# R4 cross-check by knapsack DP. Run from scripts/: /usr/bin/python3 main/may308_dp.py 7 11 13 17  (needs numpy; ~20 s)
# Expected output: expected_output/may308_dp.txt
#!/usr/bin/env python3
"""Independent cross-check of may308.py: counts May E1 monomials by a 2-D knapsack DP
(exterior generators 0/1, polynomial generators unbounded), exact Python ints, and
re-checks every vanishing / list-size claim at prime p."""
import sys, itertools
import numpy as np
def run(p):
    u = lambda i: (p**i-1)//(p-1); q = 2*p-2
    wmax = u(4) + (u(4)-1)//(q-1) + 2
    smax = q*((u(4)-1)//(q-1)) + 3
    H = [p**j*u(i) for i in range(1, 8) for j in range(8) if p**j*u(i) <= wmax]
    B = [p**(j+1)*u(i) for i in range(1, 8) for j in range(8) if p**(j+1)*u(i) <= wmax]
    dp = np.zeros((smax+1, wmax+1), dtype=object); dp[0, 0] = 1
    for c in H:   # exterior: 0/1, iterate s descending
        for s in range(smax, 0, -1):
            dp[s, c:] = dp[s, c:] + dp[s-1, :wmax+1-c]
    for c in B:   # polynomial: unbounded, s ascending
        for s in range(2, smax+1):
            dp[s, c:] = dp[s, c:] + dp[s-2, :wmax+1-c]
    cnt = lambda s, w: int(dp[s, w]) if 0 <= s <= smax and 0 <= w <= wmax else None
    bad = []
    for n in range(1, 4):
        k = 1
        while (q-1)*k <= u(n)-1:
            if cnt(q*k+1, u(n)+k): bad.append(('v%d' % n, k))
            k += 1
    for n in range(4):
        for e in itertools.product(range(3), repeat=n+1):
            rho = sum(e)
            if rho < 2: continue
            U = sum(ei*u(i) for i, ei in enumerate(e)); k = 1
            while (q-1)*k <= U+rho-1:
                if cnt(q*k+1-rho, U+k): bad.append(('prod', n, e, k))
                k += 1
    nz = {}
    k = 1
    while (q-1)*k <= u(4)-1:
        c = cnt(q*k+1, u(4)+k)
        if c: nz[k] = c
        k += 1
    W = 2*p*p+4*p+1
    small = [(j, d, cnt(6+q*j, W+d+j)) for d in (0, 1) for j in range(0 if d else 1, (W+d-6)//(q-1)+1) if cnt(6+q*j, W+d+j)]
    print(f"p={p}: lower-obstruction nonzero bidegrees: {bad[:6]}{'...' if len(bad)>6 else ''} ({len(bad)});"
          f" v4-target sizes {nz} (paper: {{1:2, 2:2, {p-1}:3, {p}:2}}); (5,W) size {cnt(5,W)}; nonzero small stems {small}")
for p in map(int, sys.argv[1:]): run(p)
