# R4 positive control at p=1009 (sampled k). Run from scripts/main/: /usr/bin/python3 spot1009.py  (~3 s)
# Expected output: expected_output/spot1009.txt
import random, may308 as M
p=1009; u,q=M.setup(p)
random.seed(1)
ks=[1,2,3,p-2,p-1,p,p+1,p+2,p+3]+random.sample(range(4,(u(4)-1)//(q-1)+1),150)
for k in ks:
    ms=M.monomials(p,q*k+1,u(4)+k)
    if ms or k in (1,2,p-1,p): print(k,[M.show(m) for m in ms])
