# R14 primary script. Env: sage -python (SageMath 9.5), run from scripts/primary/.
# Expected: expected_output/primary_exact.txt
# exact certificates: q=16 Mf <= lam f (lam<1) with paper's f (5 on t with bar t extra, else 1);
# q=8: Collatz-Wielandt lower bound M g >= mu g with mu>1 (so sum_W P(W)^2 grows; lemma false)
import sys; sys.path.insert(0,'.')
from fractions import Fraction as Fr
from kz196_spectral import build_W, pairing
import numpy as np
def mat(q):
    v=q*q+q+1; W,v,p=build_W(q,None)
    bar=pairing(v,(0,1,2),list(range(v)))
    n=v+3
    M=[[W[bar[t]][u] if u!=bar[t] else Fr(0) for u in range(n)] for t in range(n)]
    return M,bar,v
for q in [16,32]:
    M,bar,v=mat(q); n=v+3
    f=[Fr(5) if bar[t]>=v else Fr(1) for t in range(n)]
    lam=max(sum(M[t][u]*f[u] for u in range(n))/f[t] for t in range(n))
    print('q',q,'paper-f certificate max (Mf)/f =',lam,float(lam),'<1:',lam<1,' sum f=',sum(f))
q=8; M,bar,v=mat(q); n=v+3
A=np.array([[float(x) for x in r] for r in M]); w,V=np.linalg.eig(A); i=np.argmax(abs(w)); g=np.abs(np.real(V[:,i]))
g=[Fr(int(x*10**9)+1,10**9) for x in g]
mu=min(sum(M[t][u]*g[u] for u in range(n))/g[t] for t in range(n))
print('q 8 Collatz-Wielandt exact lower bound mu =',float(mu),'>1:',mu>1)
