# PROVENANCE ONLY (not part of the certificate, not run by run_all.sh).
# Heuristic search used to find the nine points: maximize min_x max_i |<x,p_i>|^2 over 9 points of CP^2
# (soft-min/soft-max objective, L-BFGS-B, Nelder-Mead hole refinement). Run as `python3 search_capcover.py 4`
# (numpy + scipy). The run with seed 4 saved audit_capcover_best_4.npy; P9.txt is that array written as
# exact hex floats (see make_P9.py). Re-running may not reproduce the same array bit-for-bit across
# numpy/scipy versions; this does not matter, because the certificate is about the fixed points in P9.txt.
# Audit: can 9 points in CP^2 have covering radius < pi/4, i.e. every unit x has max_i |<x,p_i>|^2 > 1/2 ?
# If yes, the 9 closed balls {|<x,p_i>|^2 >= t}, t>1/2, have no orthogonal pair (|<x,p>|^2+|<y,p>|^2<=1) => 9-cover of CP^2 projectors.
import numpy as np, sys
from scipy.optimize import minimize
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def rnd(n): 
    z=rng.normal(size=(n,3))+1j*rng.normal(size=(n,3)); return z/np.linalg.norm(z,axis=1,keepdims=True)
def unpack(v,k): 
    z=(v[:3*k]+1j*v[3*k:]).reshape(k,3); return z/np.linalg.norm(z,axis=1,keepdims=True)
def pack(z): z=z.reshape(-1); return np.concatenate([z.real,z.imag])
def maxov(x,P): return (np.abs(x@P.conj().T)**2).max(-1)
def holes(P,ns=20000,keep=200):
    X=rnd(ns); m=maxov(X,P); idx=np.argsort(m)[:keep]; H=[]
    for x0 in X[idx]:
        r=minimize(lambda v: maxov(unpack(v,1),P)[0], pack(x0), method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-12,'maxiter':4000})
        H.append(unpack(r.x,1)[0])
    H=np.array(H); return H, maxov(H,P)
k=9; P=rnd(k); S=rnd(3000)
best=0
for it in range(40):
    beta=200.
    def obj(v):
        Q=unpack(v,k); ov=np.abs(S@Q.conj().T)**2
        sm=np.log(np.exp(beta*ov).sum(1))/beta     # soft max over points
        return -(-np.log(np.exp(-beta*sm).sum())/beta)  # soft min over samples
    r=minimize(obj,pack(P),method='L-BFGS-B',options={'maxiter':500}); P=unpack(r.x,k)
    H,mh=holes(P,ns=20000,keep=60)
    cov=mh.min(); print(it,'cover value min_x max_i ov =',cov,flush=True)
    if cov>best: best=cov; np.save('audit_capcover_best_%s.npy'%(sys.argv[1] if len(sys.argv)>1 else 0),P)
    S=np.vstack([S,H[np.argsort(mh)[:40]]])
print('BEST',best)
