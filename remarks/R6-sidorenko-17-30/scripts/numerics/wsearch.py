"""Numerical search for W with t(H,W) < (E W)^30, H = incidence graph of the V=7 complex.
Rectangular (oriented) kernels W: [k] x [l] with weights mu, nu; also report the transpose
orientation and the symmetric crossed product t(W)t(W^T).  Multi-start Adam + LBFGS.
Run with any Python 3 environment that has PyTorch installed."""
import torch, sys, itertools, math, random
torch.set_num_threads(3); torch.set_default_dtype(torch.float64)
FACES=[(0,1,2),(0,1,6),(0,2,3),(0,3,4),(0,4,5),(0,5,6),(1,2,5),(1,5,6),(2,3,4),(2,4,5)]
L='abcdefg'
eq=','.join(L[a]+L[b]+L[c] for a,b,c in FACES)+','+','.join(L)+'->'
def logt(W,mu,nu):  # points ~ mu on rows, faces ~ nu on cols
    g=torch.einsum('ay,by,cy,y->abc',W,W,W,nu)
    t=torch.einsum(eq,*([g]*10),*([mu]*7))
    return torch.log(t)
def obj(p,k,l):
    W=torch.nn.functional.softplus(p[:k*l].view(k,l))+1e-9
    mu=torch.softmax(p[k*l:k*l+k],0); nu=torch.softmax(p[k*l+k:],0)
    m=mu@W@nu
    return logt(W,mu,nu)-30*torch.log(m), (logt(W.T.contiguous(),nu,mu)-30*torch.log(m))
def run(k,l,seed,scale):
    torch.manual_seed(seed)
    p=(torch.randn(k*l+k+l)*scale).requires_grad_()
    opt=torch.optim.Adam([p],lr=0.05)
    for it in range(1500):
        opt.zero_grad(); f,_=obj(p,k,l); f.backward(); opt.step()
    opt=torch.optim.LBFGS([p],max_iter=500,line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); f,_=obj(p,k,l); f.backward(); return f
    try: opt.step(cl)
    except Exception as e: pass
    f,fT=obj(p,k,l); return f.item(),fT.item(),p.detach()
best=(1e9,None)
sizes=[(k,l) for k in range(2,9) for l in range(2,9)]
random.seed(1)
for rep in range(int(sys.argv[1]) if len(sys.argv)>1 else 3):
    for (k,l) in sizes:
        for scale in (1.0,3.0,6.0):
            seed=random.randrange(10**9)
            f,fT,p=run(k,l,seed,scale)
            if f<best[0]: best=(f,(k,l,seed,scale))
            print('k=%d l=%d scale=%.0f seed=%d  min log(t/p^30)=%.3e  transpose=%.3e  sum=%.3e'%(k,l,scale,seed,f,fT,f+fT),flush=True)
            if f< -1e-9: torch.save(p,'neg_%d_%d_%d.pt'%(k,l,seed)); print('NEGATIVE FOUND',flush=True)
print('BEST',best)
