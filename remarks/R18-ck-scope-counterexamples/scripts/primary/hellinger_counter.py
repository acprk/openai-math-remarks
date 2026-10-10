# Universal / p-biased Hellinger analogue is false: n=1, f = biased dictator on Omega={0,1}, mu=(1-w,w)
# LHS = sqrt(1-m^2) - E sqrt(1-(R_lam f)^2), RHS = 1 - sqrt(1-lam^2)
import numpy as np
Hh=lambda t: np.sqrt(1-t*t)
for w in [0.25,0.3,0.4,0.45]:
    m=2*w-1
    out=[]
    for lam in [0.05,0.1,0.3,0.5,0.7,0.9,0.99]:
        up=m+lam*(1-m); dn=m+lam*(-1-m)
        lhs=Hh(m)-(w*Hh(up)+(1-w)*Hh(dn)); rhs=1-Hh(lam)
        out.append((lam,round(lhs,6),round(rhs,6),'VIOL' if lhs>rhs+1e-15 else 'ok'))
    print(f"w={w} m={m}:",out)
