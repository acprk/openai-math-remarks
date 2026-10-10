import sys,json,time
k=int(sys.argv[1]); B=int(sys.argv[2])
J=json.load(open("mat_%d_%d.json"%(k,B))); cols=J["cols"]; V={int(b):v for b,v in J["V"].items()}
keys=sorted({(g,q) for b in cols for g,q,v in V[b]}); ri={kk:i for i,kk in enumerate(keys)}
M=matrix(ZZ,len(keys),len(cols))
for j,b in enumerate(cols):
    for g,q,v in V[b]: M[ri[(g,q)],j]=v
t=time.time(); Kq=M.change_ring(QQ).right_kernel_matrix(); print("QQ kernel",Kq.nrows(),time.time()-t,flush=True)
den=lcm([x.denominator() for x in Kq.list()]); Kz=(den*Kq).change_ring(ZZ)
t=time.time(); Kz=Kz.LLL(); print("LLL1",time.time()-t,flush=True)
# saturate: kernel lattice = (Q-span) ∩ Z^n ; use saturation
t=time.time(); Ks=Kz.saturation(); print("saturation",time.time()-t,flush=True)
save(Ks,"K_%d_%d"%(k,B))
