import sys,json,itertools
k=int(sys.argv[1]); B=int(sys.argv[2]); d=3**k; D=3*d
J=json.load(open("mat_%d_%d.json"%(k,B))); cols=J["cols"]
K=load("K_%d_%d.sobj"%(k,B))
def chainlen(c):
    l=0
    while c: a,r=divmod(D,c); l+=1; c=r
    return l
cost=[4+sum(2*(chainlen(c)-1) for c in (3*b,b,d-b,d+b)) for b in cols]
Wd=diagonal_matrix(ZZ,cost)
Kl=(K*Wd).LLL()
if len(sys.argv)>3: Kl=Kl.BKZ(block_size=int(sys.argv[3]))
rows=[vector(ZZ,[r[j]//cost[j] for j in range(len(cols))]) for r in Kl.rows()]
def s0of(e): return sum(abs(e[j])*cost[j] for j in range(len(cols)))+2*abs(sum(e))
best=None
cands=list(rows)+[rows[i]+s*rows[j] for i,j in itertools.combinations(range(min(len(rows),40)),2) for s in (1,-1)]
for e in cands:
    if sum(e)!=0:
        s=s0of(e)
        if best is None or s<best[0]: best=(s,e)
s,e=best
if sum(e)<0: e=-e
json.dump({"d":int(d),"relation":[{"b":int(cols[j]),"e":int(e[j])} for j in range(len(cols)) if e[j]]},open("relK_%d_%d.json"%(k,B),"w"))
print("approx s0",s)
