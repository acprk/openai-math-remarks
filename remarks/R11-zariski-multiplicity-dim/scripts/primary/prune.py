import sys,time
k=int(sys.argv[1]); B=int(sys.argv[2]); d=3**k
N=3*d+10
spf=list(range(N))
for i in range(2,int(N**.5)+1):
    if spf[i]==i:
        for j in range(i*i,N,i):
            if spf[j]==j: spf[j]=i
def fac(n):
    out={}
    while n>1:
        p=spf[n]; out[p]=out.get(p,0)+1; n//=p
    return out
V={}
for b in range(1,(d+1)//2):
    if b%3==0: continue
    ints=(d-b,b,2*d+b,2*d-b,3*d-b,d+b)
    fs=[fac(n) for n in ints]
    if any(max(f) > B for f in fs if f): continue
    out={}
    for grp,idx in ((0,((0,1),(1,-1))),(1,((2,1),(3,1),(4,-1),(5,-1)))):
        for i,s in idx:
            for p,e in fs[i].items(): out[(grp,p)]=out.get((grp,p),0)+s*e
    V[b]={kk:v for kk,v in out.items() if v}
alive=set(V)
cnt={}
for b in alive:
    for kk in V[b]: cnt[kk]=cnt.get(kk,0)+1
stack=[b for b in alive if any(cnt[kk]==1 for kk in V[b])]
while stack:
    b=stack.pop()
    if b not in alive: continue
    alive.discard(b)
    for kk in V[b]:
        cnt[kk]-=1
        if cnt[kk]==1:
            for b2 in alive:
                if kk in V[b2]: stack.append(b2); break
rows=sorted({kk for b in alive for kk in V[b]})
print(k,B,"cols",len(V),"alive",len(alive),"rows",len(rows),"excess",len(alive)-len(rows))
import json
json.dump({"d":d,"cols":sorted(alive),"V":{str(b):[[g,p,v] for (g,p),v in V[b].items()] for b in alive}},open("mat_%d_%d.json"%(k,B),"w"))
