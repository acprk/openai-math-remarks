from collections import Counter
def vec(cx,a,b):
    n=Counter(('p',a[i]) for i in range(cx.V)); n.update(('f',b[j]) for j in range(cx.F))
    n.update(('c',a[i],b[j]) for (i,j) in cx.corners); return n
if __name__=='__main__':
    from cplx import load
    cx=load('cands.txt')[1]
    idv=vec(cx,tuple(range(6)),tuple(range(8)))
    for f,g in [(((0,1,2,3,0,5),(0,1,2,3,4,0,1,2)),((4,1,2,3,4,5),(5,6,7,3,4,5,6,7))),
                (((0,1,1,5,0,5),(0,1,1,4,4,0,1,1)),((4,2,2,3,4,3),(5,7,7,3,3,5,7,7)))]:
        s=vec(cx,*f); s.update(vec(cx,*g))
        print(all(s[k]==2*idv[k] for k in set(s)|set(idv)))
