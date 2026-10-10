# read plantri -a on stdin, output WL-discrete triangulations as: V F f0a f0b f0c f1a ...
import sys
exec(open('sid_complex_search.py').read().split("if __name__ == '__main__':")[0])
tot=0;k=0
for line in sys.stdin:
    if not line.strip(): continue
    n,F=faces_from_plantri(line); tot+=1
    if wl_discrete(n,F):
        k+=1; print(n,len(F),' '.join('%d %d %d'%f for f in F))
print('# total',tot,'wl-discrete',k,file=sys.stderr)
