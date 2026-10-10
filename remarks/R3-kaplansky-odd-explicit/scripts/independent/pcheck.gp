\\ R3 independent check (PARI/GP; tested with the version printed below).
\\ Run:  gp -q independent/pcheck.gp < /dev/null
\\ Checks: m = C(86,41); p = least prime > m with p = 1 (mod 4); APRCL primality;
\\ an explicit Pocklington certificate (full factorisation of p-1, each factor re-proved);
\\ p does NOT divide (m!)^2+1 (so the paper's definition of p must be replaced);
\\ the numerical side conditions of 03-modular.tex (h/p^3 > 1/48, h/M < 1/4, (p^3-1)/M < 1/4)
\\ and t/m < 1/96.
print("gp version: ", version());
m=binomial(86,41); p=6054252965430758779982953;
print("m = ", m, "   p - m = ", p-m);
q=nextprime(m+1); while(q%4!=1, q=nextprime(q+1));
print("least prime > m with p = 1 mod 4: ", q, "   equals p: ", q==p);
print("APRCL isprime(p,2): ", isprime(p,2), "   p mod 4 = ", p%4);
F=factor(p-1); print("p-1 = ", F);
ok=1; for(i=1,#F~, ok = ok && isprime(F[i,1],2)); print("all factors of p-1 proved prime (APRCL): ", ok);
\\ Pocklington/Lucas: for every prime factor r of p-1 find a with a^(p-1)=1 and gcd(a^((p-1)/r)-1,p)=1
cert=vector(#F~); for(i=1,#F~, r=F[i,1]; a=2; while(!(Mod(a,p)^(p-1)==1 && gcd(lift(Mod(a,p)^((p-1)/r))-1,p)==1), a++); cert[i]=[r,a]);
print("Pocklington witnesses [r, a]: ", cert);
\\ with p-1 completely factored, these witnesses prove p prime (Lucas-Lehmer-Pocklington test)
print("(m!)^2+1 = 0 mod p ?  ", (Mod(-1,p)/Mod(192!,p))^2+1==0, "   (p = m+193, Wilson: m! = -1/192! mod p)");
ell=(p-1)/2; h=binomial(ell+2,3); M=p^4-1;
print("h/p^3 > 1/48: ", h/p^3>1/48, "   h/M < 1/4: ", h/M<1/4, "   (p^3-1)/M < 1/4: ", (p^3-1)/M<1/4, "   p > 600: ", p>600);
t=sum(j=0,28,binomial(86,j)); print("t = ", t, "   t/m < 1/96: ", t/m<1/96, "   (t/m)/(1/96) = ", (96*t/m)*1.);
