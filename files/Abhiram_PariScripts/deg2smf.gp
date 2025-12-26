L(s,D) = lfun(lfuncreate(D),s);


h(r,N) = my(D = (-1)^r*N); if(D%4 == 0 || D%4 == 1, (-1)^(ceil(r/2))*(r-1)!*N^(r-1/2)*2^(1-r)*Pi^(-r)*L(D,r), 0);

sqfac(n) = my(l = divisors(n), t= listcreate(), len = length(l)); for(i=1,len, if(n%l[i]^2==0, listput(t,l[i]))); return(Vec(t));

T(r,f,X)= sumdiv(f,d,moebius(d)*kronecker(X,d)*d^(r-1)*sigma(f/d,2*r-1));


 H(r,N) = my(D=(-1)^r*N, f = coredisc(D,1)[2], X = coredisc(D,1)[1],t); if(D==0, t=zeta(1-2*r), if(N > 0, if(D%4 == 0 || D%4 ==1, t= L(1-r,X)*T(r,f,X),t=0), t=0)); return(t);

div(x) = if(x==0,[1],divisors(x));

comdiv(V) = my(l = length(V), k1 = div(V[1])); for(i=2,l, k1 = setintersect(k1, div(V[i]))); return(Vec(k1));

nonzero(V) = my(k=listcreate(), len = length(V),t); if(vecprod(V)!=0, t= V, for(i=1,len, if(V[i]!= 0, listput(k,V[i]));); t=Vec(k); ); return(t);

posdisc(n) = my(k= listcreate()); for(a=0,n, for(b=0,n, for(c=0,n, if(4*a*c -b^2 >= 0, listput(k,[a,b,c]))))); return(Vec(k));

C(a,b,c,k) = my(abc=comdiv(nonzero([a,b,c])),r=k-1,D=4*a*c-b^2,t);t=vecsum(apply(x->x^(k-1)*H(r,D/x^2),abc));return(bestappr(-2*k*t/(bernfrac(k)*zeta(3-2*k))));

SiegelEisensteinCoeffs(k,ord)= my(K=posdisc(ord),len=length(K));for(i=2,len,x=K[i][1];y=K[i][2];z=K[i][3];print(K[i],": ",C(x,y,z,k)));

kill(y);

qpyser(k,ord) = my(K=posdisc(ord),len=length(K), t = 0); for(i=2,len,n=K[i][1]; l=K[i][2]; m = K[i][3]; t = t+ C(n,l,m,k)*q^(n)*y^(l)*p^(m)); return(t)

print("Printing coefficients up to order ", lim, " for weight k=", wt," Siegel EEisenstein series" );

SiegelEisensteinCoeffs(wt,lim)

print("Coefficients up to order ", lim, " for weight k=", wt," Siegel Eisenstein series" );