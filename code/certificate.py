"""EXACT integer verification of the A-group certificate (no floating point).
Items: type i (n copies), P(X>=k) = floor(R_i/max(k,p_i))/2^d for 1<=k<=p_i+e, 0 beyond. Lattice unit: N = concentrated-group total.
Run:  python certificate.py instance.json
"""
import json, sys
from fractions import Fraction as Fr
try:
    from gmpy2 import mpz   # optional, only for speed (GMP multiplication)
except ImportError:
    mpz=int
def tl(R,p,e,W,k): return W if k==0 else (R//max(k,p) if k<=p+e else 0)
def check(N,e,n,d,P,R,T0,verbose=True):
    W=1<<d; m=len(P); Dn=W**(n*m); s=d*n*m+20; X=1<<s
    c=[[tl(r,p,e,W,k)-tl(r,p,e,W,k+1) for k in range(p+e+1)] for p,r in zip(P,R)]
    assert all(min(x)>=0 and sum(x)==W for x in c)
    # SREV
    srev=sum(Fr(n*max(k*tl(r,p,e,W,k) for k in range(p+e+1)),W) for p,r in zip(P,R))
    EZ=sum(Fr(n*sum(k*x for k,x in enumerate(ci)),W) for ci in c)
    def _pack(co):
        if len(co)==1: return mpz(co[0])
        h=len(co)//2; return _pack(co[:h])+(_pack(co[h:])<<(s*h))
    def powmod(a,k,T):
        M=(mpz(1)<<(s*T))-1; r=mpz(1)
        while k:
            if k&1: r=(r*a)&M
            k>>=1
            if k: a=(a*a)&M
        return r
    M=(mpz(1)<<(s*T0))-1; F=1
    for ci in c: F=(F*powmod(_pack(ci[:T0]),n,T0))&M
    C=0;S=0;worst=Fr(0); mask=X-1; ok=True
    for t in range(T0):
        # P(Z>=t): D-C
        if t>=1:
            lhs=(N+t)*(Dn-C); 
            if lhs>N*Dn: ok=False; worst=max(worst,Fr(lhs,N*Dn))
        C+= int((F>>(s*t))&mask); S+=C     # C = C_{t+1}, S=sum_{1..t+1} C
    lhs=(N+T0)*(Dn-C); ok&=lhs<=N*Dn
    low=(N+T0)*C-S    # sum_{z<T0}(N+z)a_z
    tail=Fr(N)+EZ-Fr(low,Dn)    # E[(N+Z)1{Z>=T0}]
    ok2= tail<=N
    # revenue
    g=[];h=[]
    for ci,p in zip(c,P):
        gg=[sum(ci[:p+1])]+ci[p+1:]; hh=[p*x for x in ci[p:]]
        g.append(_pack(gg[:e]));h.append(_pack(hh[:e]))
    Me=(mpz(1)<<(s*e))-1
    pw=[powmod(x,n,e) for x in g]; pw1=[powmod(x,n-1,e) for x in g]
    G=1
    for x in pw: G=G*x&Me
    H=0
    for i in range(m):
        q=h[i]*pw1[i]*n&Me
        for j in range(m):
            if j!=i: q=q*pw[j]&Me
        H+=q
    H&=Me
    g_low=sum(int((G>>(s*k))&mask) for k in range(e)); h_low=sum(int((H>>(s*k))&mask) for k in range(e))
    Epay=sum(Fr(n*p*tl(r,p,e,W,p),W) for p,r in zip(P,R))
    rev=e*(1-Fr(g_low,Dn))+Epay-Fr(h_low,Dn)
    ok3= srev<=N
    if verbose:
        print('BREV prefix ok',ok,float(worst),'tail ok',ok2,float(tail/N),'SREV/N',float(srev/N),ok3,'rev/N',float(rev/N),'ratio sup', float((N+rev)/N) if ok3 else None)
    return ok,ok2,ok3,rev/N,worst
if __name__=='__main__':
    I=json.load(open(sys.argv[1])); r=check(I['N'],I['e'],I['n'],I['d'],I['P'],I['R'],I['T0']); print('1+kappa =',float(1+r[3])); print('exceeds 2.190270:', 1+r[3]>Fr(2190270,10**6)+Fr(1,10**6))
