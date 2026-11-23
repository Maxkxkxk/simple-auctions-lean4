"""EXACT (fractions): iid lower bound.  n iid items, Pr[X>=k] = A/(n*max(k,a)) for 1<=k<=b (equal revenue on [a,b]);
mechanism: entry fee e, then item price p=a.  Run: python iid_lower_bound.py"""
from fractions import Fraction as Fr
import evaluator as E
def check(n,a,b,A,p,e):
    tl=lambda k: Fr(1) if k==0 else (Fr(A,n*max(k,a)) if k<=b else Fr(0))
    f=[tl(k)-tl(k+1) for k in range(b+1)]; assert min(f)>=0 and sum(f)==1
    den=1
    from math import lcm
    for x in f: den=lcm(den,x.denominator)
    c=[int(x*den) for x in f]
    def mul(x,y):
        z=[0]*(len(x)+len(y)-1)
        for i,u in enumerate(x):
            if u:
                for j,v in enumerate(y): z[i+j]+=u*v
        return z
    def power(x,k):
        r=[1]
        while k:
            if k&1: r=mul(r,x)
            k>>=1
            if k: x=mul(x,x)
        return r
    Z=power(c,n); D=den**n; assert sum(Z)==D
    B=0; t=0
    for v in range(len(Z)-1,-1,-1):
        t+=Z[v]; B=max(B,v*t)
    B=Fr(B,D); S=n*max(k*tl(k) for k in range(b+1))
    g=[sum(c[:p+1])]+c[p+1:]; W=power(g,n)[:e]; W1=power(g,n-1)
    loss=sum(c[k]*p*sum(W1[:max(e-(k-p),0)]) for k in range(p,b+1))*n
    R=e*(1-Fr(sum(W),D))+n*p*sum(f[p:])-Fr(loss,D)
    return S,B,R
if __name__=='__main__':
    # cross-check with brute force on a tiny instance
    S,B,R=check(4,1,3,Fr(8,5),1,2)
    F={Fr(0):Fr(3,5),Fr(1):Fr(1,5),Fr(2):Fr(1,15),Fr(3):Fr(2,15)}
    assert S==E.srev([F]*4) and B==E.brev([F]*4)
    from itertools import combinations
    menu=[(tuple(int(i in s) for i in range(4)),2+len(s)) for k in range(1,5) for s in combinations(range(4),k)]
    assert R==E.menu_revenue([F]*4,menu)
    for n,a,b,A,p,e in ((8,1,3,Fr(16,5),1,2),(20,3,9,9,3,6),(40,4,13,12,4,9),(100,4,13,12,4,9)):
        S,B,R=check(n,a,b,A,p,e); M=max(S,B)
        print('n=%d: SREV=%.6f BREV=%.6f REV>=%.6f ratio>=%.6f  exact ratio > 1.1877: %s, > 1.30: %s'%(n,S,B,R,R/M,R/M>Fr(11877,10000),R/M>Fr(13,10)))
