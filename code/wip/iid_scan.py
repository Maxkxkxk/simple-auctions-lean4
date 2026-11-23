"""NUMERICAL: iid items, ER law on [a, b] lattice (tail r/k), n items; entry-fee tariff (p,e). ratio R/max(S,B)."""
import numpy as np, itertools
def powlaw(f,n,L):
    g=np.zeros(L); g[:len(f)]=f; return np.maximum(np.fft.irfft(np.fft.rfft(g)**n,L),0)
def ratio(T,n,ps,es):
    # T[k]=Pr[X>=k], k=0..K
    K=len(T)-1; f=T-np.append(T[1:],0); L=1<<int(np.ceil(np.log2(K*n+2)))
    L=min(L,1<<14)
    Z=powlaw(f,n,L); tail=np.cumsum(Z[::-1])[::-1]; B=(np.arange(L)*tail).max(); S=n*(np.arange(K+1)*T).max()
    best=(0,)
    for p in ps:
        g=np.zeros(K-p+1); g[0]=f[:p+1].sum(); g[1:]=f[p+1:]; W1=powlaw(g,n-1,L); W=powlaw(g,n,L); A=lambda x:1-W[:max(x,0)].sum()
        c1=np.concatenate([[0],np.cumsum(W1)])
        for e in es:
            R=e*A(e)+n*p*sum(f[k]*(1-c1[max(e-(k-p),0)]) for k in range(p,K+1))
            if R>best[0]: best=(R,p,e)
    return best[0]/max(S,B),S,B,best
if __name__=='__main__':
    res=[]
    for n in (10,20,40,100,400):
        for a in (1,2,3,4,6,8):
            for b in (a+2,a+4,a+6,a+9,a+14,a+20,a+30):
                for r in (.3,.45,.6,.8,1.0):
                    for nscale in (1.,2.,3.,4.,6.):
                        # tail r*a... Pr[X>=k]=q a/k on [a,b], total rate n q = nscale
                        q=nscale/n
                        if q>1: continue
                        k=np.arange(b+1); T=np.where(k==0,1,q*a/np.maximum(k,a)); 
                        x=ratio(T,n,range(0,a+2),range(1,2*b,max(1,b//6))); res.append((x[0],n,a,b,nscale,x[1:]))
                    break
    res.sort(key=lambda t:-t[0]); 
    for r in res[:8]: print(r)
    for n in (10,20,40,100,400): print(n,max(r for r in res if r[1]==n)[:5])
