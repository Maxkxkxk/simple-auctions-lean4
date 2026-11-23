"""SLP search for an A-seed in the Poisson (many-small-items) limit, two-part-tariff mechanism.
Normalised: concentrated group B sold at 1.  Lattice spacing h=1/N.
Maximise REV_A - SREV_A  s.t. P(Z_A>=t) <= 1/(1+t); ratio = 2 + REV_A - SREV_A (if SREV_A<=1).
"""
import numpy as np, sys, time
from scipy.optimize import linprog
import scipy.sparse as sp

def cpois(nu, L):
    f = np.zeros(L); f[:len(nu)] = nu
    s = f.sum()
    g = np.fft.irfft(np.exp(np.fft.rfft(f) - s), L)
    return np.maximum(g, 0)

class Model:
    def __init__(s, N, e, prices, K, L=None):
        s.N=N; s.e=e; s.P=list(prices); s.K=K
        s.L = L or 1<<int(np.ceil(np.log2(10*K)))
        s.k=np.arange(K+1)
    def evaluate(s, lam, grad=False):
        N,e,P,K,L=s.N,s.e,s.P,s.K,s.L
        nu=lam.sum(0)
        Z=cpois(nu,L); ZT=np.cumsum(Z[::-1])[::-1]   # ZT[t]=P(Z>=t)
        om=np.zeros(K+1)
        for i,p in enumerate(P):
            om[:K+1-p]+=lam[i,p:]
        w0=om[0]; om2=om.copy(); om2[0]=0
        W=cpois(om2,L); WT=np.cumsum(W[::-1])[::-1]
        A=np.ones(L+e+K+2); A[e+K+1:e+K+1+L]=WT; A[e+K+1+L:]=0  # A(x) at index x+e+K+1
        off=e+K+1
        pay=np.zeros(K+1)
        for i,p in enumerate(P): pay[:K+1-p]+=p*lam[i,p:]
        def H(x):
            return float(np.dot(pay, A[off+x-np.arange(K+1)]))
        R=e*A[off+e]+H(e)
        out=dict(R=R/N, ZT=ZT, A=A, off=off, srev=0)
        if grad:
            Hs=np.array([H(e-d) for d in range(0,e+1)])  # surplus capped at e
            As=np.array([A[off+e-d] for d in range(0,e+1)])
            G=np.zeros_like(lam)
            for i,p in enumerate(P):
                d=np.clip(s.k-p,0,e)
                G[i]=e*(As[d]-As[0])+Hs[d]-Hs[0]+np.where(s.k>=p,p*As[d],0)
            out['G']=G/N
        return out
    def srev_types(s,lam):
        T=np.cumsum(lam[:,::-1],1)[:,::-1]
        return (T*s.k/s.N).max(1)
    def value(s,lam):
        o=s.evaluate(lam); sr=s.srev_types(lam).sum()
        t=np.arange(1,s.L//2)
        viol=((1+t/s.N)*o['ZT'][1:s.L//2]).max()
        return o['R']-sr, sr, viol

def slp(m, lam, iters=60, tr=0.3, verbose=True, TM=None):
    P,K,N,L=m.P,m.K,m.N,m.L
    n=len(P)
    TM=TM or min(L//2, 4*K)
    # variables lam (n*(K)) for k=1..K  and r (n)
    nv=n*K+n
    rows=[];cols=[];vals=[]
    # sum_{l>=k} lam_l * k/N - r <=0
    # dense lower part
    Tri=np.triu(np.ones((K,K)))*(np.arange(1,K+1)[:,None]/N)
    Asrev=sp.hstack([sp.kron(sp.eye(n),sp.csr_matrix(Tri)), sp.kron(-sp.eye(n),np.ones((K,1)))]).tocsr()
    rows_S=sp.csr_matrix(np.concatenate([np.zeros(n*K),np.ones(n)])[None])
    best=None; f,sr,vi=m.value(lam)
    merit=lambda f,vi: f-20*max(0,vi-1)
    cur=merit(f,vi)
    for it in range(iters):
        o=m.evaluate(lam,True)
        ZT=o['ZT']; t=np.arange(1,TM)
        be=getattr(m,'beta',1.0); t=t[t/N>=1-be-1e-9]; act=np.where((be+t/N)*ZT[t]>0.80)[0]
        ta=t[act]
        kk=np.arange(1,K+1)
        ext=np.concatenate([np.ones(K+1),ZT])
        D=ext[K+1+ta[:,None]-kk[None,:]]-ZT[ta][:,None]   # dP(Z>=t)/dnu_k
        D*= (be+ta/N)[:,None]
        Ab=sp.hstack([sp.csr_matrix(np.tile(D,(1,n))), sp.csr_matrix((len(ta),n))]).tocsr()
        bb=1-(be+ta/N)*ZT[ta]
        x0=lam[:,1:].ravel()
        # in terms of delta
        c=np.concatenate([-o['G'][:,1:].ravel(), np.ones(n)])
        r0=m.srev_types(lam)
        ok=False
        for tries in range(8):
            step=tr*(x0+3e-3/ (kk[None,:].repeat(n,0).ravel()/N)/N )
            lo=np.maximum(-x0,-step); hi=step
            bounds=np.concatenate([np.stack([lo,hi],1), np.stack([np.zeros(n),np.ones(n)],1)])
            Aub=sp.vstack([Asrev,Ab,rows_S]).tocsr()
            bub=np.concatenate([-Asrev[:,:n*K]@x0, np.maximum(bb,-1),[1.0]])
            res=linprog(c,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs-ds' if n*K<1500 else 'highs-ipm',options={'time_limit':15.0})
            if res.x is None: tr*=0.5; continue
            lam2=lam.copy(); lam2[:,1:]=np.maximum(x0+res.x[:n*K],0).reshape(n,K)
            f2,sr2,vi2=m.value(lam2)
            if merit(f2,vi2)>cur-1e-12:
                ok=True;break
            tr*=0.5
        if not ok: break
        gain=merit(f2,vi2)-cur
        lam,f,sr,vi,cur=lam2,f2,sr2,vi2,merit(f2,vi2)
        tr=min(tr*1.3,1.0)
        if verbose: print(it,'ratio %.6f srev %.5f viol-1 %.2e tr %.3f'%(2+f,sr,vi-1,tr),flush=True)
        if gain<1e-7 and it>3: break
    return lam

def fix(m,lam):
    """scale down to exact feasibility"""
    lo,hi=0.,1.
    for _ in range(50):
        mid=(lo+hi)/2
        if m.value(mid*lam)[2]<=1: lo=mid
        else: hi=mid
    return lo*lam

def init(m):
    lam=np.zeros((len(m.P),m.K+1))
    for i,p in enumerate(m.P):
        a=max(p,1 if p else max(1,m.N//30)); b=min(m.K,p+m.e)
        if a>=b: continue
        T=np.zeros(m.K+2); k=np.arange(a,b+1); T[a:b+1]=0.6/len(m.P)*m.N/k
        T[:a]=T[a]
        lam[i]=T[:-1]-T[1:]
        lam[i,0]=0
    return fix(m,lam)

if __name__=='__main__':
    N=int(sys.argv[1]); e=int(sys.argv[2]); K=int(sys.argv[3]); P=[int(x) for x in sys.argv[4].split(',')]
    m=Model(N,e,P,K)
    lam=init(m); print('init',m.value(lam))
    t=time.time(); lam=slp(m,lam); lam=fix(m,lam)
    f,sr,vi=m.value(lam); print('FINAL ratio',2+f,'srev',sr,vi,'time',time.time()-t)
    print('r',dict(zip(P,m.srev_types(lam).round(4))))
    np.save('wip/lam_%d_%d_%d.npy'%(N,e,K),lam)
