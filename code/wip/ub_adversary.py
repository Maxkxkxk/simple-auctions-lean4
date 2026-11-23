"""NUMERICAL: adversary against the GH master bound S + E[sum Y - max Y] in the Poisson limit.
Types a (cap c=a): tail T_a(u)=r_a/a for u<=a, T_a nonincreasing, T_a(u)<=r_a/u.  Drift beta (tiny items).
Constraint: (beta+t)P(Z>=t)<=1, SREV<=1.  Maximise  1(S, filled) + beta + int phi(N(t))dt, N(t)=sum_a T_a(a+t).
With --laplace only the (weaker) constraints sum_j int lam e^{-lam u}Fbar <= Psi(lam) are imposed (what the proof uses)."""
import sys, numpy as np, scipy.sparse as sp
sys.path.insert(0,'.')
from slp_seed import cpois
from scipy.optimize import linprog
def run(N,A,K,beta,iters=60,tr=.3,verbose=False,seed=0):
    n=len(A); A=list(A); L=1<<int(np.ceil(np.log2(8*K)))
    k=np.arange(K+1)
    # variables mass lam[i,k], k>=A[i]; T=cumsum
    rng=np.random.default_rng(seed)
    r=np.array([.6/n]*n)*(0.5+rng.random(n)); lam=np.zeros((n,K+1))
    for i,a in enumerate(A):
        T=np.zeros(K+2); kk=np.arange(a,K+1); T[a:K+1]=r[i]*N/kk; lam[i]=np.where(k>=a,T[:-1]-T[1:],0)
    def ev(l):
        Z=cpois(l.sum(0),L); ZT=np.cumsum(Z[::-1])[::-1]
        T=np.cumsum(l[:,::-1],1)[:,::-1]
        Nt=np.zeros(K+1)
        for i,a in enumerate(A): Nt[:K+1-a]+=T[i,a:]
        # value at lattice: Y=k-a; Fbar_Y(t) on (j,j+1) equals T[a+j+1]
        obj=(Nt[1:]-1+np.exp(-Nt[1:])).sum()/N
        r=np.array([T[i,a]*a/N for i,a in enumerate(A)])
        viol=(T*k/N-r[:,None]).max()
        t=np.arange(max(1,int(round((1-beta)*N))),L//2); bv=((beta+t/N)*ZT[t]).max()
        return obj,r.sum(),viol,bv,ZT,Nt
    def scale(l):
        lo,hi=0,1
        for _ in range(40):
            m=(lo+hi)/2
            if ev(m*l)[3]<=1 and ev(m*l)[1]<=1: lo=m
            else: hi=m
        return lo*l
    lam=scale(lam); o=ev(lam); merit=lambda o:o[0]-20*max(0,o[3]-1); cur=merit(o)
    # structure: T_i(k)k/N <= r_i = T_i(a)a/N
    blocks=[]
    Tri=np.triu(np.ones((K,K)))
    for i,a in enumerate(A):
        M=Tri*(np.arange(1,K+1)[:,None]/N) - np.ones((K,1))*(np.where(np.arange(1,K+1)>=a,a/N,0))[None,:]
        blocks.append(sp.csr_matrix(M[a:]))
    As=sp.block_diag(blocks).tocsr(); rcoef=np.concatenate([np.where(np.arange(1,K+1)>=a,a/N,0) for a in A])
    for it in range(iters):
        ZT,Nt=o[4],o[5]; g=(1-np.exp(-Nt))/N
        G=np.zeros((n,K))
        cg=np.concatenate([[0],np.cumsum(g[1:])])
        for i,a in enumerate(A):
            kk=np.arange(1,K+1); G[i]=np.where(kk>=a,cg[np.maximum(kk-a,0)],0)
        TM=min(L//2,4*K); t=np.arange(max(1,int(round((1-beta)*N))),TM); t=t[(beta+t/N)*ZT[t]>.75]
        ext=np.concatenate([np.ones(K+1),ZT]); kk=np.arange(1,K+1)
        D=(ext[K+1+t[:,None]-kk[None,:]]-ZT[t][:,None])*(beta+t/N)[:,None]
        Aub=sp.vstack([As,sp.csr_matrix(np.tile(D,(1,n))),sp.csr_matrix(rcoef[None])]).tocsr()
        x0=lam[:,1:].ravel(); b=np.concatenate([np.maximum(-(As@x0),0),1-(beta+t/N)*ZT[t],[1-rcoef@x0]])
        ok=False
        for _ in range(7):
            st=tr*(x0+2e-3/(np.tile(kk,n)/N)/N); mask=np.concatenate([kk>=a for a in A])
            up=np.where(mask,st,0); lo=np.maximum(-x0,-st)
            res=linprog(-G.ravel(),A_ub=Aub,b_ub=b,bounds=np.stack([lo,up],1),method='highs-ipm' if n*K>1500 else 'highs-ds',options={'time_limit':20.})
            if res.x is None: tr*=.5; continue
            l2=lam.copy(); l2[:,1:]=np.maximum(x0+res.x,0).reshape(n,K); o2=ev(l2)
            if merit(o2)>cur-1e-12: ok=True;break
            tr*=.5
        if not ok: break
        gain=merit(o2)-cur; lam,o,cur=l2,o2,merit(o2); tr=min(1,tr*1.3)
        if verbose: print(it,1+beta+o[0],o[1:4],flush=True)
        if gain<1e-6 and it>4: break
    lam=scale(lam);o=ev(lam)
    return 1+beta+o[0],o[1:4],[(a,round(float(lam[i].sum()*a/N),4)) for i,a in enumerate(A) if lam[i].sum()*a/N>1e-3],lam
if __name__=='__main__':
    from multiprocessing import Pool
    jobs=[(10,(1,2,3,5,10,20,30,40),80,1.),(10,(1,2,4,40),120,1.),(10,(1,2,3,5,8,12,20),60,1.),(10,(1,2,4,8,16,32),160,1.),(10,(1,2,3,5,10,20,40),80,.8),(10,(1,2,3,5,10,20,40),80,.5),(20,(1,2,3,5,8,12,20,40,60),160,1.)]
    with Pool(7) as p:
        for j,r in zip(jobs,p.starmap(run,jobs)): print(j[0],j[1:],round(r[0],5),r[1],r[2],flush=True)
