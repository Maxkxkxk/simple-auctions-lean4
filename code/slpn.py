"""Like slp_seed.Model but type i consists of n[i] iid items (n=inf -> Poisson limit)."""
from slp_seed import *
class ModelN(Model):
    def __init__(s,N,e,prices,n,K,beta=1.0):
        super().__init__(N,e,prices,K); s.n=list(n); s.beta=beta; s.lmax=s.n
    def _f(s,vec,tot,n,L):  # fft of single item law / poisson; returns (F_all, F_minus_one, F_minus_two)
        g=np.zeros(L); g[:len(vec)]=vec; g0=g.sum()
        fh=np.fft.rfft(g)
        if np.isinf(n):
            F=np.exp(fh-g0); return F,F,F
        one=1-g0/n+fh/n
        return one**n, one**(n-1), (one**(n-2) if n>=2 else None)
    def evaluate(s,lam,grad=False):
        N,e,P,K,L=s.N,s.e,s.P,s.K,s.L; m=len(P)
        inv=lambda F: np.maximum(np.fft.irfft(F,L),0)
        tail=lambda x: np.cumsum(x[::-1])[::-1]
        ZF=[s._f(lam[i],0,s.n[i],L) for i in range(m)]
        def prod(T,ex):  # ex: dict i->number of excluded copies
            r=np.ones(L//2+1,complex)
            for i in range(m): r=r*T[i][ex.get(i,0)]
            return r
        out=dict(ZT=tail(inv(prod(ZF,{}))))
        WF=[];PF=[]
        for i,p in enumerate(P):
            om=np.zeros(K+1); om[:K+1-p]=lam[i,p:]; pj=p*om.copy(); om[0]=0
            WF.append(s._f(om,0,s.n[i],L)); g=np.zeros(L); g[:K+1]=pj; PF.append(np.fft.rfft(g))
        off=e+K+1
        def ext(x):
            A=np.zeros(L+off); A[off:]=x; A[:off]=x[0]; return A
        def revs(ex):   # returns array over d=0..e : (P(W'>=e-d), E[pay';W'>=e-d])
            A=ext(tail(inv(prod(WF,ex))))
            B=0
            for j in range(m):
                c=(s.n[j]-ex.get(j,0))
                if c<=0: continue
                e2=dict(ex); e2[j]=e2.get(j,0)+1
                B=B+prod(WF,e2)*PF[j]*(1.0 if np.isinf(s.n[j]) else c/s.n[j])
            Bt=ext(tail(np.fft.irfft(B,L))) if not np.isscalar(B) else np.zeros(L+off)
            d=np.arange(e+1)
            return A[off+e-d],Bt[off+e-d]
        A,Bv=revs({})
        out['R']=(e*A[0]+Bv[0])/N
        if grad:
            G=np.zeros_like(lam); DZ=[]
            kk=s.k
            for i,p in enumerate(P):
                A,Bv=revs({i:1})
                d=np.clip(kk-p,0,e)
                G[i]=e*(A[d]-A[0])+Bv[d]-Bv[0]+np.where(kk>=p,p*A[d],0)
                DZ.append(tail(inv(prod(ZF,{i:1}))))
            out['G']=G/N; out['DZ']=DZ
        return out
def slpn(m,lam,iters=100,tr=0.3,verbose=False,srcap=1.0):
    P,K,N,L=m.P,m.K,m.N,m.L; n=len(P); be=m.beta
    TM=min(L//2,6*K)
    Tri=np.triu(np.ones((K,K)))*(np.arange(1,K+1)[:,None]/N)
    Asrev=sp.hstack([sp.kron(sp.eye(n),sp.csr_matrix(Tri)), sp.kron(-sp.eye(n),np.ones((K,1)))]).tocsr()
    rs=[np.concatenate([np.zeros(n*K),np.ones(n)])]; 
    for i in range(n):
        if not np.isinf(m.n[i]):
            v=np.zeros(n*K+n); v[i*K:(i+1)*K]=1; rs.append(v)
    rs=sp.csr_matrix(np.array(rs)); nfin=[i for i in range(n) if not np.isinf(m.n[i])]
    merit=lambda f,vi: f-20*max(0,vi-1)
    f,sr,vi=m.value(lam); cur=merit(f,vi); small=0
    for it in range(iters):
        o=m.evaluate(lam,True); ZT=o['ZT']; t=np.arange(1,TM); t=t[t/N>=1-be-1e-9]
        ta=t[(be+t/N)*ZT[t]>0.75]; kk=np.arange(1,K+1); blocks=[]
        for i in range(n):
            ext=np.concatenate([np.ones(K+1),o['DZ'][i]])
            blocks.append((ext[K+1+ta[:,None]-kk[None,:]]-o['DZ'][i][ta][:,None])*(be+ta/N)[:,None])
        Ab=sp.hstack([sp.csr_matrix(np.concatenate(blocks,1)),sp.csr_matrix((len(ta),n))])
        bb=1-(be+ta/N)*ZT[ta]
        x0=lam[:,1:].ravel(); c=np.concatenate([-o['G'][:,1:].ravel(),np.ones(n)]); ok=False
        rb=[srcap]+[max(m.n[i]-lam[i,1:].sum(),0) for i in nfin]
        Aub=sp.vstack([Asrev,Ab,rs]).tocsr(); bub=np.concatenate([-Asrev[:,:n*K]@x0,np.maximum(bb,-1),rb])
        for tries in range(8):
            step=tr*(x0+3e-3/(np.tile(kk,n)/N)/N)
            b=np.concatenate([np.stack([np.maximum(-x0,-step),step],1),np.stack([np.zeros(n),np.ones(n)],1)])
            res=linprog(c,A_ub=Aub,b_ub=bub,bounds=b,method='highs-ds' if n*K<1500 else 'highs-ipm',options={'time_limit':15.0})
            if res.x is None: tr*=0.5; continue
            l2=lam.copy(); l2[:,1:]=np.maximum(x0+res.x[:n*K],0).reshape(n,K)
            f2,sr2,vi2=m.value(l2)
            if merit(f2,vi2)>cur-1e-12: ok=True;break
            tr*=0.5
        if not ok: break
        gain=merit(f2,vi2)-cur; lam,f,sr,vi,cur=l2,f2,sr2,vi2,merit(f2,vi2); tr=min(tr*1.3,1.0)
        if verbose: print(it,'%.6f %.4f %.1e %.3f'%(1+be+f,sr,vi-1,tr),flush=True)
        small=small+1 if gain<2e-7 else 0
        if small>=3: break
    return lam
if __name__=='__main__':
    I=np.inf
    for n in ([I,I,I],[I,I,1],[I,I,2],[I,I,4],[I,1,I],[I,4,I],[1,I,I],[2,I,I],[4,I,I]):
        m=ModelN(20,32,(0,4,64),n,96); lam=fix(m,slpn(m,init(m),150)); f,sr,vi=m.value(lam); print(n,2+f,sr,m.srev_types(lam).round(3),lam[:,1:].sum(1).round(3),flush=True)
