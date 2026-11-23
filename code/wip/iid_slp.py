"""NUMERICAL: iid items in the Poisson limit, single price p + entry fee e.  max R s.t. SREV<=1, BREV<=1."""
import sys; sys.path.insert(0,'.')
from slp_seed import Model
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
def val(m,lam):
    o=m.evaluate(lam,True); t=np.arange(1,m.L//2); b=(t/m.N*o['ZT'][t]).max(); s=m.srev_types(lam).sum(); return o['R'],s,b,o
def run(N,e,p,K,iters=80,tr=.3,seed=0):
    m=Model(N,e,(p,),K); k=np.arange(K+1); rng=np.random.default_rng(seed)
    a=max(p,1); T=np.where(k>=1,.5*N/np.maximum(k,a),0)*(1+.1*rng.random()); lam=(T-np.append(T[1:],0))[None,:].copy(); lam[0,0]=0
    def scale(l):
        lo,hi=0,3.
        for _ in range(40):
            x=(lo+hi)/2; v=val(m,x*l)
            if v[1]<=1 and v[2]<=1: lo=x
            else: hi=x
        return lo*l
    lam=scale(lam); R,s,b,o=val(m,lam); merit=lambda R,s,b:R-20*max(0,b-1)-20*max(0,s-1); cur=merit(R,s,b)
    Tri=np.triu(np.ones((K,K)))*(np.arange(1,K+1)[:,None]/N)
    for it in range(iters):
        ZT=o['ZT'];t=np.arange(1,4*K);t=t[t/N*ZT[t]>.8]; kk=np.arange(1,K+1); ext=np.concatenate([np.ones(K+1),ZT])
        D=(ext[K+1+t[:,None]-kk[None,:]]-ZT[t][:,None])*(t/N)[:,None]
        x0=lam[0,1:]; ok=False
        for _ in range(8):
            st=tr*(x0+2e-3/(kk/N)/N)
            r=linprog(-o['G'][0,1:],A_ub=np.vstack([Tri,D]),b_ub=np.concatenate([1-Tri@x0,1-t/N*ZT[t]]),bounds=np.stack([np.maximum(-x0,-st),st],1),method='highs-ds',options={'time_limit':10.})
            if r.x is None: tr*=.5; continue
            l2=lam.copy(); l2[0,1:]=np.maximum(x0+r.x,0); R2,s2,b2,o2=val(m,l2)
            if merit(R2,s2,b2)>cur-1e-12: ok=True;break
            tr*=.5
        if not ok: break
        g=merit(R2,s2,b2)-cur; lam,o,cur=l2,o2,merit(R2,s2,b2); tr=min(1,tr*1.3)
        if g<1e-7 and it>5: break
    lam=scale(lam); R,s,b,o=val(m,lam); return R/max(s,b),(N,e,p,K),s,b,lam[0]
def job(a): 
    try: return run(*a)
    except Exception as ex: return (0,a,0,0,None)
if __name__=='__main__':
    from multiprocessing import Pool
    jobs=[(20,e,p,max(p+e,40)+20) for e in (6,10,14,18,22) for p in (0,2,4,6,8,12,16)]
    with Pool(11) as P: res=sorted(P.map(job,jobs,chunksize=1),key=lambda r:-r[0])
    for r in res[:8]: print(r[:4])
    lam=res[0][4]; T=np.cumsum(lam[::-1])[::-1]; print([round(float(T[k]*k/20),3) for k in range(1,len(T))])
