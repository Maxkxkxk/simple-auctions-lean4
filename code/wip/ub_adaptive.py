"""NUMERICAL: fully profile-adaptive scheme: V(r)=1+min_w[ sum a s L(1/s) + sum_c r_c( int (chord_c(t)-w(c+t))^+/(c+t)dt - avg w ) ]
 chord = psi(U_r(t)) (ideal t-dependent Jensen weights) or GH weights. sup over profiles r (sum r<=1, rest at infinity)."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, brentq
from ub_slack import Lf,phi
C=np.array([.02,.04,.07,.12,.2,.35,.6,1,1.8,3.2,6.,12.]); n=len(C)
s=np.geomspace(.01,200,34); t=np.geomspace(1e-3,400,220); dt=np.gradient(t); t0=1e-3
psi=lambda x: phi(x)/np.maximum(x,1e-12)
def V(r,ideal=True):
    Uc=1/(C[:,None]+t[None,:]); U=(r[:,None]*Uc).sum(0)
    if ideal: ch=np.repeat(psi(U)[None],n,0)
    else:
        f=lambda l:(r/np.sqrt(6*l*l+4*l*C)).sum()-1
        l=brentq(f,1e-5,1); h=1/np.sqrt(6*l*l+4*l*C); ch=np.array([psi(Uc[i]/h[i]) for i in range(n)])
    R=[];nv=len(s); nz=n*len(t)
    E=[np.exp(-(C[i]+t)[:,None]/s[None,:]) for i in range(n)]
    A=sp.hstack([sp.csr_matrix(-np.vstack(E)),-sp.eye(nz)]).tocsr(); b=-ch.ravel()
    A=sp.vstack([A,sp.csr_matrix(np.concatenate([-np.ones(nv),np.zeros(nz)])[None])]).tocsr(); b=np.append(b,-1)
    wt=(r[:,None]*Uc*dt[None,:]).ravel(); wt=wt+0
    cost=np.concatenate([s*Lf(1/s)-sum(r[i]*s*(1-np.exp(-C[i]/s))/C[i] for i in range(n)),wt])
    res=linprog(cost,A_ub=A,b_ub=b,bounds=[(0,None)]*(nv+nz),method='highs',options={'time_limit':60.})
    x=res.x; z=x[nv:].reshape(n,-1); a=x[:nv]
    # gradient wrt r: direct + via U
    avg=np.array([(a*s*(1-np.exp(-C[i]/s))/C[i]).sum() for i in range(n)])
    H=(z*Uc*dt[None]).sum(1)-avg
    if ideal:
        dpsi=(psi(U+1e-5)-psi(U-1e-5))/2e-5; on=(z>1e-9)
        M=(r[:,None]*Uc*on*dpsi[None]*dt[None]).sum(0)
        H=H+(M[None]*Uc).sum(1)
    return 2*0+1+res.fun,H,a
def maximise(ideal,seed=0,iters=25):
    rng=np.random.default_rng(seed); r=rng.dirichlet(np.ones(n+1))[:n]*.8+1e-3; step=.15; v,H,a=V(r,ideal)
    for it in range(iters):
        g=H-0.0
        rr=r*np.exp(step*(g-(r*g).sum())/ (abs(g).max()+1e-9)); 
        if rr.sum()>1: rr/=rr.sum()
        v2,H2,a2=V(rr,ideal)
        if v2>v: r,v,H,a=rr,v2,H2,a2; step*=1.3
        else: step*=.5
        if step<2e-3: break
    return v,np.round(r,3),[(round(x,3),round(y,3)) for x,y in zip(s,a) if y>1e-3]
def job(x): return (x,)+maximise(*x)
if __name__=='__main__':
    from multiprocessing import Pool
    with Pool(8) as p:
        for o in p.imap_unordered(job,[(True,0),(True,1),(True,2),(True,3),(False,0),(False,1),(False,2),(False,3)]): print(o,flush=True)
