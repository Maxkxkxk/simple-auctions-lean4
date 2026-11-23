"""NUMERICAL: fully profile-adaptive with sequential splitting (sorted by cap, slopes [phi(U_{>=c})-phi(U_{>c})]/U_c)."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from ub_slack import Lf,phi
C=np.array([.02,.04,.07,.12,.2,.35,.6,1,1.8,3.2,6.,12.]); n=len(C)
s=np.geomspace(.01,200,34); t=np.geomspace(1e-3,400,220); dt=np.gradient(t)
Uc=1/(C[:,None]+t[None,:])
def slopes(r):
    X=r[:,None]*Uc; cum=np.cumsum(X[::-1],0)[::-1]; nxt=np.vstack([cum[1:],np.zeros((1,len(t)))])
    return np.where(X>1e-12,(phi(cum)-phi(nxt))/np.maximum(X,1e-12),1-np.exp(-cum))
def V(r):
    ch=slopes(r); nv=len(s); nz=n*len(t)
    E=[np.exp(-(C[i]+t)[:,None]/s[None,:]) for i in range(n)]
    A=sp.vstack([sp.hstack([sp.csr_matrix(-np.vstack(E)),-sp.eye(nz)]),sp.csr_matrix(np.concatenate([-np.ones(nv),np.zeros(nz)])[None])]).tocsr(); b=np.append(-ch.ravel(),-1)
    cost=np.concatenate([s*Lf(1/s)-sum(r[i]*s*(1-np.exp(-C[i]/s))/C[i] for i in range(n)),(r[:,None]*Uc*dt[None,:]).ravel()])
    res=linprog(cost,A_ub=A,b_ub=b,bounds=[(0,None)]*(nv+nz),method='highs',options={'time_limit':60.})
    return 2*0+1+res.fun,res.x[:nv]
def maximise(seed,iters=70):
    rng=np.random.default_rng(seed); r=rng.dirichlet(np.ones(n+1))[:n]*.9+1e-3; v,a=V(r); step=.5
    for it in range(iters):
        rr=r*np.exp(step*rng.normal(size=n)*(rng.random(n)<.4)); 
        if rr.sum()>1: rr/=rr.sum()
        v2,a2=V(rr)
        if v2>v: r,v,a=rr,v2,a2
        else: step=max(.05,step*.97)
    return v,np.round(r,3),[(round(x,3),round(y,3)) for x,y in zip(s,a) if y>1e-3]
if __name__=='__main__':
    from multiprocessing import Pool
    with Pool(8) as p:
        for o in p.imap_unordered(maximise,range(8)): print(o,flush=True)
