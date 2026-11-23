"""NUMERICAL: V(l) with t-dependent Jensen weights h_c(t)=(1-mu(t)) eta(c/l)/l + mu(t) zeta(c,t)."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from ub_slack import Lf,phi
def V(l,mu,zeta,eta,nu=240,full=False):
    cs=np.concatenate([l*np.geomspace(1e-3,300,22)]); s=np.geomspace(.004,120,30); n=len(s)
    R=[];C=[];D=[];B=[];row=0;off=n+2
    for c in cs:
        u=np.geomspace(c,600,nu); du=np.gradient(u); t=u-c; h=(1-mu(t))*eta(c/l)/l+mu(t)*zeta(c,t); ph=h*phi(1/(h*u))
        E=np.exp(-u[:,None]/s[None,:])/u[:,None]; rr=np.arange(row,row+nu)
        R+=[np.repeat(rr,n),rr];C+=[np.tile(np.arange(n),nu),off+np.arange(nu)];D+=[(-E).ravel(),-np.ones(nu)];B.append(-ph);row+=nu
        R+=[np.full(n,row),np.array([row,row]),np.full(nu,row)];C+=[np.arange(n),np.array([n,n+1]),off+np.arange(nu)]
        D+=[-s*(1-np.exp(-c/s))/c,np.array([eta(c/l)/l,-1.]),du];B.append(np.zeros(1));row+=1;off+=nu
    R.append(np.full(n,row));C.append(np.arange(n));D.append(-np.ones(n));B.append(np.array([-1.]));row+=1
    A=sp.csr_matrix((np.concatenate(D),(np.concatenate(R),np.concatenate(C))),shape=(row,off))
    cost=np.zeros(off);cost[:n]=s*Lf(1/s);cost[n]=-1;cost[n+1]=1
    r=linprog(cost,A_ub=A,b_ub=np.concatenate(B),bounds=[(0,None)]*n+[(None,None),(0,None)]+[(0,None)]*(off-n-2),method='highs')
    if full: return 1+r.fun,r.x[n],r.x[n+1],[(round(x,3),round(y,3)) for x,y in zip(s,r.x[:n]) if y>1e-3]
    return 1+r.fun
E0=lambda t:1/np.sqrt(5.4+4.1*t)
def job(a):
    l,kind,tau=a
    if kind=='none': return a,V(l,lambda t:0*t,lambda c,t:1+0*t,E0)
    mu=lambda t:t/(t+tau)
    if kind=='one': return a,V(l,mu,lambda c,t:1+0*t,E0)
    if kind=='one2': return a,V(l,lambda t:t*t/(t*t+tau*tau),lambda c,t:1+0*t,E0)
    if kind=='tc': return a,V(l,mu,lambda c,t:t/(c+t),E0)
    if kind=='half': return a,V(l,lambda t:.5*t/(t+tau),lambda c,t:1+0*t,E0)
if __name__=='__main__':
    from multiprocessing import Pool
    ls=(.12,.15,.18,.22); jobs=[(l,'none',0) for l in ls]+[(l,k,tau) for k in ('one','one2','tc','half') for tau in (.3,1.,3.) for l in ls]
    with Pool(11) as p: res=dict(p.map(job,jobs,chunksize=1))
    for k,tau in sorted({(j[1],j[2]) for j in jobs}): print(k,tau,[round(res[(l,k,tau)],4) for l in ls])
