"""NUMERICAL: two-parameter variant. Units B=1,S=1. For given l: profile satisfies sum r eta(c/l)=l, sum r<=1.
Jensen weights h_c=(1-mu)*eta(c/l)/l+mu.  V(l)=min_{w,kap} sum a s L(1/s) - kap + sup_c[ H_c + kap*eta(c/l)/l ] +1"""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from ub_slack import Lf,phi
def V(l,mu,eta,nu=260,S=1.):
    cs=np.concatenate([l*np.geomspace(1e-3,400,22)]); s=np.geomspace(.004,120,30); n=len(s)
    R=[];C=[];D=[];B=[];row=0;off=n+2
    for c in cs:
        u=np.geomspace(c,600,nu); du=np.gradient(u); h=(1-mu)*eta(c/l)/l+mu; ph=h*phi(1/(h*u))
        E=np.exp(-u[:,None]/s[None,:])/u[:,None]; rr=np.arange(row,row+nu)
        R+=[np.repeat(rr,n),rr];C+=[np.tile(np.arange(n),nu),off+np.arange(nu)];D+=[(-E).ravel(),-np.ones(nu)];B.append(-ph);row+=nu
        R+=[np.full(n,row),np.array([row,row]),np.full(nu,row)];C+=[np.arange(n),np.array([n,n+1]),off+np.arange(nu)]
        D+=[-s*(1-np.exp(-c/s))/c,np.array([eta(c/l)/l,-1.]),du];B.append(np.zeros(1));row+=1;off+=nu
    R.append(np.full(n,row));C.append(np.arange(n));D.append(-np.ones(n));B.append(np.array([-1.]));row+=1
    A=sp.csr_matrix((np.concatenate(D),(np.concatenate(R),np.concatenate(C))),shape=(row,off))
    cost=np.zeros(off);cost[:n]=s*Lf(1/s);cost[n]=-1;cost[n+1]=1
    r=linprog(cost,A_ub=A,b_ub=np.concatenate(B),bounds=[(0,None)]*n+[(None,None),(0,None)]+[(0,None)]*(off-n-2),method='highs')
    return 1+r.fun
def job(a):
    l,mu,k=a; e=(lambda t:1/np.sqrt(5.5+4*t)) if k==0 else ((lambda t:1/(2.4+1.0*t)) if k==1 else (lambda t:1/(2.4+t**.75*1.3)))
    return a,V(l,mu,e)
if __name__=='__main__':
    from multiprocessing import Pool
    jobs=[(l,mu,k) for k in (0,1,2) for l in (.04,.08,.12,.16,.22,.3) for mu in (0,.1,.25,.4)]
    with Pool(10) as p: res=p.map(job,jobs,chunksize=1)
    for k in (0,1,2):
        for mu in (0,.1,.25,.4): print(k,mu,[round(v,4) for (l,m,kk),v in res if m==mu and kk==k])
