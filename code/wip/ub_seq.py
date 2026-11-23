"""NUMERICAL: sequential (sorted by cap) splitting.  phi(sum x) <= sum_c x_c * avg_{R in (R_>,R_>=]} (1-e^{-R/(c+t)}).
bound = 1 + sum a s L(1/s) + int_0^1 h(R)dR,  h(R)=sup_c [ int_c^oo (1-e^{-R/u}-w(u))^+ du/u - (1/c)int_0^c w ].  LP in w."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from ub_slack import Lf
def solve(nR=16,nu=220,s=np.geomspace(.01,300,40),cs=None,umax=3000,verbose=True):
    x,wq=np.polynomial.legendre.leggauss(nR); Rs=(x+1)/2; wq=wq/2
    u=np.geomspace(1e-3,umax,nu); lu=np.log(u); dl=np.gradient(lu)
    cs=u[::4] if cs is None else cs
    n=len(s); m=nu
    E=np.exp(-u[:,None]/s[None,:])
    nv=n+nR+nR*m
    R=[];C=[];D=[];B=[];row=0
    for i,r in enumerate(Rs):
        off=n+nR+i*m; rr=np.arange(row,row+m)
        R+=[np.repeat(rr,n),rr];C+=[np.tile(np.arange(n),m),off+np.arange(m)];D+=[(-E).ravel(),-np.ones(m)];B.append(-(1-np.exp(-r/u)));row+=m
        for c in u[::4]:
            mask=np.where(u>=c*0.999)[0]
            R+=[np.full(n,row),np.array([row]),np.full(len(mask),row)];C+=[np.arange(n),np.array([n+i]),off+mask];D+=[-s*(1-np.exp(-c/s))/c,np.array([-1.]),dl[mask]];B.append(np.zeros(1));row+=1
    R.append(np.full(n,row));C.append(np.arange(n));D.append(-np.ones(n));B.append(np.array([-1.]));row+=1
    A=sp.csr_matrix((np.concatenate(D),(np.concatenate(R),np.concatenate(C))),shape=(row,nv))
    cost=np.zeros(nv);cost[:n]=s*Lf(1/s);cost[n:n+nR]=wq
    res=linprog(cost,A_ub=A,b_ub=np.concatenate(B),bounds=[(0,None)]*n+[(None,None)]*nR+[(0,None)]*(nR*m),method='highs')
    a=res.x[:n]; h=res.x[n:n+nR]
    return 1+res.fun,[(round(p,3),round(q,4)) for p,q in zip(s,a) if q>1e-5],list(zip(Rs.round(3),h.round(4)))
if __name__=='__main__':
    for nR,nu in ((10,160),(16,260)):
        r=solve(nR,nu); print(r[0]); print(r[1]); print(r[2])
