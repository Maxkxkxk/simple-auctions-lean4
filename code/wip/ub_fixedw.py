"""NUMERICAL: fixed kernel w (B=1 scale), Jensen weights lam_c=r_c (uniform).  LP: min sum a s L(1/s)+g,
 f_w(c)=int_c^oo (phi(1/u)-w(u)/u)^+ du - (1/c)int_0^c w <= g for all c."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from ub_slack import Lf,phi
def solve(s,cs,nu=700,umax=3000):
    u=np.unique(np.concatenate([np.geomspace(1e-4,umax,nu),cs])); n=len(s); m=len(u)
    du=np.gradient(u)
    E=np.exp(-u[:,None]/s[None,:])
    # vars a(n), z(m)>=phi-w/u, g
    A1=sp.hstack([sp.csr_matrix(-E/u[:,None]),-sp.eye(m),sp.csr_matrix((m,1))]); b1=-phi(1/u)
    rows=[];
    for c in cs:
        rows.append(np.concatenate([-s*(1-np.exp(-c/s))/c, np.where(u>=c,du,0),[-1]]))
    A=sp.vstack([A1,sp.csr_matrix(np.array(rows))]).tocsr()
    # c->0: w(0)>=1 needed
    A=sp.vstack([A,sp.csr_matrix(np.concatenate([-np.ones(n),np.zeros(m+1)])[None])]).tocsr()
    b=np.concatenate([b1,np.zeros(len(cs)),[-1.]])
    cost=np.concatenate([s*Lf(1/s),np.zeros(m),[1]])
    r=linprog(cost,A_ub=A,b_ub=b,bounds=[(0,None)]*(n+m)+[(0,None)],method='highs')
    a=r.x[:n]; return 1+r.fun,r.x[-1],[(round(x,4),round(y,4)) for x,y in zip(s,a) if y>1e-5]
if __name__=='__main__':
    cs=np.geomspace(2e-4,200,80)
    print(solve(np.geomspace(.01,1000,60),cs))
    for s in (.5,1,1.5,2,3,4): print(s,solve(np.array([s]),cs)[:2])
