"""NUMERICAL: for a point-mass cap profile (all items have cap c, S=1, B=1) best dual kernel w (LP over exponentials)
G(c,w)= int_c^oo (phi(1/u)-w(u)/u)^+ du - (1/c) int_0^c w + sum a_i s_i L(1/s_i); bound 1+G.  Also real adversary for one type."""
import sys; sys.path.insert(0,'wip'); sys.path.insert(0,'.')
import numpy as np
from scipy.optimize import linprog
from ub_slack import Lf,phi
def pointLP(c,s=np.geomspace(.01,300,41),S=1.):
    u=np.geomspace(c,3000,1500); du=np.gradient(u)
    E=np.exp(-u[:,None]/s[None,:]); ph=phi(S/u)
    n=len(s); m=len(u)
    cost=np.concatenate([s*Lf(1/s)-S*s*(1-np.exp(-c/s))/c, du])
    A=np.hstack([-E*S/u[:,None],-np.eye(m)]); b=-ph
    r=linprog(cost,A_ub=A,b_ub=b,bounds=[(0,None)]*(n+m),method='highs-ds')
    a=r.x[:n]; return S+r.fun,[(round(x,3),round(y,3)) for x,y in zip(s,a) if y>1e-4]
if __name__=='__main__':
    for c in (.003,.01,.03,.1,.2,.3,.5,.7,1,1.5,2,3,5,10): print(c,pointLP(c))
    from ub_adversary import run
    for a in (1,2,3,5,7,10,15,20,30): 
        r=run(10,(a,),120,1.); print('real adversary one type c=',a/10,round(r[0],4),'S_A',round(r[1][0],3))
