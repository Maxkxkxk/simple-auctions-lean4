import sys; sys.path.insert(0,'wip')
import numpy as np
from ub_lp import lp,besteta
import ub_lp
from ub_slack import Lf,phi
# find adversary: at fixed optimal (a,kappa,eta) list regions where psi(alpha/v)>w(v), and g(t) equalisation; dual multipliers via re-solve
import scipy.sparse as sp
from scipy.optimize import linprog
ts=np.concatenate([[.002,.02,.06],np.geomspace(.1,60,24),[90,140,200]]); s=np.array([.5,1.4,20,60,400.]); 
eta=1/np.sqrt(5.4+4.1*ts)
for it in range(6):
    a,k,g,b,f=lp(eta,ts,s,nu=260); e2,g2=besteta(a,s,k,ts); eta=np.exp(.4*np.log(eta)+.6*np.log(e2))
print('bound',1+f,a,k,g,b)
# dual of g rows: rebuild LP quickly with marginals
def lp_marg(eta,ts,s,nu=260):
    n=len(s); ys=np.geomspace(.05,400,160); R=[];C=[];D=[];B=[];row=0;off=n+3;gr=[]
    for j,t in enumerate(ts):
        v=np.geomspace(t,3000,nu); dv=np.gradient(v); e=eta[j]; ph=e*phi(1/(e*v)); E=np.exp(-v[:,None]/s[None,:])/v[:,None]; rr=np.arange(row,row+nu)
        R+=[np.repeat(rr,n),rr];C+=[np.tile(np.arange(n),nu),off+np.arange(nu)];D+=[(-E).ravel(),-np.ones(nu)];B.append(-ph);row+=nu
        R+=[np.full(n,row),np.array([row]),np.array([row]),np.full(nu,row)];C+=[np.arange(n),np.array([n]),np.array([n+1]),off+np.arange(nu)]
        D+=[-s*(1-np.exp(-t/s))/t,np.array([e]),np.array([-1.]),dv];B.append(np.array([0.]));gr.append(row);row+=1;off+=nu
    y0=row
    for y in ys:
        R+=[np.full(n,row),np.array([row,row])];C+=[np.arange(n),np.array([n,n+2])];D+=[s*Lf(y/s),np.array([-1.,-y])];B.append(np.array([0.]));row+=1
    R+=[np.full(n,row)];C+=[np.arange(n)];D+=[-np.ones(n)];B.append(np.array([-1.]));dr=row;row+=1
    A=sp.csr_matrix((np.concatenate(D),(np.concatenate(R),np.concatenate(C))),shape=(row,off)); c=np.zeros(off);c[n+1]=1;c[n+2]=1
    res=linprog(c,A_ub=A,b_ub=np.concatenate(B),bounds=[(0,None)]*n+[(None,None),(None,None),(0,None)]+[(0,None)]*(off-n-3),method='highs')
    m=-res.ineqlin.marginals
    return res.fun,m[gr],[(round(y,3),round(q,4)) for y,q in zip(ys,m[y0:y0+160]) if q>1e-4],m[dr],[ (t,np.round(m[g0-nu:g0]/np.maximum(m[g0],1e-12)/np.gradient(np.geomspace(t,3000,nu)),2)[::12]) for t,g0 in zip(ts,gr) if m[g0]>1e-3]
f,mg,my,dr,on=lp_marg(eta,ts,s)
print('profile r_t:',[(round(t,3),round(x,4)) for t,x in zip(ts,mg) if x>1e-3],'sum',mg.sum()); print('y=B/l mass:',my,'drift',dr)
for t,o in on: print('c/l=',round(t,3),'on-fraction along v (geom grid t..3000):',list(o))
