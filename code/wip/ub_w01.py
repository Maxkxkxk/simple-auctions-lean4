"""NUMERICAL: best 2-/3-exponential kernel with w(0)=1 (so R bounded at 0), alpha(t) free."""
import sys; sys.path.insert(0,'wip')
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from ub_slack import kappa, phi
TS=np.concatenate([[1e-3,.02,.05],np.geomspace(.1,60,26),[90,140,200]])
def gs(a,s,k,full=False):
    out=[];al=[]
    for t in TS:
        v=np.concatenate([np.geomspace(max(t,1e-4),5000,2500)]); wv=sum(x*np.exp(-v/y) for x,y in zip(a,s)); lv=np.log(v); av=sum(x*y*(1-np.exp(-t/y)) for x,y in zip(a,s))/t
        def f(q):
            A=np.exp(q); x=A/v; h=np.maximum(1-np.where(x<1e-5,1-x/2,(1-np.exp(-x))/np.maximum(x,1e-300))-wv,0); return np.sum((h[1:]+h[:-1])/2*np.diff(lv))+k/A-av
        r=minimize_scalar(f,bounds=(np.log(1.2),np.log(60)),method='bounded'); out.append(r.fun); al.append(np.exp(r.x))
    return (np.array(out),np.array(al)) if full else max(out)
def f2(x):
    a=x[0];a=(a,1-a);s=(1.,np.exp(x[1]));b=x[2]; k=kappa(a,s,b); return 1+b+gs(a,s,k)
def f3(x):
    a1,a2=x[0],x[1];a=(a1,a2,1-a1-a2)
    if min(a)<0: return 9
    s=(1.,np.exp(x[2]),np.exp(x[3]));b=x[4];k=kappa(a,s,b); return 1+b+gs(a,s,k)
if __name__=='__main__':
    r=minimize(f2,[.775,np.log(26),1.36],method='Nelder-Mead',options=dict(maxfev=90,xatol=1e-3,fatol=1e-6)); x=r.x
    a=(x[0],1-x[0]);s=(1.,np.exp(x[1]));k=kappa(a,s,x[2]); print('2exp',r.fun,a,s,'b',x[2],'kappa',k); g,al=gs(a,s,k,True)
    print(np.round(TS,3));print(np.round(g,4));print(np.round(al,3))
    r=minimize(f3,[.6,.2,.5,np.log(26),1.36],method='Nelder-Mead',options=dict(maxfev=160,xatol=1e-3,fatol=1e-6)); print('3exp',r.fun,r.x)
