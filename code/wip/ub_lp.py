"""NUMERICAL: limit of the one-parameter (l) scheme: w=sum a_i e^{-u/s_i} (LP over many scales), eta(t)=1/alpha(t) free, b free.
 min b + g : kappa >= sum a s L(y/s) - b y  for all y;  g >= int_t (eta phi(1/(eta v)) - w/v)^+ dv + kappa*eta(t) - avg_0^t w."""
import sys; sys.path.insert(0,'wip')
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, minimize_scalar
from ub_slack import Lf,phi,total
def lp(eta,ts,s,nu=500):
    n=len(s); ys=np.geomspace(.05,400,160); ng=len(ts)
    cols=n+3  # a, kappa, g, b  then z blocks
    blocks=[];rows=[];rhs=[];Z=[]
    off=n+3; R=[];C=[];D=[];B=[];row=0
    for j,t in enumerate(ts):
        v=np.geomspace(t,3000,nu); dv=np.gradient(v); e=eta[j]; ph=e*phi(1/(e*v))
        E=np.exp(-v[:,None]/s[None,:])/v[:,None]
        # z>= ph - E a
        for m in range(nu):
            pass
        rr=np.arange(row,row+nu)
        R+= [np.repeat(rr,n),rr]; C+=[np.tile(np.arange(n),nu),off+np.arange(nu)]; D+=[(-E).ravel(),-np.ones(nu)]; B.append(-ph); row+=nu
        # g row
        R+=[np.full(n,row),np.array([row]),np.array([row]),np.full(nu,row)]; C+=[np.arange(n),np.array([n]),np.array([n+1]),off+np.arange(nu)]
        D+=[-s*(1-np.exp(-t/s))/t,np.array([e]),np.array([-1.]),dv]; B.append(np.array([0.])); row+=1
        off+=nu
    for y in ys:
        R+=[np.full(n,row),np.array([row,row])]; C+=[np.arange(n),np.array([n,n+2])]; D+=[s*Lf(y/s),np.array([-1.,-y])]; B.append(np.array([0.])); row+=1
    # w(0)>=1 (drift)
    R+=[np.full(n,row)];C+=[np.arange(n)];D+=[-np.ones(n)];B.append(np.array([-1.]));row+=1
    A=sp.csr_matrix((np.concatenate(D),(np.concatenate(R),np.concatenate(C))),shape=(row,off))
    c=np.zeros(off); c[n+1]=1;c[n+2]=1
    res=linprog(c,A_ub=A,b_ub=np.concatenate(B),bounds=[(0,None)]*n+[(None,None),(None,None),(0,None)]+[(0,None)]*(off-n-3),method='highs')
    x=res.x; return x[:n],x[n],x[n+1],x[n+2],res.fun
def besteta(a,s,k,ts):
    out=[];g=[]
    for t in ts:
        v=np.geomspace(t,4000,3000); lv=np.log(v); wv=(a[None,:]*np.exp(-v[:,None]/s[None,:])).sum(1); av=(a*s*(1-np.exp(-t/s))).sum()/t
        def f(le):
            e=np.exp(le); h=np.maximum(e*phi(1/(e*v))-wv/v,0)*v; return np.sum((h[1:]+h[:-1])/2*np.diff(lv))+k*e-av
        r=minimize_scalar(f,bounds=(np.log(.01),np.log(3)),method='bounded'); out.append(np.exp(r.x)); g.append(r.fun)
    return np.array(out),np.array(g)
if __name__=='__main__':
    ts=np.geomspace(1e-3,200,26); s=np.geomspace(.05,600,32); eta=1/np.sqrt(6+4*ts)
    for it in range(8):
        a,k,g,b,f=lp(eta,ts,s); e2,g2=besteta(a,s,k,ts)
        print(it,'bound',1+f,'b',round(b,4),'g',round(g,4),'kappa',round(k,4),'g after eta',g2.max().round(4),[(round(x,2),round(y,4)) for x,y in zip(s,a) if y>1e-4],flush=True)
        eta=np.exp(.5*np.log(eta)+.5*np.log(e2))
    print(np.round(ts,3)); print(np.round(1/eta**2,3)); print(np.round(g2,4))
