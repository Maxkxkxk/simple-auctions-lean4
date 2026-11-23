import sys; sys.path.insert(0,'.'); sys.path.insert(0,'wip')
import numpy as np, math, time
from ub_certificate import *
from ub_slack import Lf,kappa
def geom(a,b,r):
    g=[a]
    while g[-1]<b: g.append(max(g[-1]+1,int(g[-1]*r)))
    return g
def design(AN,SN,bt,grid,nlev=26,tmax=260):
    a=[x/1000 for x in AN]; s=[x/10 for x in SN]
    k=kappa(a,s,bt)
    # float alpha(t) optimum
    from scipy.optimize import minimize_scalar
    ts=np.concatenate([[0.002],np.geomspace(.15,tmax,nlev)]); lev=[]
    for i,t in enumerate(ts):
        tc=math.sqrt(t*(ts[i+1] if i+1<len(ts) else t)) if i else .1
        v=np.geomspace(tc,20000,3000); wv=sum(x*np.exp(-v/y) for x,y in zip(a,s)); lv=np.log(v)
        def f(q):
            A=np.exp(q); x=A/v; h=np.maximum(1-np.where(x<1e-5,1-x/2,(1-np.exp(-x))/np.maximum(x,1e-300))-wv,0); return np.sum((h[1:]+h[:-1])/2*np.diff(lv))+k/A
        r=minimize_scalar(f,bounds=(np.log(1.2),np.log(70)),method='bounded'); lev.append((0 if i==0 else int(t*D),int(round(np.exp(r.x)*D))))
    return k,lev
if __name__=='__main__':
    for AN,SN,bt in (([743,257],[10,225],1.3382),([728,197,75],[10,150,600],1.3522),([739,224,37],[10,180,1000],1.3549),([272,504,199,25],[5,14,220,1500],1.3674)):
        grid=geom(2,20000*D,1.006); t=time.time()
        k,lev=design(AN,SN,bt,grid); g,prof=run(AN,SN,lev,grid,int(k*SC)); 
        print(AN,SN,'cells',len(grid),'kappa',k,'g',g/SC,'TOTAL',1+bt+g/SC,round(time.time()-t,1),flush=True); print([ (p[0],p[1],round(p[2],4)) for p in prof])
