import sys; sys.path.insert(0,'.'); sys.path.insert(0,'wip')
import numpy as np, math, json, time
from ub_certificate import *
from ub_slack import kappa
from scipy.optimize import minimize_scalar
def geom(a,b,r):
    g=[a]
    while g[-1]<b: g.append(max(g[-1]+1,int(g[-1]*r)))
    return g[:-1]+[b]
def mkgrid(fine):
    g=geom(2,100,1.2)
    for a,b,r in ((100,300,1.06),(300,1000,1+5*fine),(1000,4000,1+2*fine),(4000,600000,1+fine),(600000,3000000,1+2*fine),(3000000,3*10**8,1.05)): g=g[:-1]+geom(a,b,r)
    return g
def alphas(AN,SN,k,nlev,tmax=230):
    a=[x/1000 for x in AN]; s=[x/10 for x in SN]
    ts=np.concatenate([[0,.015,.035,.06,.09],np.geomspace(.125,tmax,nlev)]); lev=[]
    for i,t in enumerate(ts):
        tc=math.sqrt(t*(ts[i+1] if i+1<len(ts) else t)) if i else .004
        v=np.geomspace(tc,30000,4000); wv=sum(x*np.exp(-v/y) for x,y in zip(a,s)); lv=np.log(v)
        def f(q):
            A=np.exp(q); x=A/v; h=np.maximum(1-np.where(x<1e-5,1-x/2,(1-np.exp(-x))/np.maximum(x,1e-300))-wv,0); return np.sum((h[1:]+h[:-1])/2*np.diff(lv))+k/A
        r=minimize_scalar(f,bounds=(np.log(1.2),np.log(80)),method='bounded'); lev.append([int(t*D),int(round(np.exp(r.x)*D))])
    for i in range(1,len(lev)): lev[i][1]=max(lev[i][1],lev[i-1][1])
    return lev
def total(AN,SN,ystar,grid,nlev,xr=1.02,full=False):
    tans=[]
    for A,S in zip(AN,SN):
        # y0 = ystar/s = ystar*10/S
        y0n,y0d=int(ystar*1000)*10,S*1000; y0=y0n/y0d
        xmax=1+min(400,14/y0); xs=[1000]
        while xs[-1]<1000*xmax: xs.append(max(xs[-1]+1,int(xs[-1]*xr)))
        bn,kn=best_tan(y0n,y0d,xs); tans.append(dict(y0=(y0n,y0d),xs=xs,b=bn,k=kn))
    b=cdiv(sum(A*t['b'] for A,t in zip(AN,tans)),1000); kap=cdiv(sum(A*S*t['k'] for A,S,t in zip(AN,SN,tans))*SC,10**10)
    tab=table(AN,SN,grid); L=alphas(AN,SN,kap/SC,nlev)
    gs=[level(al,ts,(L[i+1][0] if i+1<len(L) else None),kap,tab,grid[-1],AN,SN) for i,(ts,al) in enumerate(L)]
    g=max(gs); tot=1+b/1e6+g/SC
    if full: return tot,b,kap,g,L,tans,gs
    return tot
if __name__=='__main__':
    grid=mkgrid(.008); print(len(grid))
    for AN,SN,ys in (([272,504,199,25],[5,14,220,1500],5.9),([272,504,199,25],[5,14,220,1500],6.0),([272,504,199,25],[5,14,220,1500],6.1),([257,510,167,66],[5,14,180,600],6.),([200,576,197,27],[4,13,220,1000],6.)):
        t=time.perf_counter(); r=total(AN,SN,ys,grid,26,full=True); print(AN,SN,ys,round(r[0],6),'b',r[1],'kap',r[2]/SC,'g',r[3]/SC,[len(x['xs']) for x in r[5]],round(time.perf_counter()-t,1),flush=True)
