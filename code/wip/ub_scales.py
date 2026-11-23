import sys; sys.path.insert(0,'wip')
import numpy as np, itertools, time
from ub_lp import lp,besteta
from multiprocessing import Pool
ts=np.concatenate([[.002,.05],np.geomspace(.2,40,14),[70,110,180]])
def ev(sc):
    s=np.array(sc,float); eta=1/np.sqrt(5.4+4.1*ts)
    for it in range(4):
        a,k,g,b,f=lp(eta,ts,s,nu=260); e2,g2=besteta(a,s,k,ts); eta=np.exp(.4*np.log(eta)+.6*np.log(e2))
    a,k,g,b,f=lp(eta,ts,s,nu=400)
    return 1+f,sc,np.round(a,4),round(b,4),round(k,4),round(g,4)
if __name__=='__main__':
    cand=[(1.,s2) for s2 in (18,22,26)]
    cand+=[(s0,1.3,s2) for s0 in (.3,.5) for s2 in (18,23,28)]+[(1.,10.,40.),(1.,15,60),(1.,18,100),(1.,20,150),(1.2,20,60),(.9,5,30)]
    cand+=[(.5,1.4,22,150),(.5,1.4,18,60),(.4,1.3,22,100),(.5,1.5,20,300),(.5,1.4,20,60,400)]
    with Pool(12) as p:
        for r in p.imap_unordered(ev,cand): print(r,flush=True)
