import sys; sys.path.insert(0,'wip')
from ub_tdep import V
import numpy as np, itertools
from multiprocessing import Pool
def job(a):
    l,m,tau,(p,q)=a
    return a,V(l,lambda t:m*t/(t+tau),lambda c,t:1+0*t,lambda x:1/np.sqrt(p+q*x),nu=200)
if __name__=='__main__':
    ls=(.145,.16,.175,.19); jobs=[(l,m,tau,s) for m,tau in ((.3,1),(.5,2),(.6,1.5),(.7,3),(.8,2),(1,4),(1,8)) for s in ((5.4,4.1),(5.4,6.),(5.4,9.),(7.,5.)) for l in ls]
    with Pool(11) as p: res=dict(p.map(job,jobs,chunksize=1))
    out=[]
    for k in {(j[1],j[2],j[3]) for j in jobs}: out.append((max(res[(l,)+k] for l in ls),k,[round(res[(l,)+k],4) for l in ls]))
    for o in sorted(out)[:8]: print(o)
    print('best per l',[(l,round(min(v for j,v in res.items() if j[0]==l),4)) for l in ls])
