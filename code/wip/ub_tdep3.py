import sys; sys.path.insert(0,'wip')
from ub_tdep import V
import numpy as np
from multiprocessing import Pool
SH={'s9':lambda x:1/np.sqrt(5.4+9*x),'s14':lambda x:1/np.sqrt(5.4+14*x),'s20':lambda x:1/np.sqrt(5.4+20*x),'l1':lambda x:1/(2.3+x),'l.5':lambda x:1/(2.3+.5*x),'p75':lambda x:1/(2.3+1.3*x**.75),'p75b':lambda x:1/(2.3+2*x**.75),'q':lambda x:1/np.sqrt(5.4+6*x+.15*x*x),'q2':lambda x:1/np.sqrt(5.4+9*x+.5*x*x)}
def job(a):
    l,m,tau,k=a
    return a,V(l,lambda t:m*t/(t+tau),lambda c,t:1+0*t,SH[k],nu=180)
if __name__=='__main__':
    ls=(.08,.105,.125,.145,.17); jobs=[(l,m,tau,k) for m,tau in ((.8,2),(.9,1.),(1,2)) for k in SH for l in ls]
    with Pool(11) as p: res=dict(p.map(job,jobs,chunksize=1))
    out=[]
    for k in {(j[1],j[2],j[3]) for j in jobs}: out.append((max(res[(l,)+k] for l in ls),k,[round(res[(l,)+k],4) for l in ls]))
    for o in sorted(out)[:14]: print(o)
