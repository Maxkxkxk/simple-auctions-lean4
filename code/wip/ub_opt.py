"""NUMERICAL (float, no claim): optimise a_i, s_i, b in GH's Claim 4.7 scheme. Run: cd code; PYTHONPATH=wip python wip/ub_opt.py"""
from ub_slack import *
def f2(x):
    a=(np.exp(x[0]),np.exp(x[1])); s=(1.,np.exp(x[2])); b=x[3]
    return total(a,s,b) if 1<b<1.8 else 9
def f3(x):
    a=tuple(np.exp(x[:3])); s=(1.,np.exp(x[3]),np.exp(x[4])); b=x[5]
    return total(a,s,b) if 1<b<1.8 else 9
import time
t=time.time()
r=minimize(f2,[np.log(.75),np.log(.25),np.log(22.),4/3],method='Nelder-Mead',options=dict(maxfev=170,xatol=1e-3,fatol=1e-6)); x=r.x
print('2 exponentials:',r.fun,'a',np.exp(x[:2]),'s2',np.exp(x[2]),'b',x[3],round(time.time()-t),flush=True)
x0=[x[0]-.2,x[1],-1.2,np.log(4.),x[2],x[3]]
r=minimize(f3,x0,method='Nelder-Mead',options=dict(maxfev=260,xatol=1e-3,fatol=1e-6)); x=r.x
print('3 exponentials:',r.fun,'a',np.exp(x[:3]),'s',np.exp(x[3:5]),'b',x[5],round(time.time()-t),flush=True)
