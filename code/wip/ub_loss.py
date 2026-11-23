"""NUMERICAL: where does the GH-type relaxation lose, on the adversary's instance?"""
import sys; sys.path.insert(0,'wip'); sys.path.insert(0,'.')
import numpy as np
from ub_slack import Lf,phi,kappa
from ub_adversary import run
from scipy.optimize import brentq, minimize
N=20;A=(1,2,3,5,8,12,20,40,60);K=160
v,o,prof,lam=run(N,A,K,1.)
print('adversary bound',v,prof)
T=np.cumsum(lam[:,::-1],1)[:,::-1]; k=np.arange(K+1); c=np.array(A)/N; r=T[range(len(A)),A]*c
h=1/N
def Y(i): # T_c(c+t) on cells t in (j,j+1]/N , u at cell right end
    return np.concatenate([T[i,A[i]+1:],np.zeros(A[i])])
Ys=np.array([Y(i) for i in range(len(A))]); t=(np.arange(Ys.shape[1])+.5)*h
Nt=Ys.sum(0); PhiT=(phi(Nt)).sum()*h; print('Phi true',PhiT, 'S_A',r.sum())
for nm,a,s,al in (('GH',(.75,.25),(1.,22.),lambda q:np.sqrt(6+4*q)),):
    l=brentq(lambda l:(r/(l*al(c/l))).sum()-1,1e-3,2); hc=1/(l*al(c/l)); print('l',l)
    PJ=sum((r[i]*hc[i]*phi(Ys[i]/(r[i]*hc[i]))).sum()*h for i in range(len(A)) if r[i]>1e-9); print('Phi after Jensen(GH weights)',PJ)
U=sum(r[i]/(c[i]+t) for i in range(len(A)))
hid=[1/((c[i]+t)*U) for i in range(len(A))]
PJ2=sum((r[i]*hid[i]*phi(Ys[i]/(r[i]*hid[i]))).sum()*h for i in range(len(A)) if r[i]>1e-9); print('Phi after Jensen(ideal U weights)',PJ2)
PhiU=phi(U).sum()*h; print('phi(U) untruncated on grid',PhiU)
print('t<->U,N:',[(round(x,2),round(a,2),round(b,2)) for x,a,b in zip(t[::8],U[::8],Nt[::8])])
for i in range(len(A)):
    ratio=T[i]*k/N/max(r[i],1e-12); print(A[i]/N,round(r[i],4),[round(x,2) for x in ratio[A[i]::10]][:14])
