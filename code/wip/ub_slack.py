"""NUMERICAL (float, no claim): slack in GH Claim 4.7 (arXiv 2609.35711v3).
REV <= S + sum_j lam_j int_{c_j} phi(Fbar_j/lam_j); w(u)=sum a_i exp(-u/s_i) (units of l), 1/lam_j = l*alpha(c_j/l)/r_j.
=> REV <= (1+g) S + b B with  kappa(b)=sup_y sum a_i s_i L(y/s_i) - b y,
   g = sup_t [ int_t^oo (phi(al/v)/al - w(v)/v)^+ dv + kappa/al - (1/t) int_0^t w ], al=alpha(t) free."""
import numpy as np, sys
from scipy.special import exp1
from scipy.optimize import minimize, minimize_scalar
def Lf(y):   # -ln int_1^oo e^{-yx}/x^2 dx ,  E_2(y)=e^{-y}-y E_1(y)
    y=np.asarray(y,float)
    return -np.log(np.where(y>30, np.exp(-np.minimum(y,600))/(y+2), np.exp(-y)-y*exp1(np.maximum(y,1e-300))))
V=np.concatenate([np.geomspace(1e-5,400,6001)]); LV=np.log(V)
def phi(x): return np.where(x<1e-4,x*x/2-x**3/6,x-1+np.exp(-np.minimum(x,700)))
def kappa(a,s,b):
    y=np.geomspace(.01,400,4000); f=sum(ai*si*Lf(y/si) for ai,si in zip(a,s))-b*y
    i=f.argmax(); r=minimize_scalar(lambda q:-(sum(ai*si*Lf(q/si) for ai,si in zip(a,s))-b*q),bounds=(y[max(i-1,0)],y[min(i+1,3999)]),method='bounded',options=dict(xatol=1e-10)); return -r.fun
def G(t,al,a,s,k):
    v=np.geomspace(max(t,1e-6),4000,4000); wv=sum(ai*np.exp(-v/si) for ai,si in zip(a,s))
    h=np.maximum(phi(al/v)/al-wv/v,0); lv=np.log(v)
    I=np.sum((h*v)[1:]*np.diff(lv)+(h*v)[:-1]*np.diff(lv))/2
    av=sum(ai*(si*(1-np.exp(-t/si))/t if t>1e-9 else 1) for ai,si in zip(a,s))
    return I+k/al-av
TS=np.array([1e-4,.02,.05,.1,.2,.3,.45,.6,.8,1,1.25,1.5,1.8,2.2,2.6,3,4,6,10,20,60])
def total(a,s,b,alpha=None,full=False):
    k=kappa(a,s,b); g=[]
    for t in TS:
        if alpha is not None: g.append(G(t,alpha(t),a,s,k))
        else: g.append(minimize_scalar(lambda la:G(t,np.exp(la),a,s,k),bounds=(np.log(.3),np.log(60)),method='bounded').fun)
    m=max(g); 
    return (1+m+b,m,k,np.round(g,4)) if full else 1+m+b
if __name__=='__main__':
    for S2 in (22,2*2**.5):
        f=lambda y:.75*Lf(y)+.25*S2*Lf(y/S2)-4/3*y
        print('scale',S2,'value at y=6:',f(6.0))
    a,s=(.75,.25),(1.,22.)
    print('GH parameters, alpha=sqrt(6+4t):',total(a,s,4/3,lambda t:np.sqrt(6+4*t),True))
    print('GH parameters, best alpha(t):  ',total(a,s,4/3,None,True))
