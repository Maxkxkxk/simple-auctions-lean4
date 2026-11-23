import numpy as np
def levy(N,K):
    # Lomax tail 1/(1+t): lattice law p_k=P(t in ((k-1)h... ) use Z=k if t in [kh,(k+1)h)
    h=1/N; k=np.arange(K+1); p=1/(1+k*h)-1/(1+(k+1)*h)
    nu=np.zeros(K+1)
    for n in range(1,K+1):
        s=n*p[n]-np.dot(np.arange(1,n)*nu[1:n],p[n-1:0:-1]); nu[n]=s/(n*p[0])
    return nu,p
for N,K in ((50,4000),(100,6000),(200,12000)):
    nu,p=levy(N,K); Lam=-np.log(p[0]); tailbeyond=Lam-nu.sum()
    M=Lam-np.cumsum(nu)  # N*(t) for t in (k h,(k+1)h)
    phi=M-1+np.exp(-M); t=np.arange(K+1)/N
    I=phi.sum()/N; print(N,'min nu',nu[1:].min(),'int phi(N*) up to',K/N,'=',I,'+tail~',0.5*(N/K),'=>',1+I+.5*N/K, 'tN* at',[(round(x,1),round(x*M[int(x*N)],4)) for x in (.05,.1,.3,1,3,10,30)])
