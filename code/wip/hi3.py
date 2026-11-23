from slp_seed import *
from multiprocessing import Pool
def up(lam,f):  # refine grid by factor f : spread tail ER-like
    n,K1=lam.shape; K=K1-1; out=np.zeros((n,K*f+1))
    for i in range(n):
        T=np.cumsum(lam[i,::-1])[::-1]; T=np.append(T,0)
        for k in range(1,K+1):
            m=T[k]-T[k+1]
            if m<=0: continue
            # split mass on k*f..k*f+f-1 ER-like, only if below top
            js=np.arange(k*f,min(k*f+f,K*f+1)); w=1/js-1/(js+1.0)
            if T[k+1]<=1e-9 or m>0.5*T[k]: w=np.zeros(len(js)); w[0]=1   # atom (top)
            out[i,js]+=m*w/w.sum()
    return out
def run(a):
    N,e,P,f,t0=a; t=time.time()
    K=max(P)+e
    m=Model(N,e,P,K); lam=fix(m,slp(m,init(m),iters=200,verbose=False))
    for _ in range(f):
        N*=2; e*=2; P=tuple(2*p for p in P); K*=2
        lam=up(lam,2); m=Model(N,e,P,K); lam=fix(m,lam); lam=fix(m,slp(m,lam,iters=60,tr=0.1,verbose=False))
    v=m.value(lam); tag='N%d_e%d_P%s'%(N,e,'-'.join(map(str,P))); np.save('wip/hi_'+tag+'.npy',lam)
    print('%.6f sr %.4f'%(2+v[0],v[1]),tag,m.srev_types(lam).round(4),round(time.time()-t),flush=True)
    return 0
if __name__=='__main__':
    jobs=[]
    for e in (30,32,34):
        for P in ((1,4,8,72),(0,4,68),(0,2,4,8,70),(1,3,6,10,72),(0,2,5,12,72),(0,1,2,4,8,16,72)): jobs.append((20,e,P,3,0))
    with Pool(18) as p: p.map(run,jobs,chunksize=1)
