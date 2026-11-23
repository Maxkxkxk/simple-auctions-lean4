"""Greedy local search over (e, prices) at N=80; deterministic seeds. Usage: python local80.py <seconds>"""
from slp_seed import *
from hi import up
from multiprocessing import Pool
import random
def adapt(lam,P0,e0,P,e):
    K=max(P)+e; out=np.zeros((len(P),K+1))
    for i,p in enumerate(P):
        j=min(range(len(P0)),key=lambda j:abs(P0[j]-p)); sh=p-P0[j] if p>100 else 0
        for k in range(1,lam.shape[1]):
            if lam[j,k]>0:
                kk=min(max(1,k+sh),p+e); out[i,kk]+=lam[j,k]*(1 if len(P)==len(P0) else .5)
    return out
def run(a):
    N,lam,P0,e0,P,e=a
    try:
        m=Model(N,e,P,max(P)+e); l=fix(m,adapt(lam,P0,e0,P,e)); l=fix(m,slp(m,l,iters=30,tr=0.2,verbose=False)); return 2+m.value(l)[0],P,e,l
    except Exception as ex: return 0,P,e,None
if __name__=='__main__':
    T0=time.time(); lim=float(sys.argv[1]); N=20;e=32;P=(1,4,8,72)
    m=Model(N,e,P,104); lam=fix(m,slp(m,init(m),200,verbose=False))
    for _ in range(2):
        N*=2;e*=2;P=tuple(2*p for p in P); lam=up(lam,2); m=Model(N,e,P,max(P)+e); lam=fix(m,slp(m,fix(m,lam),iters=60,tr=.1,verbose=False))
    best=2+m.value(lam)[0]; print(best,P,e,flush=True); rnd=random.Random(1)
    with Pool(14) as pool:
        while time.time()-T0<lim:
            jobs=[]
            for _ in range(14):
                Q=list(P); e2=e
                for _ in range(rnd.choice([1,1,2])):
                    u=rnd.random()
                    if u<.25: e2+=rnd.choice([-6,-4,-2,2,4,6])
                    elif u<.9:
                        i=rnd.randrange(len(Q)); Q[i]=max(0,Q[i]+(rnd.choice([-20,-10,-5,5,10,20]) if Q[i]>100 else rnd.choice([-4,-2,-1,1,2,4])))
                    elif len(Q)<6: Q.append(rnd.choice([6,10,24,45,rnd.randint(200,400)]))
                jobs.append((N,lam,P,e,tuple(sorted(set(Q))),e2))
            r=max(pool.map(run,jobs,chunksize=1),key=lambda r:r[0])
            if r[0]>best+1e-6:
                best,P,e,lam=r; print('%.6f'%best,P,e,round(time.time()-T0),flush=True); np.savez('wip/local80.npz',lam=lam,P=P,e=e,N=N)
