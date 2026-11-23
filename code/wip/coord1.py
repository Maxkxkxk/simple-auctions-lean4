from slp_seed import *
from multiprocessing import Pool
def run(a):
    N,e,P=a; K=max(P)+e
    m=Model(N,e,P,K)
    try: lam=fix(m,slp(m,init(m),iters=200,verbose=False)); f,sr,vi=m.value(lam)
    except Exception as ex: return 0,a,None
    return 2+f,a,m.srev_types(lam).round(4)
if __name__=='__main__':
    N=int(sys.argv[1]); e=int(sys.argv[2]); P=tuple(int(x) for x in sys.argv[3].split(',')); step=int(sys.argv[4])
    best=run((N,e,P)); print(best,flush=True)
    with Pool(8) as pool:
        while step>=1:
            cand=set()
            for d in (-2*step,-step,step,2*step):
                cand.add((N,e+d,P))
                for i in range(len(P)):
                    if P[i]+d>=0: cand.add((N,e,tuple(sorted(set(P[:i]+(P[i]+d,)+P[i+1:])))))
            import random
            rs=pool.map(run,sorted(cand))
            r=max(rs,key=lambda r:r[0])
            if r[0]>best[0]+1e-7: best=r; _,(N,e,P),_=r; print(step,best,flush=True)
            else: step//=2
