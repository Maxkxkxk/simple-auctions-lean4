from slp_seed import *
from multiprocessing import Pool
def run(a):
    N,e,K,P=a
    m=Model(N,e,P,K); lam=fix(m,slp(m,init(m),iters=150,verbose=False)); f,sr,vi=m.value(lam)
    return 2+f,sr,a,m.srev_types(lam).round(3)
if __name__=='__main__':
    N=20
    jobs=[(N,e,p3+e,(0,p2,p3)) for e in (24,30,36) for p2 in (3,4,6) for p3 in (40,60,80,100)]
    jobs+=[(40,60,192,(0,8,132)),(40,60,192,(0,8,30,132)),(40,60,192,(0,4,8,20,100,132))]
    with Pool(12) as p:
        for r in sorted(p.imap_unordered(run,jobs),key=lambda r:r[0]): print(r,flush=True)
