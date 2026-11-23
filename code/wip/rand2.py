from slp2d import *
from multiprocessing import Pool
import random
def run(seed):
    rnd=random.Random(seed); np.random.seed(seed); N=5; I=INF
    e1=rnd.choice([4,5,6,7,8,9,10,12,15]); mode=seed%4
    e0=I;e01=I
    if mode==0:   # single tariff, random prices
        types=[(I,p) for p in sorted(rnd.sample(range(0,30),5))]
    elif mode==1: # two nested tariffs; q<=p+e1
        e0=rnd.choice([0,1,2,3,4,6]); types=[]
        for p in rnd.sample(range(0,26),5): types.append((rnd.choice([I,p+e1-e0, p+rnd.randint(0,e1)]),p))
    elif mode==2: # clusters mixed
        e0=rnd.choice([2,4,6,8,10,16]); e01=rnd.randint(max(e0,e1),e0+e1)
        types=[(I,p) for p in rnd.sample(range(0,22),3)]+[(p,I) for p in rnd.sample(range(0,22),3)]
    else: # bundle size-like: two tariffs on separate price vectors arbitrary
        e0=rnd.choice([0,3,6,9,12]); types=[]
        for p in rnd.sample(range(0,22),5): types.append((max(0,p+rnd.randint(-6,8)),p))
    K=min(45,max(min(a,b) for a,b in types)+16)
    try:
        m=Model2(N,types,e0,e1,e01,K,L2=256)
        lam=init2(m,w=[rnd.random()**2+.01 for _ in types],caps=[min(p1+e1 if p1<I else I,p0+min(e0,20) if p0<I else I)+rnd.randint(0,3) for p0,p1 in types])
        lam=fix(m,slp(m,lam,iters=60,verbose=False)); f,sr,vi=m.value(lam)
    except Exception as ex: return 0,seed,str(ex)
    return 2+f,seed,mode,types,e0,e1,e01,sr,list(m.srev_types(lam).round(3))
if __name__=='__main__':
    res=[]
    with Pool(12) as p:
        for r in p.imap_unordered(run,range(int(sys.argv[1]),int(sys.argv[2]))):
            res.append(r)
            if r[0]>2.12: print(r,flush=True)
    res.sort(key=lambda r:-r[0])
    for md in range(4):
        for r in [r for r in res if len(r)>3 and r[2]==md][:5]: print(r)
    print([r for r in res if len(r)==3][:3], len(res))
