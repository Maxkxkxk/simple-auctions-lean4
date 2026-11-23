from slp2d import *
from multiprocessing import Pool
def run(a):
    N,types,e0,e1,e01,K,be=a
    t0=time.time()
    m=Model2(N,types,e0,e1,e01,K,beta=be); 
    caps=[min(p1+e1 if p1<INF else INF, p0+(e0 or e1) if p0<INF else INF) for p0,p1 in types]
    lam=fix(m,slp(m,init2(m,caps=caps),iters=80,verbose=False)); f,sr,vi=m.value(lam)
    return 1+be+f,round(sr,4),a,list(m.srev_types(lam).round(3)),round(time.time()-t0)
if __name__=='__main__':
    N=10; I=INF
    jobs=[(N,[(I,0),(I,2),(I,32)],I,15,I,60,1.0)]   # sanity: 2.16707
    for be in (.9,.8): jobs.append((N,[(I,0),(I,2),(I,32)],I,15,I,60,be))
    for q in (33,36,40,44):
        jobs.append((N,[(I,0),(I,2),(q,32)],0,15,I,60,1.0))
    for q2 in (8,12):
        for q in (36,40): jobs.append((N,[(I,0),(q2,2),(q,32)],0,15,I,60,1.0))
    # mixed bundling between clusters
    jobs.append((N,[(I,0),(I,2),(0,I)],30,15,I,60,1.0))
    jobs.append((N,[(I,0),(I,2),(I,32),(0,I),(20,I)],20,15,35,60,1.0))
    jobs.append((N,[(I,0),(I,2),(0,I),(2,I)],15,15,28,60,1.0))
    jobs.append((N,[(I,0),(I,2),(I,32),(0,I),(2,I)],15,15,28,60,1.0))
    jobs.append((N,[(I,0),(I,2),(I,32),(0,I)],30,15,43,60,1.0))
    # upsell: w0 lower fee, higher prices
    for d in (3,5,8):
        jobs.append((N,[(d,0),(2+d,2),(32+d,32)],15-d,15,I,60,1.0))
        jobs.append((N,[(I,0),(2+d,2),(32+d,32)],15-d,15,I,60,1.0))
    with Pool(14) as p:
        for r in p.imap_unordered(run,jobs): print(r,flush=True)
