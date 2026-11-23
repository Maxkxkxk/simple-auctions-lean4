from slp2d import *
from multiprocessing import Pool
I=INF
def run(a):
    N,e1,e0,e01,types,tag=a
    t=time.time(); K=max(min(x,y) for x,y in types)+max(e1,e0 if e0<I else 0)
    m=Model2(N,types,e0,e1,e01,K,L2=256)
    caps=[(p1+e1) if p1<I else p0+e0 for p0,p1 in types]
    w=[1 if p1<I else .15 for p0,p1 in types]
    lam=fix(m,slp(m,init2(m,w=w,caps=caps,tot=.8),iters=100,verbose=False)); f,sr,vi=m.value(lam)
    print('%.6f'%(2+f),round(sr,3),tag,e1,e0,e01,list(m.srev_types(lam).round(4)),round(time.time()-t),flush=True)
if __name__=='__main__':
    base=[(I,1),(I,2),(I,18)]; N=10; jobs=[]
    jobs.append((N,16,I,I,base,'base'))
    for e0 in (1,2,3,4,6): 
        jobs.append((N,16,e0,I if False else e0+16,base+[(0,I),(1,I)],'small'))
    for e0 in (24,32,48): jobs.append((N,16,e0,e0+16,base+[(0,I),(6,I)],'large'))
    for e0 in (24,32): jobs.append((N,16,e0,e0+10,base+[(0,I),(6,I),(20,I)],'largemix'))
    b36=[(I,1),(I,2),(I,36)]; jobs.append((N,16,I,I,b36,'b36'))
    for e0 in (1,2,4,8,24,32): jobs.append((N,16,e0,e0+16,b36+[(0,I),(1,I),(12,I)],'b36+'))
    for q in (40,44,48): jobs.append((N,16,8,I,[(I,1),(I,2),(q,36)],'nested'))
    for d in (8,12,14,17,18,20): jobs.append((N,16,d,I,[(I,1),(I,2),(40,36),(20,I)],'nest2'))
    with Pool(12) as p: p.map(run,jobs,chunksize=1)
