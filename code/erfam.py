"""ER family: type p has tail r_p/v on [p,p+e] (lattice h=1/N), atom at p+e.  SLP over r only."""
from slp_seed import *
class ER:
    def __init__(s,N,e,P):
        s.N,s.e,s.P=N,e,list(P); s.K=max(P)+e; s.m=Model(N,e,P,s.K)
        s.B=np.zeros((len(P),s.K+1))
        for i,p in enumerate(P):
            k=np.arange(p,p+e+1); T=N/k; s.B[i,p:p+e+1]=T-np.append(T[1:],0)
    def lam(s,r): return s.B*np.asarray(r)[:,None]
    def val(s,r):
        f,sr,vi=s.m.value(s.lam(r)); return f,sr,vi
    def opt(s,r,iters=40,tr=.1,TMf=6):
        N=s.N;n=len(s.P);TM=int(TMf*N)
        merit=lambda f,vi:f-20*max(0,vi-1)
        f,sr,vi=s.val(r);cur=merit(f,vi);small=0
        for it in range(iters):
            o=s.m.evaluate(s.lam(r),True);ZT=o['ZT'];t=np.arange(1,TM);t=t[(1+t/N)*ZT[t]>.8]
            F=np.fft.rfft(np.diff(-np.append(ZT,0)))
            D=np.zeros((len(t),n))
            for i in range(n):
                b=np.zeros(s.m.L);b[:s.K+1]=s.B[i]
                g=np.fft.irfft(F*(np.fft.rfft(b)-b.sum()),s.m.L); D[:,i]=-(np.cumsum(g)[t-1])
            D*=(1+t/N)[:,None]; bb=1-(1+t/N)*ZT[t]
            c=-((o['G']*s.B).sum(1)-1); ok=False
            for _ in range(8):
                st=tr*(r+.01)
                res=linprog(c,A_ub=np.vstack([D,np.ones(n)]),b_ub=np.append(bb,1-r.sum()),bounds=np.stack([np.maximum(-r,-st),st],1),method='highs-ds',options={'time_limit':15.})
                if res.x is None: tr*=.5;continue
                r2=np.maximum(r+res.x,0);f2,sr2,vi2=s.val(r2)
                if merit(f2,vi2)>cur-1e-12: ok=True;break
                tr*=.5
            if not ok:break
            gain=merit(f2,vi2)-cur;r,cur=r2,merit(f2,vi2);tr=min(tr*1.4,1.)
            small=small+1 if gain<1e-7 else 0
            if small>=3:break
        lo,hi=0.,1.
        for _ in range(40):
            mid=(lo+hi)/2
            if s.val(mid*r)[2]<=1: lo=mid
            else: hi=mid
        return lo*r,2+s.val(lo*r)[0]
def job(a):
    N,e,P=a; t=time.time(); E=ER(N,e,P); r=np.array([.3/ sum(1 for q in P if q<N) if p<N else .6/sum(1 for q in P if q>=N) for p in P]); r=r*.7
    r,v=E.opt(r,80); return v,N,e,[(p,round(float(x),4)) for p,x in zip(P,r) if x>1e-4],round(time.time()-t,1)
if __name__=='__main__':
    from multiprocessing import Pool
    N=int(sys.argv[1]); jobs=[]
    cheap=sorted({max(1,round(N*x)) for x in (.02,.03,.04,.05,.06,.08,.1,.12,.15,.18,.22,.27,.33,.4,.5,.6,.8,1.0,1.3)})
    for e in (1.4,1.5,1.6,1.7,1.8):
        for q in (3.2,3.4,3.6,3.8): jobs.append((N,round(e*N),cheap+[round(q*N)]))
        jobs.append((N,round(e*N),cheap+[round(x*N) for x in (1.8,2.4,2.8,3.2,3.4,3.6,3.8,4.2,5)]))
    with Pool(8) as p:
        for r in sorted(p.imap_unordered(job,jobs,chunksize=1),key=lambda r:-r[0])[:10]: print(r)
