"""Pick few ER types and finite n; SLP over r with exact finite-n laws (float)."""
from erfam import *
from slpn import ModelN
class ERN(ER):
    def __init__(s,N,e,P,n):
        super().__init__(N,e,P); s.m=ModelN(N,e,P,[n]*len(P),s.K); s.n=n
    def opt(s,r,iters=40,tr=.1,TMf=6.2):
        N=s.N;n=len(s.P);TM=int(TMf*N); merit=lambda f,vi:f-20*max(0,vi-1)
        f,sr,vi=s.val(r);cur=merit(f,vi);small=0
        for it in range(iters):
            o=s.m.evaluate(s.lam(r),True);ZT=o['ZT'];t=np.arange(1,TM);t=t[(1+t/N)*ZT[t]>.8]
            D=np.zeros((len(t),n))
            for i in range(n):
                F=np.fft.rfft(np.diff(-np.append(o['DZ'][i],0)))
                b=np.zeros(s.m.L);b[:s.K+1]=s.B[i]
                g=np.fft.irfft(F*(np.fft.rfft(b)-b.sum()),s.m.L); D[:,i]=-(np.cumsum(g)[t-1])
            D*=(1+t/N)[:,None]; bb=1-(1+t/N)*ZT[t]; c=-((o['G']*s.B).sum(1)-1); ok=False
            for _ in range(8):
                st=tr*(r+.003)
                res=linprog(c,A_ub=np.vstack([D,np.ones(n)]),b_ub=np.append(bb,1-r.sum()),bounds=np.stack([np.maximum(-r,-st),st],1),method='highs-ds',options={'time_limit':15.})
                if res.x is None: tr*=.5;continue
                r2=np.maximum(r+res.x,0);f2,sr2,vi2=s.val(r2)
                if merit(f2,vi2)>cur-1e-12: ok=True;break
                tr*=.5
            if not ok:break
            gain=merit(f2,vi2)-cur;r,cur=r2,merit(f2,vi2);tr=min(tr*1.4,1.); small=small+1 if gain<2e-7 else 0
            if small>=3:break
        lo,hi=0.,1.
        for _ in range(30):
            mid=(lo+hi)/2
            if s.val(mid*r)[2]<=1: lo=mid
            else: hi=mid
        return lo*r,2+s.val(lo*r)[0]
def job2(a):
    N,e,P,n=a; t=time.time(); E=ER(N,e,P); r0=np.array([.3/ sum(1 for q in P if q<N) if p<N else .6/sum(1 for q in P if q>=N) for p in P])*.7
    r,v=E.opt(r0,80); vf=None
    if n: EN=ERN(N,e,P,n); r,vf=EN.opt(r,30,.02)
    return v,vf,N,e,n,P,[float(x) for x in r],round(time.time()-t,1)
if __name__=='__main__':
    from multiprocessing import Pool
    import random,json
    N=int(sys.argv[1]); nt=int(sys.argv[2]); rnd=random.Random(5); secs=float(sys.argv[3]); T0=time.time()
    P=tuple(sorted({max(1,round(N*x)) for x in np.geomspace(.02,.32,nt-1)}))+(round(3.2*N),); e=round(1.44*N); best=(0,)
    with Pool(10) as pool:
        cand=[(N,e,P,0)]
        while time.time()-T0<secs:
            for r in pool.imap_unordered(job2,cand,chunksize=1):
                if r[0]>best[0]+1e-7: best=r
            P,e=best[5],best[3]; cand=[]
            for _ in range(20):
                Q=list(P);e2=e
                for _ in range(rnd.choice([1,1,2,3])):
                    if rnd.random()<.15: e2+=round(N*rnd.choice([-.04,-.02,-.01,.01,.02,.04]))
                    else:
                        i=rnd.randrange(len(Q)); s=max(1,round(Q[i]*rnd.choice([.03,.07,.15,.3]))); Q[i]=max(1,Q[i]+rnd.choice([-s,s]))
                cand.append((N,e2,tuple(sorted(set(Q))),0))
    print(best[0],best[2:6],[round(x,4) for x in best[6]])
    for n in (512,1024,4096):
        r=job2(best[2:4]+(best[5],n)); print(n,r[0],r[1],r[-1],flush=True); json.dump(r,open('wip/design_N%d_t%d_n%d.json'%(N,nt,n),'w'))
