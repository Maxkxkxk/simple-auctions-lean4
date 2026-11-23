"""SLP search, 2D mechanism class: u = max(0, w0-e0, w1-e1, w0+w1-e01)
 w_i = sum_j (v_j - p_{i,j})^+ ; types = list of (p0,p1) (INF allowed). Everything on lattice 1/N.
 beta: size of concentrated group (sold at beta)."""
import numpy as np, sys, time
from slp_seed import cpois, slp, fix
INF=10**6
class Model2:
    def __init__(s,N,types,e0,e1,e01,K,kmin=1,beta=1.0,L2=None):
        s.N=N; s.P=types; s.K=K; s.kmin=kmin; s.beta=beta
        s.e=dict(a=e0,b=e1,c=e01)
        s.L=1<<int(np.ceil(np.log2(8*K))); s.k=np.arange(K+1)
        s.L2=L2 or 1<<int(np.ceil(np.log2(2.5*K)))
        L2=s.L2
        w0=np.arange(L2)[:,None]; w1=np.arange(L2)[None,:]
        # priority: c, b, a, none under ties
        us=np.stack([np.zeros((L2,L2)), (w0-e0)+0.0*w1, (w1-e1)+0.0*w0, (w0+w1-e01)+0.0])
        us=us+np.array([0,1e-6*min(e0,1e4)+1e-7,1e-6*min(e1,1e4)+2e-7,1e-6*min(e01,1e4)+3e-7])[:,None,None]
        reg=us.argmax(0)
        s.masks=[(reg==i).astype(float) for i in (1,2,3)]
        s.fm=[np.conj(np.fft.rfft2(m)) for m in s.masks]  # for correlation
        s.es=[e0 if e0<INF else 0,e1 if e1<INF else 0,e01 if e01<INF else 0]
        s.jump=[]; s.own=[]
        for (p0,p1) in types:
            a=np.maximum(s.k-p0,0) if p0<INF else np.zeros(K+1,int)
            b=np.maximum(s.k-p1,0) if p1<INF else np.zeros(K+1,int)
            o0=np.where(s.k>=p0,p0,0) if p0<INF else np.zeros(K+1)
            o1=np.where(s.k>=p1,p1,0) if p1<INF else np.zeros(K+1)
            s.jump.append((np.minimum(a,L2-1),np.minimum(b,L2-1))); s.own.append((o0,o1,o0+o1))
    def evaluate(s,lam,grad=False):
        N,K,L2=s.N,s.K,s.L2
        Z=cpois(lam.sum(0),s.L); ZT=np.cumsum(Z[::-1])[::-1]
        J=np.zeros((L2,L2)); PJ=[np.zeros((L2,L2)) for _ in range(3)]
        for i in range(len(s.P)):
            a,b=s.jump[i]
            np.add.at(J,(a,b),lam[i])
            for r in range(3): np.add.at(PJ[r],(a,b),lam[i]*s.own[i][r])
        fJ=np.fft.rfft2(J); fpi=np.exp(fJ-fJ[0,0].real)
        F=[None]*3; R=0; C=[None]*3; Rs=0
        for r in range(3):
            fp=fpi*(s.es[r]+np.fft.rfft2(PJ[r]))
            X=np.fft.irfft2(np.conj(fp)*np.conj(s.fm[r]),(L2,L2))
            Rs=Rs+X
            C[r]=np.fft.irfft2(np.conj(fpi)*np.conj(s.fm[r]),(L2,L2))
        out=dict(R=Rs[0,0]/N,ZT=ZT)
        if grad:
            G=np.zeros_like(lam)
            for i in range(len(s.P)):
                a,b=s.jump[i]
                G[i]=Rs[a,b]-Rs[0,0]
                for r in range(3): G[i]+=s.own[i][r]*C[r][a,b]
            out['G']=G/N
        return out
    def srev_types(s,lam):
        T=np.cumsum(lam[:,::-1],1)[:,::-1]
        return (T*s.k/s.N).max(1)
    def value(s,lam):
        o=s.evaluate(lam); sr=s.srev_types(lam).sum()
        t=np.arange(max(1,int(np.ceil((1-s.beta)*s.N-1e-9))),s.L//2)
        viol=((s.beta+t/s.N)*o['ZT'][t]).max()
        return o['R']-sr, sr, viol
def init2(m,w=None,caps=None,tot=0.6):
    n=len(m.P); lam=np.zeros((n,m.K+1)); w=w or [1.0]*n
    for i,(p0,p1) in enumerate(m.P):
        p=min(p0,p1); a=max(p,max(1,m.N//30)); b=min(m.K,caps[i] if caps else m.K)
        if a>=b: continue
        T=np.zeros(m.K+2); k=np.arange(a,b+1); T[a:b+1]=w[i]*tot/sum(w)*m.N/k; T[:a]=T[a]
        lam[i]=T[:-1]-T[1:]; lam[i,0]=0
    return fix(m,lam)
