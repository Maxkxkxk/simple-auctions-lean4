"""Optimal (randomised) mechanism LP on a truncated Poisson state space (<=M arrivals), exchangeable ids inside a type.
Compares with best two-part tariff on the same truncated law."""
from slp_seed import *
from itertools import combinations_with_replacement
from math import factorial
from collections import Counter
def build(N,P,lam,M,thr=1e-4):
    jumps=[(i,k) for i in range(len(P)) for k in range(1,lam.shape[1]) if lam[i,k]>thr]
    rate={j:lam[j] for j in jumps}; Lam=sum(rate.values())
    states=[];pr=[]
    for m in range(M+1):
        for c in combinations_with_replacement(jumps,m):
            q=np.exp(-Lam)
            for j,n in Counter(c).items(): q*=rate[j]**n/factorial(n)
            states.append(c);pr.append(q)
    return jumps,states,np.array(pr)
def solve(N,P,lam,M,rounds=40,thr=1e-4,verbose=True):
    jumps,states,pr=build(N,P,lam,M,thr); nT=len(P); S=len(states)
    V=np.zeros((S,nT*M)); 
    for s,c in enumerate(states):
        for i in range(nT):
            vs=sorted([k for (t,k) in c if t==i],reverse=True)
            V[s,i*M:i*M+len(vs)]=vs
    V/=N
    xid=-np.ones((S,nT*M),int); nx=0
    for s in range(S):
        for q in range(nT*M):
            if V[s,q]>0: xid[s,q]=nx; nx+=1
    # vars: u_s (S), x (nx).  max sum pr (V.x - u)
    c=np.zeros(S+nx); c[:S]=pr
    for s in range(S):
        for q in range(nT*M):
            if xid[s,q]>=0: c[S+xid[s,q]]-=pr[s]*V[s,q]
    rows=set()
    for s in range(S):  # initial: drop-one & zero
        rows.add((s,0))
    from scipy.sparse import csr_matrix
    consts=[]
    # x sorted inside a type is wlog up to ties
    import scipy.sparse as sp
    cons=list(rows)
    idx={c:i for i,c in enumerate(states)}
    for s,cst in enumerate(states):
        for j in range(len(cst)):
            cons.append((s,idx[cst[:j]+cst[j+1:]]))
            for jj in jumps:
                if jj[0]==cst[j][0] and jj!=cst[j] and abs(jj[1]-cst[j][1])<=2:
                    cons.append((s,idx[tuple(sorted(cst[:j]+(jj,)+cst[j+1:]))]))
    cons=list(dict.fromkeys(cons))
    t0=time.time()
    for rd in range(rounds):
        R=[];C=[];D=[]
        for r,(s,s2) in enumerate(cons):  # u_s >= u_s2 + x_s2.(V_s-V_s2)  -> -u_s+u_s2+x_s2.(..)<=0
            if s==s2: continue
            R+= [r,r];C+=[s,s2];D+=[-1,1]
            for q in np.nonzero(V[s2])[0]:
                d=V[s,q]-V[s2,q]
                if d!=0: R.append(r);C.append(S+xid[s2,q]);D.append(d)
        nr=len(cons)
        for ss in range(S):
            for q in range(nT*M-1):
                if (q+1)%M and xid[ss,q+1]>=0: R+=[nr,nr];C+=[S+xid[ss,q],S+xid[ss,q+1]];D+=[-1,1];nr+=1
        A=csr_matrix((D,(R,C)),shape=(nr,S+nx))
        res=linprog(c,A_ub=A,b_ub=np.zeros(nr),bounds=[(0,None)]*S+[(0,1)]*nx,method='highs-ipm')
        u=res.x[:S]; X=np.zeros((S,nT*M)); msk=xid>=0; X[msk]=res.x[S:][xid[msk]]
        val=V@X.T - ((V*X).sum(1)-u)[None,:]   # utility of s reporting s2
        viol=val-u[:,None]
        # note: if x not sorted inside type best matching could be higher; check by sorting
        Xs=X.reshape(S,nT,M); Xs=-np.sort(-Xs,axis=2); val2=V@Xs.reshape(S,-1).T-((V*X).sum(1)-u)[None,:]
        top=np.argsort(-viol,axis=1)[:,:4]; new=0; have=set(cons)
        for s in range(S):
            for s2 in top[s]:
                if viol[s,s2]>1e-6 and (s,s2) not in have: cons.append((s,s2));new+=1
        if verbose: print('round',rd,'obj',-res.fun,'maxviol',viol.max(),'sortviol',(val2-u[:,None]).max(),'ncons',len(cons),round(time.time()-t0),flush=True)
        if new==0: break
    return jumps,states,pr,V,X,u,-res.fun
def tariff(N,P,states,pr,V,M,e,ps):
    off=np.repeat(np.array(ps),M)
    sur=np.maximum(V-off,0).sum(1); buy=((V>=off-1e-12)&(V>0))
    return (pr*((sur>=e-1e-12)*(e+(buy*off).sum(1)))).sum()
if __name__=='__main__':
    N=5; P=(0,1,16); e=8; M=int(sys.argv[1]); K=24
    m=Model(N,e,P,K); lam=fix(m,slp(m,init(m),300,verbose=False)); print(m.value(lam)); print(m.srev_types(lam)); 
    for i in range(3): print(P[i],{k:round(float(lam[i,k]),4) for k in range(1,K+1) if lam[i,k]>1e-4})
    jumps,states,pr,V,X,u,obj=solve(N,P,lam,M)
    print('OPT',obj,'mass',pr.sum(),'tariff',tariff(N,P,states,pr,V,M,e/N,[p/N for p in P]))
    best=max((tariff(N,P,states,pr,V,M,ee/N,(a/N,b/N,cc/N)),ee,a,b,cc) for ee in range(1,13) for a in range(0,3) for b in range(0,4) for cc in range(10,22)); print('best tariff',best)
    np.savez('wip/optmech.npz',V=V,X=X,u=u,pr=pr)
    o=np.argsort(-pr)[:150]
    for s in o: print(states[s],'%.4f'%pr[s],'u %.3f pay %.3f'%(u[s],(V[s]*X[s]).sum()-u[s]),'x',X[s][V[s]>0].round(2))
