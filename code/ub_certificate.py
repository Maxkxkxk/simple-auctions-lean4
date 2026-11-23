"""EXACT integer certificate for REV <= (1+g) SREV + b BREV  (mirrors lean/SimpleAuctions/CertKernel.lean).
Scheme of Giannakopoulos-Hahn (arXiv 2609.35711v3, App. C.4.1) with a 4-exponential kernel and step-function
Jensen weights.  Every rounding is in the safe direction; no floating point.   Run: python ub_certificate.py
w(v)=sum AN/1000 e^{-v/(SN/10)}, lengths in units 1/D."""
import json
SC=10**12; D=10000
cdiv=lambda a,b:(a+b-1)//b
def lg(n): return n.bit_length()-1 if n>0 else 0
def kk(a,b): return lg(a//b)+7
def exp_lo(a,b):
    k=kk(a,b); x=cdiv(a*SC,b<<k)
    y=max(0,(SC+x*x//(2*SC))-(x+cdiv(x*x*x,6*SC*SC)))
    for _ in range(k): y=y*y//SC
    return y
def exp_up(a,b):
    k=kk(a,b); x=a*SC//(b<<k)
    y=max(0,SC+cdiv(x*x,2*SC)-x)
    for _ in range(k): y=cdiv(y*y,SC)
    return y
def psi_up(al,v): return SC+cdiv(v*exp_up(al,v),al)-SC*v//al
def wt(A,S,a,b,sign):
    n=2*D*S; d=(b-a)*10
    if d>n: return 0
    return A*exp_lo((a+b)*10,n)*(n+sign*d)//(1000*n)
def wA(AN,SN,a,b): return sum(wt(A,S,a,b,1) for A,S in zip(AN,SN))
def wB(AN,SN,a,b): return sum(wt(A,S,a,b,-1) for A,S in zip(AN,SN))
def avg_lo(AN,SN,b): return sum(A*S*max(0,SC-exp_up(b*10,D*S))*D//(10000*b) for A,S in zip(AN,SN))
def w1(AN,SN): return sum(cdiv(A*10*SC,S*1000) for A,S in zip(AN,SN))
def table(AN,SN,grid):
    """descending cells (a,wa,wb,av), starting below grid[-1]"""
    return [(grid[i],wA(AN,SN,grid[i],grid[i+1]),wB(AN,SN,grid[i],grid[i+1]),avg_lo(AN,SN,grid[i+1])) for i in range(len(grid)-2,-1,-1)]
def level(al,ts,te,kap,tab,Vmax,AN,SN):
    """smallest admissible g (SC-scaled) for this level; te=None for the last level"""
    G=cdiv(SC*al,2*Vmax); ka=cdiv(kap*D,al); pb=psi_up(al,Vmax); b=Vmax; g=0
    if te is None: g=ka+G
    for a,wa,wb,av in tab:
        if b<=ts: return g
        pa=psi_up(al,a); la=max(0,pa-wa); lb=max(0,pb-wb)
        G+=cdiv((b-a)*(2*la*b+la*a+lb*b+2*lb*a),6*a*b); pb=pa
        if te is None or a<te: g=max(g,ka+G-av)
        b=a
    if b<=ts: return g
    assert ts==0
    return max(g,ka+G+cdiv(w1(AN,SN)*b,D)-avg_lo(AN,SN,b))
def tangent(y0n,y0d,xs,Dx=1000):
    """Q=sum om, M2=sum om (p+q) with om<=SC delta e^{-y0 m}, delta=Dx(q-p)/(pq), m=(p+q)/(2Dx)"""
    Q=M=0
    for p,q in zip(xs,xs[1:]):
        om=exp_lo(y0n*(p+q),y0d*2*Dx)*Dx*(q-p)//(p*q); Q+=om; M+=om*(p+q)
    return Q,M
def tan_ok(y0n,y0d,xs,bn,kn,Dx=1000,UN=10**6):
    """check with b=bn/UN, kap=kn/UN: M<=bQ and exp(-(kap+b y0)) <= Q-y0(bQ-M)   (all scaled)"""
    Q,M=tangent(y0n,y0d,xs)                       # M/(2Dx)
    if M*UN>bn*Q*2*Dx: return False
    # exp(-(kn*y0d+bn*y0n)/(UN*y0d)) * SC <= Q - y0 (bQ-M)
    E=exp_up(kn*y0d+bn*y0n,UN*y0d)
    # multiply by UN*y0d*2Dx
    return E*UN*y0d*2*Dx + y0n*(bn*Q*2*Dx-M*UN) <= Q*UN*y0d*2*Dx
def best_tan(y0n,y0d,xs,Dx=1000,UN=10**6):
    Q,M=tangent(y0n,y0d,xs); bn=cdiv(M*UN,Q*2*Dx)
    lo,hi=-5*UN,5*UN
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid*y0d+bn*y0n>=0 and tan_ok(y0n,y0d,xs,bn,mid): hi=mid
        else: lo=mid
    return bn,hi
if __name__=='__main__':
    I=json.load(open('../data/ub_certificate.json'))
    AN,SN=I['AN'],I['SN']; tab=table(AN,SN,I['grid']); Vmax=I['grid'][-1]
    assert sum(AN)==1000
    assert all(tan_ok(t['y0'][0],t['y0'][1],t['xs'],t['b'],t['k']) for t in I['tan'])
    assert sum(A*S*t['k'] for A,S,t in zip(AN,SN,I['tan']))*SC<=I['kap']*10**4*10**6      # kappa (SC scaled)
    assert sum(A*t['b'] for A,t in zip(AN,I['tan']))<=I['b']*1000                         # b in 1e-6
    L=I['levels']; gs=[level(al,ts,(L[i+1][0] if i+1<len(L) else None),I['kap'],tab,Vmax,AN,SN) for i,(ts,al) in enumerate(L)]
    assert max(gs)<=I['g'] and all(L[i][1]<=L[i+1][1] and L[i][0]<L[i+1][0] for i in range(len(L)-1))
    from fractions import Fraction as Fr
    print('verified: REV <= (1+g) SREV + b BREV, g=',Fr(I['g'],SC),'b=',Fr(I['b'],10**6),'  c_ub=',float(1+Fr(I['g'],SC)+Fr(I['b'],10**6)))
