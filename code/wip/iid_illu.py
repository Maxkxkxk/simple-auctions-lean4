import evaluator as E
from fractions import Fraction as Fr
from itertools import product
F={Fr(0):Fr(3,5),Fr(1):Fr(1,5),Fr(2):Fr(1,15),Fr(3):Fr(2,15)}
for m in (3,4,5,6,8):
    D=[F]*m; S=E.srev(D); B=E.brev(D); best=0
    for p in (0,1,2,3):
        # surplus/payment joint law for symmetric price p
        J={(Fr(0),Fr(0)):Fr(1)}
        for _ in range(m):
            K={}
            for (w,q),pr in J.items():
                for v,pv in F.items():
                    k=(w+max(v-p,0),q+(p if v>=p else 0)); K[k]=K.get(k,Fr(0))+pr*pv
            J=K
        for e in range(0,10):
            r=sum(pr*(e+q) for (w,q),pr in J.items() if w>=e) if e else sum(pr*q for (w,q),pr in J.items())
            if r>best: best=r;arg=(p,e)
    M=max(S,B); print(m,float(S),float(B),float(best),arg,float(best/M),float((Fr(11187,10000)*S+Fr(1367299,10**6)*B)/M))
