"""Round the float design wip/r400.json to integers R_i and verify exactly; writes ../data/instance_A.json"""
import json, certificate as C
N,e,n,d,P,T0=400,580,256,24,(8,21,65,132,1280),2000
r=json.load(open('wip/r400.json')); R0=[x*N*2**d/n for x in r]
for sc in (1.0,.9999,.9998,.9996,.999):
    R=[int(x*sc) for x in R0[:-1]]+[int(R0[-1])]
    print(sc,R,flush=True); o=C.check(N,e,n,d,P,R,T0)
    if o[0] and o[1] and o[2]:
        json.dump(dict(N=N,e=e,n=n,d=d,P=list(P),R=R,T0=T0),open('../data/instance_A.json','w')); break
