"""Float design step (deterministic): weights r for N=400, e=580, prices (8,21,65,132,1280), n=256.
Run:  python design_r400.py   (writes wip/r400.json; then python make_instance.py rounds + verifies exactly)"""
import json
from finite_design import *
N,e,n,P=400,580,256,(8,21,65,132,1280)
E=ER(N,e,P); r,v=E.opt(np.array([.04,.04,.08,.1,.7]),80); print('Poisson limit (float)',v)
r,v=ERN(N,e,P,n).opt(r,40,.02); print('n=256 (float)',v)
json.dump(list(map(float,r)),open('wip/r400.json','w'))
