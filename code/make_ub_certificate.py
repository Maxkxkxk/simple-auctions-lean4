"""Build ../data/ub_certificate.json (float design of alpha levels, then exact integer data).  Run: PYTHONPATH=wip python make_ub_certificate.py"""
import sys; sys.path.insert(0,'.'); sys.path.insert(0,'wip')
import json
from ub_tune2 import *
A,S,y=[266,509,203,22],[5,14,220,2093],5.72
if __name__=='__main__':
    for fine,nlev,xr in ((.006,24,1.02),(.005,26,1.02),(.004,30,1.015)):
        g=mkgrid(fine); r=total(A,S,y,g,nlev,xr=xr,full=True)
        print(fine,nlev,len(g),len(r[4]),r[0],[len(t['xs']) for t in r[5]],'psi evals',sum(sum(1 for x in g if x>=l[0]) for l in r[4]))
    tot,b,kap,gg,L,tans,gs=r
    json.dump(dict(AN=A,SN=S,grid=g,levels=L,tan=tans,b=b,kap=kap,g=gg),open('../data/ub_certificate.json','w'))
