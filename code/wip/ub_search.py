import sys; sys.path.insert(0,'.'); sys.path.insert(0,'wip')
import random, time
from ub_tune2 import *
from multiprocessing import Pool
grid=mkgrid(.006)
def f(x):
    try: return total(list(x[0]),list(x[1]),x[2],grid,24,xr=1.02),x
    except Exception as ex: return 9,x
if __name__=='__main__':
    rnd=random.Random(3); T0=time.time(); lim=float(sys.argv[1])
    best=min((f((a,s,y)) for a,s,y in (((272,504,199,25),(5,14,220,1500),5.6),((248,521,186,45),(5,14,200,800),5.7),((257,510,167,66),(5,14,180,600),5.7),((728,197,60,15),(10,150,500,2000),5.8))),key=lambda r:r[0]); print(best,flush=True)
    with Pool(12) as p:
        while time.time()-T0<lim:
            c=[]
            for _ in range(48):
                a,s,y=list(best[1][0]),list(best[1][1]),best[1][2]
                for _ in range(rnd.choice([1,2,3])):
                    u=rnd.random()
                    if u<.5:
                        i,j=rnd.sample(range(4),2); d=rnd.choice([1,2,4,8,16]); 
                        if a[i]>d: a[i]-=d;a[j]+=d
                    elif u<.8:
                        i=rnd.randrange(4); s[i]=max(1,s[i]+rnd.choice([-1,1])*max(1,round(s[i]*rnd.choice([.02,.06,.15]))))
                    else: y=round(y+rnd.choice([-.2,-.1,-.05,-.02,.02,.05,.1,.2]),3)
                c.append((tuple(a),tuple(s),y))
            r=min(p.map(f,c,chunksize=2),key=lambda r:r[0])
            if r[0]<best[0]-1e-7: best=r; print(round(best[0],6),best[1],round(time.time()-T0),flush=True)
