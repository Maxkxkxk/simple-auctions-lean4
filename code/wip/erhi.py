from erfam import *
from multiprocessing import Pool
if __name__=='__main__':
    N=int(sys.argv[1]); jobs=[]
    cheap=sorted({max(1,round(N*x)) for x in list(np.geomspace(max(.01,1./N),.6,int(sys.argv[2])))})
    for e in (1.44,1.5,1.56):
        for q in ((3.1,3.2,3.3,3.4),(3.2,),(3.3,),(3.0,3.2,4,5)): jobs.append((N,round(e*N),cheap+[round(x*N) for x in q]))
    with Pool(8) as p:
        for r in sorted(p.imap_unordered(job,jobs,chunksize=1),key=lambda r:-r[0])[:4]: print(r)
