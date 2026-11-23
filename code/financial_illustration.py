"""Phase-3 illustration (STYLISED numbers, not market data).  Run: python financial_illustration.py
One strategic bidder (e.g. a dominant primary dealer) with additive, independent private surplus values for the tranches
of an issue.  Exact rational SREV, BREV, best menu found, and the proven envelope
   max(SREV,BREV) <= REV <= 1.1187 SREV + 1.367299 BREV    (Lean: upper_bound_two)."""
import json
from fractions import Fraction as Fr
from itertools import product, combinations
import evaluator as E
I=json.load(open('../data/bond_illustration.json'))
D=[{Fr(v):Fr(p) for v,p in t['dist']} for t in I['tranches']]
for F in D: E.check_dist(F)
m=len(D); S=E.srev(D); B=E.brev(D)
ub=Fr(11187,10000)*S+Fr(1367299,1000000)*B
types=list(product(*[list(F.items()) for F in D]))
def revenue(menu):
    R=Fr(0)
    for c in types:
        pr=Fr(1)
        for _,p in c: pr*=p
        v=[x for x,_ in c]; R+=pr*max([(Fr(0),Fr(0))]+[(sum(v[i] for i in s)-p,p) for s,p in menu])[1]
    return R
# entry-fee (two-part tariff) menus with an optional separately priced core: prices on the support grid
best=(Fr(0),None)
subs=[s for k in range(1,m+1) for s in combinations(range(m),k)]
grid=[sorted(F) for F in D]
import itertools
fees=[Fr(x) for x in I['fees']]
for e in fees:
    for pr in itertools.product(*grid):
        menu=[(s,e+sum(pr[i] for i in s)) for s in subs]
        r=revenue(menu)
        if r>best[0]: best=(r,('entry fee',e,[str(x) for x in pr]))
print('tranches:',[t['name'] for t in I['tranches']])
print('SREV = %s = %.6f'%(S,S)); print('BREV = %s = %.6f'%(B,B))
print('best entry-fee menu: %s = %.6f  (%s)'%(best[0],best[0],best[1]))
print('proven upper envelope 1.1187 SREV + 1.367299 BREV = %.6f'%ub)
M=max(S,B); print('best menu / max(SREV,BREV) = %.6f ;  envelope / max = %.6f (worst case over all instances: 2.486)'%(best[0]/M,ub/M))
# second stylised case: 40 identical lots (proved in Lean: iid_lower_bound)
from iid_lower_bound import check
S,B,R=check(40,4,13,12,4,9); M=max(S,B)
print('\n40 identical lots, Pr[surplus>=k]=3/(10 max(k,4)), k<=13 bp:')
print('lot-by-lot pricing %.4f | one package %.4f | participation fee 9 + 4 per lot %.4f  -> +%.1f%% over the better simple format'%(S,B,R,float(100*(R/M-1))))
