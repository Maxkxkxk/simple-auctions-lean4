# LITERATURE (checked 2026-11-23 on arXiv)

## Lower bounds
**Giannakopoulos & Hahn, arXiv:2609.35711v3 (v3: 3 Nov 2026), Theorem D.1: REV > 2.19 max{SREV,BREV}.**
- Exact published constant: `E + S_A > 2 S_A + 22/25`, `S_A = 3/10+19/20+27/8 = 4.625`, i.e. `9.25+0.88 over 4.625 = 2.190270...`.
- Structure: group A = three groups of N_a = 1000 iid items, each *equal revenue* (ER):
  free items (price 0, values (5+k)/40, k=0..275, N_a*tail 3/(10 v)), cheap (price 1, values (40+k)/40, k=0..280, 19/(20 v)), expensive (price 15, values (600+k)/40, k=0..280, 27/(8 v)). Entry fee 7.
  Group B = N_b^2 iid ER items whose sum concentrates at S_A while SREV_B = S_A/H_{N_b} -> 0; sold as a bundle at about S_A.
  Key condition: `sup_p p Pr[S_A + Z_A >= p] = S_A`, so adding A to B does not increase bundle revenue; then BREV -> S_A and SREV -> S_A.
- The proof's rounding is loose: our directed-rounding evaluation of *their* instance gives >= 2.190949 (`code/baselines.py`).
- They state that AI-assisted exploration found slightly better, unpublished instances, so 2.19 "should be improvable". No numbers are given.

**Huang, arXiv:2610.24973v2 (v2: 3 Nov 2026), Thm 5.5: > 2.11.**
- "Completing a seed" lemma: if `c >= SREV_seed` and `REV(c + S_seed) = c`, the ratio `2 + (REV_seed - SREV_seed)/c` is approached (concentrated group at c + seed + rare expensive item filling SREV up to c).
- Two-item seed: 2.02 exactly (reproduced: 101/50). Compound-Poisson seed, c = 1/4, price 1/2: proved 2.11, actual value about 2.1265 (reproduced numerically).

**Rubinstein (ITCS 2016): 2.** Concentrated group + rare items (limit; finite exact instances in `code/baselines.py`).

## Upper bounds
- **GH v3 Prop 4.6: REV <= 7/6 SREV + 4/3 BREV <= 2.5 max.** Mechanics: master inequality (8): for caps c_j, Y_j=(X_j-c_j)^+, Z_j=min(X_j,c_j), `REV <= sum_j max{r_j, E Z_j} + E[sum Y_j - max_j Y_j]`; caps chosen with `E Z_j = r_j`; then Jensen separates items; tails bounded with `Pr[X_j>u] <= r_j/u` and `sum_j Pr[X_j>u] <= -ln(1-B/u)`. (Correction of an earlier version of this file: the values 2.4996..2.4999 belong to Claim 4.11, the iid bound `REV <= 2.5 BREV`, not to Prop 4.6.) Claim 4.7, which proves Prop 4.6, has visible slack: weights `lam_j = r_j/sqrt(6 l^2 + 4 l c_j)`, `w(u) = 3/4 e^{-u/l} + 1/4 e^{-u/(22 l)}`; they bound a maximum by 1.92 that is at most 1.914124, and show `R+K <= 0.165` where 1/6 is needed. Their own printed 0.165 in place of 1/6 would already read `1.165 SREV + 4/3 BREV`, i.e. about 2.4984 (not checked by us).
- GH iid: 2.15 max. GH v3 (p. 3) call 1.1877 [Ma & Simchi-Levi, Thm 3.4, two items] "the currently best known lower-bound" for iid items. Our exact 40-item instance gives 1.3168 (see `code/iid_lower_bound.py`).
- GH say (App. C.4) that AI agents found slightly better analyses than 2.5 which they did not include; so an improved constant in their scheme is expected by the authors.
- Huang: `REV <= SREV + pi^2/6 BREV` (2.645); regular items: `REV <= 2 BREV`.
- Cai et al. arXiv:2609.28873v1: 3.52.
- Newer cs.GT listing up to 18 Nov 2026 scanned by keyword: arXiv:2611.06101 ("Optimality of Separate Selling with Finite-Support Values", 4 Nov 2026) characterises when separate selling is exactly optimal (abstract read; it states no bound for c*). No newer bound for this constant found.

## Lean infrastructure
- EconCSLib (arXiv:2606.13306; github.com/nikhgarg/EconCSLib, Lean v4.30.0-rc2, about 986k lines). We use a small standalone development on Mathlib v4.36.0 instead: the statement only needs finite sums and suprema, and the dependency would add a huge build and an older toolchain. Definitions are in `lean/SimpleAuctions/Basic.lean`.
