# DEAD ENDS (all numbers float; search evidence, not proofs)

Setting: concentrated group B normalised to 1, group A in the Poisson (many small items) limit, ratio = 2 + REV_A - SREV_A, subject to (1+t) Pr[Z_A >= t] <= 1 and SREV_A <= 1. Lattice 1/N.

1. **Smaller concentrated group** (`beta < 1`, so A may contain small sure-ish items): N=10: 2.167 (beta=1) vs 2.0876 (0.9) vs 2.005 (0.8). `beta=1` is a kink optimum: shrinking B by x loses x directly and A gains less than x.
2. **Selling expensive items without entry fee** (second tariff with e0=0, N=10, code/wip/scan2.py): 2.107 to 2.157, always below 2.167. Buyers escape the entry fee.
3. **Two nested tariffs / mixed bundling of two clusters / arbitrary price vectors** (`code/wip/rand2.py`, N=5, random starts up to seed 558 of 1200; the run was stopped early; only results above 2.12 were logged, 227 of them): best per mode: single 2.1554, nested 2.1528, two clusters 2.1516, arbitrary price vectors 2.1516. Good results put (almost) zero weight on the second tariff. `code/wip/clus.py` (N=10): adding small- or large-scale clusters to the optimum gives weight 0, or at most 0.001 with a lower ratio.
4. **Few large items instead of many small ones** (`code/slpn.py`, n = 1, 2, 4 copies per type, N=20): 2.1665 to 2.1784 vs 2.1812 in the Poisson limit. Caveat: fixed prices (0,4,64)/20.
5. **Non-tariff menus / lotteries for the current A** (`code/wip/optmech.py`, N=5, at most 2 arrivals): the LP over all randomised mechanisms returns exactly the best two-part tariff (0.491833). With at most 3 arrivals two runs were cut off (last LP value 0.8056, still decreasing, after 5 rounds): undecided.
6. **Convex pricing g(W) of the surplus**: revenue is linear in g over convex functions with slope in [0,1], so a single entry fee (extreme point) is always optimal (this one is a proof).

## Upper-bound (Track 2) structural tests, Run 2 (all float unless marked proof)

7. **Other caps (theta-scaling, item-adaptive): dead, proof.** In GH's master inequality, raising `c_j` beyond `E min(X_j,c_j)=r_j` costs `F̄_j(c_j)` per unit and gains at most `F̄_j(c_j)`; lowering it saves nothing on `max{r_j, E Z_j}` and only enlarges the excess. So `E min(X_j,c_j)=r_j` is optimal for that inequality.
8. **Other coordinate routings: no gain, proof for equal-revenue items.** For items with tail r/v on [a,∞) a line contributes `r + r∫(1-w)(v-a)/v²`, minimised pointwise by routing to the largest excess.
9. **Value of the master inequality itself:** an adversary (SLP over instances, `wip/ub_adversary.py`) reaches 2.325 (N=20), i.e. the inequality cannot give less than about 2.33 to 2.35. Most of the loss up to 2.5 is in the later relaxation, not in the master inequality.
10. **No caps:** about 2.65. **Fixed kernel + uniform Jensen weights:** 2.600 (`wip/ub_fixedw.py`).
11. **One-parameter GH scheme, any number of exponentials, optimal weights:** 2.482 to 2.484 (`wip/ub_lp.py`). More exponentials than 4 gain less than 0.001. This is the scheme we proved (2.486).
12. **Two normalisations (t-dependent Jensen weights, own kernel per l-range):** about 2.461 (`wip/ub_tdep*.py`). Needs a 2D certificate; not done.
13. **Kernel fully adapted to the profile (r_c):** about 2.432 (`wip/ub_adaptive.py`). With the sorted/telescoping split of φ instead of the chord bound (`wip/ub_adaptive_seq.py`): again about 2.43. So the Laplace + Lagrangian family is numerically stuck near 2.43: the remaining gap to 2.34 is a genuine duality gap (convex maximisation).
14. **Sorting items by cap and telescoping with R/u instead of U:** 2.87 (`wip/ub_seq.py`), useless because R/u overestimates the expensive items at small u.

Technical traps:
- Lean: never let huge closed `Nat` constants (2^30744, generating functions) appear in unification goals; keep them `irreducible` and prove everything over the abstract `Data` structure.
- Python big ints are too slow for 60-Mbit products; use `gmpy2`.
- `PYTHONPATH` needs `;` on Windows (`.:wip` is silently wrong).
- In this Mathlib, `.neg/.add` of `HasDerivAt` speak about `-f`, `f+g`; use `.fun_neg/.fun_add`. `fun_prop` inside `ContinuousOn.div` with powers can time out: name the pieces.
- Structural recursion with a literal like 10^12 inside `Nat.sub` makes the kernel unfold it; pass such constants as arguments.
- On this Windows machine use Pool workers only in background jobs with logs; killing the parent leaves idle orphan workers.
