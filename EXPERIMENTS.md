# EXPERIMENTS (Run 1). Float unless marked EXACT.

| # | Script | Setting | Best ratio |
|---|---|---|---|
| 1 | slp_seed.py / wip/scan1.py | free lambda per price, 3 to 6 prices, lattice N=10 / 20 / 40 | 2.1671 / 2.1817 / 2.1880 |
| 2 | wip/hi3.py | refinement to N=160, e=256, prices 8,32,64,576 | 2.192760 |
| 3 | observation | optimal lambda of #2 is (almost) exactly equal-revenue on [p, p+e] for every price p | - |
| 4 | erfam.py / wip/erhi.py | ER family, many prices, N=100 / 200 / 400 / 800 | 2.19316 / 2.19450 / 2.19512 / 2.19543 |
| 5 | finite_design.py | few prices: N=160 with 4 / 5 / 6 / 8 types | 2.19373 / 2.19399 / 2.19403 / 2.19410 |
| 6 | finite_design.py | N=400, 5 types (8,21,65,132,1280), e=580 | 2.19491 (Poisson), 2.19482 (n=512) |
| 7 | finite_design.py | N=480, 6 types | 2.19510 (Poisson) |
| 8 | make_instance.py + certificate.py | **EXACT**: N=400, n=256, d=24, R=(1067039,1128540,2129082,2782033,19107703) | 1+kappa = 2.194649949... |
| 9 | baselines.py | **EXACT / rigorous rounding** for Huang two-item (101/50) and GH (>= 2.190949) | - |
| 10 | wip/ub_slack.py, wip/ub_opt.py | Track 2, float only: GH Claim 4.7 scheme, own / re-optimised parameters (upper bounds) | 2.4904 / 2.4819 |

Extrapolation of #4: differences 0.00134, 0.00062, 0.00031 halve with N, so the limit is about 2.1957.

## Run 2 (Track 2 and iid)
| # | Script | What | Value |
|---|---|---|---|
| 11 | wip/ub_adversary.py | adversary against GH master bound (lower estimate of its value) | 2.325 |
| 12 | wip/ub_lp.py, ub_scales.py | one-parameter scheme, LP over kernels | 2.4837 (w(0)>=1), 4 exponentials 2.4844 |
| 13 | wip/ub_tdep2.py, ub_tdep3.py | two normalisations | 2.461 |
| 14 | wip/ub_adaptive.py, ub_adaptive_seq.py | profile-adaptive kernel (floor of the family) | 2.432 / 2.43 |
| 15 | ub_certificate.py | **EXACT integer certificate** g=0.118477349, b=1.367299 | 2.485776349 |
| 16 | lean `upper_bound` | **Lean** (g rounded to 0.1187) | 2.486 |
| 17 | iid_lower_bound.py | **EXACT**, iid equal-revenue + entry fee, n = 8 / 20 / 40 / 100 | 1.2360 / 1.3046 / 1.3168 / 1.3249 |
| 18 | wip/iid_slp.py | iid, Poisson limit, one entry fee + one price | 1.3292 (N=20) |
