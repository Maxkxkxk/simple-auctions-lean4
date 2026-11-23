# CHANGELOG

## Run 1, 2026-11-23, Phases 0 to 2 (Track 1)
- Read arXiv 2609.35711v3 and 2610.24973v2 (lower-bound sections and the 2.5 proof outline); abstracts of 2609.28873, 2606.13306, 2606.10169. Keyword scan of cs.GT to 18 Nov 2026.
- Built exact evaluator, baselines, search code (SLP), exact integer certificate.
- Found the equal-revenue-on-[p,p+e] structure; instance with exact sup ratio 2.194649949.
- Lean 4: `SimpleAuctions.lower_bound`, constant 2.1946, axioms propext / Classical.choice / Quot.sound only.
- Installed: elan 4.4.0, Lean v4.36.0, Mathlib v4.36.0 cache (`lean/.lake`, several GB).
- Not done: Track 2, financial bridge; M=3 mechanism LP unfinished; random multi-tariff search stopped early.

## Run 1 (continued), 2026-11-23
- A `leanchecker` call was interrupted by a session end; re-run module by module: all ten modules exit 0.
- Corrected wrong notes about GH Prop 4.6 (numbers of Claim 4.11 had been attributed to it). Their Claim 4.7 has slack; float lead for Track 2: 2.4819 (`code/wip/ub_slack.py`, `ub_opt.py`), no proof.
- `baselines.py`: GH baseline now also checks the bundle bound for t >= 2000.
- Files: LITERATURE, CONJECTURES, DEAD_ENDS, EXPERIMENTS, STATE, lean/AUDIT.md, code/baselines.py, code/wip/ub_*.py.

## Run 2, 2026-11-23, Track 2 (Phases 1 and 2)
- Re-derived GH's master inequality, layer-cake, Jensen and Laplace steps for finite distributions; structural tests logged in DEAD_ENDS 7 to 14 and EXPERIMENTS 11 to 18.
- Exact integer certificate 2.485776; Lean: `upper_bound` (2.486), `upper_bound_two`, `main_bracket`, `cStar_mem`; axioms propext / Classical.choice / Quot.sound. leanchecker exit 0 for the 13 new modules.
- Side result: exact iid lower bound 1.3249 (not in Lean; draft in `lean/wip/IID.lean` does not finish compiling).
- Not done: constant below 2.45, lottery test with 3 arrivals, LaTeX.
- Run 1 files untouched (`lower_bound`, `KernelRun` not changed).

## Run 2 (finish & polish), 2026-11-23
- Completed the unfinished Lean 4 i.i.d. lower bound by refactoring into modular `SimpleAuctions/IIDKernel.lean` (Horner packing, binary squaring, mod-`X^9` truncation, shifting loop: builds in 1.2 s) and `SimpleAuctions/IID.lean` (abstract `IIDData` record avoiding large-integer unfolding: builds in 21 s, zero warnings). Proved `SimpleAuctions.iid_lower_bound` with strictly stronger constant `> 1.316` (40 i.i.d. items), exported in `SimpleAuctions.lean` and audited in `Audit.lean` (`[propext, Classical.choice, Quot.sound]` only, `leanchecker` replay exit 0). Removed obsolete `lean/wip/IID.lean`.
- Produced arXiv-ready LaTeX paper (`paper/paper.tex`), GitHub Actions verification CI (`.github/workflows/verify.yml`), root `README.md`, and `.gitignore`.
- Cleaned up 9 verified dead-parent Python `multiprocessing.spawn` orphan processes left over from Run 1's `Pool` workers.
