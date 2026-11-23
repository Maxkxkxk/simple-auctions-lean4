# Lean audit (2026-11-23)

Main theorem (`SimpleAuctions/Main.lean`):

```lean
theorem SimpleAuctions.lower_bound : ∃ (m : ℕ) (D : Fin m → FinDist),
    (2.1946 : ℝ) * max (SREV D) (BREV D) < REV D
```

`lake env lean Audit.lean` output:

```
'SimpleAuctions.lower_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
```

No `sorry`, `admit`, custom `axiom` or `native_decide` (so no `Lean.ofReduceBool`). The two heavy facts in `KernelRun.lean` use `decide +kernel`: they are evaluated by the Lean kernel itself, using its built-in GMP arithmetic for `Nat` literals. Trust base: Lean kernel incl. its GMP acceleration.

Replay check: `lake env leanchecker SimpleAuctions.<module>` exits with 0 for all ten modules (Basic, Reindex, Kernel, KernelRun, Spec, GroupAData, Instance, RealA, GroupB, Main; about 11 minutes in total, 2026-11-23). It re-sends every declaration of the module to the Lean kernel. It is the same kernel implementation, so not a second independent checker; Mathlib itself was not replayed.

## Anti-trivialisation checklist
- (a) **Type spaces non-empty.** `FinDist` requires `sum prob = 1`, so `size >= 1`; every profile space is non-empty.
- (b) **Independent product.** `prob D ω = ∏ i, (D i).prob (ω i)`; `prob_sum` proves total mass 1.
- (c) **All prices.** `FinDist.rev` and `BREV` are suprema over `{p : ℝ // 0 ≤ p}`. Boundedness is proven (`FinDist.rev_bdd`, `BREV_bdd`), so the real `⨆` is never the junk value 0.
- (d) **All types/deviations.** `Mechanism.ic` and `.ir` quantify over all nonnegative vectors `v, v' : ι → ℝ`, not only the support. Allocations lie in [0,1] (lotteries allowed). Payments are arbitrary reals.
- (e) **REV not junk.** `REV_bdd` bounds every mechanism's revenue by expected welfare; `le_REV` is proven from it. `Mechanism` is inhabited.
- (f) **Non-vacuous.** The statement is an unconditional existential. `Data` is inhabited by `Cert.theData` (concrete constants), all fields proven.
- (g) **No hidden hypothesis.** `lower_bound` has no arguments.
- (h) **Trivial case excluded.** `max(SREV,BREV) = 0` would make the statement weaker only if `REV > 0`; for a nonnegative finite distribution `BREV = 0` forces all values to be 0, hence `REV = 0` by the welfare bound, so the strict inequality cannot come from a degenerate instance. (Argument in prose, not formalised.)

## File map
- `Basic.lean` definitions; `Reindex.lean` item names do not matter (transport to `Fin m`).
- `Kernel.lean` numerals and computations; `KernelRun.lean` kernel evaluation; `Spec.lean` meaning of the computations (Kronecker substitution, digits, equal-revenue tails).
- `GroupAData.lean` abstract analysis of group A; `Instance.lean` concrete instance; `RealA.lean` real-valued statements (bundle bound, SREV, revenue).
- `GroupB.lean` concentrated group; `Main.lean` union, mechanism (IC, IR), limit argument, theorem.

# Upper bound (Run 2, 2026-11-23)

```lean
theorem SimpleAuctions.upper_bound : ∀ (m : ℕ) (D : Fin m → FinDist),
    REV D ≤ (2486 / 1000 : ℝ) * max (SREV D) (BREV D)
theorem SimpleAuctions.upper_bound_two (D : ι → FinDist) :
    REV D ≤ 1.1187 * SREV D + 1.367299 * BREV D
theorem SimpleAuctions.main_bracket : (lower_bound statement) ∧ (upper_bound statement)
theorem SimpleAuctions.cStar_mem : 2.1946 ≤ cStar ∧ cStar ≤ 2486 / 1000
```

`lake env lean Audit.lean`: all five theorems depend on `[propext, Classical.choice, Quot.sound]` only. No `sorry`, `admit`, `axiom`, `native_decide` in `SimpleAuctions/` (checked by grep). `lean/wip/` is not part of the build.

`leanchecker` exit 0 for Expect, Duality, LayerCake, Laplace, Abstract, AnalyticCert, CertKernel, CertData, CertRun, CertSound, CertLoop, UpperBound, Bracket (same kernel, Mathlib not replayed).

Anti-trivialisation, upper bound:
- Same definitions as the lower bound (`Basic.lean`, unchanged since Run 1). `REV` is a supremum over *all* IC+IR mechanisms with allocations in [0,1] and real payments; `Mechanism` is inhabited (`Reindex.lean`), and `REV_bdd` shows the supremum is not the junk value 0. So the statement bounds every single mechanism (`UB.rev_bound`).
- Quantifiers: all `m` (including 0), all `D`. Degenerate cases `SREV = 0` and `BREV = 0` are proved separately (`rev_le_zero`, `BREV_pos`).
- `SREV`, `BREV` are suprema over all real prices `p ≥ 0` (bounded families, `rev_bdd`, `BREV_bdd`).
- The certificate structure `UCert` is inhabited by `theCert`; all its fields are proved. The integer facts (`CertRun.lean`) are checked by `decide +kernel`; `CertSound.lean` and `CertLoop.lean` prove what they mean.
- Limitation: distributions have finite support. Extension to general distributions is not formalised.

# i.i.d. lower bound (Run 2 finish, 2026-11-23)

```lean
theorem SimpleAuctions.iid_lower_bound : ∃ (m : ℕ) (F : FinDist),
    (1.316 : ℝ) * max (SREV fun _ : Fin m => F) (BREV fun _ : Fin m => F) <
      REV fun _ : Fin m => F
```

`lake env lean Audit.lean` output (all six exported theorems):
- `'SimpleAuctions.lower_bound' depends on axioms: [propext, Classical.choice, Quot.sound]`
- `'SimpleAuctions.upper_bound' depends on axioms: [propext, Classical.choice, Quot.sound]`
- `'SimpleAuctions.upper_bound_two' depends on axioms: [propext, Classical.choice, Quot.sound]`
- `'SimpleAuctions.main_bracket' depends on axioms: [propext, Classical.choice, Quot.sound]`
- `'SimpleAuctions.cStar_mem' depends on axioms: [propext, Classical.choice, Quot.sound]`
- `'SimpleAuctions.iid_lower_bound' depends on axioms: [propext, Classical.choice, Quot.sound]`

Zero `sorry`, `admit`, custom `axiom`, `native_decide`, `unsafe`, or `implemented_by` across all 25 modules in `SimpleAuctions/`, `SimpleAuctions.lean`, and `Audit.lean` (verified by grep). `lean/wip/` is empty.

`leanchecker` replay: `lake env leanchecker SimpleAuctions.IIDKernel` and `lake env leanchecker SimpleAuctions.IID` both exit with 0 (bringing all 25 modules to exit 0 under `leanchecker`; same kernel implementation, Mathlib not replayed).

Anti-trivialisation, i.i.d. lower bound:
- Uses the exact same `FinDist`, `SREV`, `BREV`, `Mechanism`, `REV` definitions from `Basic.lean`.
- `m = 40` identical copies of `F` (supported on `{0, 4, …, 13}`, total probability `1`).
- `SREV ≤ 12` and `BREV ≤ 12.151` are proven globally over all real prices `p ≥ 0`.
- `mech` reuses `menuMech` (entry fee `9`, item price `4`, empty Group B), inheriting full IC and IR over all nonnegative real valuation vectors.
- Kernel checks (`IIDKernel.lean`) use `decide +kernel` only (Horner packing, binary squaring, modular truncation mod `X^9`, shifting loop; ~1.2 s). All algebraic reasoning in `IID.lean` runs over the abstract `IIDData` record so large integers are never unfolded during elaboration.
