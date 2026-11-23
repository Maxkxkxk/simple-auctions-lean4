# Machine-Verified Bounds for Simple Multi-Item Auctions (Lean 4)

[![Lean 4 & Exact Python Verification](https://github.com/Maxkxkxk/simple-auctions-lean4/actions/workflows/verify.yml/badge.svg)](https://github.com/Maxkxkxk/simple-auctions-lean4/actions/workflows/verify.yml)

**Author:** Max Sievert (`sievert-max@gmx.de`)  
**Date:** November 23, 2026

Formal **Lean 4** (zero `sorry`, zero `native_decide`, standard axioms `[propext, Classical.choice, Quot.sound]` only) and exact rational Python proofs tightening both sides of the universal approximation ratio $c^*$ between optimal multi-item mechanisms ($\mathrm{REV}$) and the better of separate selling ($\mathrm{SREV}$) and grand bundling ($\mathrm{BREV}$) for a single additive buyer with independent finitely supported item valuations.

## Main Verified Theorems

| Result | Lean 4 Theorem | Value | Prior Published Bound |
| :--- | :--- | :--- | :--- |
| **Universal upper bound (two-parameter)** | `SimpleAuctions.upper_bound_two` | $\mathrm{REV}(D) \le 1.1187\,\mathrm{SREV}(D) + 1.367299\,\mathrm{BREV}(D)$ | $\frac{7}{6}\mathrm{SREV} + \frac{4}{3}\mathrm{BREV}$ |
| **Universal upper bound** | `SimpleAuctions.upper_bound` | $\forall m, D,\; \mathrm{REV}(D) \le 2.486 \max\{\mathrm{SREV}(D), \mathrm{BREV}(D)\}$ | $2.50$ (Giannakopoulos–Hahn, `arXiv:2609.35711v3`) |
| **General lower bound** | `SimpleAuctions.lower_bound` | $\exists m, D,\; \mathrm{REV}(D) > 2.1946 \max\{\mathrm{SREV}(D), \mathrm{BREV}(D)\}$ | $2.19027$ (`> 2.19`, `arXiv:2609.35711v3`) |
| **Combined bracket** | `SimpleAuctions.main_bracket`, `SimpleAuctions.cStar_mem` | $2.1946 \le c^* \le 2.486$ | $[2.19,\; 2.50]$ |
| **i.i.d. lower bound (40 items)** | `SimpleAuctions.iid_lower_bound` | $\exists m, F,\; \mathrm{REV}(F^{\otimes m}) > 1.316 \max\{\mathrm{SREV}, \mathrm{BREV}\}$ | $1.1877$ (Ma–Simchi-Levi; cited in `arXiv:2609.35711v3`) |
| **i.i.d. lower bound (100 items, exact Python)** | `code/iid_lower_bound.py` | $\ge 1.324907$ (exact `fractions.Fraction`) | $1.1877$ |

## Reproducibility

### 1. Verify all Lean 4 Proofs (Lean `v4.36.0`, Mathlib `v4.36.0`)

```bash
cd lean
lake exe cache get
lake build
lake env lean Audit.lean
```

Expected output of `lake env lean Audit.lean`: every exported theorem depends **only** on `[propext, Classical.choice, Quot.sound]` (no `sorryAx`, no `Lean.ofReduceBool`). Optional full kernel replay: `lake env leanchecker SimpleAuctions.<Module>`.

### 2. Run Exact Rational / Integer Python Certificates

```bash
cd code
python evaluator.py
python ub_certificate.py
python iid_lower_bound.py
python financial_illustration.py
# Optional (fast with `pip install gmpy2`):
python certificate.py ../data/instance_A.json
python baselines.py
```

## Repository Structure

- `lean/SimpleAuctions/` — Complete Lean 4 formalization (25 modules, all passing `leanchecker`):
  - `Basic.lean`, `Reindex.lean` — Definitions of `FinDist`, `SREV`, `BREV`, `Mechanism` (IC/IR on $\mathbb{R}_{\ge 0}^m$), `REV`.
  - `Kernel.lean`, `KernelRun.lean`, `Spec.lean`, `GroupAData.lean`, `Instance.lean`, `RealA.lean`, `GroupB.lean`, `Main.lean` — General lower bound (`> 2.1946`).
  - `Expect.lean`, `Duality.lean`, `LayerCake.lean`, `Laplace.lean`, `Abstract.lean`, `AnalyticCert.lean`, `CertKernel.lean`, `CertData.lean`, `CertRun.lean`, `CertSound.lean`, `CertLoop.lean`, `UpperBound.lean`, `Bracket.lean` — Universal upper bound (`≤ 2.486`) and combined bracket `cStar_mem`.
  - `IIDKernel.lean`, `IID.lean` — 40-item i.i.d. lower bound (`> 1.316`).
- `lean/Audit.lean`, `lean/AUDIT.md` — Axiom inspection and anti-trivialization audit.
- `code/` — Exact rational evaluators, integer certificates, baseline reproductions, and stylised financial-auction illustration.
- `paper/` — Working paper draft (`working_paper.md`) and arXiv-ready LaTeX source (`paper.tex`).
- `CITATION.cff`, `LICENSE` — Authorship, citation metadata, and licensing.

## How to Cite

```bibtex
@misc{sievert2026simpleauctions,
  author       = {Max Sievert},
  title        = {Machine-Verified Two-Sided Bounds for Simple Multi-Item Auctions: $2.1946 \le c^* \le 2.486$ and an i.i.d.\ Lower Bound above $1.316$},
  year         = {2026},
  month        = {November},
  howpublished = {\url{https://github.com/Maxkxkxk/simple-auctions-lean4}}
}
```

## Credit, License & Scope

The core mathematical proof frameworks build on Giannakopoulos & Hahn (`arXiv:2609.35711v3`), Huang (`arXiv:2610.24973v2`), Cai–Devanur–Weinberg, and Rubinstein. All formal theorems are proven for arbitrary finite independent product distributions `D : Fin m → FinDist`. See `paper/paper.tex` and `lean/AUDIT.md` for full attribution, methodology (human-directed AI-assisted discovery + Lean 4 verification), and scope. Code and Lean proofs are MIT-licensed; the manuscript in `paper/` is licensed under CC BY 4.0 (Copyright © 2026 Max Sievert).
