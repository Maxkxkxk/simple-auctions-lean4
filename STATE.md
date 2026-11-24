# STATE (Run 2 — Finalized & Audited, 2026-11-23)

## Human Summary (Einfaches Deutsch)

- **Alle drei Hauptresultate sind jetzt 100 % in Lean 4 bewiesen und geprüft:**
  1. **Allgemeine Untergrenze (`lower_bound`):** `REV > 2,1946 · max(SREV, BREV)` (veröffentlicht war `2,19027`).
  2. **Allgemeine Obergrenze (`upper_bound`, `upper_bound_two`):** Für *jede* endliche unabhängige Verteilung und *jeden* Mechanismus (auch Lotterien) gilt `REV ≤ 1,1187·SREV + 1,367299·BREV ≤ 2,486 · max(SREV, BREV)` (veröffentlicht war `2,50`). Zusammen ergibt das **`2,1946 ≤ c* ≤ 2,486`** (`main_bracket`, `cStar_mem`).
  3. **i.i.d.-Untergrenze (`iid_lower_bound`) — jetzt ebenfalls fertig in Lean 4:** Für 40 identisch verteilte Güter gilt **`REV > 1,316 · max(SREV, BREV)`** (und mit 100 Gütern exakt rational in Python **`≥ 1,3249`**). Giannakopoulos & Hahn zitieren `1,1877` als bisher beste bekannte Schranke – das ist ein Sprung um über `+0,128` (Lean) bzw. `+0,137` (Python).
- **Saubere Verifikation ohne Tricks:** Alle 25 Lean-Module bauen sauber (`lake build`), alle 6 exportierten Theoreme hängen laut `Audit.lean` **ausschließlich** von den 3 Standard-Axiomen `[propext, Classical.choice, Quot.sound]` ab (kein `sorry`, kein `native_decide`), und alle 25 Module bestehen `leanchecker`.
- **Veröffentlichungsfertig vorbereitet:** `paper/paper.tex` (arXiv-/Overleaf-taugliches LaTeX), `.github/workflows/verify.yml` (automatischer Lean-4- & Python-CI-Check), `README.md` und `.gitignore` liegen bereit. Die 9 verwaisten Python-Prozesse aus Run 1 wurden bereinigt.
- **Ehrliche Einordnung & Grenzen:** Die Schranken `2,1946` und `2,486` verbessern den Beweisrahmen von Giannakopoulos & Hahn (2026); das stärkste inhaltlich neue Ergebnis neben dem vollständigen Maschinenbeweis ist der große Sprung bei i.i.d.-Gütern (`1,1877 → > 1,316` in Lean / `1,3249` exakt). Alle formalen Sätze gelten für Verteilungen mit endlich vielen Werten; die Finanzauktions-Illustration (`code/financial_illustration.py`) ist stilisiert (keine Marktdaten).

## Live Scoreboard

| Quantity | Value | Status |
|---|---|---|
| SOTA lower bound (published) | `2 + 0.88/4.625 = 2.190270...` | `arXiv:2609.35711v3`, Thm D.1 (their instance re-evaluated rigorously: `>= 2.190949`) |
| **Our lower bound** | **`> 2.1946`** (`2.194649949...`) | **Lean 4 proven** (`lower_bound`) + exact integer check (`code/certificate.py`) |
| SOTA upper bound (published) | `2.5` (`7/6 SREV + 4/3 BREV`) | `arXiv:2609.35711v3`, Prop 4.6 |
| **Our upper bound** | **`<= 2.486`** (`1.1187 SREV + 1.367299 BREV`) | **Lean 4 proven** (`upper_bound`, `upper_bound_two`) + exact integer check (`2.485776`) |
| **Combined bracket** | **`2.1946 <= c* <= 2.486`** | **Lean 4 proven** (`main_bracket`, `cStar_mem`) |
| SOTA i.i.d. lower bound (cited in GH v3) | `1.1877` | Ma & Simchi-Levi (2 items) |
| **Our i.i.d. lower bound (40 items)** | **`> 1.316`** (`1.316884...`) | **Lean 4 proven** (`iid_lower_bound`) + exact rational check |
| **Our i.i.d. lower bound (100 items)** | **`>= 1.324907`** | Exact rational Python (`code/iid_lower_bound.py`) |

## Current Phase & Active Track

**All core tracks (Phases 0, 1, 2, and 3) are complete, verified, and publicly archived:**
- Public GitHub repository: https://github.com/Maxkxkxk/simple-auctions-lean4 (green GitHub Actions Lean 4 & Python CI).
- Permanent CERN Zenodo DOI: https://doi.org/10.5281/zenodo.24482558 (`Max Sievert`, Release `v1.0.1`).
- Track 1 (`lower_bound`), Track 2 (`upper_bound`, `upper_bound_two`, `main_bracket`, `cStar_mem`), and i.i.d. Track (`iid_lower_bound`) are all machine-verified in Lean 4.
- Release & publication assets ready: `paper/paper.tex`, `paper/paper.pdf`, `paper/working_paper.md`, `README.md`, `CITATION.cff`, `LICENSE`, `.gitignore`, `.github/workflows/verify.yml`.

## Verified Files & Re-Run Commands

- `cd lean && lake build` (Lean `v4.36.0`, Mathlib `v4.36.0`; run `lake exe cache get` first on a new machine; `elan` is in `%USERPROFILE%\.elan\bin`). Heavy kernel runs: `KernelRun` (~2 min), `CertRun` (~35 s), `IIDKernel` (~1.2 s).
- `cd lean && lake env lean Audit.lean` prints statements and axioms for all 6 exported theorems.
- `cd lean && lake env leanchecker SimpleAuctions.<Module>`: exit 0 for all 25 modules (see `lean/AUDIT.md`).
- `cd code && python evaluator.py && python ub_certificate.py && python iid_lower_bound.py && python financial_illustration.py`
- `cd code && python baselines.py`, `python certificate.py ../data/instance_A.json` (fast with `pip install gmpy2`).

## Immediate Next Tasks (Optional / Publication)

1. **Publish & claim priority:** Fill in author name/email in `paper/paper.tex`, push repo to GitHub (verify green GitHub Actions CI), mint a Zenodo DOI for release `v1.0`, upload PDF to SSRN, and email Giannakopoulos–Hahn, Yichen Huang (Harvard), and Nikhil Garg (`EconCSLib`) requesting feedback and a `cs.GT` arXiv endorsement.
2. **Optional future research (Track 2 below `2.47`):** Two-normalisation 2D certificate (`code/wip/ub_tdep3.py` float lead `2.461`).
3. **Optional future research (i.i.d. 100 items in Lean):** Upgrade `IIDKernel.lean` from `n = 40` (`> 1.316`) to `n = 100` (`> 1.324`) using the same Horner + binary squaring architecture if desired.
