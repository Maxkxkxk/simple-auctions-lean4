# Machine-verified bounds $2.1946 \le c^* \le 2.486$ and an i.i.d. lower bound above $1.316$

*Working paper draft, 2026-11-23 (LaTeX source in `paper/paper.tex`). Verified claims only, unless marked "numerical".*

## Setting and results
One additive buyer, $m$ items, independent finitely supported values. `REV` is the supremum over all IC+IR mechanisms (lotteries allowed), `SREV` and `BREV` are the optimal separate and grand-bundle posted-price revenues over all real prices. Let $c^*$ be the infimum of all $c$ with $\mathrm{REV} \le c \cdot \max\{\mathrm{SREV}, \mathrm{BREV}\}$ on all finite instances.

- **Theorem 1** (Lean `SimpleAuctions.lower_bound`): There is an instance with $\mathrm{REV} > 2.1946 \cdot \max\{\mathrm{SREV}, \mathrm{BREV}\}$.
- **Theorem 2** (Lean `SimpleAuctions.upper_bound_two`, `upper_bound`): For every instance, $\mathrm{REV} \le 1.1187\,\mathrm{SREV} + 1.367299\,\mathrm{BREV}$, hence $\mathrm{REV} \le 2.486 \cdot \max\{\mathrm{SREV}, \mathrm{BREV}\}$.
- **Corollary 3** (Lean `SimpleAuctions.cStar_mem`, `main_bracket`): $2.1946 \le c^* \le 2.486$.
- **Theorem 4** (Lean `SimpleAuctions.iid_lower_bound`): For $40$ i.i.d. items with $\Pr[X \ge k] = \frac{3}{10\max(k,4)}$ ($1 \le k \le 13$) and deterministic entry fee $9$ + item price $4$, $\mathrm{REV}(F^{\otimes 40}) > 1.316 \cdot \max\{\mathrm{SREV}, \mathrm{BREV}\}$.
- **Proposition 5** (exact rational computation, `code/iid_lower_bound.py`): Scaling the same i.i.d. construction to $100$ items yields $\ge 1.324907$.

Published bounds (Nov 2026): $[2.19, 2.50]$ and, for i.i.d. items, $1.1877$ as the best known lower bound (Giannakopoulos–Hahn, `arXiv:2609.35711v3`, and Ma & Simchi-Levi).

## Theorem 1: construction
Unit $N = 400$. *Group A*: five types, $256$ i.i.d. items each, prices $P = (8, 21, 65, 132, 1280)$, entry fee $e = 580$, $R = (1067039, 1128540, 2129082, 2782033, 19107703)$, tails $\Pr[X \ge k] = \lfloor R_i/\max(k,P_i)\rfloor / 2^{24}$ for $1 \le k \le P_i + e$. *Group B*: $L^2$ i.i.d. items on $\{j\delta : j \le L\}$, $\Pr[Y \ge j\delta] = 1/j$, $\delta = N/(L^2 H_L)$. *Mechanism*: sell $B$ as a bundle at $N(1 - 10^{-6})$; independently, pay $e$ and buy any $A$-items at prices $P_i$. Proof: equal revenue gives $\mathrm{SREV} \le N + N/H_L$; $(N+t)\Pr[Z_A \ge t] \le N$ for all $t$ gives $\mathrm{BREV} \le N + o(1)$; the $A$-revenue is at least $1.19464\,N$ by a kernel computation with packed generating functions.

## Theorem 2: proof
1. **Master inequality** (Giannakopoulos–Hahn): Caps $c_j$ with $\mathbb{E}\min(X_j, c_j) = r_j = \mathrm{REV}(X_j)$, $Y_j = (X_j - c_j)^+$. Routing the utility term to the item with the largest excess and bounding each coordinate line by threshold allocations gives $\mathrm{REV} \le \mathrm{SREV} + \mathbb{E}[\sum Y_j - \max Y_j]$. Our formal proof (`Duality.lean`) is elementary (finite sums, one induction over the support).
2. **Layer cake and Jensen** (`LayerCake.lean`, `Abstract.lean`): $\mathbb{E}[\sum Y_j - \max Y_j] = \int (\sum x_j - 1 + \prod(1 - x_j))\,dt \le \int \sum x_j \psi(\alpha_j/v_j)\,dt$, with $x_j = \Pr[X_j \ge c_j + t]$, $v_j = (c_j+t)/\ell$, $\psi(z) = 1 - (1-e^{-z})/z$, weights $r_j/(\ell\alpha_j)$, $\alpha_j = \alpha(c_j/\ell)$. Taking $\alpha$ to be a monotone *step function* gives the normalizing scale $\ell$ as a monotone fixed point without continuity.
3. **Bundle revenue** (`Laplace.lean`): With $w(v) = \sum_{i=1}^4 a_i e^{-v/s_i}$, $\sum_j \mathbb{E}[1 - e^{-X_j/\sigma}] \le \kappa_i + b_i \mathrm{BREV}/\sigma$, from $\mathbb{E} e^{-\sum X/\sigma} \ge \int_1^\infty e^{-yx}x^{-2}\,dx$ ($y = \mathrm{BREV}/\sigma$) and certified linear tangents.
4. **Certificate** (`AnalyticCert.lean`–`CertLoop.lean`): The 1D inequality in $t = c/\ell$ is checked by directed integer arithmetic over $1846$ cells and $35$ steps of $\alpha$ (`CertRun.lean`); soundness is proved formally. Parameters: $a = (0.266, 0.509, 0.203, 0.022)$, $s = (0.5, 1.4, 22, 209.3)$.

## Theorem 4: i.i.d. proof (`IIDKernel.lean`, `IID.lean`)
$40$ i.i.d. items on $\{0, 4, \dots, 13\}$ with weights $c_N / 3603600$ giving $\Pr[X \ge k] = \frac{3}{10\max(k,4)}$ ($1 \le k \le 13$). $\mathrm{SREV} \le 12$; packed generating polynomial $g_z^{40}$ at base $X = 2^{900}$ with a shifting loop certifies $\mathrm{BREV} \le 12.151$; modular surplus generating function mod $X^9$ certifies deterministic entry-fee revenue $\ge 16$; $16 / 12.151 > 1.316$.

## What limits the method (numerical)
The master inequality itself allows about $2.33$–$2.35$. The one-parameter scheme plateaus at $2.482$–$2.484$; a two-normalisation 2D variant reaches about $2.461$; a fully profile-adapted kernel plateaus near $2.43$. Constants below $\approx 2.43$ require a stronger treatment of the bundle constraint.

## Illustration: multi-tranche issues (stylised, not market data)
Consider a single strategic bidder (e.g., dominant primary dealer) with independent private surplus across tranches (`code/financial_illustration.py`). Theorem 2 guarantees that the better of tranche-by-tranche pricing and one package captures at least $1/2.486 = 40.2\%$ of any mechanism; Theorem 1 shows this worst-case guarantee cannot exceed $1/2.1946 = 45.6\%$. For $40$ identical lots with Theorem 4's law, separate pricing yields $12.0000$, one package $12.1509$, and participation fee $9$ + price $4$/lot yields $16.0014$ ($+31.7\%$ over the better simple format). For $4$ heterogeneous tranches (`2y, 5y, 10y, 30y`), separate pricing ($7.50$) already dominates bundling ($4.50$) and all tested entry-fee menus, with Theorem 2 envelope $14.5431$. Caveats: single bidder, additive independent values, no common-value component; real bond tranches are strongly correlated.

## Credit
Master inequality, caps, Jensen weights, Laplace argument and entry-fee constructions: Giannakopoulos–Hahn. Core/tail duality: Cai–Devanur–Weinberg. Seed/completion view: Huang. Concentrated group: Rubinstein. Ours: the equal-revenue-on-$[p, p+e]$ parameterization and constants in Theorem 1; the 4-exponential kernel, step-function weights, and constants in Theorem 2; Theorem 4 / Proposition 5 (i.i.d. lower bound); and the complete Lean 4 formalization of all three theorems.

## Reproducibility
`cd lean && lake build && lake env lean Audit.lean`; `cd code && python evaluator.py && python ub_certificate.py && python iid_lower_bound.py && python financial_illustration.py`.
