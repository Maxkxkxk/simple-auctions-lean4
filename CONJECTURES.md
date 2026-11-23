# CONJECTURES

- **C1 (PROVEN, Lean `lower_bound`):** $c^* > 2.1946$.
- **C1' (PROVEN, Lean `upper_bound`):** $c^* \le 2.486$; sharper: $\mathrm{REV} \le 1.1187\,\mathrm{SREV} + 1.367299\,\mathrm{BREV}$ (`upper_bound_two`).
- **C1'' (PROVEN, Lean `cStar_mem`):** $2.1946 \le c^* \le 2.486$, with $c^*$ the infimum of valid constants over finite instances.
- **C2 (OPEN, numerical):** the equal-revenue entry-fee family has supremum about $2.1957$.
- **C3 (OPEN):** within "concentrated group + independent group A" instances, a single two-part tariff on A is optimal for the worst-case A.
- **C4 (OPEN, numerical):** a two-normalisation version of the scheme gives about $2.461$; the Laplace/Lagrangian family cannot go below about $2.43$.
- **C5 (PROVEN, Lean `iid_lower_bound` + exact Python):** for i.i.d. items, $\mathrm{REV}(F^{\otimes 40}) > 1.316 \max\{\mathrm{SREV}, \mathrm{BREV}\}$ (Lean 4: `SimpleAuctions.iid_lower_bound`); exact rational computation gives $\ge 1.324907$ for $100$ i.i.d. items (`code/iid_lower_bound.py`), improving the $1.1877$ bound of Ma & Simchi-Levi cited in `arXiv:2609.35711v3`. Numerical limit of "one entry fee + one price" found so far: about $1.33$.
