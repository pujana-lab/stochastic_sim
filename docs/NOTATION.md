# Model notation

Symbols used in code comments, docstrings, docs and plans. Use these names;
do not import ad-hoc symbols from sketches or notebooks (e.g. `W`, `T_A`, `φ`,
`f_N`, `r_F`, `d_F`, `K_0`), map them to the names below.

There are two layers, and each keeps its own notation.

## Mean-field / ODE

Used by the report (`main.pdf`, eq. 4),
[`scaling_calculator.py`](../src/gillespie/application/scaling_calculator.py) and
[`omega_convergence.py`](../notebooks/omega_convergence.py).

- Populations: `C` (tumour = `mutated`), `I` (active immune = `immune`),
  `E` (`exhausted`). Wild type (`base`) is not part of the 3-population model.
- Dimensional parameters: `λ_C, λ_I, μ_C, μ_E, θ, β, ω, K_C, K_I`; initial
  values `C0, I0, E0`.
- Dimensionless (rates / `λ_C`, populations / `K`): `c, i, e`, `τ = λ_C t`;
  groups `γ_C = μ_C/λ_C`, `σ = θ K_I/λ_C`, `ρ = λ_I/λ_C`, `β ← β K_C/λ_C`,
  `ω ← ω K_C/λ_C`, `μ_E ← μ_E/λ_C` (same names as `DimensionlessParameters`).

## Gillespie

Used by [`README.md`](../README.md) and the poster.

- Counts `N_<type>` with the code type name: `N_base`, `N_mutated`,
  `N_immune`, `N_exhausted`.
- Effective rates per clone: `r_B` (birth), `r_D` (death), `r_M` (mutation),
  `r_E` (exhaustion/conversion), written as `r_B(<type>) = ...`.
- Base forms: `r_B = λ · N · g(Y_X, K_X)`, `r_D = μ · N`,
  `r_M = ν · N · (1 + instability)`, `r_E = ω · N`.
- Interactions: `θ_I` (killing by immune), `β` (immune activation by mutated),
  `ω` (exhaustion by mutated). System size is `OMEGA` (`Ω`); `K` inputs are
  fractions of `OMEGA`.
- Code names follow the symbols: `lambda0`, `mu`, `nu`, `omega_exhaust`, `K`,
  `theta_I`, `beta`.

## Crowding (both layers)

`g(Y, K_X) = 1 − Y / K_X`, where `Y` is the competing population of type X
(the sum returned by `crowding_numerator()`) and `K_X` is its carrying
capacity:

- unscaled in the mean-field equations (`K_C`, `K_I`);
- scaled in the Gillespie (`actual_K`: × `OMEGA`, inflated by
  `AdaptedCrowding`, time-dependent if `decay`).

The code clamps it, `max(0, g)`. Do not use the older `C(t)` / `N_crowd` /
`K_t` names; write the tumour as `N_mutated` (Gillespie) or `C` (mean-field).

## Naming new terms

- Own rates take the type subscript (`λ_C`, `μ_E` → `λ_F`, `μ_F`); capacities
  are `K_<type>`.
- Interactions are Greek letters: no subscript when driven by the tumour
  (`β`, `ω`, `κ`), subscript = the driving type otherwise (`θ_I`, `ω_F`, `α_F`).
- Every new parameter gets a dimensional symbol, a dimensionless group and a
  code name, and its `OMEGA` scaling is stated (bilinear ÷ `OMEGA`,
  capacities and half-saturations × `OMEGA`, saturating rates unscaled).

## Planned fibroblast module

Not implemented; see [`PLAN_fibroblastos.md`](./PLAN_fibroblastos.md).
`F_N` (`fibroblast`) and `F_T` (`caf`); `λ_F, μ_F, K_F, κ, C_½, ω_F, α_F`;
code `kappa`, `C_half`, `omega_F`, `alpha_F`.
