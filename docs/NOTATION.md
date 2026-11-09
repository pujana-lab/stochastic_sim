# Model notation

The model is written in two layers, and each keeps its own symbols:

- **Mean-field (ODE)**: the report (`main.pdf`, eq. 4) and the
  equations of [`PLAN_fibroblastos.md`](./PLAN_fibroblastos.md) §2.
- **Gillespie**: [`README.md`](../README.md), the poster and the docstrings
  in `src/gillespie/`.

Use these names in code comments, docstrings, docs and plans. Symbols from
sketches and notebooks are mapped in the last section; do not copy them.

## Populations

| Code type | Mean-field | Dimensionless | Gillespie count |
|---|---|---|---|
| `base` | — | — | `N_base` |
| `mutated` | `C` | `c = C / K_C` | `N_mutated` |
| `immune` | `I` | `i = I / K_I` | `N_immune` |
| `exhausted` | `E` | `e = E / K_I` | `N_exhausted` |
| `fibroblast` | `F_N` | `f_N = F_N / K_F` | `N_fibroblast` |
| `caf` | `F_T` | `f_T = F_T / K_F` | `N_caf` |

- The wild type (`base`) only exists in the Gillespie.
- Initial values: `C0, I0, E0, F_N0, F_T0`. Time: `t`; dimensionless
  `τ = λ_C t`. In the Gillespie layer `τ` is instead the time step (the SSA
  waiting time in `README.md`, the leap in `PLAN_tau_leaping.md`).
- Lowercase letters are always the dimensionless variables. In particular
  `f_N`, `f_T` are `F_N / K_F`, `F_T / K_F`, never the dimensional populations.

## Parameters

| Meaning | Mean-field | Dimensionless group | Code | `OMEGA` scaling |
|---|---|---|---|---|
| tumour birth | `λ_C` | time unit (`τ = λ_C t`) | `cell_parameters["mutated"].lambda0` | — |
| tumour death | `μ_C` | `γ_C = μ_C / λ_C` | `cell_parameters["mutated"].mu` | — |
| tumour capacity | `K_C` | unit of `c` | `cell_parameters["mutated"].K` | fraction of `OMEGA`, × `OMEGA` |
| immune birth | `λ_I` | `ρ = λ_I / λ_C` | `cell_parameters["immune"].lambda0` | — |
| immune capacity | `K_I` | unit of `i`, `e` | `cell_parameters["immune"].K` | fraction of `OMEGA`, × `OMEGA` |
| exhausted death | `μ_E` | `μ_E ← μ_E / λ_C` | `cell_parameters["exhausted"].mu` | — |
| killing of tumour by immune | `θ` | `σ = θ K_I / λ_C` | `theta_I` | ÷ `OMEGA` |
| immune activation by tumour | `β` | `β ← β K_C / λ_C` | `beta` | ÷ `OMEGA` |
| immune exhaustion by tumour | `ω` | `ω ← ω K_C / λ_C` | `cell_parameters["immune"].omega_exhaust` | ÷ `OMEGA` |
| fibroblast birth | `λ_F` | `ρ_F = λ_F / λ_C` | `cell_parameters["fibroblast"].lambda0` | — |
| fibroblast capacity | `K_F` | `k = K_C / K_F` | `cell_parameters["fibroblast"].K` | fraction of `OMEGA`, × `OMEGA` |
| CAF death | `μ_F` | `μ_F ← μ_F / λ_C` | `cell_parameters["caf"].mu` | — |
| activation `F_N → F_T` by tumour | `κ` | `κ ← κ / λ_C` | `kappa` | unscaled (saturating) |
| tumour size for half-maximal activation | `C_½` | `c_½ = C_½ / K_C` | `C_half` | fraction of `OMEGA`, × `OMEGA` |
| immune exhaustion by CAFs | `ω_F` | `ω_F ← ω_F K_F / λ_C` | `omega_F` | ÷ `OMEGA` |
| CAF support of tumour birth | `α_F` | `α_F ← α_F K_F` | `alpha_F` | ÷ `OMEGA` |
| crowding weight: share of type T in the crowding numerator of X | `w_T` (in the factor of X) | `w_T` (already dimensionless) | `cell_parameters[X].crowding_pool[T]` | unscaled (default 1) |

- `←` means the dimensionless group reuses the name of the dimensional
  parameter.
- The killing rate is `θ` in the mean-field layer and `θ_I` in the Gillespie
  (code `theta_I`). It is the same parameter.
- "÷ `OMEGA`" and "× `OMEGA`" are applied once by `SimulationConfig` when
  `scale=True`; config files write unscaled values.
- The effective Gillespie birth rate is `lambda0 · (1 + fitness_gain)`.
- Gillespie only: mutation rate `ν` (`nu`), `instability`, `K_min` (`Kmin`),
  and the `base` rates.

## Gillespie rates

- Effective rates per clone: `r_B` (birth), `r_D` (death), `r_M` (mutation),
  `r_E` (conversion: exhaustion or fibroblast activation), written
  `r_B(<type>) = ...`.
- Base forms: `r_B = λ · N · g(Y_X, K_X)`, `r_D = μ · N`,
  `r_M = ν · N · (1 + instability)`, `r_E = ω · N`.
- In the per-type formulas, a cell parameter carries the code type name:
  `λ_mutated` (= `λ_C`), `μ_mutated` (= `μ_C`), `λ_immune` (= `λ_I`),
  `ω_immune` (= `ω`), `λ_fibroblast` (= `λ_F`), `μ_caf` (= `μ_F`).
  Interaction parameters keep their symbol: `θ_I, β, κ, C_½, ω_F, α_F`.

## Crowding

The mean-field equations write the logistic factor explicitly, with unscaled
capacities. The Gillespie writes it as `g(Y_X, K_X) = 1 − Y_X / K_X`, where
`Y_X` is the (weighted) population that competes with type X (`crowding_numerator()`)
and `K_X` is `actual_K`: × `OMEGA`, inflated by `AdaptedCrowding`, and
time-dependent if `decay`. The code clamps it, `max(0, g)`. `g` is Gillespie
notation only; do not use it in the ODEs.

| Type X | Mean-field factor | Gillespie `Y_X` |
|---|---|---|
| `base` | — | `N_base + N_mutated` |
| `mutated` | `1 − C / K_C` | `N_mutated` |
| `immune` | `1 − (I + E) / K_I` | `N_immune + N_exhausted` |
| `fibroblast` | `1 − (w_N F_N + w_T F_T + w_C C) / K_F` | `w_N N_fibroblast + w_T N_caf + w_C N_mutated` |

The table writes the weights only where they are biologically relevant (the fibroblast
niche). In general `Y_X = Σ_T w_T N_T` over the crowding pool of X (`cell_parameters[X].crowding_pool`,
the types that compete), and `w_T` is the share that one cell of T takes in X's crowding
numerator. Every member has its own weight, X itself included (`w_N` for F_N in the
fibroblast niche): no type is the reference. By default every weight is 1, so the factors reduce to the plain sums; a config
changes only the weights it names, and weight 0 removes a type from the pool. The
weight belongs to the niche X: if two niches weight the same type differently, write
`w_{X,T}`. Dimensionless: `1 − w_N f_N − w_T f_T − k w_C c`.

Do not use the older `C(t)`, `N_crowd` or `K_t`.

## Naming new terms

- Own rates take the type subscript: `λ_C`, `μ_E` → `λ_F`, `μ_F`. Capacities
  are `K_<type>`.
- Interactions are Greek letters, with no subscript when the tumour drives
  them (`β`, `ω`, `κ`) and with the driving type as subscript otherwise
  (`θ_I`, `ω_F`, `α_F`).
- Every new parameter gets a row in the table above: mean-field symbol,
  dimensionless group, code name and `OMEGA` scaling (bilinear terms
  ÷ `OMEGA`; capacities and half-saturations × `OMEGA`; saturating rates
  unscaled).

## Symbols from sketches and notebooks

| Sketch | Use |
|---|---|
| `T_A` | `I` |
| `f_N`, `f_T` (dimensional) | `F_N`, `F_T` |
| `r_F` | `λ_F` |
| `d_F` | `μ_F` |
| `d_E` | `μ_E` |
| `K` (fibroblast equation) | `K_F` |
| `K₀` / `K_0` | `C_½` |
| `κ` (fibroblast equations) | `κ` |
| `κ` (equation of `E`) | `ω` |
| `κ̂` | `ω_F` |

Other sketch symbols (e.g. `W`, `φ`) have no counterpart in the model yet.
Name them with the rules above when they enter it.
