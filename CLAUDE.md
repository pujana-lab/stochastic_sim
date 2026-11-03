# CLAUDE.md

The canonical agent guidance for this repo lives in `.github/` (shared with
GitHub Copilot). It is imported here so both tools follow the same rules; edit
the `.github/` sources, not copies.

@.github/copilot-instructions.md

## Claude Code mapping of `.github/`

| `.github/` source | Claude Code equivalent |
| --- | --- |
| `copilot-instructions.md` | imported above |
| `instructions/gillespie-python.instructions.md` | `.claude/rules/gillespie-python.md` (path-scoped) |
| `agents/*.agent.md` | `.claude/agents/*.md` subagents |
| `skills/*/SKILL.md` | `.claude/skills/*` (symlinks to `.github/skills/*`) |

When an agent or instruction file in `.github/` changes, update its
`.claude/` counterpart in the same change.

## Environment notes

- Use the project virtualenv: `.venv/bin/python -m pytest tests/gillespie -q`.
  (The `makefile` still references `venv/`, which does not exist.)
- The root-level `debug.csv`, `history.csv` and `rate_his.csv` are multi-GB
  generated outputs; never read them whole and never commit them.
- `CLONE_TYPE_REFACTORING.md`, `TISSUE_STATE_REFACTOR.md` and
  `docs/gillespie_full_diagnosis_roadmap.md` are historical notes, not current
  contracts; verify against `src/gillespie/` and tests.

## Model notation

<!-- Pending: move this section to `.github/copilot-instructions.md` so Copilot shares it (that file could not be edited from this session). -->

Use the notation of the existing documentation in code comments, docstrings,
docs and plans. Do not import ad-hoc symbols from sketches or notebooks (e.g.
`W`, `T_A`, `φ`, `f_N`, `r_F`, `d_F`, `K_0`); map them to the names below.

There are two layers, and each keeps its own notation.

**Mean-field / ODE** (report `main.pdf` eq. 4, `application/scaling_calculator.py`,
`notebooks/omega_convergence.py`):

- Populations: `C` (tumour = `mutated`), `I` (active immune = `immune`),
  `E` (`exhausted`). Wild type (`base`) is not part of the 3-population model.
- Dimensional parameters: `λ_C, λ_I, μ_C, μ_E, θ, β, ω, K_C, K_I`; initial
  values `C0, I0, E0`.
- Dimensionless (rates / `λ_C`, populations / `K`): `c, i, e`, `τ = λ_C t`;
  groups `γ_C = μ_C/λ_C`, `σ = θ K_I/λ_C`, `ρ = λ_I/λ_C`, `β ← β K_C/λ_C`,
  `ω ← ω K_C/λ_C`, `μ_E ← μ_E/λ_C` (same names as `DimensionlessParameters`).

**Gillespie** (`README.md`, poster):

- Counts `N_<type>` with the code type name: `N_base`, `N_mutated`,
  `N_immune`, `N_exhausted`.
- Effective rates per clone: `r_B` (birth), `r_D` (death), `r_M` (mutation),
  `r_E` (exhaustion/conversion), written as `r_B(<type>) = ...`.
- Base forms: `r_B = λ · N · g(Y_X, K_X)`, `r_D = μ · N`,
  `r_M = ν · N · (1 + instability)`, `r_E = ω · N`.
- Crowding (both layers): `g(Y, K_X) = 1 − Y / K_X`, where `Y` is the
  competing population of type X (the sum returned by `crowding_numerator()`)
  and `K_X` is its carrying capacity: unscaled in the mean-field equations
  (`K_C`, `K_I`), scaled in the Gillespie (`actual_K`: × `OMEGA`, inflated by
  `AdaptedCrowding`, time-dependent if `decay`). The code clamps it,
  `max(0, g)`. Do not use the README's older `C(t)` / `N_crowd` / `K_t`
  names; write the tumour as `N_mutated` (Gillespie) or `C` (mean-field).
- Interactions: `θ_I` (killing by immune), `β` (immune activation by mutated),
  `ω` (exhaustion by mutated). System size is `OMEGA` (`Ω`); `K` inputs are
  fractions of `OMEGA`.
- Code names follow the symbols: `lambda0`, `mu`, `nu`, `omega_exhaust`, `K`,
  `theta_I`, `beta`.

**Naming new terms**

- Own rates take the type subscript (`λ_C`, `μ_E` → `λ_F`, `μ_F`); capacities
  are `K_<type>`.
- Interactions are Greek letters: no subscript when driven by the tumour
  (`β`, `ω`, `κ`), subscript = the driving type otherwise (`θ_I`, `ω_F`, `α_F`).
- Every new parameter gets a dimensional symbol, a dimensionless group and a
  code name, and its `OMEGA` scaling is stated (bilinear ÷ `OMEGA`,
  capacities and half-saturations × `OMEGA`, saturating rates unscaled).

**Planned fibroblast module** (`docs/PLAN_fibroblastos.md`, not implemented):
`F_N` (`fibroblast`) and `F_T` (`caf`); `λ_F, μ_F, K_F, κ, C_½, ω_F, α_F`;
code `kappa`, `C_half`, `omega_F`, `alpha_F`.
