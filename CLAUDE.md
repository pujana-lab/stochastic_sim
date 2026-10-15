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
