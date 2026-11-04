from __future__ import annotations

import json
from pathlib import Path

from src.gillespie.config_report import format_config
from src.gillespie.simulation_config import SimulationConfig


def save_run_config(path: Path, config: SimulationConfig) -> tuple[Path, Path]:
    """Write the conditions of a run: input values (unscaled) and used values (after OMEGA scaling).

    <path>.txt: readable table, one row per parameter.
    <path>.json: {"input": ..., "used": ...}. "input" is a valid JSON config:
    pass it to --config (or load_config_json) to reproduce the run.
    Shared by the CLI, scripts/run_batch_simulations.py and param_sweep.py.
    """
    txt, js = path.with_suffix(".txt"), path.with_suffix(".json")
    txt.parent.mkdir(parents=True, exist_ok=True)
    txt.write_text(format_config(config) + "\n", encoding="utf-8")
    payload = {"input": config.input_values, "used": config.as_dict()}
    js.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return txt, js
