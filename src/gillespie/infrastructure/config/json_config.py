from __future__ import annotations

import json
from pathlib import Path


def load_config_json(path: Path) -> dict:
    """Load the canonical nested Gillespie configuration."""
    with path.open(encoding="utf-8") as config_file:
        config = json.load(config_file)

    if "cell_types" in config:
        raise ValueError("Use the nested 'cell_parameters' configuration key")
    if "cell_parameters" in config and not isinstance(config["cell_parameters"], dict):
        raise TypeError("'cell_parameters' must be an object keyed by cell type")

    # Calculator provenance is useful in the file but is not a
    # SimulationConfig field.
    config.pop("_calculator", None)

    # JSON has no tuples: clone ids (default_id) come back as lists.
    for params in config.get("cell_parameters", {}).values():
        if isinstance(params, dict) and isinstance(params.get("default_id"), list):
            params["default_id"] = tuple(params["default_id"])
    return config
