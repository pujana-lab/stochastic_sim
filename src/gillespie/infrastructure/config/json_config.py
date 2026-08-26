from __future__ import annotations

from pathlib import Path
import json as _json

from src.gillespie.simulation_config import CellTypeConfig


_CELL_TYPE_FIELD_MAP = {
    "N": "N",
    "K": "K",
    "lambda_0": "lambda0",
    "lambda0": "lambda0",
    "mu": "mu",
    "nu": "nu",
    "omega_exhaust": "omega_exhaust",
    "next_mutation": "next_mutation",
    "fitness_gain": "fitness_gain",
}

_ALLOWED_CELL_TYPES = {"base", "immune", "mutated", "exhausted"}


def _coerce_cell_type_config(raw: dict | CellTypeConfig) -> CellTypeConfig:
    if isinstance(raw, CellTypeConfig):
        return raw
    if not isinstance(raw, dict):
        raise TypeError(f"Cell type config must be a dict or CellTypeConfig, got {type(raw).__name__}")

    mapped = {}
    for key, value in raw.items():
        if key in _CELL_TYPE_FIELD_MAP:
            mapped[_CELL_TYPE_FIELD_MAP[key]] = value
    print(f"Mapped cell type config: {mapped}")
    return CellTypeConfig(**mapped)


def flatten_cell_types(raw: dict) -> dict:
    result: dict = {}
    cell_map = raw.get("cell_parameters") if "cell_parameters" in raw else raw.get("cell_types")

    if cell_map is not None:
        converted: dict[str, CellTypeConfig] = {}
        for cell_name, spec in cell_map.items():
            if cell_name not in _ALLOWED_CELL_TYPES:
                raise ValueError(f"Unknown cell type '{cell_name}' in config")
            converted[cell_name] = _coerce_cell_type_config(spec)
        result["cell_parameters"] = converted

    for k, v in raw.items():
        if k not in {"cell_types", "cell_parameters"}:
            result[k] = v
    return result


def load_config_json(path: Path) -> dict:
    with open(path) as f:
        return _json.load(f)
