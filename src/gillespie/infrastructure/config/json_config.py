from pathlib import Path
import json as _json


_CELL_TYPE_MAP: dict[str, dict[str, str]] = {
    "base": {
        "N": "N0",
        "lambda_0": "lambda0",
        "mu": "mu0",
        "nu": "nu0",
        "K": "K0",
    },
    "immune": {
        "N": "N_immune",
        "lambda_0": "lambda_Immune",
        "mu": "mu_Immune",
        "K": "K_immune",
    },
    "mutant": {
        "N": "N_mutant",
        "K": "K_mutant",
    },
    "exhausted": {
        "N": "N_exhausted",
        "mu": "mu_Exhausted",
    },
}

def flatten_cell_types(raw: dict) -> dict:
    result = {}
    for ct_name, ct_params in raw.get("cell_types", {}).items():
        mapping = _CELL_TYPE_MAP.get(ct_name)
        if mapping is None:
            raise ValueError(f"Unknown cell type '{ct_name}' in config")
        for inner_key, flat_key in mapping.items():
            if inner_key in ct_params:
                result[flat_key] = ct_params[inner_key]
    for k, v in raw.items():
        if k != "cell_types":
            result[k] = v
    return result


def load_config_json(path: Path) -> dict:
    with open(path) as f:
        return _json.load(f)
