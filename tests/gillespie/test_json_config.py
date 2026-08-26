import json
import pytest
from pathlib import Path

from src.gillespie.infrastructure.config.json_config import load_config_json, flatten_cell_types
from src.gillespie.simulation_config import SimulationConfig


def test_load_empty_json(tmp_path):
    path = tmp_path / "config.json"
    path.write_text("{}")
    assert load_config_json(path) == {}


def test_load_simple_values(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"N0": 100, "lambda0": 0.5}))
    result = load_config_json(path)
    assert result == {"N0": 100, "lambda0": 0.5}


def test_load_full_config(tmp_path):
    path = tmp_path / "config.json"
    data = {"N0": 50, "T_max": 500, "seed": 42, "use_logistic": True}
    path.write_text(json.dumps(data))
    assert load_config_json(path) == data


def test_load_raises_on_missing_file(tmp_path):
    path = tmp_path / "does_not_exist.json"
    with pytest.raises(FileNotFoundError):
        load_config_json(path)


def test_load_raises_on_invalid_json(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{invalid")
    with pytest.raises(json.JSONDecodeError):
        load_config_json(path)


# ── flatten_cell_types ─────────────────────────────────────────────────────────

def test_flatten_empty_dict():
    assert flatten_cell_types({}) == {}


def test_flatten_passes_through_non_cell_type_keys():
    assert flatten_cell_types({"T_max": 500, "seed": 42}) == {"T_max": 500, "seed": 42}


def test_flatten_cell_types_base():
    raw = {"cell_types": {"base": {"N": 100, "lambda_0": 0.005, "mu": 0.002, "nu": 0.0002, "K": 50}}}
    result = flatten_cell_types(raw)
    assert result["cell_parameters"]["base"].N == 100
    assert result["cell_parameters"]["base"].lambda0 == 0.005
    assert result["cell_parameters"]["base"].mu == 0.002
    assert result["cell_parameters"]["base"].nu == 0.0002
    assert result["cell_parameters"]["base"].K == 50


def test_flatten_cell_types_all_types():
    raw = {
        "cell_types": {
            "base": {"N": 10, "lambda_0": 0.1, "mu": 0.01, "nu": 0.001, "K": 100},
            "immune": {"N": 5, "lambda_0": 0.2, "omega_exhaust": 0.02, "K": 50},
            "mutated": {"N": 3, "K": 200, "fitness_gain": 0.4},
            "exhausted": {"N": 1, "mu": 0.05},
        },
        "T_max": 1000,
    }
    result = flatten_cell_types(raw)
    assert result["cell_parameters"]["base"].N == 10
    assert result["cell_parameters"]["immune"].N == 5
    assert result["cell_parameters"]["mutated"].N == 3
    assert result["cell_parameters"]["exhausted"].N == 1
    assert result["cell_parameters"]["base"].lambda0 == 0.1
    assert result["cell_parameters"]["immune"].lambda0 == 0.2
    assert result["cell_parameters"]["base"].mu == 0.01
    assert result["cell_parameters"]["immune"].omega_exhaust == 0.02
    config = SimulationConfig(**result)
    assert config.cell_parameters["immune"].omega_exhaust == 0.02 / 4000.0
    assert result["cell_parameters"]["exhausted"].mu == 0.05
    assert result["cell_parameters"]["base"].nu == 0.001
    assert result["cell_parameters"]["base"].K == 100
    assert result["cell_parameters"]["immune"].K == 50
    assert result["cell_parameters"]["mutated"].K == 200
    assert result["cell_parameters"]["mutated"].fitness_gain == 0.4
    assert result["T_max"] == 1000


def test_flatten_unknown_cell_type_raises():
    raw = {"cell_types": {"unknown_type": {"N": 0}}}
    with pytest.raises(ValueError, match="Unknown cell type"):
        flatten_cell_types(raw)


def test_flatten_unknown_inner_key_ignored():
    raw = {"cell_types": {"base": {"N": 10, "nonexistent": 999}}}
    result = flatten_cell_types(raw)
    assert result["cell_parameters"]["base"].N == 10
    assert "nonexistent" not in result["cell_parameters"]["base"].__dict__


def test_flatten_combined_with_global_keys():
    raw = {
        "cell_types": {"base": {"N": 10, "lambda_0": 0.1}},
        "T_max": 500,
        "fitness_gain": 0.3,
    }
    result = flatten_cell_types(raw)
    assert result["cell_parameters"]["base"].N == 10
    assert result["cell_parameters"]["base"].lambda0 == 0.1
    assert result["T_max"] == 500
    assert result["fitness_gain"] == 0.3


def test_flatten_current_nested_schema_builds_valid_simulation_config():
    raw = {
        "cell_types": {
            "base": {"N": 4000, "K": 1, "lambda_0": 0.005, "mu": 0.002, "nu": 0.00002, "next_mutation": "mutated"},
            "immune": {"N": 2000, "K": 0.5, "lambda_0": 0.005, "omega_exhaust": 0.003, "mu": 0.0},
            "mutated": {"N": 1, "K": 2, "fitness_gain": 0.2},
        },
        "T_max": 1000,
        "use_logistic": True,
    }
    result = flatten_cell_types(raw)
    config = SimulationConfig(**result)
    assert config.cell_parameters["base"].N == 4000
    assert config.cell_parameters["base"].next_mutation == "mutated"
    assert config.cell_parameters["mutated"].fitness_gain == 0.2
    assert config.T_max == 1000
