import json

import pytest

from src.gillespie.infrastructure.config.json_config import load_config_json


def test_load_empty_json(tmp_path):
    path = tmp_path / "config.json"
    path.write_text("{}")
    assert load_config_json(path) == {}


def test_load_nested_config(tmp_path):
    path = tmp_path / "config.json"
    data = {
        "cell_parameters": {"base": {"N": 100, "lambda0": 0.005}},
        "T_max": 500,
    }
    path.write_text(json.dumps(data))
    assert load_config_json(path) == data


def test_calculator_metadata_is_removed(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"_calculator": {}, "T_max": 500}))
    assert load_config_json(path) == {"T_max": 500}


def test_legacy_cell_types_are_rejected(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"cell_types": {"base": {"N": 100}}}))
    with pytest.raises(ValueError, match="cell_parameters"):
        load_config_json(path)


def test_invalid_cell_parameters_are_rejected(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"cell_parameters": []}))
    with pytest.raises(TypeError, match="cell_parameters"):
        load_config_json(path)


def test_load_raises_on_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config_json(tmp_path / "does_not_exist.json")


def test_load_raises_on_invalid_json(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{invalid")
    with pytest.raises(json.JSONDecodeError):
        load_config_json(path)
