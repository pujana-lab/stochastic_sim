import json
import pytest
from pathlib import Path

from src.gillespie.infrastructure.config.json_config import load_config_json


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
