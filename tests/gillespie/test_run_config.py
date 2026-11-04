"""Run conditions: input values (before OMEGA scaling) next to the values actually used."""
import json

import pytest

from src.gillespie.config_report import config_rows, format_config
from src.gillespie.infrastructure.config.json_config import load_config_json
from src.gillespie.infrastructure.run_config_output import save_run_config
from src.gillespie.simulation_config import SimulationConfig


def _rows(cfg):
    return {p: (i, u, changed) for p, i, u, changed in config_rows(cfg)}


def test_input_values_are_unscaled_and_as_dict_is_scaled():
    cfg = SimulationConfig(OMEGA=4000)
    assert cfg.input_values["beta"] == 0.0004
    assert cfg.as_dict()["beta"] == pytest.approx(0.0004 / 4000)
    assert cfg.input_values["cell_parameters"]["mutated"]["K"] == 2
    assert cfg.as_dict()["cell_parameters"]["mutated"]["K"] == 8000


def test_input_values_include_json_and_defaults():
    cfg = SimulationConfig(OMEGA=4000, theta_I=0.001, cell_parameters={"base": {"N": 10}})
    assert cfg.input_values["theta_I"] == 0.001                       # given
    assert cfg.input_values["cell_parameters"]["base"]["N"] == 10     # given
    assert cfg.input_values["cell_parameters"]["base"]["mu"] == 0.002  # default


def test_rows_mark_exactly_the_scaled_values():
    rows = _rows(SimulationConfig(OMEGA=4000))
    changed = {p for p, (_, _, c) in rows.items() if c}
    assert {"theta_I", "beta", "C_half", "base.K", "mutated.K", "mutated.lambda0"} <= changed
    assert not rows["T_max"][2] and not rows["base.mu"][2]


def test_nothing_is_marked_without_scaling():
    assert not any(c for _, _, _, c in config_rows(SimulationConfig(scale=False)))


def test_format_config_has_input_and_used_columns():
    text = format_config(SimulationConfig(OMEGA=4000))
    header = text.splitlines()[0].split()
    assert header == ["parameter", "input", "used"]
    beta = next(line for line in text.splitlines() if line.startswith("beta "))
    assert beta.split() == ["beta", "0.0004", "1e-07", "*"]


def test_saved_table_has_input_and_used(tmp_path):
    txt, _ = save_run_config(tmp_path / "run_config", SimulationConfig(OMEGA=4000))
    assert txt.read_text() == format_config(SimulationConfig(OMEGA=4000)) + "\n"


def test_saved_input_reproduces_the_config(tmp_path):
    cfg = SimulationConfig(OMEGA=4000, theta_I=0.001,
                           cell_parameters={"fibroblast": {"crowding_pool": {"mutated": 0.5}}})
    _, js = save_run_config(tmp_path / "run_config", cfg)
    saved = json.loads(js.read_text())
    assert set(saved) == {"input", "used"}
    (tmp_path / "input.json").write_text(json.dumps(saved["input"]))
    again = SimulationConfig(**load_config_json(tmp_path / "input.json"))
    assert again.as_dict() == cfg.as_dict()


def test_json_default_id_lists_become_tuples(tmp_path):
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps({"cell_parameters": {"mutated": {"default_id": [-3]}}}))
    assert load_config_json(path)["cell_parameters"]["mutated"]["default_id"] == (-3,)
