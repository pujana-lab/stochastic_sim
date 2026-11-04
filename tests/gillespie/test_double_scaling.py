"""SimulationConfig must never scale by OMEGA twice (it overwrites its fields when it scales)."""
import dataclasses

import pytest

from src.gillespie.simulation_config import SimulationConfig


def test_scaled_config_marks_its_cell_parameters():
    cfg = SimulationConfig(OMEGA=4000)
    assert all(c.scaled for c in cfg.cell_parameters.values())


def test_unscaled_config_is_not_marked():
    cfg = SimulationConfig(scale=False)
    assert not any(c.scaled for c in cfg.cell_parameters.values())


def test_replace_on_a_scaled_config_raises():
    """Without the guard this gave beta / OMEGA**2 and K * OMEGA**2 silently."""
    cfg = SimulationConfig(OMEGA=4000)
    with pytest.raises(ValueError, match="already scaled"):
        dataclasses.replace(cfg, T_max=500)


def test_reusing_scaled_cell_parameters_raises():
    cfg = SimulationConfig(OMEGA=4000)
    with pytest.raises(ValueError, match="already scaled"):
        SimulationConfig(OMEGA=4000, cell_parameters=cfg.cell_parameters)


def test_replace_on_an_unscaled_config_is_allowed():
    cfg = SimulationConfig(scale=False)
    assert dataclasses.replace(cfg, T_max=500).T_max == 500


def test_building_from_unscaled_values_twice_gives_the_same_config():
    """The supported way to vary a parameter: rebuild from the unscaled inputs."""
    raw = {"OMEGA": 4000, "theta_I": 0.0005, "cell_parameters": {"base": {"K": 1}}}
    a, b = SimulationConfig(**raw), SimulationConfig(**raw)
    assert a.theta_I == b.theta_I == pytest.approx(0.0005 / 4000)
    assert a.cell_parameters["base"].K == b.cell_parameters["base"].K == 4000
