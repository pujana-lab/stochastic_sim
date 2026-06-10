import pytest

from src.gillespie.application.config_service import merge_and_build, build_config_from_cli
from src.gillespie.simulation_config import SimulationConfig


class TestMergeAndBuild:
    def test_empty_dicts_uses_simulation_config_defaults(self):
        config = merge_and_build()
        assert isinstance(config, SimulationConfig)
        assert config.N0 == 100

    def test_json_only(self):
        config = merge_and_build(json_data={"N0": 50, "T_max": 500})
        assert config.N0 == 50
        assert config.T_max == 500
        assert config.lambda0 == 0.005

    def test_cli_overrides_json(self):
        config = merge_and_build(
            json_data={"N0": 50, "T_max": 500},
            cli_overrides={"N0": 999},
        )
        assert config.N0 == 999
        assert config.T_max == 500

    def test_cli_only(self):
        config = merge_and_build(cli_overrides={"N0": 777, "lambda0": 0.99})
        assert config.N0 == 777
        assert config.lambda0 == 0.99

    def test_cli_overrides_none_values(self):
        config = merge_and_build(
            json_data={"seed": None},
            cli_overrides={"seed": 42},
        )
        assert config.seed == 42
