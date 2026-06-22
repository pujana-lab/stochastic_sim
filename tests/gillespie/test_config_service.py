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

    def test_grouped_cell_types_format(self):
        config = merge_and_build(json_data={
            "cell_types": {
                "base": {"N": 10, "lambda_0": 0.1, "mu": 0.01, "nu": 0.001, "K": 50},
                "immune": {"N": 3, "lambda_0": 0.2, "mu": 0.02, "K": 25},
                "mutant": {"N": 1, "K": 500},
                "exhausted": {"N": 0, "mu": 0.005},
            },
            "T_max": 100,
        })
        assert config.N0 == 10
        assert config.N_immune == 3
        assert config.N_mutant == 1
        assert config.N_exhausted == 0
        assert config.lambda0 == 0.1
        assert config.lambda_Immune == 0.2
        assert config.mu0 == 0.01
        assert config.mu_Immune == 0.02
        assert config.mu_Exhausted == 0.005
        assert config.nu0 == 0.001
        assert config.K0 == 50
        assert config.K_immune == 25
        assert config.K_mutant == 500
        assert config.T_max == 100
