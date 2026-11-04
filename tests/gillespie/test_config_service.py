import pytest

from src.gillespie.application.config_service import merge_and_build
from src.gillespie.simulation_config import SimulationConfig


class TestMergeAndBuild:
    def test_empty_dicts_uses_simulation_config_defaults(self):
        config = merge_and_build()
        assert isinstance(config, SimulationConfig)
        assert config.cell_parameters["base"].N == 4000

    def test_json_only(self):
        config = merge_and_build(
            json_data={
                "cell_parameters": {"base": {"N": 50}},
                "T_max": 500,
            }
        )
        assert config.cell_parameters["base"].N == 50
        assert config.T_max == 500

    def test_json_cell_parameters_preserve_values_and_fill_type_defaults(self):
        config = merge_and_build(
            json_data={"cell_parameters": {"base": {"N": 17, "K": 3}}},
            cli_overrides={"scale": False},
        )

        base = config.cell_parameters["base"]
        assert base.N == 17
        assert base.K == 3
        assert base.lambda0 == 0.005
        assert base.mu == 0.002
        assert base.nu == 0.00002
        assert base.next_mutation == "mutated"
        assert base.default_id == ()

        assert config.cell_parameters["immune"] == SimulationConfig(
            scale=False
        ).cell_parameters["immune"]
        assert set(config.cell_parameters) == {"base", "immune", "mutated", "exhausted", "fibroblast", "caf"}

    def test_cli_overrides_json(self):
        config = merge_and_build(
            json_data={"T_max": 500},
            cli_overrides={"T_max": 999},
        )
        assert config.T_max == 999

    def test_nested_cell_parameters_are_preserved(self):
        config = merge_and_build(
            json_data={
                "cell_parameters": {
                    "base": {"N": 10, "lambda0": 0.1, "mu": 0.01, "nu": 0.001, "K": 50},
                    "immune": {"N": 3, "lambda0": 0.2, "mu": 0.02, "K": 25},
                    "mutated": {"N": 1, "K": 500},
                    "exhausted": {"N": 0, "mu": 0.005},
                },
                "T_max": 100,
            }
        )
        assert config.cell_parameters["base"].N == 10
        assert config.cell_parameters["immune"].N == 3
        assert config.cell_parameters["mutated"].N == 1
        assert config.cell_parameters["exhausted"].N == 0
        assert config.cell_parameters["immune"].lambda0 == 0.2
        assert config.cell_parameters["exhausted"].mu == 0.005
        assert config.cell_parameters["immune"].K == 25 * config.OMEGA
        assert config.T_max == 100

    def test_legacy_cell_types_are_rejected(self):
        with pytest.raises(ValueError, match="cell_parameters"):
            merge_and_build(json_data={"cell_types": {"base": {"N": 10}}})

    def test_new_cell_type_without_default_id_is_rejected(self):
        with pytest.raises(ValueError, match="default_id"):
            merge_and_build(json_data={"cell_parameters": {"new_type": {"N": 10}}})

    def test_duplicate_default_id_is_rejected(self):
        with pytest.raises(ValueError, match="default_id"):
            merge_and_build(json_data={"cell_parameters": {"immune": {"default_id": (-2,)}}})
