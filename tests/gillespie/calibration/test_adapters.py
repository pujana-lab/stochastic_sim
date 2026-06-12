from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gillespie.calibration.adapters.csv_output import CsvOutputAdapter
from src.gillespie.calibration.adapters.csv_reference import CsvReferenceAdapter


class TestCsvOutputAdapter:
    @pytest.fixture
    def tmp_out(self):
        with tempfile.TemporaryDirectory() as d:
            yield Path(d)

    def test_save_generation(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        particles = [{"lambda0": 0.05, "mu0": 0.02}, {"lambda0": 0.08, "mu0": 0.03}]
        weights = np.array([0.6, 0.4])
        distances = np.array([10.0, 20.0])
        adapter.save_generation(0, particles, weights, distances)

        path = tmp_out / "gen_00.csv"
        assert path.exists()
        df = pd.read_csv(path)
        assert list(df.columns) == ["lambda0", "mu0", "weight", "distance"]
        assert len(df) == 2
        assert df["weight"].iloc[0] == 0.6

    def test_save_epsilon_schedule(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        adapter.save_epsilon_schedule([100.0, 50.0, 25.0])

        path = tmp_out / "epsilon_schedule.csv"
        assert path.exists()
        df = pd.read_csv(path)
        assert list(df.columns) == ["generation", "epsilon"]
        assert df["epsilon"].iloc[2] == 25.0

    def test_save_diagnostics(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        diag = [
            {"gen": 0, "ess": 200.0, "acceptance_rate": 0.5},
            {"gen": 1, "ess": 180.0, "acceptance_rate": 0.3},
        ]
        adapter.save_generation_diagnostics(diag)

        path = tmp_out / "generation_diagnostics.csv"
        assert path.exists()
        df = pd.read_csv(path)
        assert len(df) == 2

    def test_save_manifest(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        adapter.save_manifest(
            calibration_config={"n_particles": 50, "n_generations": 3},
            priors={"lambda0": type("obj", (object,), {"param_name": "lambda0", "lo": 0.001, "hi": 0.1, "distribution": "uniform"})()},
        )

        path = tmp_out / "run_manifest.json"
        assert path.exists()
        with open(path) as f:
            manifest = json.load(f)
        assert "timestamp" in manifest
        assert manifest["calibration_config"]["n_particles"] == 50
        assert "lambda0" in manifest["priors"]

    def test_find_last_generation(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        assert adapter.find_last_generation() is None

        pd.DataFrame({"x": [1]}).to_csv(tmp_out / "gen_00.csv", index=False)
        pd.DataFrame({"x": [1]}).to_csv(tmp_out / "gen_02.csv", index=False)
        assert adapter.find_last_generation() == 2

    def test_load_generation_roundtrip(self, tmp_out):
        adapter = CsvOutputAdapter(str(tmp_out))
        particles = [{"lambda0": 0.05, "mu0": 0.02}, {"lambda0": 0.08, "mu0": 0.03}]
        weights = np.array([0.6, 0.4])
        distances = np.array([10.0, 20.0])
        adapter.save_generation(0, particles, weights, distances)

        loaded_p, loaded_w, loaded_d = adapter.load_generation(0)
        assert loaded_p == particles
        assert np.allclose(loaded_w, weights)
        assert np.allclose(loaded_d, distances)


class TestCsvReferenceAdapter:
    def test_load_homeostasis_csv(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            csv_path = tmp / "reference_test.csv"
            with open(csv_path, "w") as f:
                f.write("time,type,mean_N,std_N,n_replicates\n")
                f.write("0,base,50.0,5.0,50\n")
                f.write("10,base,55.0,6.0,50\n")
                f.write("0,mutated,0.0,0.0,50\n")
                f.write("10,mutated,5.0,1.0,50\n")

            adapter = CsvReferenceAdapter(str(csv_path))
            data = adapter.load()
            assert len(data.stat_names) == 2
            assert "base" in data.cell_types
            assert "mutated" in data.cell_types
            assert len(data.ref_stats) == 4  # 2 types × 2 stats
            assert len(data.precision_weights) == 4
            assert np.all(data.precision_weights == 1.0)  # no weights file

    def test_with_weights(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            csv_path = tmp / "ref.csv"
            w_path = tmp / "weights.json"
            with open(csv_path, "w") as f:
                f.write("time,type,mean_N,std_N,n_replicates\n")
                f.write("0,base,50.0,5.0,50\n")
                f.write("10,base,55.0,6.0,50\n")
            with open(w_path, "w") as f:
                json.dump({"mean_base": 0.04, "final_base": 0.02}, f)

            adapter = CsvReferenceAdapter(str(csv_path), weights_path=str(w_path))
            data = adapter.load()
            assert np.any(data.precision_weights != 1.0)

    def test_missing_columns_raises(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            csv_path = tmp / "bad.csv"
            with open(csv_path, "w") as f:
                f.write("time,value\n0,10\n")
            adapter = CsvReferenceAdapter(str(csv_path))
            with pytest.raises(ValueError, match="missing required columns"):
                adapter.load()
