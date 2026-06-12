from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pytest

from src.gillespie.calibration.adapters.csv_output import CsvOutputAdapter
from src.gillespie.calibration.application.abc_smc_service import AbcSmcService
from src.gillespie.calibration.domain.calibration_config import (
    CalibrationConfig,
    Prior,
)
from src.gillespie.calibration.ports.calibration_ports import ReferenceData


class MockSimulator:
    def __init__(self, true_params: Optional[Dict[str, float]] = None):
        self.true_params = true_params or {"lambda0": 0.05, "mu0": 0.02}
        self.call_count = 0

    def run(self, params: Dict[str, float]) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        self.call_count += 1
        times = np.array([0.0, 10.0, 20.0])
        pop = {
            "base": np.array([50.0, 55.0, 60.0]),
            "mutated": np.array([0.0, 5.0, 10.0]),
        }
        return times, pop


class MockReference:
    def __init__(self):
        self.ref_stats = np.array([55.0, 60.0, 5.0, 10.0])  # mean_base, final_base, mean_mut, final_mut
        self.precision_weights = np.array([1.0, 1.0, 1.0, 1.0])
        self.stat_names = ["mean", "final"]
        self.cell_types = ["base", "mutated"]

    def load(self) -> ReferenceData:
        return ReferenceData(
            ref_stats=self.ref_stats.copy(),
            precision_weights=self.precision_weights.copy(),
            stat_names=self.stat_names,
            cell_types=self.cell_types,
        )


class MockOutput:
    def __init__(self):
        self.generations = []
        self.eps_history = []
        self.diagnostics = []
        self.manifest = None

    def save_generation(self, gen, particles, weights, distances):
        self.generations.append((gen, particles, weights, distances))

    def save_epsilon_schedule(self, eps_history):
        self.eps_history = eps_history

    def save_generation_diagnostics(self, diagnostics):
        self.diagnostics = diagnostics

    def save_manifest(self, calibration_config, priors):
        self.manifest = {"config": calibration_config, "priors": priors}

    def find_last_generation(self):
        return None

    def load_generation(self, gen):
        raise NotImplementedError


class TestAbcSmcService:
    def test_initialization(self):
        config = CalibrationConfig(n_particles=10, n_generations=1)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        sim = MockSimulator()
        ref = MockReference()
        out = MockOutput()
        rng = np.random.default_rng(42)

        service = AbcSmcService(
            config=config, priors=priors, simulator=sim, reference=ref, output=out, rng=rng
        )
        assert service.param_names == ["lambda0", "mu0"]
        assert len(service.ref_stats) == 4

    def test_sample_prior_respects_bounds(self):
        config = CalibrationConfig(n_particles=10, n_generations=1)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        ref = MockReference()
        out = MockOutput()
        rng = np.random.default_rng(42)
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=rng
        )

        for _ in range(500):
            theta = service._sample_prior()
            if theta is None:
                continue
            assert 0.001 <= theta["lambda0"] <= 0.1
            assert 0.001 <= theta["mu0"] <= 0.1

    def test_sample_prior_constraint(self):
        def constraint(theta):
            return theta["mu0"] < theta["lambda0"]

        config = CalibrationConfig(n_particles=10, n_generations=1)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1, constraint=constraint),
            "mu0": Prior("mu0", 0.001, 0.1, constraint=constraint),
        }
        ref = MockReference()
        out = MockOutput()
        rng = np.random.default_rng(42)
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=rng
        )

        for _ in range(500):
            theta = service._sample_prior()
            if theta is None:
                continue
            assert theta["mu0"] < theta["lambda0"]

    def test_run_abc_smc_two_generations(self):
        config = CalibrationConfig(
            n_particles=5, n_generations=2, alpha=0.6, epsilon_0=1e8,
            max_attempts=500, seed=42, n_workers=1,
        )
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        ref = MockReference()
        out = MockOutput()
        rng = np.random.default_rng(42)
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=rng,
        )

        result = service.run_abc_smc()
        assert result.n_generations() == 2
        assert len(result.particles[0]) == 5
        assert len(result.particles[1]) == 5
        assert len(result.epsilon_history) == 2
        assert len(result.diagnostics) == 2
        assert result.diagnostics[0]["acceptance_rate"] > 0

    def test_weights_normalize_to_one(self):
        config = CalibrationConfig(n_particles=5, n_generations=2, seed=42)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        ref = MockReference()
        out = MockOutput()
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=np.random.default_rng(42),
        )

        result = service.run_abc_smc()
        for g in range(result.n_generations()):
            assert abs(np.sum(result.weights[g]) - 1.0) < 1e-10
            assert np.all(result.weights[g] > 0)

    def test_epsilon_decreases(self):
        config = CalibrationConfig(n_particles=5, n_generations=3, alpha=0.5, seed=42)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        ref = MockReference()
        out = MockOutput()
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=np.random.default_rng(42),
        )

        result = service.run_abc_smc()
        for g in range(1, result.n_generations()):
            assert result.epsilon_history[g] <= result.epsilon_history[g - 1]

    def test_prior_sensitivity(self):
        config = CalibrationConfig(n_particles=10, n_generations=1, seed=42)
        priors = {
            "lambda0": Prior("lambda0", 0.001, 0.1),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        ref = MockReference()
        out = MockOutput()
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=np.random.default_rng(42),
        )

        sens = service.prior_sensitivity(
            width_factors=[0.5, 1.0, 1.5], n_samples=50
        )
        assert 0.5 in sens
        assert 1.0 in sens
        assert 1.5 in sens
        for factor, info in sens.items():
            assert info["n_samples"] <= 50
            assert "acceptance_rate" in info
            assert "mean_distance" in info


class TestIntegrationCsvRoundtrip:
    def test_csv_output_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            out = CsvOutputAdapter(str(tmp))

            particles = [{"lambda0": 0.05}, {"lambda0": 0.08}]
            weights = np.array([0.6, 0.4])
            distances = np.array([10.0, 20.0])
            out.save_generation(0, particles, weights, distances)

            out.save_epsilon_schedule([100.0, 50.0])
            out.save_generation_diagnostics([
                {"gen": 0, "ess": 200.0},
            ])
            out.save_manifest(
                calibration_config={"n_particles": 2},
                priors={"lambda0": type("obj", (object,), {"param_name": "lambda0", "lo": 0.001, "hi": 0.1, "distribution": "uniform"})()},
            )

            assert (tmp / "gen_00.csv").exists()
            assert (tmp / "epsilon_schedule.csv").exists()
            assert (tmp / "generation_diagnostics.csv").exists()
            assert (tmp / "run_manifest.json").exists()

            loaded_p, loaded_w, loaded_d = out.load_generation(0)
            assert loaded_p == particles
            assert np.allclose(loaded_w, weights)
            assert np.allclose(loaded_d, distances)


class TestManualPriorRejectCount:
    def test_prior_reject_counted(self):
        priors = {
            "lambda0": Prior("lambda0", 0.5, 0.6),
            "mu0": Prior("mu0", 0.001, 0.1),
        }
        config = CalibrationConfig(
            n_particles=5, n_generations=1, epsilon_0=1e8, seed=42,
        )
        ref = MockReference()
        out = MockOutput()
        service = AbcSmcService(
            config=config, priors=priors, simulator=MockSimulator(),
            reference=ref, output=out, rng=np.random.default_rng(42),
        )

        result = service.run_abc_smc()
        assert result.diagnostics[0]["prior_reject"] >= 0
