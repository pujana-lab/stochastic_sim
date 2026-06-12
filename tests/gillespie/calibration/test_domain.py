from __future__ import annotations

import numpy as np
import pytest

from src.gillespie.calibration.domain.abc_result import AbcResult
from src.gillespie.calibration.domain.calibration_config import (
    CalibrationConfig,
    Prior,
)
from src.gillespie.calibration.domain.distance_metric import weighted_sse
from src.gillespie.calibration.domain.summary_statistic import (
    compute_summary_stats,
)


class TestPrior:
    def test_uniform_sample_within_bounds(self):
        p = Prior("lambda0", 0.001, 0.1)
        rng = np.random.default_rng(42)
        vals = [p.sample(rng) for _ in range(1000)]
        assert all(p.lo <= v <= p.hi for v in vals)
        assert 0.03 < np.mean(vals) < 0.07

    def test_uniform_pdf(self):
        p = Prior("lambda0", 0.0, 1.0)
        assert p.pdf(0.5) == 1.0
        assert p.pdf(-0.1) == 0.0
        assert p.pdf(1.1) == 0.0

    def test_lognorm_sample_positive(self):
        p = Prior("nu0", 1e-6, 1e-2, distribution="lognorm")
        rng = np.random.default_rng(42)
        vals = [p.sample(rng) for _ in range(1000)]
        assert all(v > 0 for v in vals)
        assert all(p.lo <= v <= p.hi for v in vals)

    def test_truncnorm_sample_within_bounds(self):
        p = Prior("theta_I", 0.0, 0.01, distribution="truncnorm")
        rng = np.random.default_rng(42)
        vals = [p.sample(rng) for _ in range(1000)]
        assert all(p.lo <= v <= p.hi for v in vals)

    def test_invalid_distribution_raises(self):
        with pytest.raises(ValueError, match="Invalid distribution"):
            Prior("x", 0, 1, distribution="invalid_dist")

    def test_lo_gte_hi_raises(self):
        with pytest.raises(ValueError, match="lo.*must be < hi"):
            Prior("x", 1.0, 0.5)

    def test_pdf_sum_uniform(self):
        p = Prior("x", 0.0, 2.0)
        assert abs(p.pdf(0.5) - 0.5) < 1e-10

    def test_constraint_validated_on_proposal(self):
        def constraint(theta):
            return theta["high_delta"] > theta["low_delta"]

        p_high = Prior("high_delta", 0.0, 1.0, constraint=constraint)
        p_low = Prior("low_delta", 0.0, 1.0, constraint=constraint)
        rng = np.random.default_rng(42)
        valid = False
        for _ in range(500):
            theta = {
                "high_delta": p_high.sample(rng),
                "low_delta": p_low.sample(rng),
            }
            if constraint(theta):
                valid = True
                break
        assert valid


class TestCalibrationConfig:
    def test_default_values(self):
        c = CalibrationConfig()
        assert c.n_particles == 200
        assert c.n_generations == 5
        assert c.alpha == 0.6
        assert c.epsilon_0 == 1e6
        assert c.epsilon_final == 1.0
        assert c.max_attempts == 100000
        assert c.n_workers == 4
        assert c.n_reps == 1
        assert c.seed == 42
        assert c.output_dir == "calibration_results"

    def test_custom_config(self):
        c = CalibrationConfig(
            n_particles=50, n_generations=3, alpha=0.5, n_workers=2, output_dir="test_out"
        )
        assert c.n_particles == 50
        assert c.output_dir == "test_out"

    def test_config_is_frozen(self):
        c = CalibrationConfig()
        with pytest.raises(Exception):
            c.n_particles = 100


class TestSummaryStatistic:
    def test_default_stats(self):
        times = np.array([0, 10, 20, 30])
        pop = {
            "base": np.array([50, 55, 60, 65]),
            "mutated": np.array([0, 5, 10, 15]),
        }
        stats = compute_summary_stats(pop, times)
        assert len(stats) == 4  # 2 types × 2 stats (mean, final)
        assert np.isclose(stats[0], np.mean(pop["base"]))  # mean_base
        assert np.isclose(stats[1], 65.0)  # final_base
        assert np.isclose(stats[3], 15.0)  # final_mutated

    def test_single_stat(self):
        times = np.array([0, 10])
        pop = {"base": np.array([100, 200])}
        stats = compute_summary_stats(pop, times, stat_names=["final"])
        assert len(stats) == 1
        assert stats[0] == 200.0

    def test_custom_cell_types(self):
        times = np.array([0, 10])
        pop = {
            "base": np.array([100, 200]),
            "immune": np.array([50, 30]),
            "exhausted": np.array([0, 0]),
        }
        stats = compute_summary_stats(pop, times, cell_types=["base", "immune"])
        assert len(stats) == 4  # 2 types × 2 default stats (mean, final)

    def test_unknown_stat_raises(self):
        times = np.array([0, 10])
        pop = {"base": np.array([1, 2])}
        with pytest.raises(ValueError, match="Unknown summary statistic"):
            compute_summary_stats(pop, times, stat_names=["invalid_stat"])


class TestDistanceMetric:
    def test_uniform_weights(self):
        sim = np.array([10.0, 20.0])
        ref = np.array([12.0, 18.0])
        w = np.array([1.0, 1.0])
        d = weighted_sse(sim, ref, w)
        assert abs(d - ((10-12)**2 + (20-18)**2)) < 1e-10

    def test_precision_weights(self):
        sim = np.array([10.0])
        ref = np.array([12.0])
        w = np.array([0.5])
        d = weighted_sse(sim, ref, w)
        assert abs(d - 0.5 * (10-12)**2) < 1e-10

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="Shape mismatch"):
            weighted_sse(np.array([1.0, 2.0]), np.array([1.0]), np.array([1.0]))

    def test_identical_stats(self):
        sim = np.array([5.0, 10.0, 15.0])
        ref = np.array([5.0, 10.0, 15.0])
        w = np.array([1.0, 2.0, 3.0])
        assert weighted_sse(sim, ref, w) == 0.0


class TestAbcResult:
    def test_empty_result(self):
        r = AbcResult()
        assert r.n_generations() == 0

    def test_add_generation(self):
        r = AbcResult()
        r.particles.append([{"lambda0": 0.05}, {"lambda0": 0.08}])
        r.weights.append(np.array([0.5, 0.5]))
        r.distances.append(np.array([10.0, 20.0]))
        r.diagnostics.append({"gen": 0})
        r.epsilon_history.append(100.0)
        assert r.n_generations() == 1
        assert r.ess(0) == pytest.approx(2.0)  # 1/(0.5^2 + 0.5^2)

    def test_summary(self):
        r = AbcResult()
        r.particles.append([{"x": 1.0}])
        r.weights.append(np.array([1.0]))
        r.distances.append(np.array([5.0]))
        r.diagnostics.append({"gen": 0, "attempts": 10})
        r.epsilon_history.append(100.0)
        s = r.summary(0)
        assert s["n_particles"] == 1
        assert s["mean_distance"] == 5.0
        assert s["generation"] == 0

    def test_all_summaries(self):
        r = AbcResult()
        for g in range(3):
            r.particles.append([{"x": float(g)}])
            r.weights.append(np.array([1.0]))
            r.distances.append(np.array([float(g * 10)]))
            r.diagnostics.append({"gen": g})
            r.epsilon_history.append(100.0)
        summaries = r.all_summaries()
        assert len(summaries) == 3
        assert summaries[1]["generation"] == 1
