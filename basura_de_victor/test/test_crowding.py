import pytest
from src.simulation_config import SimulationConfig
from src.clone import Clone
from src.crowding_strategy import SimpleCrowding, AdaptedCrowding


# ── fixtures ──────────────────────────────────────────────────────────────────

def make_clone(birth_rate=0.5, death_rate=0.2, N=10) -> Clone:
    return Clone(
        clone_id=(),
        N=N,
        birth_rate=birth_rate,
        death_rate=death_rate,
        mutation_rate=0.0,
    )


def simple_config(**kwargs) -> SimulationConfig:
    defaults = dict(use_logistic=True, use_logistic_adapted=False, K0=100.0, decline=0.0, Kmin=1.0)
    defaults.update(kwargs)
    return SimulationConfig(**defaults)


# ── SimpleCrowding ─────────────────────────────────────────────────────────────

class TestSimpleCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = SimpleCrowding(cfg)
        assert strategy.crowding(make_clone(), t=0.0, total_N=50) == 1.0

    def test_empty_population_returns_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        assert strategy.crowding(make_clone(), t=0.0, total_N=0) == 1.0

    def test_at_capacity_returns_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        result = strategy.crowding(make_clone(), t=0.0, total_N=100)
        assert result == pytest.approx(0.0)

    def test_above_capacity_clamped_to_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        result = strategy.crowding(make_clone(), t=0.0, total_N=150)
        assert result == 0.0

    def test_crowding_between_zero_and_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        result = strategy.crowding(make_clone(), t=0.0, total_N=50)
        
        assert 0.0 <= result <= 1.0

    def test_K_declines_with_time(self):
        cfg = simple_config(K0=100.0, decline=5.0, Kmin=1.0)
        strategy = SimpleCrowding(cfg)
        c0 = strategy.crowding(make_clone(), t=0.0, total_N=50)
        c10 = strategy.crowding(make_clone(), t=10.0, total_N=50)
        # K decreases → crowding factor decreases (more crowded)
        assert c10 <= c0

    def test_K_never_below_Kmin(self):
        cfg = simple_config(K0=100.0, decline=1000.0, Kmin=10.0)
        strategy = SimpleCrowding(cfg)
        # at t=1000, K is clamped to Kmin=10; N=5 < K → crowding > 0
        result = strategy.crowding(make_clone(), t=1000.0, total_N=5)
        assert result >= 0.0


# ── AdaptedCrowding ────────────────────────────────────────────────────────────

class TestAdaptedCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = AdaptedCrowding(cfg)
        assert strategy.crowding(make_clone(), t=0.0, total_N=50) == 1.0

    def test_result_is_non_negative(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        result = strategy.crowding(make_clone(birth_rate=0.5, death_rate=0.2), t=0.0, total_N=200)
        assert result >= 0.0

    def test_result_at_most_one(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        result = strategy.crowding(make_clone(birth_rate=0.5, death_rate=0.2), t=0.0, total_N=1)
        assert result <= 1.0

    def test_higher_fitness_yields_higher_effective_birth_rate(self):
        """Clone with higher birth_rate has higher effective birth rate despite lower crowding."""
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        c_low  = make_clone(birth_rate=0.3, death_rate=0.2)
        c_high = make_clone(birth_rate=0.8, death_rate=0.2)
        rb_low  = c_low.birth_rate_effective(strategy.crowding(c_low,  t=0.0, total_N=50))
        rb_high = c_high.birth_rate_effective(strategy.crowding(c_high, t=0.0, total_N=50))
        assert rb_high > rb_low