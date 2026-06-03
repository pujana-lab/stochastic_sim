import pytest
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone, WildTypeClone
from src.gillespie.crowding_strategy import SimpleCrowding, AdaptedCrowding
from src.gillespie.tissue_state import TissueState


# ── fixtures ──────────────────────────────────────────────────────────────────

def make_clone(birth_rate=0.5, death_rate=0.2, N=10, K=100.0) -> WildTypeClone:
    clone = WildTypeClone(
        clone_id=(),
        N=N,
        config=SimulationConfig(lambda0=birth_rate, mu0=death_rate)
    )
    clone.K = K
    return clone


def make_tissue_state(N: int = 0) -> TissueState:
    """Create a TissueState with N base (WildType) cells."""
    if N > 0:
        clone = WildTypeClone(clone_id=(0,), N=N, config=SimulationConfig())
        return TissueState(clones={(0,): clone})
    return TissueState(clones={})


def simple_config(**kwargs) -> SimulationConfig:
    defaults = dict(use_logistic=True, use_logistic_adapted=False, K0=100.0, decline=0.0, Kmin=1.0)
    defaults.update(kwargs)
    return SimulationConfig(**defaults)


# ── SimpleCrowding ─────────────────────────────────────────────────────────────

class TestSimpleCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state()
        assert strategy.crowding(make_clone(), t=0.0, tissue_state=ts) == 1.0

    def test_empty_population_returns_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=0)  # empty population → crowding_numerator = 0
        assert strategy.crowding(make_clone(), t=0.0, tissue_state=ts) == 1.0

    def test_at_capacity_returns_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=100)  # at capacity
        result = strategy.crowding(make_clone(), t=0.0, tissue_state=ts)
        assert result == pytest.approx(0.0)

    def test_above_capacity_clamped_to_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=150)  # above capacity
        result = strategy.crowding(make_clone(), t=0.0, tissue_state=ts)
        assert result == 0.0

    def test_crowding_between_zero_and_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=50)  # between 0 and K
        result = strategy.crowding(make_clone(), t=0.0, tissue_state=ts)
        assert 0.0 <= result <= 1.0

    def test_K_declines_with_time(self):
        cfg = simple_config(K0=100.0, decline=5.0, Kmin=1.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=10)
        c0 = strategy.crowding(make_clone(), t=0.0, tissue_state=ts)
        c10 = strategy.crowding(make_clone(), t=10.0, tissue_state=ts)
        # K decreases → crowding factor decreases (more crowded)
        assert c10 <= c0

    def test_K_never_below_Kmin(self):
        cfg = simple_config(K0=100.0, decline=1000.0, Kmin=10.0)
        strategy = SimpleCrowding(cfg)
        ts = make_tissue_state(N=5)
        # at t=1000, K is clamped to K_min; result must be non-negative
        result = strategy.crowding(make_clone(), t=1000.0, tissue_state=ts)
        assert result >= 0.0


# ── AdaptedCrowding ────────────────────────────────────────────────────────────

class TestAdaptedCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = AdaptedCrowding(cfg)
        ts = make_tissue_state()
        assert strategy.crowding(make_clone(), t=0.0, tissue_state=ts) == 1.0

    def test_result_is_non_negative(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        ts = make_tissue_state(N=50)
        result = strategy.crowding(make_clone(birth_rate=0.5, death_rate=0.2), t=0.0, tissue_state=ts)
        assert result >= 0.0

    def test_result_at_most_one(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        ts = make_tissue_state(N=50)
        result = strategy.crowding(make_clone(birth_rate=0.5, death_rate=0.2), t=0.0, tissue_state=ts)
        assert result <= 1.0

    def test_higher_fitness_yields_higher_effective_birth_rate(self):
        """Clone with higher birth_rate has higher effective birth rate."""
        c_low  = make_clone(birth_rate=0.3, death_rate=0.2)
        c_high = make_clone(birth_rate=0.8, death_rate=0.2)
        ts = make_tissue_state(N=10)
        rb_low  = c_low.birth_rate_effective(tissue_state=ts, crowding=1.0)
        rb_high = c_high.birth_rate_effective(tissue_state=ts, crowding=1.0)
        assert rb_high > rb_low
