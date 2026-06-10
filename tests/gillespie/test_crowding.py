import pytest
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone, WildTypeClone
from src.gillespie.crowding_strategy import SimpleCrowding, AdaptedCrowding
from src.gillespie.tissue_state import TissueState


# ── fixtures ──────────────────────────────────────────────────────────────────

def make_clone(birth_rate=0.5, death_rate=0.2, N=10) -> Clone:
    cfg = SimulationConfig(lambda0=birth_rate, mu0=death_rate)
    clone = Clone(clone_id=(), N=N, config=cfg)
    clone.K = cfg.K0
    return clone


def make_wt_clone(N=10, birth_rate=0.5, death_rate=0.2, **kw) -> WildTypeClone:
    cfg = SimulationConfig(lambda0=birth_rate, mu0=death_rate, K0=kw.pop("K0", 100), **kw)
    return WildTypeClone(clone_id=(), N=N, config=cfg)


def make_tissue_state(clone) -> TissueState:
    return TissueState(clones={clone.clone_id: clone})


def simple_config(**kwargs) -> SimulationConfig:
    defaults = dict(use_logistic=True, use_logistic_adapted=False, K0=100.0, decline=0.0, Kmin=1.0)
    defaults.update(kwargs)
    return SimulationConfig(**defaults)


# ── SimpleCrowding ─────────────────────────────────────────────────────────────

class TestSimpleCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = SimpleCrowding(cfg)
        c = make_clone()
        assert strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c)) == 1.0

    def test_empty_population_returns_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        c = make_clone(N=0)
        assert strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c)) == 1.0

    def test_at_capacity_returns_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        c = make_wt_clone(N=100)
        result = strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c))
        assert result == pytest.approx(0.0)

    def test_above_capacity_clamped_to_zero(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        c = make_wt_clone(N=200)
        result = strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c))
        assert result == 0.0

    def test_crowding_between_zero_and_one(self):
        cfg = simple_config(K0=100.0)
        strategy = SimpleCrowding(cfg)
        c = make_wt_clone(N=50)
        result = strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c))
        assert 0.0 <= result <= 1.0

    def test_K_declines_with_time(self):
        cfg = simple_config(K0=100.0, decline=5.0, Kmin=1.0)
        strategy = SimpleCrowding(cfg)
        c = make_clone()
        ts = make_tissue_state(c)
        c0 = strategy.crowding(c, t=0.0, tissue_state=ts)
        c10 = strategy.crowding(c, t=10.0, tissue_state=ts)
        assert c10 <= c0

    def test_K_never_below_Kmin(self):
        cfg = simple_config(K0=100.0, decline=1000.0, Kmin=10.0)
        strategy = SimpleCrowding(cfg)
        c = make_clone()
        result = strategy.crowding(c, t=1000.0, tissue_state=make_tissue_state(c))
        assert result >= 0.0


# ── AdaptedCrowding ────────────────────────────────────────────────────────────

class TestAdaptedCrowding:
    def test_no_logistic_returns_one(self):
        cfg = simple_config(use_logistic=False)
        strategy = AdaptedCrowding(cfg)
        c = make_clone()
        assert strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c)) == 1.0

    def test_result_is_non_negative(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        c = make_clone(birth_rate=0.5, death_rate=0.2)
        result = strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c))
        assert result >= 0.0

    def test_result_at_most_one(self):
        cfg = simple_config(K0=100.0)
        strategy = AdaptedCrowding(cfg)
        c = make_clone(birth_rate=0.5, death_rate=0.2)
        result = strategy.crowding(c, t=0.0, tissue_state=make_tissue_state(c))
        assert result <= 1.0

    def test_higher_fitness_yields_higher_effective_birth_rate(self):
        """Clone with higher birth_rate has higher effective birth rate despite lower crowding."""
        cfg = simple_config(K0=100.0)
        c_low  = make_clone(birth_rate=0.3, death_rate=0.2, N=10)
        c_high = make_clone(birth_rate=0.8, death_rate=0.2, N=10)
        ts_low = make_tissue_state(c_low)
        ts_high = make_tissue_state(c_high)
        rb_low  = c_low.birth_rate_effective(tissue_state=ts_low, crowding=0.5)
        rb_high = c_high.birth_rate_effective(tissue_state=ts_high, crowding=0.5)
        assert rb_high > rb_low