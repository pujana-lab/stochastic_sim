import pytest
from src.gillespie.clone import Clone, WildTypeClone
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory
from src.gillespie.tissue_state import TissueState

# ── fixtures ──────────────────────────────────────────────────────────────────

def make_clone(**kwargs) -> Clone:
    defaults = dict(
        clone_id=(),
        N=10,
        birth_rate=0.5,
        death_rate=0.2,
        mutation_rate=0.01,
        instability=0.0,
        buildup=0.0,
    )
    defaults.update(kwargs)
    
    config = SimulationConfig(
        lambda0=defaults["birth_rate"], 
        mu0=defaults["death_rate"], 
        nu0=defaults["mutation_rate"]
    )
    clone = Clone(clone_id=defaults["clone_id"], N=defaults["N"], config=config)
    clone.instability = defaults["instability"]
    clone.buildup = defaults["buildup"]
    return clone


def make_wt_clone_with_state(birth_rate=0.5, death_rate=0.2, mutation_rate=0.01,
                               instability=0.0, N=10):
    """Create a WildTypeClone and a TissueState containing it."""
    config = SimulationConfig(lambda0=birth_rate, mu0=death_rate, nu0=mutation_rate)
    clone = WildTypeClone(clone_id=(), N=N, config=config)
    clone.mutation_rate = mutation_rate
    clone.instability = instability
    ts = TissueState(clones={(): clone})
    return clone, ts


def make_clone_factory():
    config = SimulationConfig(
        N0=10,
        lambda0=0.5,
        mu0=0.2,
        nu0=0.01,
        instability_0=0.0,
        buildup_0=0.0,
        d1_0=0.0,
        d2_0=0.0,
        fitness_gain=0.1,
        mutation_instability_jump=0.05,
        mutation_buildup_gain=0.02,
        seed=42
    )
    return CloneFactory(config=config)

# ── is_alive ──────────────────────────────────────────────────────────────────

def test_is_alive_when_N_positive():
    assert make_clone(N=5).is_alive()



def test_not_alive_when_N_zero():
    assert not make_clone(N=0).is_alive()


# ── divide / die ──────────────────────────────────────────────────────────────

def test_divide_increments_N():
    c = make_clone(N=5)
    c.divide()
    assert c.N == 6

def test_kill_decrements_N():
    c = make_clone(N=5)
    c.kill()
    assert c.N == 4

def test_kill_does_not_go_below_zero():
    c = make_clone(N=0)
    c.kill()
    assert c.N == 0


# ── effective rates ───────────────────────────────────────────────────────────

def test_birth_rate_effective_proportional_to_N():
    c, ts = make_wt_clone_with_state(birth_rate=0.5, N=10)
    # crowding=1.0 → effective = 0.5 * 10 * 1.0 = 5.0
    assert c.birth_rate_effective(tissue_state=ts, crowding=1.0) == pytest.approx(5.0)

def test_birth_rate_effective_zero_crowding():
    c, ts = make_wt_clone_with_state(birth_rate=0.5, N=10)
    # crowding=0.0 → effective = 0.5 * 10 * 0.0 = 0.0
    assert c.birth_rate_effective(tissue_state=ts, crowding=0.0) == pytest.approx(0.0)

def test_death_rate_effective_proportional_to_N():
    c, ts = make_wt_clone_with_state(death_rate=0.2, N=10)
    # effective = 0.2 * 10 = 2.0
    assert c.death_rate_effective(tissue_state=ts) == pytest.approx(2.0)

def test_mutation_rate_effective_includes_multiplier():
    c, ts = make_wt_clone_with_state(mutation_rate=0.01, instability=1.0, N=10)
    # multiplier = 1 + 1.0 = 2.0; effective = 0.01 * 10 * 2.0 = 0.2
    assert c.mutation_rate_effective(tissue_state=ts) == pytest.approx(0.2)


# ── mutation_multiplier ───────────────────────────────────────────────────────

def test_mutation_multiplier_no_instability():
    # instability=0 → multiplier = 1 + 0 = 1
    assert make_clone(instability=0.0).mutation_multiplier() == 1.0

def test_mutation_multiplier_with_instability():
    # instability=0.5 → multiplier = 1 + 0.5 = 1.5
    assert make_clone(instability=0.5).mutation_multiplier() == pytest.approx(1.5)

# ── next_child_id ───────────────────────────────────────────────────────────
def test_next_child_id_increments():
    c = make_clone(clone_id=())
    assert c.next_child_id() == (1,)
    assert c.next_child_id() == (2,)
    assert c.next_child_id() == (3,)

def test_next_child_id_does_not_modify_clone_id():
    c = make_clone(clone_id=(42,))
    assert c.next_child_id() == (42, 1)
    assert c.clone_id == (42,)

def test_next_child_id_increments_children_count():
    c = make_clone()
    assert c.children_count == 0
    assert c.next_child_id() == (1,)
    assert c.children_count == 1
    assert c.next_child_id() == (2,)
    assert c.children_count == 2
    assert c.next_child_id() == (3,)


# ── advance_instability ───────────────────────────────────────────────────────

def test_advance_instability_increases_instability():
    c = make_clone(instability=0.0, buildup=0.01)
    c.advance_instability(dt=10.0, base_buildup=0.0)
    assert c.instability == pytest.approx(0.1)

def test_advance_instability_uses_base_buildup():
    c = make_clone(instability=0.0, buildup=0.0)
    c.advance_instability(dt=5.0, base_buildup=0.02)
    assert c.instability == pytest.approx(0.1)

def test_clone_factory_create_clone():
    factory = make_clone_factory()
    clone_id = (1,)
    clone = factory.create_clone(clone_id, clone_type="base", N=10)
    assert clone.clone_id == clone_id
    assert clone.N == 10
    assert clone.birth_rate == 0.5
    assert clone.death_rate == 0.2
    assert clone.mutation_rate == 0.01
    assert clone.instability == 0.0
    assert str(clone) == "base"
