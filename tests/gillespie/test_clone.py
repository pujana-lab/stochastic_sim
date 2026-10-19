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
        d1=0.0,
        d2=0.0,
    )
    defaults.update(kwargs)
    
    config = SimulationConfig(cell_parameters={
        "base": {
            "N": defaults["N"],
            "lambda0": defaults["birth_rate"],
            "mu": defaults["death_rate"],
            "nu": defaults["mutation_rate"],
            "K": 100,
        },
        "mutated": {"K": 100},
        "exhausted": {"K": 100},
        "immune": {"K": 100},
    }, use_logistic=False)
    clone = WildTypeClone(clone_id=defaults["clone_id"], N=defaults["N"], config=config)
    for attr in ("instability", "buildup", "d1", "d2"):
        setattr(clone, attr, defaults[attr])
    return clone

def make_tissue_state(clone: Clone) -> TissueState:
    return TissueState(clones={clone.clone_id: clone}, t=0.0)

def make_clone_factory():
    config = SimulationConfig(
        cell_parameters={
            "base": {"N": 10, "lambda0": 0.5, "mu": 0.2, "nu": 0.01},
            "mutated": {"N": 1, "fitness_gain": 0.1, "K": 100},
        },
        instability_0=0.0,
        buildup_0=0.0,
        d1_0=0.0,
        d2_0=0.0,
        mutation_instability_jump=0.05,
        mutation_buildup_gain=0.02,
        seed=42
    )
    return CloneFactory(config  = config)

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
    c = make_clone(birth_rate=0.5, N=10)
    ts = make_tissue_state(c)
    # effective = birth_rate * N * crowding = 0.5 * 10 * 1
    assert c.birth_rate_effective(tissue_state=ts) == pytest.approx(5.0)

def test_birth_rate_effective_zero_crowding():
    c = make_clone(birth_rate=0.5, N=10)
    ts = make_tissue_state(c)
    # effective = birth_rate * N * 0
    c.config.crowding_strategy.crowding = lambda *_args, **_kwargs: 0.0
    assert c.birth_rate_effective(tissue_state=ts) == pytest.approx(0.0)

def test_death_rate_effective_proportional_to_N():
    c = make_clone(death_rate=0.2, N=10)
    ts = make_tissue_state(c)
    # effective = death_rate * N = 0.2 * 10
    assert c.death_rate_effective(tissue_state=ts) == pytest.approx(2.0)

def test_mutation_rate_effective_scales_with_N_and_multiplier():
    c = make_clone(mutation_rate=0.01, instability=1.0, N=10)
    c.mutation_rate = 0.01
    ts = make_tissue_state(c)
    # effective = mutation_rate * N * (1 + instability) = 0.01 * 10 * 2.0
    assert c.mutation_rate_effective(tissue_state=ts) == pytest.approx(0.2)


# # ── mutation_multiplier ───────────────────────────────────────────────────────

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
    clone = factory.create_clone(clone_id, clone_type="mutated")
    assert clone.clone_id == clone_id
    assert clone.N == 1
    assert clone.birth_rate == pytest.approx(0.0055)
    assert clone.death_rate == 0.002
    assert clone.mutation_rate == 0.0
    assert clone.instability == 0.0
    assert clone.get_type() == "mutated"