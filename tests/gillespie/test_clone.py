import pytest
from src.gillespie.clone import Clone


# ── fixtures ──────────────────────────────────────────────────────────────────

def make_clone(**kwargs) -> Clone:
    defaults = dict(
        clone_id=(),
        N=10,
        birth_rate=0.5,
        death_rate=0.2,
        mutation_rate=0.01,
        instability=0.0,
    )
    defaults.update(kwargs)
    return Clone(**defaults)


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
    # crowding=1 → effective = 0.5 * 10 * 1
    assert c.birth_rate_effective(crowding=1.0) == pytest.approx(5.0)

def test_birth_rate_effective_zero_crowding():
    c = make_clone(birth_rate=0.5, N=10)
    # crowding=0 → effective = 0.5 * 10 * 0
    assert c.birth_rate_effective(crowding=0.0) == pytest.approx(0.0)

def test_death_rate_effective_proportional_to_N():
    c = make_clone(death_rate=0.2, N=10)
    # effective = 0.2 * 10
    assert c.death_rate_effective() == pytest.approx(2.0)

def test_mutation_rate_effective_uses_multiplier():
    c = make_clone(mutation_rate=0.01, instability=1.0, N=10)
    # multiplier = 1 + instability = 2.0
    assert c.mutation_rate_effective() == pytest.approx(0.01 * 10 * 2.0)


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

# ── mutated_child ─────────────────────────────────────────────────────────────


def test_mutated_child_N_is_one():
    c = make_clone(N=5)
    child = c.mutate(fitness_gain=0.1, instability_jump=0.0, buildup_gain=0.0)
    assert child.N == 1


def test_mutated_child_decrements_parent_N():
    c = make_clone(N=5)
    c.mutate(fitness_gain=0.0, instability_jump=0.0, buildup_gain=0.0)
    assert c.N == 4
    
def test_mutated_child_parent_is_parent_clone_id():
    c = make_clone(clone_id=(3, 1), N=5)
    child = c.mutate(fitness_gain=0.0, instability_jump=0.0, buildup_gain=0.0)
    assert child.parent == (3, 1)

def test_mutated_child_has_higher_birth_rate():
    c = make_clone(birth_rate=0.5, N=5)
    child = c.mutate(fitness_gain=0.1, instability_jump=0.0, buildup_gain=0.0)
    assert child.birth_rate == pytest.approx(0.5 * 1.1)

def test_mutated_child_inherits_death_rate():
    c = make_clone(death_rate=0.3, N=5)
    child = c.mutate(fitness_gain=0.0, instability_jump=0.0, buildup_gain=0.0)
    assert child.death_rate == pytest.approx(0.3)

def test_mutated_child_raises_on_dead_clone():
    c = make_clone(N=0)
    with pytest.raises(ValueError):
        c.mutate(fitness_gain=0.1, instability_jump=0.0, buildup_gain=0.0)

def test_mutated_child_id_increments_and_validates_parent_id():
    c = make_clone(clone_id=(), N=10)
    child1 = c.mutate(fitness_gain=0.0, instability_jump=0.0, buildup_gain=0.0)
    child2 = c.mutate(fitness_gain=0.0, instability_jump=0.0, buildup_gain=0.0)
    assert child1.parent == ()
    assert child2.parent == ()
    assert child1.clone_id == (1,)
    assert child2.clone_id == (2,)


# ── advance_instability ───────────────────────────────────────────────────────

def test_advance_instability_increases_instability():
    c = make_clone(instability=0.0, buildup=0.01)
    c.advance_instability(dt=10.0, base_buildup=0.0)
    assert c.instability == pytest.approx(0.1)

def test_advance_instability_uses_base_buildup():
    c = make_clone(instability=0.0, buildup=0.0)
    c.advance_instability(dt=5.0, base_buildup=0.02)
    assert c.instability == pytest.approx(0.1)
