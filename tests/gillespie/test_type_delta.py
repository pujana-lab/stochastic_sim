"""Phase 0.1 of docs/PLAN_tau_leaping.md: the declarative stoichiometry
(EventType.consumes_parent, Clone.type_delta) must match what the SSA does."""
import contextlib
import io

import pytest

from src.gillespie.clone import MUTANT_SEED_SIZE
from src.gillespie.event_type import EventType
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


def _sim():
    """All six types alive and every interaction on, so every channel has an event."""
    cp = {"base": {"N": 200}, "mutated": {"N": 20}, "immune": {"N": 50}, "exhausted": {"N": 5, "mu": 0.01},
          "fibroblast": {"N": 80, "K": 0.5, "lambda0": 0.01}, "caf": {"N": 10, "mu": 0.003}}
    cfg = SimulationConfig(seed=1, T_max=50, save_all_steps=False, OMEGA=200, cell_parameters=cp,
                           kappa=0.05, omega_F=0.002, alpha_F=1.0)
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(cfg)


def _events():
    return [(e.kind, e.clone_id) for e in _sim()._build_rate_matrix().events]


def test_consumes_parent_is_false_only_for_birth():
    assert {k: k.consumes_parent for k in EventType} == {
        EventType.BIRTH: False, EventType.DEATH: True,
        EventType.MUTATION: True, EventType.CONVERSION: True,
    }


def test_every_event_kind_is_exercised():
    assert {kind for kind, _ in _events()} == set(EventType)


@pytest.mark.parametrize("kind, clone_id", _events())
def test_type_delta_matches_apply_event(kind, clone_id):
    sim = _sim()
    event = next(e for e in sim._build_rate_matrix().events if (e.kind, e.clone_id) == (kind, clone_id))
    clone = sim.tissue_state.clones[clone_id]
    expected = clone.type_delta(kind)

    before = dict(sim.tissue_state.pop_map)
    n_clone = clone.N
    sim._apply_event(event)
    sim.tissue_state.update_pop_map()
    after = sim.tissue_state.pop_map

    observed = {t: after.get(t, 0) - before.get(t, 0) for t in set(before) | set(after)}
    assert {t: d for t, d in observed.items() if d} == expected
    assert clone.N == n_clone - (1 if kind.consumes_parent else -1)


def test_mutation_delta_uses_seed_size():
    sim = _sim()
    base = sim.tissue_state.clones[sim.config.cell_parameters["base"].default_id]
    assert base.type_delta(EventType.MUTATION) == {"base": -1, "mutated": MUTANT_SEED_SIZE}


def test_mutation_into_same_type_sums_deltas():
    sim = _sim()
    mutated = sim.tissue_state.clones[sim.config.cell_parameters["mutated"].default_id]
    mutated.next_mutation = "mutated"   # e.g. a future mutated -> mutated subclone
    assert mutated.type_delta(EventType.MUTATION) == {"mutated": MUTANT_SEED_SIZE - 1}


@pytest.mark.parametrize("clone_type, kind, attr", [
    ("caf", EventType.CONVERSION, "converts_to"),
    ("immune", EventType.MUTATION, "next_mutation"),
])
def test_missing_target_is_rejected(clone_type, kind, attr):
    sim = _sim()
    clone = sim.tissue_state.clones[sim.config.cell_parameters[clone_type].default_id]
    with pytest.raises(ValueError, match=attr):
        clone.type_delta(kind)
