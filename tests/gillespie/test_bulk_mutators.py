"""Phase 0.5 of docs/PLAN_tau_leaping.md: bulk mutators.

`_apply_event(event, k)` must leave the same tissue as k calls of
`_apply_event(event)`, and k = 1 must stay the SSA path."""
import contextlib
import io

import pytest

from src.gillespie.clone import Clone
from src.gillespie.event_type import EventType
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

K = 3


def _sim():
    """All six types alive and every interaction on, so every channel has an event."""
    cp = {"base": {"N": 200}, "mutated": {"N": 20}, "immune": {"N": 50}, "exhausted": {"N": 5, "mu": 0.01},
          "fibroblast": {"N": 80, "K": 0.5, "lambda0": 0.01}, "caf": {"N": 10, "mu": 0.003}}
    cfg = SimulationConfig(seed=1, T_max=50, save_all_steps=False, OMEGA=200, cell_parameters=cp,
                           kappa=0.05, omega_F=0.002, alpha_F=1.0)
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(cfg)


def _event(sim, kind, clone_id):
    return next(e for e in sim._build_rate_matrix().events if (e.kind, e.clone_id) == (kind, clone_id))


def _state(sim):
    """Every clone with its N, parent and child counter: a bulk mutation must create the same clones."""
    return {cid: (c.get_type(), c.N, c.parent, c.children_count) for cid, c in sim.tissue_state.clones.items()}


def _events():
    return [(e.kind, e.clone_id) for e in _sim()._build_rate_matrix().events]


def test_every_event_kind_is_exercised():
    assert {kind for kind, _ in _events()} == set(EventType)


@pytest.mark.parametrize("kind, clone_id", _events())
def test_bulk_event_equals_repeated_single_events(kind, clone_id):
    bulk, single = _sim(), _sim()
    bulk._apply_event(_event(bulk, kind, clone_id), K)
    event = _event(single, kind, clone_id)
    for _ in range(K):
        single._apply_event(event)
    assert _state(bulk) == _state(single)


@pytest.mark.parametrize("kind, clone_id", _events())
def test_bulk_event_changes_type_totals_by_k_type_deltas(kind, clone_id):
    sim = _sim()
    delta = sim.tissue_state.clones[clone_id].type_delta(kind)
    before = dict(sim.tissue_state.pop_map)
    sim._apply_event(_event(sim, kind, clone_id), K)
    sim.tissue_state.update_pop_map()
    after = sim.tissue_state.pop_map
    observed = {t: after.get(t, 0) - before.get(t, 0) for t in set(before) | set(after)}
    assert {t: d for t, d in observed.items() if d} == {t: K * d for t, d in delta.items()}


def test_zero_count_is_a_no_op():
    sim = _sim()
    before = _state(sim)
    for kind, cid in _events():
        sim._apply_event(_event(sim, kind, cid), 0)
    assert _state(sim) == before


def test_negative_count_is_rejected():
    sim = _sim()
    kind, cid = _events()[0]
    with pytest.raises(ValueError):
        sim._apply_event(_event(sim, kind, cid), -1)


@pytest.mark.parametrize("kind", [EventType.DEATH, EventType.MUTATION, EventType.CONVERSION])
def test_consuming_more_cells_than_the_clone_has_is_rejected(kind):
    sim = _sim()
    event = next(e for e in sim._build_rate_matrix().events if e.kind == kind)
    before = _state(sim)
    n = sim.tissue_state.clones[event.clone_id].N
    with pytest.raises(ValueError):
        sim._apply_event(event, n + 1)
    assert _state(sim) == before


def test_birth_may_exceed_n():
    sim = _sim()
    event = next(e for e in sim._build_rate_matrix().events if e.kind == EventType.BIRTH)
    clone = sim.tissue_state.clones[event.clone_id]
    n = clone.N
    sim._apply_event(event, 2 * n)
    assert clone.N == 3 * n


def test_divide_and_kill_take_k():
    clone = next(iter(_sim().tissue_state.clones.values()))
    n = clone.N
    clone.divide(5)
    assert clone.N == n + 5
    clone.kill(4)
    assert clone.N == n + 1
    clone.kill(n + 10)
    assert clone.N == 0


@pytest.mark.parametrize("method", [Clone.divide, Clone.kill])
def test_divide_and_kill_reject_negative_k(method):
    clone = next(iter(_sim().tissue_state.clones.values()))
    with pytest.raises(ValueError):
        method(clone, -1)
