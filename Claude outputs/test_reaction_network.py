"""The stoichiometry in reaction_network.py must match what the SSA really does."""
import contextlib
import copy
import io

import numpy as np
import pytest

from src.gillespie.event import Event
from src.gillespie.event_type import EventType
from src.gillespie.reaction_network import ReactionNetwork, apply_counts, type_delta
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


@pytest.fixture
def sim():
    """All six types alive (fibroblasts on) and one mutant child clone."""
    config = SimulationConfig(
        seed=1, save_all_steps=False,
        cell_parameters={"mutated": {"N": 300}, "exhausted": {"N": 50},
                         "fibroblast": {"N": 400}, "caf": {"N": 40}},
        kappa=0.002, omega_F=0.5, alpha_F=0.5,
    )
    with contextlib.redirect_stdout(io.StringIO()):
        s = TumorSimulation(config)
        s._introduce_mutation(s.tissue_state.clones[()])
        s._rate_matrix = s._build_rate_matrix()
    s.network = ReactionNetwork(config)
    return s


def pop(s):
    s.tissue_state.update_pop_map()
    return dict(s.tissue_state.pop_map)


def test_type_delta_matches_apply_event(sim):
    """Fire every (type, kind) once through _apply_event and compare with type_delta."""
    for clone in list(sim.tissue_state.clones.values()):
        for kind in EventType:
            try:
                expected = type_delta(clone, kind)
            except ValueError:
                continue  # this type has no such event (e.g. no next_mutation)
            probe = copy.deepcopy(sim)
            c = probe.tissue_state.clones[clone.clone_id]
            c.N = max(c.N, 1)
            before = pop(probe)
            with contextlib.redirect_stdout(io.StringIO()):
                probe._apply_event(Event(kind=kind, clone_id=c.clone_id, rate=1.0))
            after = pop(probe)
            change = {t: after.get(t, 0) - before.get(t, 0) for t in after}
            assert {t: v for t, v in change.items() if v} == {t: v for t, v in expected.items() if v}


def test_apply_counts_once_equals_apply_event(sim):
    """R = one firing of event j gives the same clones (ids, N, children) as the SSA."""
    def state(s):
        return {cid: (c.get_type(), c.N, c.children_count) for cid, c in s.tissue_state.clones.items()}

    arrays = sim.network.to_arrays(sim._current_rate_matrix(), sim.tissue_state)
    for j, event in enumerate(arrays.events):
        ssa, leap = copy.deepcopy(sim), copy.deepcopy(sim)
        with contextlib.redirect_stdout(io.StringIO()):
            ssa._apply_event(event)
            R = np.zeros(len(arrays.events), dtype=int)
            R[j] = 1
            apply_counts(leap, leap.network.to_arrays(leap._current_rate_matrix(), leap.tissue_state), R)
        assert state(ssa) == state(leap)


def test_bulk_counts_change_pop_map_by_R_times_S(sim):
    arrays = sim.network.to_arrays(sim._current_rate_matrix(), sim.tissue_state)
    R = np.where(arrays.consumes, 1, 3)  # 1 exit per channel stays under every n_i here
    before = np.array([pop(sim).get(t, 0) for t in sim.network.types])
    with contextlib.redirect_stdout(io.StringIO()):
        apply_counts(sim, arrays, R)
    after = np.array([pop(sim).get(t, 0) for t in sim.network.types])
    assert (after - before == R @ arrays.S).all()


def test_crowding_pools(sim):
    assert sim.network.pools == {
        "base": ["base", "mutated"],
        "mutated": ["mutated"],
        "immune": ["immune", "exhausted"],
        "exhausted": [],
        "fibroblast": ["mutated", "fibroblast", "caf"],
        "caf": [],
    }
