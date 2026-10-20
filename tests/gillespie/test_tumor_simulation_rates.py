import contextlib
import io
from collections import defaultdict

import pytest

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


def _quiet_sim(**kwargs) -> TumorSimulation:
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(SimulationConfig(seed=3, **kwargs))


def _with_mutant_clones(sim: TumorSimulation, sizes=(50, 120, 300)) -> TumorSimulation:
    root = sim.tissue_state.clones[()]
    for n in sizes:
        sim.create_clone(clone_id=root.next_child_id(), clone_type="mutated", N=n, parent=())
    return sim


def _rates_by_clone(rate_matrix) -> dict:
    keys = {"BIRTH": "rb", "DEATH": "rd", "MUTATION": "rm", "EXHAUSTION": "re"}
    rates = defaultdict(lambda: {"rb": 0.0, "rd": 0.0, "rm": 0.0, "re": 0.0})
    for e in rate_matrix.events:
        rates[e.clone_id][keys[e.kind.name]] = e.rate
    return rates


def test_snapshot_reports_effective_rates_of_rate_matrix():
    sim = _with_mutant_clones(_quiet_sim(save_interval=1))
    rate_matrix = sim._build_rate_matrix()
    expected = _rates_by_clone(rate_matrix)
    snap = sim.tissue_state.snapshot()

    for cid, values in snap.items():
        for key in ("rb", "rd", "rm", "re"):
            assert values[key] == pytest.approx(expected[cid][key]), (cid, key)


def test_per_clone_rates_sum_to_type_level_formula():
    sim = _with_mutant_clones(_quiet_sim())
    rate_matrix = sim._build_rate_matrix()
    cfg, pop = sim.config, sim.tissue_state.pop_map
    params = cfg.cell_parameters["mutated"]
    K = sim.tissue_state.clones[(1,)].actual_K

    birth = sum(e.rate for e in rate_matrix.events if e.kind.name == "BIRTH" and e.clone_type == "mutated")
    death = sum(e.rate for e in rate_matrix.events if e.kind.name == "DEATH" and e.clone_type == "mutated")

    assert birth == pytest.approx(params.lambda0 * pop["mutated"] * max(0.0, 1 - pop["mutated"] / K))
    assert death == pytest.approx(params.mu * pop["mutated"] + pop["mutated"] * pop["immune"] * cfg.theta_I)


def test_snapshots_during_run_match_rates_of_their_state():
    sim = _with_mutant_clones(_quiet_sim(save_interval=1, T_max=5))
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(50):
            if not sim.step():
                break
    last = sim.history[-1]
    expected = _rates_by_clone(sim._build_rate_matrix())

    for cid, values in last.items():
        assert values["rb"] == pytest.approx(expected[cid]["rb"]), cid
        assert values["rd"] == pytest.approx(expected[cid]["rd"]), cid
    assert last[()]["t"] == sim.t


def test_dead_clone_reports_zero_rates():
    sim = _quiet_sim()
    snap = sim.tissue_state.snapshot()
    # The exhausted pool starts empty
    assert sim.tissue_state.clones[(-2,)].N == 0
    assert snap[(-2,)]["rb"] == snap[(-2,)]["rd"] == 0.0
