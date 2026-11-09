import contextlib
import io
from collections import Counter, defaultdict

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
    keys = {"BIRTH": "rb", "DEATH": "rd", "MUTATION": "rm", "CONVERSION": "re"}
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


def test_one_root_clone_per_cell_type_keyed_by_default_id():
    sim = _quiet_sim()
    roots = {cid: c.get_type() for cid, c in sim.tissue_state.clones.items()}
    assert roots == {cfg.default_id: t for t, cfg in sim.config.cell_parameters.items()}


def test_seeded_trajectory_is_unchanged():
    """Regression: a seeded run gives the trajectory recorded on 03-11-2026.

    It depends on the RNG stream and on the event order in the rate matrix (the
    order of cell_parameters; the mutated root starts alive so that its position
    matters). If a change is meant to alter the dynamics, record the new values
    here in the same commit.
    """
    config = SimulationConfig(seed=7, T_max=300, OMEGA=2000, cell_parameters={"mutated": {"N": 50}})
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        times, _, tissue_state, _ = TumorSimulation(config).run()
    assert len(times) == 2583
    final = {t: tissue_state.pop_map[t] for t in ("base", "mutated", "immune", "exhausted")}
    assert final == {"base": 2480, "mutated": 190, "immune": 1901, "exhausted": 93}
    mutant_clones = len(tissue_state.clones) - len(config.cell_parameters)  # minus one root per type
    assert mutant_clones == 26


# Parameters of configs/fibroblasts.json on 09-11-2026, copied so that tuning the
# example config does not break the regression below.
_FIBROBLAST_REGRESSION_CONFIG = dict(
    cell_parameters={
        "base": {"N": 0, "K": 1.0, "lambda0": 1.0, "mu": 0.25, "nu": 0.0, "next_mutation": "mutated"},
        "immune": {"N": 100, "K": 1.0, "lambda0": 0.1, "omega_exhaust": 0.1, "mu": 0.0},
        "mutated": {"N": 800, "K": 1.0, "lambda0": 1.0, "mu": 0.1, "fitness_gain": 0.0},
        "exhausted": {"N": 0, "mu": 0.4, "lambda0": 0.0},
        "fibroblast": {"N": 2000, "K": 2.0, "lambda0": 0.1, "mu": 0.0},
        "caf": {"N": 0, "mu": 0.05, "lambda0": 0.0},
    },
    OMEGA=1000, theta_I=0.8, beta=0.45,
    kappa=0.2, C_half=0.5, omega_F=0.05, alpha_F=0.25,
    use_logistic_adapted=False, verbose=False,
)


def test_seeded_fibroblast_trajectory_is_unchanged():
    """Regression with the fibroblast module on, recorded on 09-11-2026 before tau-leaping Phase 1.

    Covers what the test above does not: fibroblast -> caf activation (κ, C_½),
    exhaustion by caf (ω_F), caf support of tumour birth (α_F), caf death and the
    six-type event order. If a change is meant to alter the dynamics, record the
    new values here in the same commit.
    """
    config = SimulationConfig(seed=11, T_max=3, **_FIBROBLAST_REGRESSION_CONFIG)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        sim = TumorSimulation(config)
        times, _, tissue_state, _ = sim.run()
    assert len(times) == 2014
    assert times[-1] == 3
    # The waiting times depend on the exact total rate of every step, so the last
    # event time also catches rate changes too small to flip any event choice.
    assert times[-2] == pytest.approx(2.9997935745704, rel=1e-12)
    assert dict(tissue_state.pop_map) == {
        "base": 0, "mutated": 752, "immune": 248, "exhausted": 43, "fibroblast": 1397, "caf": 563,
    }
    assert Counter((e.kind.name, e.clone_type) for e in sim.events) == {
        ("BIRTH", "mutated"): 515, ("DEATH", "mutated"): 563,
        ("BIRTH", "immune"): 210, ("CONVERSION", "immune"): 62,
        ("DEATH", "exhausted"): 19,
        ("CONVERSION", "fibroblast"): 603, ("DEATH", "caf"): 40,
    }
