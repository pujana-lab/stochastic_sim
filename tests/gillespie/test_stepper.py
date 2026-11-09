"""Stepper interface (PLAN_tau_leaping.md §7, Phase 1)."""
import contextlib
import io

import pytest

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.stepper import SSAStepper, Stepper, make_stepper
from src.gillespie.tumor_simulation import TumorSimulation
from tests.gillespie.test_tumor_simulation_rates import _FIBROBLAST_REGRESSION_CONFIG


def _sim(**kwargs) -> TumorSimulation:
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(SimulationConfig(**kwargs))


def _state(sim: TumorSimulation):
    """Everything a step can change, to compare two runs exactly."""
    return (
        sim.step_count,
        list(sim.times),
        dict(sim.tissue_state.pop_map),
        [(e.kind, e.clone_id) for e in sim.events],
        len(sim.history),
        len(sim.rate_history),
        sorted(sim.tissue_state.clones),
    )


def test_stepper_is_abstract():
    with pytest.raises(TypeError):
        Stepper()


@pytest.mark.parametrize("config", [
    dict(seed=7, T_max=300, OMEGA=2000, cell_parameters={"mutated": {"N": 50}}),   # mutations on
    dict(seed=11, T_max=3, **_FIBROBLAST_REGRESSION_CONFIG),                        # fibroblasts on
], ids=["mutations", "fibroblasts"])
def test_ssa_stepper_matches_simulation_step(config):
    """SSAStepper().step(sim) advances exactly as TumorSimulation.step() does."""
    by_sim, by_stepper = _sim(**config), _sim(**config)
    stepper = SSAStepper()
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(500):
            go_sim, go_stepper = by_sim.step(), stepper.step(by_stepper)
            assert go_sim == go_stepper
            if not go_sim:
                break
    assert _state(by_stepper) == _state(by_sim)


def test_run_burst_stops_after_n_steps():
    sim = _sim(seed=3, T_max=1e9, verbose=False)
    assert SSAStepper().run_burst(sim, 25) is True
    assert sim.step_count == 25
    assert len(sim.times) == 1 + 25


def test_run_burst_stops_at_t_max():
    sim = _sim(seed=3, T_max=1.0, verbose=False)
    with contextlib.redirect_stdout(io.StringIO()):
        assert SSAStepper().run_burst(sim, 10_000) is False
    assert sim.t == 1.0
    assert sim.step_count < 10_000


def test_run_burst_stops_when_no_event_is_left():
    """All rates zero (every type at N = 0): the burst stops at once."""
    empty = {t: {"N": 0} for t in ("base", "mutated", "immune", "exhausted", "fibroblast", "caf")}
    sim = _sim(seed=3, cell_parameters=empty, verbose=False)
    assert SSAStepper().run_burst(sim, 10) is False
    assert sim.step_count == 1
    assert sim.times == [0.0]


# ── Stepper selection (Phase 1.3) ────────────────────────────────────────────

def test_default_stepper_is_ssa():
    assert SimulationConfig().stepper == "ssa"
    assert isinstance(_sim(seed=1).stepper, SSAStepper)


def test_unknown_stepper_is_rejected_by_the_config():
    with pytest.raises(ValueError, match="Unknown stepper 'sas'"):
        SimulationConfig(stepper="sas")


def test_make_stepper():
    assert isinstance(make_stepper("ssa"), SSAStepper)
    with pytest.raises(ValueError, match="known: \\['ssa'\\]"):
        make_stepper("tau")


def test_stepper_is_saved_in_the_run_config(tmp_path):
    import json
    from src.gillespie.infrastructure.run_config_output import save_run_config
    _, path = save_run_config(tmp_path / "config", SimulationConfig(seed=1))
    saved = json.loads(path.read_text())
    assert saved["input"]["stepper"] == saved["used"]["stepper"] == "ssa"
