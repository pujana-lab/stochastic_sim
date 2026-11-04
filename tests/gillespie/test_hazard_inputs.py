"""Phase 0.4 of docs/PLAN_tau_leaping.md: Clone.hazard_inputs must list exactly the type
totals that each type's rates read outside its crowding pool.

The dependencies are found numerically, not from the source: with crowding off
(use_logistic=False, so the pools have no effect) and every interaction parameter
non-zero, change one type total at a time and see which rates of each type move.
A new interaction that is not declared in hazard_inputs makes this test fail.
"""
import contextlib
import io

import pytest

from src.gillespie.clone import Clone
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

RATE_METHODS = ("birth_rate_effective", "death_rate_effective",
                "mutation_rate_effective", "conversion_rate_effective")


def _sim():
    cp = {"base": {"N": 200}, "mutated": {"N": 20}, "immune": {"N": 50}, "exhausted": {"N": 5, "mu": 0.01},
          "fibroblast": {"N": 80}, "caf": {"N": 10}}
    cfg = SimulationConfig(seed=1, OMEGA=200, save_all_steps=False, use_logistic=False, cell_parameters=cp,
                           kappa=0.05, omega_F=0.002, alpha_F=1.0)
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(cfg)


def _rates(clone, tissue_state):
    return [getattr(clone, m)(tissue_state=tissue_state) for m in RATE_METHODS]


def _measured_inputs(clone_type: str) -> set[str]:
    sim = _sim()
    ts = sim.tissue_state
    clone = ts.clones[sim.config.cell_parameters[clone_type].default_id]
    base_rates = _rates(clone, ts)
    found = set()
    for other in Clone._registry:
        original = ts.pop_map.get(other, 0)
        ts.pop_map[other] = original + 7
        if _rates(clone, ts) != pytest.approx(base_rates):
            found.add(other)
        ts.pop_map[other] = original
    return found


@pytest.mark.parametrize("clone_type", sorted(Clone._registry))
def test_hazard_inputs_match_the_rates(clone_type):
    declared = set(Clone._registry[clone_type].hazard_inputs)
    measured = _measured_inputs(clone_type)
    assert declared == measured, (
        f"{clone_type}.hazard_inputs = {sorted(declared)}, but its rates depend on {sorted(measured)}"
    )


def test_hazard_inputs_are_registered_types():
    for t, cls in Clone._registry.items():
        assert set(cls.hazard_inputs) <= set(Clone._registry), t
