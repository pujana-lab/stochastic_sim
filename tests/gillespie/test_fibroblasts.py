import contextlib
import io

import pytest

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


def _sim(cell_parameters=None, **kw):
    cp = {"base": {"N": 200}, "mutated": {"N": 20}, "immune": {"N": 50}, "exhausted": {"N": 0},
          "fibroblast": {"N": 80, "K": 0.5, "lambda0": 0.01}, "caf": {"N": 10, "mu": 0.003}}
    cp.update(cell_parameters or {})
    cfg = SimulationConfig(seed=1, T_max=50, save_all_steps=False, OMEGA=200, cell_parameters=cp, **kw)
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(cfg)


def _rate(sim, kind, clone_type):
    return sum(e.rate for e in sim._build_rate_matrix().events
               if e.kind.name == kind and e.clone_type == clone_type)


def test_fibroblast_off_by_default_keeps_4_type_model():
    cfg = SimulationConfig()
    assert cfg.cell_parameters["fibroblast"].N == 0 and cfg.cell_parameters["caf"].N == 0
    assert cfg.kappa == cfg.omega_F == cfg.alpha_F == 0.0


def test_fibroblast_logistic_birth_counts_fn_ft_and_c():
    sim = _sim()
    cfg, pop = sim.config, sim.tissue_state.pop_map
    p = cfg.cell_parameters["fibroblast"]
    K = sim.tissue_state.clones[p.default_id].actual_K
    expected = p.lambda0 * pop["fibroblast"] * max(0.0, 1 - (pop["fibroblast"] + pop["caf"] + pop["mutated"]) / K)
    assert _rate(sim, "BIRTH", "fibroblast") == pytest.approx(expected)


def test_activation_is_michaelis_menten_in_cancer():
    sim = _sim(kappa=0.05, C_half=0.1)
    cfg, pop = sim.config, sim.tissue_state.pop_map
    assert cfg.C_half == pytest.approx(0.1 * cfg.OMEGA)
    expected = cfg.kappa * pop["fibroblast"] * pop["mutated"] / (cfg.C_half + pop["mutated"])
    assert _rate(sim, "CONVERSION", "fibroblast") == pytest.approx(expected)


def test_activation_event_moves_one_fn_to_ft():
    sim = _sim(kappa=0.05)
    cp, clones = sim.config.cell_parameters, sim.tissue_state.clones
    fn, ft = clones[cp["fibroblast"].default_id], clones[cp["caf"].default_id]
    n_fn, n_ft = fn.N, ft.N
    sim._induce_conversion(fn)
    assert (fn.N, ft.N) == (n_fn - 1, n_ft + 1)


def test_conversion_without_target_is_rejected():
    sim = _sim()
    with pytest.raises(ValueError, match="converts_to"):
        sim._induce_conversion(sim.tissue_state.clones[sim.config.cell_parameters["caf"].default_id])


def test_caf_does_not_divide_and_dies_at_d_F():
    sim = _sim()
    assert _rate(sim, "BIRTH", "caf") == 0.0
    assert _rate(sim, "DEATH", "caf") == pytest.approx(0.003 * sim.tissue_state.pop_map["caf"])


def test_caf_induced_exhaustion_adds_to_cancer_induced():
    sim0, sim1 = _sim(), _sim(omega_F=0.002)
    pop = sim1.tissue_state.pop_map
    extra = sim1.config.omega_F * pop["immune"] * pop["caf"]
    assert _rate(sim1, "CONVERSION", "immune") == pytest.approx(_rate(sim0, "CONVERSION", "immune") + extra)


def test_caf_boosts_cancer_birth():
    sim0, sim1 = _sim(), _sim(alpha_F=1.0)
    pop = sim1.tissue_state.pop_map
    factor = 1 + sim1.config.alpha_F * pop["caf"]
    assert _rate(sim1, "BIRTH", "mutated") == pytest.approx(_rate(sim0, "BIRTH", "mutated") * factor)
