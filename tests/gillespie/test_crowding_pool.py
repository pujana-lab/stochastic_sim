"""Crowding pools: CellTypeConfig.crowding_pool = {T: w_T}, Y_X = Σ_T w_T · N_T."""
import contextlib
import io

import pytest

from src.gillespie.clone import Clone
from src.gillespie.simulation_config import CellTypeConfig, SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

# The table of docs/NOTATION.md, "Crowding", with every weight 1.
EXPECTED_POOLS = {
    "base": {"base": 1.0, "mutated": 1.0},
    "mutated": {"mutated": 1.0},
    "immune": {"immune": 1.0, "exhausted": 1.0},
    "exhausted": {},
    "fibroblast": {"fibroblast": 1.0, "caf": 1.0, "mutated": 1.0},
    "caf": {},
}


def _sim(fibroblast_pool=None, **kw):
    cp = {"base": {"N": 200}, "mutated": {"N": 20}, "immune": {"N": 50}, "exhausted": {"N": 5},
          "fibroblast": {"N": 80, "K": 0.5, "lambda0": 0.01}, "caf": {"N": 10}}
    if fibroblast_pool is not None:
        cp["fibroblast"]["crowding_pool"] = fibroblast_pool
    cfg = SimulationConfig(seed=1, T_max=50, save_all_steps=False, OMEGA=200, cell_parameters=cp, **kw)
    with contextlib.redirect_stdout(io.StringIO()):
        return TumorSimulation(cfg)


def _pool(**cell_parameters):
    return SimulationConfig(cell_parameters=cell_parameters).cell_parameters


def _fibroblast(sim):
    return sim.tissue_state.clones[sim.config.cell_parameters["fibroblast"].default_id]


# ── defaults ──────────────────────────────────────────────────────────────────

def test_default_pools_match_notation_with_all_weights_one():
    assert {t: cfg.crowding_pool for t, cfg in SimulationConfig().cell_parameters.items()} == EXPECTED_POOLS


def test_every_registered_type_has_a_pool_entry():
    assert set(Clone._registry) <= set(SimulationConfig().cell_parameters)


def test_crowding_numerator_sums_the_pool():
    sim = _sim()
    pop = sim.tissue_state.pop_map
    for clone in sim.tissue_state.clones.values():
        assert clone.crowding_numerator(sim.tissue_state) == sum(pop[t] for t in EXPECTED_POOLS[clone.get_type()])


@pytest.mark.parametrize("clone_type", ["exhausted", "caf"])
@pytest.mark.parametrize("decay", [False, True])
def test_non_dividing_types_have_no_crowding(clone_type, decay):
    """Replaces the old `"exhausted"` special case in CrowdingStrategy.crowding."""
    sim = _sim(decay=decay)
    clone = sim.tissue_state.clones[sim.config.cell_parameters[clone_type].default_id]
    assert sim.config.crowding_strategy.crowding(clone, sim.tissue_state) == 1.0


# ── merge with the defaults (defaults -> JSON) ───────────────────────────────

def test_json_pool_changes_only_the_weights_it_names():
    pool = _pool(fibroblast={"crowding_pool": {"mutated": 0.5}})["fibroblast"].crowding_pool
    assert pool == {"fibroblast": 1.0, "caf": 1.0, "mutated": 0.5}


def test_zero_weight_removes_a_member():
    pool = _pool(fibroblast={"crowding_pool": {"caf": 0}})["fibroblast"].crowding_pool
    assert pool == {"fibroblast": 1.0, "mutated": 1.0}


def test_json_pool_can_add_a_member():
    pool = _pool(fibroblast={"crowding_pool": {"immune": 0.3}})["fibroblast"].crowding_pool
    assert pool == {"fibroblast": 1.0, "caf": 1.0, "mutated": 1.0, "immune": 0.3}


def test_other_fields_do_not_touch_the_pool():
    assert _pool(fibroblast={"N": 10})["fibroblast"].crowding_pool == EXPECTED_POOLS["fibroblast"]


def test_a_full_celltypeconfig_is_taken_as_is():
    cfg = CellTypeConfig(default_id=(-4,), K=1, crowding_pool={"fibroblast": 2})
    assert _pool(fibroblast=cfg)["fibroblast"].crowding_pool == {"fibroblast": 2.0}


def test_pool_survives_omega_scaling():
    cfg = SimulationConfig(OMEGA=500, cell_parameters={"fibroblast": {"crowding_pool": {"mutated": 0.3}}})
    assert cfg.cell_parameters["fibroblast"].crowding_pool["mutated"] == 0.3


@pytest.mark.parametrize("pool, match", [
    ({"mutatd": 1}, "unknown cell types"),
    ({"mutated": -1.0}, ">= 0"),
])
def test_invalid_pools_are_rejected(pool, match):
    with pytest.raises(ValueError, match=match):
        _pool(fibroblast={"crowding_pool": pool})


# ── weights in the rates ──────────────────────────────────────────────────────

def test_weighted_pool_no_reference_type():
    """Every member has its own weight, the niche's own type included."""
    sim = _sim({"fibroblast": 0.8, "caf": 1.5, "mutated": 0.5})
    pop = sim.tissue_state.pop_map
    expected = 0.8 * pop["fibroblast"] + 1.5 * pop["caf"] + 0.5 * pop["mutated"]
    assert _fibroblast(sim).crowding_numerator(sim.tissue_state) == pytest.approx(expected)


def test_weights_enter_the_birth_rate():
    sim = _sim({"mutated": 0})
    fn, pop = _fibroblast(sim), sim.tissue_state.pop_map
    g = max(0.0, 1 - (pop["fibroblast"] + pop["caf"]) / fn.actual_K)
    assert fn.birth_rate_effective(sim.tissue_state) == pytest.approx(fn.birth_rate * fn.N * g)
